classdef ResolutionReportTest < matlab.unittest.TestCase
%RESOLUTIONREPORTTEST  The run-start check on the baseline against the ruler.
%
%   Layer 1 of the project's test plan. The behaviour under test is a single asymmetry: a
%   baseline finer than the discretization can resolve is an ERROR in
%   certification-driven mode, because the run would claim what it cannot support,
%   and a REPORT in budget-driven mode, because nothing is certified there.

    properties (Constant)
        Tol = 1e-12
    end

    methods (Test)

        function reportsHeadroomAndTheCoverageBound(tc)
            report = rid.crit.resolutionReport(0.2, 0.01, "certification");
            tc.verifyEqual(report.headroom, 20, "AbsTol", tc.Tol);
            tc.verifyEqual(report.bound, 0.21, "AbsTol", tc.Tol);
            tc.verifyEqual(report.finestCertifiable, 0.01);
            tc.verifyFalse(report.baselineBelowFloor);
            tc.verifyEqual(report.mode, "certification");
        end

        function statesTheFinestResolvableValueInTheMessage(tc)
            % +crit may not print, so the required run-start statement travels as a
            % string for +rhc and +io to place in the log and the run record.
            report = rid.crit.resolutionReport(0.2, 0.01, "certification");
            tc.verifyClass(report.message, "string");
            tc.verifyTrue(contains(report.message, "fillDistanceQ"));
            tc.verifyTrue(contains(report.message, "estimate"));
            tc.verifyTrue(contains(report.message, "headroom"));
        end

        function equalBaselineAndFillDistanceIsStillCertifiable(tc)
            % The check is epsReq0 < fillDistanceQ, so the boundary case passes:
            % a baseline exactly at the ruler is the finest certifiable baseline.
            report = rid.crit.resolutionReport(0.05, 0.05, "certification");
            tc.verifyEqual(report.headroom, 1, "AbsTol", tc.Tol);
            tc.verifyFalse(report.baselineBelowFloor);
        end

        function certificationModeRejectsABaselineBelowTheFloor(tc)
            % Every epsReq would be floored: the refinement structure is erased and
            % the run becomes a uniform design at the discretization limit while
            % claiming otherwise.
            tc.verifyError(@() rid.crit.resolutionReport(0.005, 0.01, "certification"), ...
                "rid:crit:uncertifiableBaseline");
        end

        function budgetModeReportsTheSameConfigurationInsteadOfFailing(tc)
            report = rid.crit.resolutionReport(0.005, 0.01, "budget");
            tc.verifyTrue(report.baselineBelowFloor);
            tc.verifyEqual(report.headroom, 0.5, "AbsTol", tc.Tol);
            tc.verifyTrue(contains(report.message, "not an error"));
        end

        function budgetModeMessageSaysNothingIsCertified(tc)
            report = rid.crit.resolutionReport(0.2, 0.01, "budget");
            tc.verifyTrue(contains(report.message, "nothing is certified"));
            tc.verifyFalse(report.baselineBelowFloor);
        end

        function rejectsInvalidArguments(tc)
            tc.verifyError(@() rid.crit.resolutionReport(0, 0.01, "budget"), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.resolutionReport(NaN, 0.01, "budget"), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.resolutionReport(0.2, 0, "budget"), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.resolutionReport(0.2, NaN, "budget"), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.resolutionReport(0.2, 0.01, "certified"), ...
                "rid:crit:unknownMode");
        end


        % --- the fill-distance estimate's standing --------------------

        function anUnconvergedEstimateIsStatedBesideTheHeadroom(tc)
            % The headroom is computed from fillDistanceQ, so an estimate that was
            % still climbing makes the headroom optimistic. The run-start
            % statement -- the thing that reaches the log and the record -- must
            % carry that, not merely the returned struct.
            report = rid.crit.resolutionReport(0.2, 0.01, "certification", ...
                "estimateConverged", false);
            tc.verifyFalse(report.fillDistanceConverged);
            tc.verifySubstring(report.message, "UNCONVERGED");
            tc.verifySubstring(report.message, "optimistic");
        end

        function aConvergedEstimateSaysNothingExtra(tc)
            report = rid.crit.resolutionReport(0.2, 0.01, "certification", ...
                "estimateConverged", true);
            tc.verifyTrue(report.fillDistanceConverged);
            tc.verifyFalse(contains(report.message, "UNCONVERGED"));
        end

        function theUnconvergedEstimateIsAReportInBothModes(tc)
            % It is not an error condition: an under-resolved probe set is worth
            % saying loudly, not worth refusing to run over.
            budget = rid.crit.resolutionReport(0.2, 0.01, "budget", ...
                "estimateConverged", false);
            tc.verifyFalse(budget.fillDistanceConverged);
            tc.verifySubstring(budget.message, "UNCONVERGED");
        end

        function theVerdictDefaultsToConvergedForTheAnalyticPath(tc)
            % The grid value is analytic: there is no probe set and nothing to
            % converge, which is the only reason the default is true.
            report = rid.crit.resolutionReport(0.2, 0.01, "budget");
            tc.verifyTrue(report.fillDistanceConverged);
        end


        % --- the bracket on the certifiable floor ---------------------

        function theProbeSetsOwnFillDistanceBracketsTheFloor(tc)
            report = rid.crit.resolutionReport(0.2, 0.01, "certification", ...
                "probeFillDistance", 0.004);
            tc.verifyEqual(report.probeFillDistance, 0.004);
            tc.verifyEqual(report.fillDistanceBracket, [0.01 0.014], ...
                "RelTol", 1e-12);
            tc.verifySubstring(report.message, "1-Lipschitz");
        end

        function withoutAProbeSetNothingIsBracketed(tc)
            % The analytic grid path has no probe set, so no bracket is claimed.
            report = rid.crit.resolutionReport(0.2, 0.01, "certification");
            tc.verifyTrue(all(isnan(report.fillDistanceBracket)));
            tc.verifyFalse(contains(report.message, "1-Lipschitz"));
        end

        function theBracketDoesNotMoveTheFloorOrTheHeadroom(tc)
            % The probe bracket is a statement, not a new threshold: the floor, the headroom
            % and the bound continue to use the measured fillDistanceQ.
            bare      = rid.crit.resolutionReport(0.2, 0.01, "certification");
            bracketed = rid.crit.resolutionReport(0.2, 0.01, "certification", ...
                "probeFillDistance", 0.5);
            tc.verifyEqual(bracketed.finestCertifiable, bare.finestCertifiable);
            tc.verifyEqual(bracketed.headroom, bare.headroom);
            tc.verifyEqual(bracketed.bound, bare.bound);
        end

        function rejectsANegativeProbeFillDistance(tc)
            tc.verifyError(@() rid.crit.resolutionReport(0.2, 0.01, "budget", ...
                "probeFillDistance", -1), "rid:crit:invalidResolution");
        end


        % --- the magnitude beside the flag --------------------

        function theTraceClimbIsQuotedBesideTheFlag(tc)
            % traceTolerance is a chosen value that fires on most runs at
            % defaults, so the boolean alone becomes wallpaper. A 2.1 %% climb and
            % a 40 %% climb are not the same news.
            report = rid.crit.resolutionReport(0.2, 0.01, "budget", ...
                "estimateConverged", false, "relativeClimb", [0.0091 0.0321]);
            tc.verifyEqual(report.relativeClimb, [0.0091 0.0321]);
            tc.verifySubstring(report.message, "0.91 %");
            tc.verifySubstring(report.message, "3.21 %");
        end

        function theClimbIsQuotedEvenWhenTheVerdictIsConverged(tc)
            % The magnitude is the readable quantity whichever way the flag fell.
            report = rid.crit.resolutionReport(0.2, 0.01, "budget", ...
                "estimateConverged", true, "relativeClimb", [0 0.0015]);
            tc.verifySubstring(report.message, "0.15 %");
            tc.verifyFalse(contains(report.message, "UNCONVERGED"));
        end

        function withoutATraceNoClimbIsQuoted(tc)
            report = rid.crit.resolutionReport(0.2, 0.01, "budget");
            tc.verifyEmpty(report.relativeClimb);
            tc.verifyFalse(contains(report.message, "Trace climb"));
        end

    end

end
