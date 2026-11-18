classdef ReachabilityReportTest < matlab.unittest.TestCase
%REACHABILITYREPORTTEST  THE VALIDATE-NEVER-CORRECT RULE: the report validates and corrects nothing.
%
%   The load-bearing pair is the first two tests. A report that always warns passes
%   "a box outside is flagged" on its own, so the silent case is what gives the
%   first one meaning -- neither test is worth anything without the other.
%
%   The volume fraction is computed BY HAND in the test from a case with partial
%   overlap in exactly one dimension, so an implementation that averages the
%   per-dimension overlaps instead of multiplying them cannot pass.

    properties (Constant)
        Tol   = 1e-12
        Names = ["u1_k" "du1_k" "y_km1"]
    end

    methods (Test)

        % --- check 1: the ROI against the admissible box ----------------------

        function aRoiLeavingTheAdmissibleBoxIsFlaggedWithTheRightFraction(tc)
            % Admissible u1_k is [-2 3]; the ROI asks for [-2 8], width 10, of which
            % [-2 3] (width 5) is inside -- exactly half. The other two dimensions
            % are fully inside, so the volume fraction outside is 1 - 0.5*1*1 = 0.5.
            %
            % Fails if: the per-dimension overlaps are averaged rather than
            % multiplied (that would give 1 - 5/6 = 0.167), or the fraction is
            % computed over the admissible box instead of the ROI box.
            cons = ReachabilityReportTest.constraints();
            roi  = [-2 8; -0.4 0.4; -10 10];
            tc.verifyWarning(@() rid.gen.reachabilityReport(roi, ...
                ReachabilityReportTest.Names, cons), "rid:gen:roiOutsideAdmissible");
            report = ReachabilityReportTest.quiet(roi, cons, "budget");

            tc.verifyTrue(report.roiOutside);
            tc.verifyEqual(report.fractionOutside, 0.5, "AbsTol", tc.Tol);
            tc.verifyEqual(report.perDimOverlap, [0.5 1 1], "AbsTol", tc.Tol);
        end

        function aRoiStrictlyInsideIsSilent(tc)
            % Without this the first test is satisfied by a report that always
            % warns. Both directions are required.
            %
            % Fails if: the comparison is inverted, or a boundary-touching ROI is
            % counted as outside.
            cons = ReachabilityReportTest.constraints();
            roi  = [-1 2; -0.2 0.2; -5 5];
            tc.verifyWarningFree(@() rid.gen.reachabilityReport(roi, ...
                ReachabilityReportTest.Names, cons));
            report = ReachabilityReportTest.quiet(roi, cons, "budget");
            tc.verifyFalse(report.roiOutside);
            tc.verifyEqual(report.fractionOutside, 0, "AbsTol", tc.Tol);
            tc.verifyEmpty(report.messages);
        end

        function aRoiTouchingTheBoundaryExactlyIsStillInside(tc)
            % The ROI equal to the admissible box is the ordinary case, not a
            % violation.
            cons = ReachabilityReportTest.constraints();
            roi  = [-2 3; -0.4 0.4; -10 10];
            tc.verifyWarningFree(@() rid.gen.reachabilityReport(roi, ...
                ReachabilityReportTest.Names, cons));
        end

        function aRoiWithNoOverlapAtAllReportsExactlyOne(tc)
            % Fails if: the overlap is allowed to go negative and the fraction
            % exceeds 1, which a missing max(0, .) would produce.
            cons   = ReachabilityReportTest.constraints();
            roi    = [20 30; -0.4 0.4; -10 10];
            report = ReachabilityReportTest.quiet(roi, cons, "budget");
            tc.verifyEqual(report.fractionOutside, 1, "AbsTol", tc.Tol);
            tc.verifyEqual(report.perDimOverlap(1), 0, "AbsTol", tc.Tol);
        end

        function theSameDeclarationWarnsInBudgetAndErrorsInCertification(tc)
            % In certification mode J = 0 is the termination condition, and an
            % uncoverable supporting point makes it unattainable, so the run cannot
            % terminate. That is a stop, not a warning.
            %
            % Fails if: the mode is ignored, or certification only warns.
            cons = ReachabilityReportTest.constraints();
            roi  = [-2 8; -0.4 0.4; -10 10];
            tc.verifyWarning(@() rid.gen.reachabilityReport(roi, ...
                ReachabilityReportTest.Names, cons, "mode", "budget"), ...
                "rid:gen:roiOutsideAdmissible");
            tc.verifyError(@() rid.gen.reachabilityReport(roi, ...
                ReachabilityReportTest.Names, cons, "mode", "certification"), ...
                "rid:gen:roiOutsideAdmissible");
        end

        % --- check 2: rate consistency ----------------------------------------

        function aRateDisagreementWarnsInBothModesAndNeverErrors(tc)
            % Pinned deliberately. The two are separate declarations by design
            %, so escalating this in certification mode would be wrong --
            % this test exists to stop a later "helpful" escalation.
            %
            % Fails if: certification mode promotes the rate warning to an error.
            cons = ReachabilityReportTest.constraints("rateLimit", ...
                       struct("u1", [-0.1 0.25]));       % != du1_k's [-0.4 0.4]
            roi  = [-1 2; -0.2 0.2; -5 5];               % ROI strictly inside
            for mode = ["budget", "certification"]
                tc.verifyWarning(@() rid.gen.reachabilityReport(roi, ...
                    ReachabilityReportTest.Names, cons, "mode", mode), ...
                    "rid:gen:rateDisagreement");
                report = ReachabilityReportTest.quiet(roi, cons, mode);
                tc.verifyEqual(numel(report.rateDisagreement), 1);
                tc.verifyEqual(report.rateDisagreement(1).channel, "u1");
                tc.verifyEqual(report.rateDisagreement(1).declared, [-0.1 0.25]);
                tc.verifyEqual(report.rateDisagreement(1).admissible, [-0.4 0.4]);
            end
        end

        function agreeingRateDeclarationsAreSilent(tc)
            % Fails if: the check fires whenever a rate limit is declared at all,
            % rather than when the two numbers differ.
            cons = ReachabilityReportTest.constraints("rateLimit", ...
                       struct("u1", [-0.4 0.4]));        % == du1_k
            roi  = [-1 2; -0.2 0.2; -5 5];
            tc.verifyWarningFree(@() rid.gen.reachabilityReport(roi, ...
                ReachabilityReportTest.Names, cons));
            report = ReachabilityReportTest.quiet(roi, cons, "budget");
            tc.verifyEmpty(report.rateDisagreement);
        end

        function anUndeclaredRateLimitIsNotComparedAtAll(tc)
            % The check needs BOTH declarations to exist. Without an explicit rate
            % limit the default IS the du1_k bound, so comparing them would
            % always agree and the check would be decorative -- but the point is
            % that nothing is reported, not that it happens to match.
            cons   = ReachabilityReportTest.constraints();   % no rateLimit given
            report = ReachabilityReportTest.quiet([-1 2; -0.2 0.2; -5 5], ...
                         cons, "budget");
            tc.verifyEmpty(report.rateDisagreement);
        end

        % --- it corrects nothing ------------------------------------------------

        function theInputsAreUnmodifiedAfterTheCall(tc)
            % "Corrects nothing" is a claim, so it is asserted. The validate-never-correct rule: a clip only
            % ever narrows, which looks safe, but it changes the normalization and
            % therefore the meaning of every coordinate in every stored run.
            %
            % Fails if: the report clips the ROI box into the admissible one.
            cons     = ReachabilityReportTest.constraints();
            roi      = [-2 8; -0.4 0.4; -10 10];
            roiCopy  = roi;
            ampCopy  = cons.amplitude;
            boxCopy  = cons.admissibleBox;
            rateCopy = cons.rate;

            ReachabilityReportTest.quiet(roi, cons, "budget");

            tc.verifyEqual(roi, roiCopy);
            tc.verifyEqual(cons.amplitude, ampCopy);
            tc.verifyEqual(cons.admissibleBox, boxCopy);
            tc.verifyEqual(cons.rate, rateCopy);
        end

        % --- names present in one box and not the other -------------------------

        function aNameInOneBoxAndNotTheOtherIsAnError(tc)
            % Reported rather than chosen silently. The declaration-validation rule puts the missing-name and
            % unknown-name checks on `bounds`, which BOTH boxes already passed, so
            % two boxes disagreeing about which coordinates exist did not come from
            % one representation. Comparing only the shared dimensions would report
            % a volume fraction over a different space than the one declared -- a
            % quietly wrong number rather than an absent one.
            %
            % Fails if: the mismatched dimensions are skipped.
            cons = ReachabilityReportTest.constraints();
            tc.verifyError(@() rid.gen.reachabilityReport( ...
                [-2 3; -0.4 0.4], ["u1_k" "du1_k"], cons), ...
                "rid:gen:dimensionNameMismatch");
            tc.verifyError(@() rid.gen.reachabilityReport( ...
                [-2 3; -0.4 0.4; -10 10; 0 1], ...
                [ReachabilityReportTest.Names, "x1"], cons), ...
                "rid:gen:dimensionNameMismatch");
        end

        function aMisSizedNameListIsAnError(tc)
            cons = ReachabilityReportTest.constraints();
            tc.verifyError(@() rid.gen.reachabilityReport( ...
                [-2 3; -0.4 0.4; -10 10], ["u1_k" "du1_k"], cons), ...
                "rid:gen:invalidRoiBox");
        end

        % --- what it must not claim ---------------------------------------------

        function theReportSaysNothingAboutDynamicReachability(tc)
            % The report compares DECLARED boxes and cannot see
            % the reachable region. An asymmetric rate limit shrinks the
            % horizon-reachable set far below the amplitude box, and this ROI is
            % dynamically unreachable in large part -- yet the report is silent,
            % correctly, because both boxes agree.
            %
            % Fails if: someone adds a dynamic-reachability claim here; the report
            % would then have to warn and this test pins that it does not.
            cons = ReachabilityReportTest.constraints("rateLimit", ...
                       struct("u1", [-0.001 0.001]));    % barely able to move
            roi  = [-2 3; -0.001 0.001; -10 10];
            report = ReachabilityReportTest.quiet(roi, cons, "budget");
            tc.verifyFalse(report.roiOutside);
            tc.verifyEqual(report.fractionOutside, 0, "AbsTol", tc.Tol);
        end

    end

    methods (Static, Access = private)

        function cons = constraints(varargin)
            %CONSTRAINTS  u1_k [-2 3], du1_k [-0.4 0.4], y_km1 [-10 10].
            cons = rid.gen.Constraints(rid.rep.Narx(2, 1), ...
                struct("u1_k", [-2 3], "du1_k", [-0.4 0.4], "y_km1", [-10 10]), ...
                varargin{:});
        end

        function report = quiet(roi, cons, mode)
            %QUIET  Run the report with warnings suppressed, to inspect its struct.
            state = warning("off", "rid:gen:roiOutsideAdmissible");
            restore1 = onCleanup(@() warning(state));
            state2 = warning("off", "rid:gen:rateDisagreement");
            restore2 = onCleanup(@() warning(state2));
            report = rid.gen.reachabilityReport(roi, ...
                ReachabilityReportTest.Names, cons, "mode", mode);
        end

    end

end
