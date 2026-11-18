classdef CoverageReportTest < matlab.unittest.TestCase
%COVERAGEREPORTTEST  rid.crit.coverageReport: the three-state label has one owner.
%
%   WHAT THIS FUNCTION IS FOR. `rid.crit.evaluate` returns dTilde and psi
%   but needs `Znorm`; `probe` returns them for a CANDIDATE. Nothing returned them
%   for the state the loop had already committed, so `+viz` formed the saturation
%   boundary itself -- criterion arithmetic in the package the architecture rules
%   forbid it in, and a second implementation of the definition that explains every
%   run's behaviour. This is that accessor, and the tests below are the contract it
%   now owns on everyone's behalf.
%
%   THE BOUNDARIES ARE ASSERTED WHERE THEY ARE SHARP. A point exactly at
%   `dTilde = 1` is RESOLVED -- `psi = 0` there and it is not in the shell -- and a
%   point exactly at `dTilde = 1 + psiMax` is SATURATED, since its psi has reached
%   the clip and its gradient is gone. Those two equalities are the whole of the
%   difference between the three states, so a test that only checks points in the
%   middle of each band checks nothing that can be got wrong.

    methods (Test)

        function theThreeStatesSplitOnTheDeclaredBoundaries(tc)
            % Fails if: either boundary lands on the wrong side, or the shell is
            % opened/closed at the wrong end. Distances are chosen so dTilde is
            % exactly 0.5, 1, 1.5, 2 and 3 against epsReq = 1 and psiMax = 1.
            state = CoverageReportTest.stateWithDistances([0.5; 1; 1.5; 2; 3], 1, 1);

            report = rid.crit.coverageReport(state);

            tc.verifyEqual(report.dTilde, [0.5; 1; 1.5; 2; 3]);
            tc.verifyEqual(report.label, ...
                ["resolved"; "resolved"; "active"; "saturated"; "saturated"], ...
                "dTilde = 1 is resolved (psi = 0, not in the shell); " + ...
                "dTilde = 1 + psiMax is saturated (psi is at the clip)");
            tc.verifyEqual(report.count, [2 1 2]);
        end

        function psiIsTheClippedPenaltyOfSectionTwo(tc)
            % Fails if: psi is returned unclipped, or clipped at the wrong value.
            state  = CoverageReportTest.stateWithDistances([0.5; 1; 1.5; 2; 3], 1, 1);
            report = rid.crit.coverageReport(state);
            tc.verifyEqual(report.psi, [0; 0; 0.5; 1; 1], "AbsTol", 0);
        end

        function theLabelAgreesWithTheLogicalMasks(tc)
            % Fails if: the string label and the masks are computed separately and
            % can disagree -- which would be this entry's own defect reintroduced
            % one level down.
            state  = CoverageReportTest.stateWithDistances([0.4; 1.2; 2.5; 0.9], 1, 1);
            report = rid.crit.coverageReport(state);

            tc.verifyEqual(report.label == "resolved",  report.isResolved);
            tc.verifyEqual(report.label == "active",    report.isActive);
            tc.verifyEqual(report.label == "saturated", report.isSaturated);
            tc.verifyEqual(report.isResolved + report.isActive + report.isSaturated, ...
                ones(4, 1), "every point is in exactly one state");
        end

        function anUncoveredPointIsSaturated(tc)
            % `dRaw` is Inf where nothing has been sampled, so dTilde is Inf and psi
            % sits at its clip. Saturated is the correct answer, and it stays correct
            % under psiMax = Inf, where Inf >= Inf.
            fresh  = rid.crit.CoverageState([0 0; 1 1], [1; 1]);
            report = rid.crit.coverageReport(fresh);
            tc.verifyEqual(report.label, ["saturated"; "saturated"]);

            noClip = rid.crit.CoverageState([0 0; 1 1], [1; 1], "psiMax", Inf);
            tc.verifyEqual(rid.crit.coverageReport(noClip).label, ...
                ["saturated"; "saturated"], "psiMax = Inf: Inf >= Inf");
        end

        function psiMaxInfLeavesOnlyResolvedAndActive(tc)
            % With no clip there is nothing to saturate against at a finite
            % distance, so the third state is empty by construction.
            state  = CoverageReportTest.stateWithDistances([0.5; 1.5; 9], 1, Inf);
            report = rid.crit.coverageReport(state);
            tc.verifyEqual(report.label, ["resolved"; "active"; "active"]);
            tc.verifyEqual(report.psi, [0; 0.5; 8], "AbsTol", 0);
        end

        function epsReqIsPerSupportingPointAndNotAScalar(tc)
            % Fails if: a single epsReq is assumed. Refinement makes epsReq vary
            % across Q, so the same raw distance is a different dTilde at
            % different supporting points -- which is the whole mechanism.
            state  = CoverageReportTest.stateWithDistances([1; 1; 1], [2; 1; 0.4], 1);
            report = rid.crit.coverageReport(state);
            tc.verifyEqual(report.dTilde, [0.5; 1; 2.5], "AbsTol", 1e-15);
            tc.verifyEqual(report.label, ["resolved"; "resolved"; "saturated"]);
        end

        function itAgreesWithTheOracleOnACommittedState(tc)
            % THE EQUIVALENCE THAT MATTERS. `rid.crit.evaluate` is the permanent
            % independent recomputation and shares no code with the cached
            % path, so this is the test that catches the accessor drifting away
            % from the definition it is supposed to own.
            rng(7, "twister");
            Q      = rand(15, 3);
            epsReq = 0.15 + 0.2 * rand(15, 1);
            Znorm  = rand(9, 3);

            state = rid.crit.CoverageState(Q, epsReq, "psiMax", 0.8);
            state.add(Znorm);

            report    = rid.crit.coverageReport(state);
            [~, info] = rid.crit.evaluate(Znorm, Q, epsReq, "psiMax", 0.8);

            tc.verifyEqual(report.dTilde, info.dTilde, "RelTol", 1e-12);
            tc.verifyEqual(report.psi,    info.psi,    "RelTol", 1e-12);
            tc.verifyEqual(sum(report.isResolved), sum(info.dTilde <= 1));
            tc.verifyEqual(sum(report.isSaturated), sum(info.dTilde >= 1 + 0.8));
        end

        function aFixtureInWhichAllThreeStatesAreOccupied(tc)
            % The seen-to-fail rule, for the oracle test above: with an empty saturated
            % set, folding it into the active one would agree by construction.
            rng(7, "twister");
            Q      = rand(15, 3);
            epsReq = 0.15 + 0.2 * rand(15, 1);
            state  = rid.crit.CoverageState(Q, epsReq, "psiMax", 0.8);
            state.add(rand(9, 3));
            tc.verifyEqual(nnz(rid.crit.coverageReport(state).count), 3);
        end

        function theCountsAreOrderedAsTheNamesSay(tc)
            state  = CoverageReportTest.stateWithDistances([0.5; 1.5; 3], 1, 1);
            report = rid.crit.coverageReport(state);
            tc.verifyEqual(report.stateName, ["resolved", "active", "saturated"]);
            tc.verifyEqual(report.count,     [1 1 1]);
            tc.verifyEqual(sum(report.count), report.nSup);
        end

        function itTakesACoverageStateAndNothingElse(tc)
            tc.verifyError(@() rid.crit.coverageReport(struct("dRaw", 1)), ...
                "MATLAB:validation:UnableToConvert");
        end

    end

    methods (Static, Access = private)

        function state = stateWithDistances(dRaw, epsReq, psiMax)
            %STATEWITHDISTANCES  A committed state whose dRaw is exactly as given.
            %
            %   Built through `add`, not by writing the cache: `dRaw` is private
            %   to the class and a test that reached around the constructor would
            %   be testing a state the loop can never produce. One supporting
            %   point per requested distance, each with its own sample placed at
            %   exactly that distance along the first coordinate.
            nSup   = numel(dRaw);
            Q      = [zeros(nSup, 1), (1:nSup).' * 100];   % far apart in dim 2
            Znorm  = [dRaw(:), (1:nSup).' * 100];
            if isscalar(epsReq)
                epsReq = repmat(epsReq, nSup, 1);
            end
            state = rid.crit.CoverageState(Q, epsReq(:), "psiMax", psiMax);
            state.add(Znorm);
        end

    end

end
