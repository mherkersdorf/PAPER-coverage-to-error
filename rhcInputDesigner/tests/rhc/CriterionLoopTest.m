classdef CriterionLoopTest < matlab.unittest.TestCase
%CRITERIONLOOPTEST  The loop driven by a criterion with none of the coverage features.
%
%   MeanSquareCriterion keeps every default of rid.crit.Criterion: it never
%   certifies, is not monotone, has no resolution, no active set and no
%   supporting points. Each test pins one place where the loop reaches a
%   coverage-specific feature through the interface rather than through the class
%   it holds. The coverage side of the same gates is pinned by RunTest and by the
%   canaries, which must stay unchanged.

    methods (Test)

        function aNonMonotoneCriterionRunsToItsBudget(tc)
            % Fails if: the invariant-3 assert is applied to a criterion that does
            % not declare monotonicity -- it would raise rid:rhc:monotonicityViolated
            % at the first rise.
            r = CriterionLoopTest.loopRun(6);

            tc.verifyEqual(r.status, "budgetReached");
            tc.verifyEqual(r.nSteps, 6);
            % Non-vacuity: the J trace must actually RISE somewhere (from J = 0
            % with nothing committed, the first commit raises it), or the run
            % would pass with or without the gate.
            tc.assertTrue(any(diff([0; r.Jtrace]) > 0), ...
                "the fixture must produce a rise in J, or the gate is not tested");
        end

        function theAbsentResolutionIsRecordedAsNaN(tc)
            % NaN encodes the absence of the concept, never a default.
            %
            % Fails if: the record copies an empty matrix, or a resolution from
            % the options, into a place where the criterion has none.
            r = CriterionLoopTest.loopRun(3);

            tc.verifyTrue(isnan(r.record.setup.epsReqInitial));
            tc.verifyTrue(isnan(r.record.problem.epsReq));
            tc.verifyTrue(isnan(r.certification.epsReqMax));
            tc.verifyTrue(isnan(r.certification.bound));
            tc.verifyFalse(r.certification.isSatisfied);
            tc.verifyTrue(all(isnan([r.steps.nActive])));
            tc.verifyEqual(r.record.problem.criterion, ...
                "MeanSquareCriterion: target [0.9 0.9 0.9]");
        end

        function aCriterionThatNeverCertifiesStopsOnlyOnAGuard(tc)
            % isSatisfied is false by default, so a certification-mode run
            % ends on its guard and never reports "certified".
            % epsReq0 = 0.75: certification mode refuses a baseline finer than the
            % fill distance of this Q (about 0.70), which is a different refusal.
            r = CriterionLoopTest.loopRun(NaN, "mode", "certification", ...
                    "maxSteps", 3, "epsReq0", 0.75);

            tc.verifyEqual(r.status, "guardTripped");
            tc.verifyEqual(r.nSteps, 3);
            tc.verifyEmpty(r.eventLog, ...
                "no reset and no update may be logged for this criterion");
        end

        function aFieldUpdateForACriterionWithoutResolutionRaises(tc)
            % Fails if: a field update is silently dropped for a criterion that
            % cannot take it.
            tc.verifyError(@() CriterionLoopTest.loopRun(3, "eagerFields", true), ...
                "rid:crit:noResolution");
        end

        function aRunOfSuchACriterionReplays(tc)
            % Invariant 5, through `freshCopy()` rather than the coverage
            % reconstruction from epsReqInitial.
            r = CriterionLoopTest.loopRun(4);
            [ok, difference] = rid.io.replay(r.record);

            tc.verifyTrue(ok, difference.reason);
            tc.verifyEqual(difference.result.Jtrace, r.Jtrace);
        end

    end

    methods (Static)

        function r = loopRun(nSteps, opts)
            %LOOPRUN  The SetupTest fixture with MeanSquareCriterion in place of the state.
            arguments
                nSteps           (1,1) double
                opts.mode        (1,1) string = "budget"
                opts.maxSteps    (1,1) double = Inf
                opts.eagerFields (1,1) logical = false
                opts.epsReq0     (1,1) double = RunTest.EpsReq0
            end
            f        = SetupTest.fixture("epsReq0", opts.epsReq0, "mode", opts.mode);
            fieldSch = f.fieldSch;
            if opts.eagerFields
                fieldSch = EagerFieldScheduler();
            end
            setup = rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, f.fields, ...
                        f.gen, f.sol, MeanSquareCriterion([0.9 0.9 0.9]), f.sch, ...
                        fieldSch, f.qInfo);
            ov = RunTest.overrides("mode", opts.mode, "maxSteps", opts.maxSteps, ...
                     "epsReq0", opts.epsReq0);
            if ~isnan(nSteps)
                ov.rhc.nSteps = nSteps;
            end
            ov.io.storeQ = true;   % replay needs Q
            r = rid.rhc.run([], setup, RunTest.history("srg", f.srg), ov);
        end

    end

end
