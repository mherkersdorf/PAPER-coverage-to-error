classdef CoverageStateWithoutQTest < matlab.unittest.TestCase
%COVERAGESTATEWITHOUTQTEST  The record-only copy of a coverage state.
%
%   WHY THE CLASS HAS TO DO THIS AT ALL. `storeQ = false` must remove `Q` from
%   the record completely, the stored `Setup` included -- `setup.Q` and
%   `setup.state` both. `Q` is not assignable from outside,
%   and the constructor refuses an empty one, so a Q-less state can only come
%   from the class that owns the invariant.
%
%   IT IS A COPY, AND THAT IS THE POINT. `rid.crit.CoverageState` is the
%   project's one intentional handle class, so stripping in place would empty
%   the live state the caller is still holding -- the run's own `result`. The
%   original must be untouched, which is the sharpest thing here to get wrong.
%
%   WHAT SURVIVES is everything a reader needs and nothing that identifies a
%   point in space: `epsReq`, `dRaw`, the active set, `nSup`, `nDim`, `psiMax`
%   and `pNorm` all stay, so `rid.crit.coverageReport` still works on a stripped
%   record. Only the coordinates go.

    methods (Test)

        function theCopyHasNoQAndTheOriginalIsUntouched(tc)
            % THE SHARP ONE. A handle class stripped in place would empty the
            % caller's own object, and the defect would show up as a run result
            % losing its supporting points after being recorded.
            state = CoverageStateWithoutQTest.committedState();
            Qwas  = state.Q;

            copy = state.withoutQ();

            tc.verifyEmpty(copy.Q, "the copy carries no supporting points");
            tc.verifyEqual(state.Q, Qwas, ...
                "the ORIGINAL must be untouched: it is a handle the caller holds");
            tc.verifyNotSameHandle(copy, state);
        end

        function everythingButTheCoordinatesSurvives(tc)
            % Fails if: the cache is rebuilt rather than carried over, which
            % would silently reset dRaw to Inf and report a run that covered
            % nothing.
            state = CoverageStateWithoutQTest.committedState();
            copy  = state.withoutQ();

            tc.verifyEqual(copy.nSup,     state.nSup, "nSup is the TRUE count");
            tc.verifyEqual(copy.nDim,     state.nDim);
            tc.verifyEqual(copy.epsReq,   state.epsReq);
            tc.verifyEqual(copy.dRaw,     state.dRaw);
            tc.verifyEqual(copy.isActive, state.isActive);
            tc.verifyEqual(copy.psiMax,   state.psiMax);
            tc.verifyEqual(copy.pNorm,    state.pNorm);
            tc.verifyEqual(copy.value(),  state.value());
        end

        function theThreeStateReportStillWorksOnTheCopy(tc)
            % The reason the cache has to survive: a record written with
            % storeQ = false must still be able to say what the run resolved.
            state = CoverageStateWithoutQTest.committedState();
            tc.verifyEqual(rid.crit.coverageReport(state.withoutQ()), ...
                rid.crit.coverageReport(state));
        end

        function theCopyIsStillTheRightShapeForItsSupportingPoints(tc)
            % nSup must not collapse to zero with Q. A record reporting nSup = 0
            % would describe a run that had no supporting points at all, which
            % is a different and much worse lie than one that omits their
            % coordinates.
            state = CoverageStateWithoutQTest.committedState();
            copy  = state.withoutQ();
            tc.verifyEqual(copy.nSup, 5);
            tc.verifySize(copy.epsReq, [5 1]);
            tc.verifySize(copy.Q, [0 2], "empty, but still nDim columns wide");
        end

    end

    methods (Static, Access = private)

        function state = committedState()
            %COMMITTEDSTATE  A state with data behind it, not a fresh one.
            %
            %   A fresh state has dRaw = Inf everywhere, so carrying the cache
            %   over and rebuilding it would give the same answer and the test
            %   could not tell them apart.
            rng(3, "twister");
            Q     = rand(5, 2);
            state = rid.crit.CoverageState(Q, 0.3 * ones(5, 1), "psiMax", 0.9);
            state.add(rand(4, 2));
        end

    end

end
