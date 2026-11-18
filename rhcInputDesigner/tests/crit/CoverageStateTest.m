classdef CoverageStateTest < matlab.unittest.TestCase
%COVERAGESTATETEST  Behaviour of the cached criterion state.
%
%   Layers 1 and 2 of the project's test plan. The randomized comparison against the oracle lives
%   in CoverageStateOracleTest; this file covers the semantics of the individual
%   operations and the invariants of the cache.
%
%   ON THE COMPLEXITY CONTRACT. The method description requires that scoring a
%   candidate cost O(nSup * nHorizon) and never O(nSup * N). That contract is
%   ENFORCED BY REVIEW, not by a test in this file. Asserting it would need either
%   timing, which is not reproducible inside runtests, or a counting double injected
%   into the distance path, which would mean changing the implementation to suit the
%   test. What is tested here is the behavioural property the contract rests on:
%   `probe` reads only the cache, so it returns an identical value whether the same
%   dRaw was produced by ten committed samples or by ten thousand.
%
%   The other invariant this file is responsible for is that `probe` does not
%   mutate: it is asserted on the handle itself, since a handle class gives no
%   protection of its own.
%
%   ON REBUILD. In the default configuration -- online, measured past values, a
%   representation built from raw signals -- `rebuild` NEVER RUNS.
%   That removes it as a performance
%   concern and turns it into a maintenance risk: code exercised only in rare
%   configurations rots unnoticed and fails first in the experiment that finally
%   needs it -- simulated past values, a representation with
%   dependsOnSurrogate = true, replay reconstructing a run, or a mid-run pastSource
%   switch. These tests are therefore the only thing keeping it working, so the
%   equivalence with incremental `add` is asserted BIT-IDENTICALLY rather than to a
%   tolerance, across several batchings, and `rebuild` also stays in the operation
%   mix of CoverageStateOracleTest.

    properties (Constant)
        Seed = 20260830
        Tol  = 1e-12
    end

    methods (Access = private)

        function s = stream(tc, offset)
            s = RandStream("twister", "Seed", tc.Seed + offset);
        end

        function d = seedDiagnostic(tc)
            d = sprintf("seed = %d", tc.Seed);
        end

    end

    methods (Test)

        function startsUncoveredAtPsiMax(tc)
            % Before the first sample every supporting point is infinitely far from
            % the distribution, so J saturates.
            state = rid.crit.CoverageState((0:0.25:1).', 0.2 * ones(5, 1), ...
                "psiMax", 1, "pNorm", 1);
            tc.verifyEqual(state.dRaw, inf(5, 1));
            tc.verifyEqual(state.nSamples, 0);
            tc.verifyEqual(state.value(), 1, "AbsTol", tc.Tol);
        end

        function addUpdatesTheCacheAsARunningMinimum(tc)
            Q     = (0:0.25:1).';
            state = rid.crit.CoverageState(Q, 0.25 * ones(5, 1));
            state.add(0.5);
            tc.verifyEqual(state.dRaw, [0.5; 0.25; 0; 0.25; 0.5], "AbsTol", tc.Tol);
            tc.verifyEqual(state.nSamples, 1);

            state.add(1);
            tc.verifyEqual(state.dRaw, [0.5; 0.25; 0; 0.25; 0], "AbsTol", tc.Tol);
            tc.verifyEqual(state.nSamples, 2);

            % A sample farther away than what is already committed changes nothing.
            cacheBefore = state.dRaw;
            state.add(50);
            tc.verifyEqual(state.dRaw, cacheBefore);
            tc.verifyEqual(state.nSamples, 3);
        end

        function addIsInvariantUnderSampleOrderAndBatching(tc)
            s      = tc.stream(1);
            Q      = rand(s, 30, 2);
            epsReq = 0.05 * ones(30, 1);
            Znorm      = rand(s, 12, 2);

            oneBatch = rid.crit.CoverageState(Q, epsReq);
            oneBatch.add(Znorm);

            oneByOne = rid.crit.CoverageState(Q, epsReq);
            for iSample = randperm(s, 12)
                oneByOne.add(Znorm(iSample, :));
            end

            tc.verifyEqual(oneByOne.dRaw, oneBatch.dRaw, ...
                "AbsTol", tc.Tol, tc.seedDiagnostic());
            tc.verifyEqual(oneByOne.value(), oneBatch.value(), ...
                "AbsTol", tc.Tol, tc.seedDiagnostic());
        end

        function addEqualsRebuildBitForBit(tc)
            % The maintenance guard for the rebuild path. Both routes reduce the
            % same set of pairwise distances to the same per-point minimum, so they
            % must agree EXACTLY -- no tolerance. A tolerance here would let the two
            % paths drift apart unnoticed, which is precisely the failure mode a
            % rarely executed code path has.
            %
            % Several batchings, because how the samples were split is the only
            % thing that differs between the two routes: one sample at a time, in
            % uneven groups, and all at once.
            s      = tc.stream(2);
            nSup   = 60;
            Q      = rand(s, nSup, 2);
            epsReq = 0.02 + 0.1 * rand(s, nSup, 1);
            Znorm      = rand(s, 37, 2);

            fromScratch = rid.crit.CoverageState(Q, epsReq);
            fromScratch.rebuild(Znorm);

            batchings = {ones(1, 37), [5 5 5 5 5 5 5 2], [1 6 30], 37};
            for iBatching = 1:numel(batchings)
                sizes       = batchings{iBatching};
                incremental = rid.crit.CoverageState(Q, epsReq);
                iFirst      = 1;
                for nBatch = sizes
                    incremental.add(Znorm(iFirst:iFirst + nBatch - 1, :));
                    iFirst = iFirst + nBatch;
                end

                what = sprintf("%s, batching %d of %d", tc.seedDiagnostic(), ...
                    iBatching, numel(batchings));
                tc.verifyTrue(isequal(incremental.dRaw, fromScratch.dRaw), ...
                    "the cache differs between incremental add and rebuild: " + what);
                tc.verifyEqual(incremental.value(), fromScratch.value(), what);
                tc.verifyEqual(incremental.nSamples, fromScratch.nSamples, what);
            end
        end

        function rebuildAfterAddsRecoversTheCacheBitForBit(tc)
            % The surrogate-update path itself: samples have been committed, the
            % coordinates then move, and the cache is recovered from the full
            % history. Its result must equal a state that only ever saw the moved
            % coordinates -- exactly, and whatever was committed before.
            s      = tc.stream(8);
            nSup   = 45;
            Q      = rand(s, nSup, 3);
            epsReq = 0.02 + 0.1 * rand(s, nSup, 1);

            used = rid.crit.CoverageState(Q, epsReq);
            used.add(rand(s, 9, 3));          % a history in the OLD coordinates
            used.add(rand(s, 4, 3));

            Zmoved = rand(s, 13, 3);          % the same history, re-derived
            used.rebuild(Zmoved);

            fresh = rid.crit.CoverageState(Q, epsReq);
            fresh.rebuild(Zmoved);

            tc.verifyTrue(isequal(used.dRaw, fresh.dRaw), ...
                "rebuild did not fully discard the pre-update cache: " + ...
                tc.seedDiagnostic());
            tc.verifyEqual(used.nSamples, 13);

            % And the recovered cache is the one the oracle computes.
            [~, oracle] = rid.crit.evaluate(Zmoved, Q, epsReq);
            tc.verifyEqual(used.dRaw, oracle.dRaw, "RelTol", tc.Tol, ...
                tc.seedDiagnostic());
        end

        function rebuildIsBitIdenticalAcrossShapes(tc)
            % The same statement over a range of sizes and dimensions, since the
            % underlying distance primitive is free to choose its algorithm by input
            % size. A configuration where the two routes disagree would be a real
            % finding, not a tolerance to widen.
            s = tc.stream(9);
            for kTrial = 1:12
                nSup = randi(s, [1 120]);
                nDim = randi(s, [1 4]);
                nAll = randi(s, [1 150]);
                Q    = rand(s, nSup, nDim);
                Znorm    = rand(s, nAll, nDim);

                incremental = rid.crit.CoverageState(Q, ones(nSup, 1));
                for iSample = 1:nAll
                    incremental.add(Znorm(iSample, :));
                end
                fromScratch = rid.crit.CoverageState(Q, ones(nSup, 1));
                fromScratch.rebuild(Znorm);

                tc.verifyTrue(isequal(incremental.dRaw, fromScratch.dRaw), ...
                    sprintf("%s, trial %d: nSup = %d, nDim = %d, nAll = %d", ...
                        tc.seedDiagnostic(), kTrial, nSup, nDim, nAll));
            end
        end

        function probeDoesNotMutateAnyState(tc)
            % The optimizer calls probe once per candidate; a mutation here would
            % make the candidates depend on the order they happened to be scored in.
            s      = tc.stream(3);
            Q      = rand(s, 30, 2);
            state  = rid.crit.CoverageState(Q, 0.05 * ones(30, 1));
            state.add(rand(s, 6, 2));

            cacheBefore    = state.dRaw;
            epsReqBefore   = state.epsReq;
            nSamplesBefore = state.nSamples;
            valueBefore    = state.value();

            for kCand = 1:10
                state.probe(rand(s, 4, 2));
            end

            tc.verifyEqual(state.dRaw, cacheBefore, tc.seedDiagnostic());
            tc.verifyEqual(state.epsReq, epsReqBefore, tc.seedDiagnostic());
            tc.verifyEqual(state.nSamples, nSamplesBefore, tc.seedDiagnostic());
            tc.verifyEqual(state.value(), valueBefore, tc.seedDiagnostic());
        end

        function probeEqualsTheValueAfterCommittingTheSameCandidate(tc)
            s      = tc.stream(4);
            Q      = rand(s, 30, 2);
            epsReq = 0.05 * ones(30, 1);
            state  = rid.crit.CoverageState(Q, epsReq);
            state.add(rand(s, 6, 2));

            Zcand = rand(s, 4, 2);
            jProbe = state.probe(Zcand);
            state.add(Zcand);
            tc.verifyEqual(state.value(), jProbe, ...
                "AbsTol", tc.Tol, tc.seedDiagnostic());
        end

        function probeDependsOnTheCacheAloneNotOnTheHistorySize(tc)
            % The behavioural form of the O(nSup * nHorizon) contract: two states
            % whose caches agree score a candidate identically, even though one
            % carries ten samples behind that cache and the other ten thousand.
            s      = tc.stream(5);
            Q      = rand(s, 40, 2);
            epsReq = 0.05 * ones(40, 1);

            % The same cache with vastly more history behind it: the extra samples
            % are copies of the first ten, so every nearest-neighbour distance is
            % unchanged and only nSamples differs.
            zFew = rand(s, 10, 2);
            few  = rid.crit.CoverageState(Q, epsReq);
            few.add(zFew);
            many = rid.crit.CoverageState(Q, epsReq);
            many.add(repmat(zFew, 1000, 1));

            tc.assertEqual(many.dRaw, few.dRaw, "AbsTol", tc.Tol);
            tc.verifyEqual(many.nSamples, 10000);
            tc.verifyEqual(few.nSamples, 10);

            Zcand = rand(s, 5, 2);
            tc.verifyEqual(many.probe(Zcand), few.probe(Zcand), ...
                "AbsTol", tc.Tol, tc.seedDiagnostic());
        end

        function setResolutionLeavesTheCacheUntouched(tc)
            % The cache holds RAW distances, which do not depend on epsReq. That is
            % what makes a resolution update O(nSup) and what keeps the cache valid
            % across it.
            s     = tc.stream(6);
            Q     = rand(s, 30, 2);
            state = rid.crit.CoverageState(Q, 0.05 * ones(30, 1));
            state.add(rand(s, 8, 2));

            cacheBefore = state.dRaw;
            state.setResolution(0.01 + 0.05 * rand(s, 30, 1));
            tc.verifyEqual(state.dRaw, cacheBefore, tc.seedDiagnostic());
        end

        function setResolutionChangesTheValue(tc)
            % ... and it does change J, so the cache being untouched is not because
            % nothing happened.
            Q     = (0:0.25:1).';
            state = rid.crit.CoverageState(Q, 0.25 * ones(5, 1), "psiMax", 1);
            state.add(0.5);
            coarse = state.value();
            state.setResolution(0.1 * ones(5, 1));
            tc.verifyGreaterThan(state.value(), coarse);
        end

        function rebuildReplacesTheCacheAndCanResetIt(tc)
            s      = tc.stream(7);
            Q      = rand(s, 20, 2);
            state  = rid.crit.CoverageState(Q, 0.05 * ones(20, 1));
            state.add(rand(s, 5, 2));

            % Coordinates move after a surrogate update: the incremental minimum no
            % longer describes the same geometry and must be discarded.
            Zmoved = rand(s, 5, 2);
            state.rebuild(Zmoved);
            reference = rid.crit.CoverageState(Q, 0.05 * ones(20, 1));
            reference.add(Zmoved);
            tc.verifyTrue(isequal(state.dRaw, reference.dRaw), tc.seedDiagnostic());
            tc.verifyEqual(state.nSamples, 5);

            state.rebuild(zeros(0, 2));
            tc.verifyEqual(state.dRaw, inf(20, 1));
            tc.verifyEqual(state.nSamples, 0);
        end

        function isSatisfiedIsTrueExactlyWhenEveryPointIsResolved(tc)
            % the loop must not test J == 0 itself, so the criterion says
            % whether it has certified itself. Certification means EVERY supporting
            % point is within its own epsReq, so one unresolved point is enough to
            % make it false.
            %
            % Fails if: isSatisfied reports on the majority of points, on the mean
            % distance, or on any single point rather than on J.
            Q      = (0:0.25:1).';
            epsReq = 0.2 * ones(5, 1);
            state  = rid.crit.CoverageState(Q, epsReq);

            % Nothing committed: every point is infinitely far away.
            tc.verifyFalse(state.isSatisfied());
            tc.verifyGreaterThan(state.value(), 0);

            % Four of five covered is NOT certification.
            state.add(Q(1:4, :));
            tc.verifyFalse(state.isSatisfied());

            state.add(Q(5, :));
            tc.verifyTrue(state.isSatisfied());
            tc.verifyEqual(state.value(), 0);
        end

        function isSatisfiedIsExactAndAgreesWithTheValueAtTheBoundary(tc)
            % The boundary is dTilde == 1 exactly, which is RESOLVED: psi is zero
            % there and the point is not in the shell. A point just beyond it is
            % not, and no tolerance may absorb the difference -- a `J < tol` test
            % would certify a run that has not finished (the method description refuses the same epsilon in the gradient's zero branch).
            %
            % Fails if: isSatisfied carries a tolerance, or uses a strict `<`
            % against epsReq where the criterion uses `<=`.
            state = rid.crit.CoverageState(0, 1);
            state.add(1);                        % dRaw == epsReq exactly
            tc.verifyTrue(state.isSatisfied());
            tc.verifyEqual(state.value(), 0);

            beyond = rid.crit.CoverageState(0, 1);
            beyond.add(1 + eps(1));
            tc.verifyFalse(beyond.isSatisfied());
            tc.verifyGreaterThan(beyond.value(), 0);
        end

        function isSatisfiedFollowsTheResolutionAndTheReset(tc)
            % Certification is held in units of epsReq, so a finer demand
            % un-certifies a state that was certified (invariant 2), and the
            % budget-mode reset of the coverage reset -- rebuild with an empty matrix -- returns
            % it to the uncovered state.
            %
            % Fails if: isSatisfied is cached at add() time rather than read from
            % the current J, which is exactly the staleness invariant 2 exists for.
            Q     = (0:0.25:1).';
            state = rid.crit.CoverageState(Q, 0.2 * ones(5, 1));
            state.add(Q);
            tc.verifyTrue(state.isSatisfied());

            state.setResolution(0.05 * ones(5, 1));
            tc.verifyTrue(state.isSatisfied(), ...
                "the samples sit ON the supporting points, so dRaw is 0");

            % Move the samples off the supporting points, then sharpen the demand
            % past the distance that remains.
            state.rebuild(Q + 0.1);
            tc.verifyFalse(state.isSatisfied());

            state.rebuild(zeros(0, 1));
            tc.verifyFalse(state.isSatisfied());
            tc.verifyEqual(state.dRaw, inf(5, 1));
        end

        function emptyArgumentsAreNoOps(tc)
            Q     = (0:0.25:1).';
            state = rid.crit.CoverageState(Q, 0.25 * ones(5, 1));
            state.add(0.5);
            before = state.dRaw;

            state.add(zeros(0, 1));
            tc.verifyEqual(state.dRaw, before);
            tc.verifyEqual(state.nSamples, 1);
            tc.verifyEqual(state.probe(zeros(0, 1)), state.value(), "AbsTol", tc.Tol);
        end

        function theSupportingPointsAndTheOptionsAreImmutable(tc)
            % Q is fixed for the duration of a run, and psiMax and pNorm set
            % the criterion's character, so neither may drift mid-run.
            state = rid.crit.CoverageState((0:0.25:1).', 0.25 * ones(5, 1));
            assign = @(name, value) CoverageStateTest.assignTo(state, name, value);
            tc.verifyError(@() assign("Q", zeros(5, 1)),      "MATLAB:class:SetProhibited");
            tc.verifyError(@() assign("psiMax", 3),           "MATLAB:class:SetProhibited");
            tc.verifyError(@() assign("pNorm", 2),            "MATLAB:class:SetProhibited");
            tc.verifyError(@() assign("dRaw", zeros(5, 1)),   "MATLAB:class:SetProhibited");
            tc.verifyError(@() assign("epsReq", ones(5, 1)),  "MATLAB:class:SetProhibited");
            tc.verifyError(@() assign("nSamples", 7),         "MATLAB:class:SetProhibited");
        end

        function rejectsInvalidConstructionAndUpdates(tc)
            Q = (0:0.25:1).';
            tc.verifyError(@() rid.crit.CoverageState(zeros(0, 1), zeros(0, 1)), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.CoverageState(Q, 0.2 * ones(4, 1)), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.CoverageState(Q, zeros(5, 1)), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.CoverageState(Q, 0.2 * ones(5, 1), "psiMax", 0), ...
                "rid:crit:invalidPenaltyClip");
            tc.verifyError(@() rid.crit.CoverageState(Q, 0.2 * ones(5, 1), "pNorm", 0.5), ...
                "rid:crit:invalidNorm");

            state = rid.crit.CoverageState(Q, 0.2 * ones(5, 1));
            tc.verifyError(@() state.setResolution(0.2 * ones(4, 1)), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() state.setResolution([0.2; NaN; 0.2; 0.2; 0.2]), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() state.add([0.5 0.5]),   "rid:crit:dimensionMismatch");
            tc.verifyError(@() state.probe([0.5 0.5]), "rid:crit:dimensionMismatch");
            tc.verifyError(@() state.rebuild([0.5 0.5]), "rid:crit:dimensionMismatch");
        end

    end

    methods (Static, Access = private)

        function assignTo(obj, name, value)
            %ASSIGNTO  Property assignment from outside the class, for the
            %   immutability tests. A nested function inside the test method would
            %   do the same; this keeps the assignment a single, obvious statement.
            obj.(name) = value;
        end

    end

end
