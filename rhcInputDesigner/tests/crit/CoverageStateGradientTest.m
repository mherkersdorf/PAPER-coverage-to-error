classdef CoverageStateGradientTest < matlab.unittest.TestCase
%COVERAGESTATEGRADIENTTEST  The cached gradient against the reference path.
%
%   The method description requires the PRODUCTION gradient to go
%   through the cache. `rid.crit.gradient` takes its value from
%   `rid.crit.evaluate`, which walks every sample -- O(nSup * nSamples), and
%   nSamples grows for the whole run -- so wiring it into the optimizer would
%   violate algorithmic invariant 1 and produce a loop that gets steadily slower
%   for no visible reason.
%
%   Four things are asserted here, and each catches something the others do not:
%
%     1. EQUIVALENCE with rid.crit.gradient over randomized sequences of
%        add / probe / setResolution / rebuild, to 1e-12. Same shape as the
%        existing oracle test, and for the same reason: the cached path's one
%        failure mode is a stale cache, which raises no error and merely degrades
%        the search.
%     2. THE HISTORY CASE. A supporting point whose nearest sample lies in the
%        history contributes EXACTLY ZERO, however deep inside the shell it sits.
%        Constructed deliberately with such points IN the shell, or the case is
%        vacuous -- and a finite-difference reference makes the same restriction
%        automatically, so nothing else exposes an implementation that gets it
%        wrong (section 3).
%     3. COST. The candidate path touches O(nSup * nHorizon) distances and not
%        O(nSup * N), counted through the seam, with a history two orders of
%        magnitude longer than the horizon so the two numbers cannot be confused.
%     4. NO MUTATION. State after probeGradient is bit-identical to before.
%
%   pNorm = Inf is absent on purpose: The method description covers
%   p >= 1 finite, and both implementations now refuse the worst-case limit rather
%   than returning the plausible-looking wrong value it would otherwise produce.

    properties (Constant)
        Seeds  = [20260905, 20260906, 314, 2718, 66]
        RelTol = 1e-12
        NumOps = 60

        % A resolved supporting point contributes zero to the cached gradient AND
        % zero to the reference, so a sequence that resolves itself compares zero
        % with zero. Measured over these five seeds the
        % median active fraction is 0.45 to 0.73. See
        % spacingScaledResolution for why the draw is scaled.
        MinActiveFraction = 0.25
    end

    methods (Test)

        % --- 1: the oracle ------------------------------------------------------

        function randomizedSequencesAgreeWithTheReferenceGradient(tc)
            % The cache is the only thing that can be wrong here, and it fails
            % silently. Recomputing from scratch after every operation is what
            % catches it.
            %
            % Fails if: dRaw is not consulted when deciding whether a candidate
            % wins, or the tie convention differs from the reference's lowest-index
            % rule over [history; candidate].
            for seed = tc.Seeds
                tc.runOneSequence(seed);
            end
        end

        function itAgreesInTheFineRegimeWhereAlmostNothingResolves(tc)
            % Every supporting point active for the whole sequence, so every one is
            % compared at every operation -- by construction, not by luck.
            %
            % Fails if: the cached gradient is wrong for points that never resolve.
            tc.runOneSequence(5152, "resolutionScale", 0.05);
        end

        function itAgreesWhenThePenaltyIsUnclipped(tc)
            % psiMax = Inf removes the upper clip, so the shell becomes
            % dTilde > 1 with no outer boundary and every unresolved point
            % contributes. A separate sequence rather than relying on the draw.
            tc.runOneSequence(4242, "psiMax", Inf);
        end

        % --- 2: the history case ------------------------------------------------

        function supportingPointsOwnedByTheHistoryContributeExactlyZero(tc)
            % Section 3 states this as a REQUIREMENT, not a corollary. The
            % construction puts three supporting points nearest to history samples
            % and squarely INSIDE the shell -- verified before the assertion, or
            % the case would be vacuous -- and one nearest to a candidate.
            %
            % Fails if: the gradient accumulates a term wherever the point is in
            % the shell, regardless of which sample owns it. A finite-difference
            % reference agrees with that mistake, so only this case exposes it.
            % Spaced so that each supporting point's nearest sample is its own
            % companion at distance 0.6, hence dTilde = 1.2, squarely inside the
            % shell (1, 2) at epsReq = 0.5 and psiMax = 1.
            Q      = [0 0; 3 0; 6 0; 9 0];
            epsReq = 0.5 * ones(4, 1);
            Zhist  = [0.6 0; 3.6 0; 6.6 0];     % owns q1, q2, q3
            Zcand  = [9.6 0];                   % owns q4 only

            state = rid.crit.CoverageState(Q, epsReq, "psiMax", 1, "pNorm", 1);
            state.add(Zhist);
            [gradCand, info] = state.probeGradient(Zcand);

            % The regime is asserted, not assumed: all four points are in the shell
            % and only the last is owned by the candidate.
            tc.verifyEqual(info.isShell, true(4, 1), ...
                "every supporting point must be in the shell, or the case is vacuous");
            tc.verifyEqual(info.fromCandidate, [false; false; false; true]);

            % The reference agrees on the value...
            [gradRef, refInfo] = rid.crit.gradient([Zhist; Zcand], Q, epsReq, ...
                "psiMax", 1, "pNorm", 1);
            tc.verifyEqual(info.J, refInfo.J, "RelTol", tc.RelTol);

            % ...and the history rows of the reference carry the three terms this
            % path correctly never forms, while the candidate row matches exactly.
            tc.verifyTrue(any(gradRef(1:3, :) ~= 0, "all"), ...
                "the reference must put nonzero rows on the history, or this " + ...
                "test asserts nothing");
            tc.verifyEqual(gradCand, gradRef(4, :), "RelTol", tc.RelTol);
        end

        % --- 3: the cost --------------------------------------------------------

        function theCandidatePathCostsTheHorizonAndNotTheHistory(tc)
            % Algorithmic invariant 1, counted rather than reviewed. With
            % nSup = 30, nHorizon = 8 and N = 2000 committed samples, the contract
            % allows 240 distances and the naive path would take 60000 -- two
            % orders of magnitude apart, so no accounting detail can confuse them.
            %
            % epsReq is far below any distance 2000 samples leave behind, so NO
            % point resolves and the active set stays full: this test measures the
            % history/horizon split alone, and the active-set saving is measured
            % separately in CoverageStateActiveSetTest. Every point is saturated
            % here, which incidentally exercises the active-set rule that saturated points
            % stay in the set -- if they left, the count would be zero.
            %
            % Fails if: probeGradient recomputes over the history, or calls
            % rid.crit.evaluate for its value.
            rs     = RandStream("twister", "Seed", 99);
            nSup   = 30;
            nHor   = 8;
            nHist  = 2000;
            Q      = rand(rs, nSup, 2);
            epsReq = 1e-4 * ones(nSup, 1);

            state = CountingCoverageState(Q, epsReq, "psiMax", 1, "pNorm", 1);
            state.add(rand(rs, nHist, 2));
            state.reset();                       % the setup is not the measurement

            Zcand = rand(rs, nHor, 2);
            state.probeGradient(Zcand);

            tc.verifyEqual(state.nDistances, nSup * nHor, ...
                "the candidate path must touch nSup*nHorizon distances exactly");
            tc.verifyLessThan(state.nDistances, nSup * nHist / 100);
        end

        function theResultDoesNotDependOnHowLongTheHistoryIs(tc)
            % The behavioural half of the same statement, and it does NOT subsume
            % the count: an implementation that computed all N distances and threw
            % them away would pass this while violating the invariant. Conversely
            % the count alone would not catch a path that read the wrong cache.
            %
            % The short state is rebuilt from the subset of the history that
            % actually attains the minima, so dRaw is bit-identical by
            % construction, and the gradients must be too.
            %
            % Fails if: the gradient reads a history sample's coordinates -- it
            % never needs them (section 3).
            rs     = RandStream("twister", "Seed", 1234);
            nSup   = 25;
            Q      = rand(rs, nSup, 2);
            epsReq = 0.05 + 0.1 * rand(rs, nSup, 1);
            Zhist  = rand(rs, 1500, 2);
            Zcand  = rand(rs, 6, 2);

            long = rid.crit.CoverageState(Q, epsReq, "psiMax", 1, "pNorm", 1);
            long.add(Zhist);

            % The at most nSup samples that own a supporting point.
            [~, owner] = min(rid.ext.pairwiseDistance(Q, Zhist), [], 2);
            short = rid.crit.CoverageState(Q, epsReq, "psiMax", 1, "pNorm", 1);
            short.rebuild(Zhist(unique(owner), :));
            tc.verifyEqual(short.dRaw, long.dRaw, ...
                "the two states must hold the same cache, or the test compares " + ...
                "two different problems");
            tc.verifyLessThan(size(unique(owner), 1), size(Zhist, 1) / 10);

            tc.verifyEqual(short.probeGradient(Zcand), long.probeGradient(Zcand));
        end

        % --- 4: no mutation -----------------------------------------------------

        function probeGradientLeavesTheStateBitIdentical(tc)
            % Loop invariant 2: scoring is side-effect free and only `add` commits.
            % A mutation here would make the candidates depend on the order they
            % were scored in, which is also what makes the candidate loop unsafe
            % to parallelize.
            %
            % Fails if: probeGradient updates dRaw the way `add` does.
            rs     = RandStream("twister", "Seed", 7);
            Q      = rand(rs, 15, 3);
            epsReq = 0.05 + 0.2 * rand(rs, 15, 1);
            state  = rid.crit.CoverageState(Q, epsReq, "psiMax", 1.5, "pNorm", 2);
            state.add(rand(rs, 40, 3));

            before = struct("dRaw", state.dRaw, "epsReq", state.epsReq, ...
                            "nSamples", state.nSamples);
            valueBefore = state.value();

            state.probeGradient(rand(rs, 5, 3));

            tc.verifyEqual(state.dRaw,     before.dRaw);
            tc.verifyEqual(state.epsReq,   before.epsReq);
            tc.verifyEqual(state.nSamples, before.nSamples);
            tc.verifyEqual(state.value(),  valueBefore);
        end

        % --- boundaries ---------------------------------------------------------

        function theWorstCaseNormIsRefusedRatherThanApproximated(tc)
            % At pNorm = Inf the criterion is a max, so the derivative is that of
            % the argmax term alone. The p-norm form would return 1/nSup at the max
            % point where 1 is correct: finite, plausible, and wrong by a factor of
            % nSup. Both paths refuse it identically.
            %
            % Fails if: either implementation evaluates (psi/J)^(pNorm-1) at Inf.
            Q      = [0 0; 1 0];
            epsReq = [0.5; 0.5];
            state  = rid.crit.CoverageState(Q, epsReq, "pNorm", Inf);
            state.add([0.7 0]);

            tc.verifyError(@() state.probeGradient([1.7 0]), ...
                "rid:crit:gradientNormUnsupported");
            tc.verifyError(@() rid.crit.gradient([0.7 0], Q, epsReq, "pNorm", Inf), ...
                "rid:crit:gradientNormUnsupported");
        end

        function anEmptyCandidateBlockGivesAZeroRowGradient(tc)
            % The first step of a run probes nothing, and a shape error there would
            % surface as a dimension mismatch deep inside the composition.
            Q     = [0 0; 1 0];
            state = rid.crit.CoverageState(Q, [0.5; 0.5]);
            state.add([0.7 0]);
            [g, info] = state.probeGradient(zeros(0, 2));
            tc.verifyEqual(size(g), [0 2]);
            tc.verifyEqual(info.J, state.value(), "RelTol", tc.RelTol);
        end

        function everyPointResolvedGivesAnExactlyZeroGradient(tc)
            % J = 0 is the certification termination point, reached by design. The
            % gradient must be finite and exactly zero there, never NaN -- and at
            % pNorm > 1 the naive form is 0/0.
            %
            % Fails if: the J == 0 branch is missing, or carries a tolerance.
            Q      = [0 0; 1 0];
            epsReq = [1; 1];
            state  = rid.crit.CoverageState(Q, epsReq, "pNorm", 2);
            state.add([0.1 0; 0.9 0]);
            [g, info] = state.probeGradient([0.5 0]);

            tc.verifyEqual(info.J, 0);
            tc.verifyEqual(g, zeros(1, 2));
            tc.verifyFalse(any(isnan(g), "all"));
        end

    end

    methods (Access = private)

        function runOneSequence(tc, seed, opts)
            arguments
                tc
                seed        (1,1) double
                opts.psiMax (1,1) double = NaN
                opts.pNorm  (1,1) double = NaN
                opts.resolutionScale (1,1) double = 1
            end
            s    = RandStream("twister", "Seed", seed);
            nSup = randi(s, [5 25]);
            nDim = randi(s, [1 3]);
            Q    = rand(s, nSup, nDim);

            psiMax = opts.psiMax;
            pNorm  = opts.pNorm;
            if isnan(psiMax), psiMax = 0.2 + 2 * rand(s); end
            if isnan(pNorm),  pNorm  = 1 + 3 * rand(s);   end

            epsReq = spacingScaledResolution(s, nSup, nDim, ...
                         opts.resolutionScale);
            state  = rid.crit.CoverageState(Q, epsReq, ...
                         "psiMax", psiMax, "pNorm", pNorm);
            Zhist  = zeros(0, nDim);
            activeFraction = zeros(tc.NumOps, 1);

            for iOp = 1:tc.NumOps
                switch randi(s, 4)
                    case 1
                        Znew  = rand(s, randi(s, [1 5]), nDim);
                        state.add(Znew);
                        Zhist = [Zhist; Znew];               %#ok<AGROW>
                    case 2
                        epsReq = spacingScaledResolution( ...
                                     s, nSup, nDim, opts.resolutionScale);
                        state.setResolution(epsReq);
                    case 3
                        if ~isempty(Zhist)
                            state.rebuild(Zhist);
                        end
                    case 4
                        % nothing: exercises a probe against an unchanged state
                end

                Zcand = rand(s, randi(s, [1 6]), nDim);
                [gradCand, info] = state.probeGradient(Zcand);

                % THE REFERENCE ORDERING IS [history; candidate]. The tie-break of
                % rid.crit.gradient is the lowest row index, so this ordering is
                % what makes "the history wins ties" the equivalent convention.
                [gradRef, refInfo] = rid.crit.gradient([Zhist; Zcand], Q, epsReq, ...
                    "psiMax", psiMax, "pNorm", pNorm);
                nHist = size(Zhist, 1);

                ctx = sprintf("seed %d, op %d", seed, iOp);
                tc.verifyEqual(info.J, refInfo.J, "RelTol", tc.RelTol, ctx);
                tc.verifyEqual(gradCand, gradRef(nHist+1:end, :), ...
                    "AbsTol", 1e-12, "RelTol", tc.RelTol, ctx);
                activeFraction(iOp) = sum(state.isActive) / state.nSup;
            end

            % Non-degeneracy: the sequence must still have unresolved points to
            % compare, or it asserts that zero equals zero.
            tc.verifyGreaterThanOrEqual(median(activeFraction), ...
                CoverageStateGradientTest.MinActiveFraction, ...
                sprintf("seed %d: median active fraction %.3f, min %.3f, %d of " + ...
                        "%d operations fully resolved", seed, ...
                        median(activeFraction), min(activeFraction), ...
                        sum(activeFraction == 0), tc.NumOps));
        end

    end

end
