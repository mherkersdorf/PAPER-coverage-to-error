classdef CoverageStateOracleTest < matlab.unittest.TestCase
%COVERAGESTATEORACLETEST  rid.crit.CoverageState against the stateless oracle.
%
%   Layer 3 of the project's test plan, and the reason rid.crit.evaluate is permanent.
%
%   The cached path has exactly one failure mode and it is SILENT. A stale cache
%   raises no error; it returns a slightly wrong J, the optimizer picks a slightly
%   worse candidate, and the run looks like a criterion that simply performs badly.
%   Nothing downstream can tell the difference. The only thing that catches it is
%   recomputing the value from scratch, independently, after every operation.
%
%   So: random sequences of add / probe / setResolution / rebuild are applied to a
%   CoverageState while a plain history of the committed samples is kept alongside.
%   After EVERY operation, state.value() is compared with rid.crit.evaluate
%   recomputed from that history, to a relative 1e-12, and every probe is compared
%   with evaluate over history plus candidate. A failure reports the seed, the
%   operation index and the operation that caused it, so the sequence can be
%   replayed exactly.
%
%   The operation mix deliberately includes the two orderings that break a cache
%   holding anything in units of epsReq: setResolution between an add and a value,
%   and setResolution between an add and a probe.
%
%   `rebuild` STAYS IN THE MIX for a different reason. In the default configuration
%   -- online, measured past values, a representation built from raw signals -- it
%   never runs in production, so
%   nothing else exercises it until a run needs simulated past values, a
%   representation with dependsOnSurrogate = true, or a replay. Dropping it here to
%   shorten the sequences would leave that path covered by CoverageStateTest alone.

    properties (Constant)
        % Several independent sequences. Each is a literal, so a failure is
        % reproducible by running this file alone.
        Seeds  = [20260830, 20260831, 20260832, 91, 4242]
        RelTol = 1e-12
        NumOps = 200

        % THE SEQUENCE MUST NOT RESOLVE ITSELF INTO VACUITY. A resolved supporting
        % point contributes zero to the cached path AND zero to the oracle, so once
        % everything is resolved the comparison is 0 == 0 and proves nothing
        % (a randomized sample must be asserted
        % non-degenerate). Measured over these five seeds the median active
        % fraction is 0.47 to 0.75, so 0.25 leaves room without being slack. The
        % draw this replaced fell to 0.06 and 0.08 on two of the five and spent 148
        % of 1000 operations comparing a fully resolved state.
        MinActiveFraction = 0.25
    end

    methods (Test)

        function randomizedOperationSequencesAgreeWithTheOracle(tc)
            for seed = tc.Seeds
                tc.runOneSequence(seed);
            end
        end

        function agreesInTheFineRegimeWhereAlmostNothingResolves(tc)
            % The fully-populated comparison, by construction rather than by luck:
            % at a twentieth of the sample spacing essentially every supporting
            % point stays active for the whole sequence, so every one of them is
            % compared at every operation. The same device as the cost case in
            % CoverageStateGradientTest.
            %
            % Fails if: the cached path is wrong for points that never resolve --
            % which the mixed sequences would increasingly stop noticing as they
            % resolve.
            tc.runOneSequence(5150, "resolutionScale", 0.05);
            tc.runOneSequence(5151, "resolutionScale", 0.05);
        end

        function agreesWhenThePenaltyIsUnclippedAndTheNormIsTheWorstCase(tc)
            % The two limits of the aggregation are separate code paths in both
            % implementations (psiMax = Inf leaves psi unbounded; pNorm = Inf takes
            % its own branch), so they get their own sequences rather than relying
            % on the random draw to hit them.
            tc.runOneSequence(777, "psiMax", Inf, "pNorm", 1);
            tc.runOneSequence(778, "psiMax", 1,   "pNorm", Inf);
            tc.runOneSequence(779, "psiMax", Inf, "pNorm", Inf);
        end

    end

    methods (Access = private)

        function runOneSequence(tc, seed, opts)
            arguments
                tc
                seed          (1,1) double
                opts.psiMax   (1,1) double = NaN     % NaN: drawn from the stream
                opts.pNorm    (1,1) double = NaN
                % Multiplier on the demanded resolution. 1 is the mixed regime;
                % a small value is the fine regime, where almost nothing resolves.
                opts.resolutionScale (1,1) double = 1
            end

            s    = RandStream("twister", "Seed", seed);
            nSup = randi(s, [5 40]);
            nDim = randi(s, [1 3]);
            Q    = rand(s, nSup, nDim);

            psiMax = opts.psiMax;
            pNorm  = opts.pNorm;
            if isnan(psiMax), psiMax = 0.2 + 2 * rand(s);   end
            if isnan(pNorm),  pNorm  = 1 + 3 * rand(s);     end

            epsReq  = spacingScaledResolution(s, nSup, nDim, opts.resolutionScale);
            state   = rid.crit.CoverageState(Q, epsReq, "psiMax", psiMax, "pNorm", pNorm);
            history = zeros(0, nDim);
            activeFraction = zeros(tc.NumOps, 1);

            % The empty state, before anything has happened.
            tc.compare(state, history, epsReq, psiMax, pNorm, seed, 0, "construct");

            for iOp = 1:tc.NumOps
                switch randi(s, [1 4])

                    case 1                                          % add
                        Znew    = rand(s, randi(s, [1 5]), nDim);
                        state.add(Znew);
                        history = [history; Znew]; %#ok<AGROW> the history is the oracle's input
                        opName  = "add";

                    case 2                                          % probe
                        Zcand  = rand(s, randi(s, [1 5]), nDim);
                        jProbe = state.probe(Zcand);
                        jOracle = rid.crit.evaluate([history; Zcand], Q, epsReq, ...
                            "psiMax", psiMax, "pNorm", pNorm);
                        tc.verifyEqual(jProbe, jOracle, "RelTol", tc.RelTol, ...
                            tc.failureDiagnostic(seed, iOp, "probe"));
                        opName = "probe";

                    case 3                                          % setResolution
                        % A refinement-field update mid-sequence: the raw cache must
                        % survive it unchanged, and J must move with it.
                        epsReq = spacingScaledResolution( ...
                                     s, nSup, nDim, opts.resolutionScale);
                        state.setResolution(epsReq);
                        opName = "setResolution";

                    case 4                                          % rebuild
                        % A surrogate update that moved every past coordinate. The
                        % history is replaced wholesale, exactly as +rhc would when
                        % rep.dependsOnSurrogate is true or pastSource is
                        % "simulated".
                        history = rand(s, randi(s, [0 25]), nDim);
                        state.rebuild(history);
                        opName = "rebuild";
                end

                tc.compare(state, history, epsReq, psiMax, pNorm, seed, iOp, opName);
                activeFraction(iOp) = sum(state.isActive) / state.nSup;
            end

            % The comparison must still have something to compare. See
            % MinActiveFraction: this is the assertion that would have failed on
            % the earlier epsReq draw.
            tc.verifyGreaterThanOrEqual(median(activeFraction), ...
                CoverageStateOracleTest.MinActiveFraction, ...
                sprintf("seed %d: median active fraction %.3f, min %.3f, %d of " + ...
                        "%d operations fully resolved -- the oracle is comparing " + ...
                        "zero with zero", seed, median(activeFraction), ...
                        min(activeFraction), sum(activeFraction == 0), tc.NumOps));

            if opts.resolutionScale >= 1
                % ...and points must actually resolve along the way, or the active
                % set is never exercised inside the oracle at all.
                tc.verifyLessThan(min(activeFraction), 0.9, ...
                    sprintf("seed %d: nothing ever resolved", seed));
            end
        end

        function compare(tc, state, history, epsReq, psiMax, pNorm, seed, iOp, opName)
            %COMPARE  state.value() against the oracle recomputed from the history.
            [jOracle, oracleInfo] = rid.crit.evaluate(history, state.Q, epsReq, ...
                "psiMax", psiMax, "pNorm", pNorm);
            jCached = state.value();

            if isinf(jOracle)
                % psiMax = Inf with nothing committed: both must be infinite, and a
                % relative tolerance is meaningless there.
                tc.verifyEqual(jCached, jOracle, ...
                    tc.failureDiagnostic(seed, iOp, opName));
            else
                tc.verifyEqual(jCached, jOracle, "RelTol", tc.RelTol, "AbsTol", 1e-14, ...
                    tc.failureDiagnostic(seed, iOp, opName));
            end

            % The cache itself, not only the value it produces: a compensating pair
            % of errors in the aggregation would otherwise hide a wrong cache.
            tc.verifyEqual(state.dRaw, oracleInfo.dRaw, "RelTol", tc.RelTol, ...
                tc.failureDiagnostic(seed, iOp, opName));
        end

        function d = failureDiagnostic(~, seed, iOp, opName)
            d = sprintf("seed = %d, operation %d (%s)", seed, iOp, opName);
        end

    end

end
