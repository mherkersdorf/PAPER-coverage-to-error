classdef CoverageStateActiveSetTest < matlab.unittest.TestCase
%COVERAGESTATEACTIVESETTEST  The active set: add maintains it, probe exploits it.
%
%   A supporting point whose psi has reached zero can never contribute
%   again while epsReq and the surrogate hold, so `probe` and `probeGradient`
%   measure distances only to the unresolved ones and a candidate evaluation costs
%   O(nActive * nHorizon).
%
%   THE SPLIT IS THE POINT, and two of the tests below exist only to pin it.
%   `probe` skipping a resolved point is exact -- `min(dRaw, anything) <= epsReq`
%   still gives psi = 0. `add` skipping one is NOT: dRaw would freeze at an upper
%   bound, which is invisible while epsReq holds and wrong the moment it moves.
%   The three-operation case is asserted directly in
%   `resolvedPointsDistanceKeepsShrinkingAcrossAResolutionChange`.
%
%   The other two mistakes this file is built to catch:
%
%     - dropping SATURATED points as well, since they also have zero gradient. The
%       criterion for leaving is monotonicity, not zero gradient: dRaw decreases,
%       so a saturated point can move INTO the shell.
%     - averaging J over the active set instead of over all nSup. Resolved points
%       contribute zero to the numerator and one to the denominator.
%
%   The randomized oracle sequences that guard the whole cache live in
%   CoverageStateOracleTest and CoverageStateGradientTest and now run against the
%   active set as well; what is here is the constructed cases those would only
%   reach by luck.

    methods (Test)

        % --- the split: what add measures and what probe measures ---------------

        function probeMeasuresOnlyUnresolvedPointsWhileAddMeasuresEveryOne(tc)
            % Counted through the distanceMatrix seam, with a KNOWN number of
            % points resolving, so the drop is an exact figure and not a trend.
            % Twenty supporting points, five resolved by construction: probe falls
            % from 20 to 15 per candidate while add stays at 20 per sample.
            %
            % Fails if: probe measures the whole set (no drop), or add exploits the
            % set (its count would fall too, and dRaw would start freezing).
            nSup   = 20;
            Q      = [(0:nSup-1).', zeros(nSup, 1)];
            epsReq = 0.3 * ones(nSup, 1);
            state  = CountingCoverageState(Q, epsReq, "psiMax", 1, "pNorm", 1);
            Zcand  = [0.4 0; 1.4 0; 2.4 0; 3.4 0];

            state.reset();
            state.probe(Zcand);
            tc.verifyEqual(state.nDistances, nSup * size(Zcand, 1), ...
                "with nothing resolved, probe must measure every point");

            % Five samples placed exactly on five supporting points resolve those
            % five and nothing else.
            state.reset();
            state.add(Q(1:5, :));
            tc.verifyEqual(state.nDistances, nSup * 5, ...
                "add must measure every point, resolved or not");
            tc.verifyEqual(sum(state.isActive), nSup - 5);

            state.reset();
            state.probe(Zcand);
            tc.verifyEqual(state.nDistances, (nSup - 5) * size(Zcand, 1), ...
                "probe must now skip the five resolved points");

            % More samples near the ALREADY resolved points: nothing new resolves,
            % and add's per-sample count is unchanged. This is what pins "add's
            % count stays flat" as a property rather than a coincidence.
            state.reset();
            state.add([0.05 0; 1.05 0; 2.05 0]);
            tc.verifyEqual(state.nDistances, nSup * 3);
            tc.verifyEqual(sum(state.isActive), nSup - 5);
        end

        function theSetEmptiesExactlyWhenTheCriterionReachesZero(tc)
            % The active-set rule's consequence: in certification mode J = 0 means every point is
            % resolved, so the design gets cheaper the closer it is to finishing.
            % At J = 0 a probe measures nothing at all.
            %
            % Fails if: membership is not recomputed by add, or uses a strict
            % inequality on the wrong side -- a point at dTilde = 1 exactly has
            % psi = 0 and is resolved.
            Q     = [0 0; 1 0; 2 0];
            state = CountingCoverageState(Q, ones(3, 1), "psiMax", 1, "pNorm", 1);

            % Distances of exactly epsReq: the boundary case, resolved.
            state.add([1 0]);
            tc.verifyEqual(state.value(), 0);
            tc.verifyFalse(any(state.isActive));

            state.reset();
            tc.verifyEqual(state.probe([0.5 0]), 0);
            tc.verifyEqual(state.nDistances, 0, ...
                "an empty active set must cost no distances at all");
        end

        % --- why add may not exploit the set ------------------------------------

        function resolvedPointsDistanceKeepsShrinkingAcrossAResolutionChange(tc)
            % THE THREE-OPERATION CASE. add at 0.9 resolves the point at epsReq = 1;
            % add at 0.5 must still be measured, or dRaw freezes at 0.9. The freeze
            % is invisible until epsReq moves, and then it is wrong: at
            % setResolution(0.6) a frozen 0.9 un-resolves the point and gives
            % J = 0.5 where the truth is 0.
            %
            % Fails if: add measures only the active set. Both branches below are
            % needed -- the first catches the freeze, the second exercises a
            % genuine un-resolution so the invalidation path is not untested.
            Q      = [0 0; 5 0];
            state  = rid.crit.CoverageState(Q, [1; 1], "psiMax", 1, "pNorm", 1);

            state.add([0.9 0]);
            tc.verifyEqual(state.isActive, [false; true], ...
                "0.9 <= epsReq = 1 resolves the first point");
            state.add([0.5 0]);
            Zall = [0.9 0; 0.5 0];

            % Still resolved at 0.6 -- but only because dRaw is 0.5 and not the
            % frozen 0.9. This is the assertion the skipping implementation fails.
            state.setResolution([0.6; 1]);
            tc.verifyEqual(state.isActive, [false; true]);
            tc.verifyEqual(state.value(), ...
                rid.crit.evaluate(Zall, Q, [0.6; 1], "psiMax", 1, "pNorm", 1), ...
                "RelTol", 1e-12);

            % Now genuinely un-resolved: 0.5 > 0.4. dTilde = 1.25, psi = 0.25, and
            % the second point is saturated at psi = 1, so J = 0.625 exactly.
            state.setResolution([0.4; 1]);
            tc.verifyEqual(state.isActive, [true; true], ...
                "a shrinking epsReq must put the point back in the active set");
            tc.verifyEqual(state.value(), 0.625, "RelTol", 1e-12);
            tc.verifyEqual(state.value(), ...
                rid.crit.evaluate(Zall, Q, [0.4; 1], "psiMax", 1, "pNorm", 1), ...
                "RelTol", 1e-12);
        end

        function aRebuildThatUnResolvesAPointIsFollowedByTheActiveSet(tc)
            % Invariant 3: a surrogate update moves historical coordinates, so dRaw
            % can come back LARGER and a resolved point can un-resolve. Here the
            % rebuild flips both points at once, in opposite directions.
            %
            % Fails if: rebuild recomputes the cache but leaves membership behind.
            Q     = [0 0; 5 0];
            state = rid.crit.CoverageState(Q, [1; 1], "psiMax", 1, "pNorm", 1);
            state.add([0.5 0]);
            tc.verifyEqual(state.isActive, [false; true]);

            Zmoved = [3 0; 5.2 0];
            state.rebuild(Zmoved);

            tc.verifyEqual(state.isActive, [true; false], ...
                "the rebuild un-resolves the first point and resolves the second");
            tc.verifyEqual(state.value(), 0.5, "RelTol", 1e-12);
            tc.verifyEqual(state.value(), ...
                rid.crit.evaluate(Zmoved, Q, [1; 1], "psiMax", 1, "pNorm", 1), ...
                "RelTol", 1e-12);
        end

        % --- saturated is not symmetric with resolved ----------------------------

        function aSaturatedPointStaysInTheSetAndContributesOnEnteringTheShell(tc)
            % THE TEMPTING MISTAKE. A point beyond 1 + psiMax has zero gradient and
            % a constant contribution, so it looks as skippable as a resolved one.
            % It is not: dRaw only decreases, so it can move into the shell.
            %
            % Fails if: the active set is "points with a nonzero gradient" rather
            % than "points that can still change". The mutant measures no distance,
            % leaves J at psiMax and returns a zero gradient -- all three
            % assertions below fail together.
            Q     = [0 0];
            state = CountingCoverageState(Q, 1, "psiMax", 1, "pNorm", 1);

            state.add([2.5 0]);                 % dTilde = 2.5 >= 1 + psiMax
            tc.verifyEqual(state.value(), 1, "the point is saturated");
            tc.verifyTrue(state.isActive, ...
                "saturated is not resolved: dRaw can still shrink");

            state.reset();
            J = state.probe([1.5 0]);
            tc.verifyEqual(state.nDistances, 1, ...
                "the saturated point must still be measured");
            tc.verifyEqual(J, 0.5, "RelTol", 1e-12, ...
                "entering the shell must lower J");

            g = state.probeGradient([1.5 0]);
            tc.verifyTrue(any(g ~= 0), ...
                "a point that has entered the shell must contribute a gradient");
        end

        % --- the bookkeeping trap ------------------------------------------------

        function theCriterionStaysAMeanOverEverySupportingPoint(tc)
            % The active-set rule's bookkeeping trap: resolved points contribute zero to the
            % numerator and ONE to the denominator. Two of these four points are
            % resolved and two sit at psi = 0.5 exactly, so J = 0.25.
            %
            % Fails if: J is averaged over the active set, which gives 0.5 --
            % double, and plausible.
            Q      = [0 0; 10 0; 20 0; 30 0];
            epsReq = ones(4, 1);
            Zall   = [0 0; 10 0; 21.5 0; 31.5 0];
            state  = rid.crit.CoverageState(Q, epsReq, "psiMax", 1, "pNorm", 1);
            state.add(Zall);

            tc.verifyEqual(sum(state.isActive), 2, ...
                "exactly two points must remain unresolved, or the case is vacuous");
            tc.verifyEqual(state.value(), 0.25, "RelTol", 1e-12);
            tc.verifyEqual(state.value(), ...
                rid.crit.evaluate(Zall, Q, epsReq, "psiMax", 1, "pNorm", 1), ...
                "RelTol", 1e-12);
        end

        function probeAndProbeGradientLeaveMembershipUntouched(tc)
            % Loop invariant 2 extended to the new state: scoring is side-effect
            % free, and membership is state.
            %
            % Fails if: probe calls refreshActive -- candidates would then resolve
            % points they were never committed to.
            Q      = [0 0; 5 0];
            state  = rid.crit.CoverageState(Q, [1; 1], "psiMax", 1, "pNorm", 1);
            state.add([2 0]);
            before = state.isActive;

            state.probe([0.1 0; 5.05 0]);
            state.probeGradient([0.1 0; 5.05 0]);

            tc.verifyEqual(state.isActive, before, ...
                "a candidate must not resolve anything");
        end

    end

end
