classdef RebuildRuleTest < matlab.unittest.TestCase
%REBUILDRULETEST  The four-row rebuild table, asserted as a table.
%
%   The method description states when re-deriving the historical
%   regressor points is required after a surrogate update. The rule the code must
%   implement is
%
%       rebuild  <=>  surrogate changed  AND  coordinates depend on it
%
%   and "coordinates depend on it" has TWO independent causes: pastSource is
%   "simulated", or the representation reconstructs its coordinates from surrogate
%   parameters. The specification is explicit that the code reads those two flags
%   and NEVER the online/offline mode -- which is why no test in this file mentions
%   a mode.
%
%   The two are deliberately not collapsed into one flag. Collapsing them would
%   make rows 3 and 4 of the table indistinguishable, and
%   `requiresRebuildOnSurrogateUpdate` could then be implemented as
%   `pastSource == "simulated"` with every other test still passing --
%   surrogateDependentCoordinatesForceARebuild is what fails.

    methods (Test)

        function narxCoordinatesNeverDependOnTheSurrogate(tc)
            % "Uses the surrogate" is not "depends on the surrogate": a NARX
            % representation uses it to predict the future output of a candidate,
            % but its coordinates are the raw signals. The flag is false in BOTH
            % pastSource settings -- it is a property of the coordinates, not of
            % where the output came from.
            for src = ["measured", "simulated"]
                rep = rid.rep.Narx(2, 1, "pastSource", src);
                tc.verifyFalse(rep.dependsOnSurrogate, ...
                    sprintf("pastSource = %s", src));
            end
        end

        function measuredNarxNeverNeedsARebuild(tc)
            % Row 2: online, measured, built from raw signals. This is the default
            % configuration, in which rebuild never runs.
            rep = rid.rep.Narx(2, 1, "pastSource", "measured");
            tc.verifyFalse(rep.requiresRebuildOnSurrogateUpdate());
        end

        function simulatedPastForcesARebuild(tc)
            % Row 3: the past outputs are themselves simulated, so refitting moves
            % every historical point.
            rep = rid.rep.Narx(2, 1, "pastSource", "simulated");
            tc.verifyTrue(rep.requiresRebuildOnSurrogateUpdate());
        end

        function surrogateDependentCoordinatesForceARebuild(tc)
            % Row 4, the one that is easy to overlook: the output was MEASURED, yet
            % the coordinates are reconstructed using surrogate parameters -- an
            % estimated state, an OBF pole -- so refitting still moves them.
            rep = SurrogateDependentRep("measured");
            tc.verifyTrue(rep.dependsOnSurrogate);
            tc.verifyTrue(rep.requiresRebuildOnSurrogateUpdate());
        end

        function bothCausesTogetherStillRequireARebuild(tc)
            rep = SurrogateDependentRep("simulated");
            tc.verifyTrue(rep.requiresRebuildOnSurrogateUpdate());
        end

        function simulatedStatesForceARebuildAndMeasuredOnesDoNot(tc)
            % Row 3 and row 2 again, for the second real representation. The
            % flag is FALSE in both settings: that simulated
            % states move with the surrogate is what pastSource already says,
            % and `dependsOnSurrogate` is reserved for states ESTIMATED from
            % measured outputs by an observer, which is not built.
            for src = ["measured", "simulated"]
                rep = rid.rep.StateSpace(2, "pastSource", src);
                tc.verifyFalse(rep.dependsOnSurrogate, src);
                tc.verifyEqual(rep.requiresRebuildOnSurrogateUpdate(), ...
                    src == "simulated", src);
            end
        end

        function theRuleIsExactlyTheDisjunctionOfTheTwoFlags(tc)
            % Exhaustive over the 2x2 of (dependsOnSurrogate, pastSource), with
            % both real representations in it. SurrogateDependentRep stays for
            % row 4: no real representation can exercise it until estimated
            % states exist.
            cases = { rid.rep.Narx(2, 1, "pastSource", "measured"),        false
                      rid.rep.Narx(2, 1, "pastSource", "simulated"),       true
                      rid.rep.StateSpace(2, "pastSource", "measured"),     false
                      rid.rep.StateSpace(2, "pastSource", "simulated"),    true
                      SurrogateDependentRep("measured"),             true
                      SurrogateDependentRep("simulated"),            true };
            for iCase = 1:size(cases, 1)
                rep = cases{iCase, 1};
                tc.verifyEqual(rep.requiresRebuildOnSurrogateUpdate(), ...
                    cases{iCase, 2}, sprintf("case %d", iCase));
                tc.verifyEqual(rep.requiresRebuildOnSurrogateUpdate(), ...
                    rep.dependsOnSurrogate || rep.pastSource == "simulated");
            end
        end

        function anInvalidPastSourceIsRejectedAtConstruction(tc)
            tc.verifyError(@() rid.rep.Narx(2, 1, "pastSource", "online"), ...
                "rid:rep:invalidPastSource");
        end

    end

    methods (Test)   % the rule EXERCISED, with a real representation

        function surrogateUpdateMovesSimulatedPointsAndRebuildRepairsTheCache(tc)
            % ROW 3, END TO END. Every test above asserts what the flags SAY;
            % this one asserts what they are for -- that after a surrogate
            % update the cached nearest-neighbour distances are wrong until
            % `rebuild` runs, and right afterwards. Until a state-space representation existed, no real
            % representation could move a point under a surrogate update at all,
            % so the consequence of getting the rule wrong was untested.
            %
            % THE TWO SURROGATES SHARE A STATE BASIS AND A REST STATE (both are
            % unit-gain second-order systems in position/velocity coordinates),
            % and differ only in bandwidth. That is the honest version of a
            % refit: the same physical states, moved.
            f = RebuildRuleTest.movingPointsFixture();

            % The INPUT column cannot move -- it is recorded data, not a
            % prediction -- and the first state row cannot either, because x(1)
            % is the run's initial state.
            tc.verifyEqual(f.ZB(:, 1), f.ZA(:, 1), ...
                "The input coordinate is recorded data and cannot move.");
            tc.verifyEqual(f.ZB(1, 2:3), f.ZA(1, 2:3), ...
                "Row 1 is x0, which is not simulated.");

            % ... and every LATER state row must actually move, or the fixture
            % would not show the defect it exists to show.
            moved = vecnorm(f.ZB(2:end, 2:3) - f.ZA(2:end, 2:3), 2, 2);
            tc.verifyGreaterThan(min(moved), 1e-3, ...
                "The two surrogates give nearly the same states; the fixture " + ...
                "cannot distinguish a rebuild from a skipped one.");

            % The state the run carries: built under srgA, then srgB is fitted.
            state = rid.crit.CoverageState(f.Q, f.epsReq);
            state.add(f.nrm.toNormalized(f.ZA));
            % A NON-DEGENERATE mix of resolved and unresolved points, or the
            % active-set comparison below proves nothing.
            tc.verifyTrue(any(state.isActive) && ~all(state.isActive), ...
                "epsReq resolves all or none of Q; the active set is vacuous.");

            % The oracle: a state that never saw the old coordinates at all.
            fresh = rid.crit.CoverageState(f.Q, f.epsReq);
            fresh.add(f.nrm.toNormalized(f.ZB));

            % WITHOUT the rebuild the cache is a fact about coordinates that no
            % longer exist. This is the assertion that makes the rest bite.
            tc.verifyFalse(isequal(state.dRaw, fresh.dRaw), ...
                "The stale cache already equals the fresh one, so the test " + ...
                "would pass with rebuild removed.");
            tc.verifyNotEqual(state.value(), fresh.value());

            state.rebuild(f.nrm.toNormalized(f.ZB));
            tc.verifyEqual(state.dRaw,    fresh.dRaw);
            tc.verifyEqual(state.value(), fresh.value());
            tc.verifyEqual(state.isActive, fresh.isActive);
        end

        function measuredStatePointsDoNotMoveWhenTheSurrogateIsRefitted(tc)
            % ROW 2 for the same representation: the states were MEASURED, so a
            % refit cannot move them and no rebuild is owed. The two surrogates
            % of the fixture above give bit-identical points here -- which is
            % also how this test would catch a `build` that simulated anyway and
            % ignored hist.x.
            f   = RebuildRuleTest.movingPointsFixture();
            rep = rid.rep.StateSpace(2, "pastSource", "measured");

            hist   = f.hist;
            hist.x = f.ZA(:, 2:3);          % as if the process had reported them

            tc.verifyEqual(rep.build(hist, f.srgB), rep.build(hist, f.srgA));
            tc.verifyFalse(rep.requiresRebuildOnSurrogateUpdate());
        end

    end

    methods (Static, Access = private)

        function f = movingPointsFixture()
            %MOVINGPOINTSFIXTURE  Two unit-gain second-order models, one history.
            %
            %   xdot = [0 1; -w0^2, -2*zeta*w0] x + [0; w0^2] u,  y = x1
            %
            %   Unit gain, so both models rest at [u, 0] and the state basis --
            %   position and velocity -- is the same one. Only the bandwidth
            %   differs, which is what a refit changes.
            zeta = 0.3;
            tS   = 0.05;
            mk   = @(w0) rid.surrogate.LinearStateSpace.fromContinuous( ...
                       [0 1; -w0^2, -2*zeta*w0], [0; w0^2], [1 0], tS);
            f.srgA = mk(2);
            f.srgB = mk(3);

            % Eight held levels, so the states are genuinely excited rather than
            % averaged away by a white input. The magnitude is bounded away from
            % zero: a segment at u = 0 would drive neither model anywhere and
            % would weaken the "the points moved" assertion for five samples.
            rs     = RandStream("twister", "Seed", 20260917);
            level  = (0.3 + 0.5 * rand(rs, 8, 1)) .* sign(randn(rs, 8, 1));
            u      = repelem(level, 5, 1);

            % Rest at the ORIGIN rather than under u(1): the rest state of u(1)
            % is a fixed point for both models, so the first segment would move
            % neither and the rows would agree by construction.
            f.hist = struct("u", u, "x0", [0 0], "tS", tS);

            rep  = rid.rep.StateSpace(2, "pastSource", "simulated");
            f.ZA = rep.build(f.hist, f.srgA);
            f.ZB = rep.build(f.hist, f.srgB);

            roi  = struct("u1_k", [-1 1], "x1_k", [-1.5 1.5], "x2_k", [-5 5]);
            [box, dimName] = rep.bounds(roi);
            f.nrm = rid.roi.Normalization(box, dimName);

            % A small deterministic grid, interior to the normalized box: the
            % rebuild rule is about the cache, not about how Q was generated,
            % and a grid keeps the fixture free of the supporting-point
            % generators and their toolbox calls.
            g          = linspace(0.1, 0.9, 4);
            [g1, g2, g3] = ndgrid(g, g, g);
            f.Q        = [g1(:), g2(:), g3(:)];
            f.epsReq   = 0.15 * ones(size(f.Q, 1), 1);
        end

    end

end
