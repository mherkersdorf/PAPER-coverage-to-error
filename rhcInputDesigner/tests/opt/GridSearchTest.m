classdef GridSearchTest < matlab.unittest.TestCase
%GRIDSEARCHTEST  The exhaustive oracle's own behaviour.
%
%   `rid.opt.GridSearch` is the correctness oracle of the method description, so its agreement with MultiStart is asserted in
%   OptimizerContractTest. What is asserted here is what makes it usable AS an
%   oracle: that it enumerates the declared feasible set and nothing else, that it
%   refuses a problem it should not be asked to enumerate, and that it reports the
%   absence of starts rather than inventing values for them.

    methods (Test)

        function itScoresEveryFeasibleGridPointAndNoInfeasibleOne(tc)
            % Loop invariant 4: no infeasible candidate is ever scored. The rate
            % band couples the two variables, so it cannot shape the axes and must
            % act as a filter -- and the filter must run BEFORE scoring.
            %
            % Fails if: the linear constraints are applied after scoring, or not
            % at all. With a rate box of [-1 1] over a [0 4] amplitude box, 9 of
            % the 25 grid points survive; without the filter, 25 would be scored.
            gen  = OptTestProblem.generator(2, [0 4], [-1 1]);
            hist = wellFormedHistory(2);
            prob = OptTestProblem(@(t) sum(t.^2), @(t) 2 * t, 0);

            [best, diag] = rid.opt.GridSearch("nLevels", 5).solve( ...
                gen, hist, @(t) prob.withGradient(t));

            [lb, ub, A, b] = gen.feasibleSet(hist);
            tc.verifyLessThan(diag.startInfo.nFeasible, diag.startInfo.nGridPoints);
            tc.verifyEqual(numel(diag.scores), diag.startInfo.nFeasible);
            % One further call is the single question asked at the optimum for J
            % and the penalty; the hot path is the rest.
            tc.verifyEqual(diag.nEvaluations, diag.startInfo.nFeasible + 1);
            tc.verifyEqual(prob.count, diag.nEvaluations);
            tc.verifyTrue(all(best.theta >= lb) && all(best.theta <= ub));
            tc.verifyTrue(all(A * best.theta <= b + 1e-12));
        end

        function aVariablePinnedByTheRateLimitContributesOneLevelNotFive(tc)
            % feasibleSet intersects the first variable's box with u(k-1) + rate.
            % With the actuator at the bottom of the amplitude box and a rate
            % limit that cannot open, that intersection is the single point 0 --
            % and five identical levels there would multiply the work by five
            % while changing nothing. The declared rate box is still a proper one:
            % rid.gen.Constraints rejects a degenerate [0, 0] declaration, so the
            % pinning has to come from the intersection, which is where it comes
            % from in a real run.
            %
            % Fails if: the axis is always linspace(lb, ub, nLevels) -- the grid
            % would then have 25 points rather than 5.
            gen  = OptTestProblem.generator(2, [0 4], [-0.5 0]);
            hist = wellFormedHistory(0);
            prob = OptTestProblem(@(t) sum((t - 9).^2), @(t) 2 * (t - 9), 0);

            [lb, ub] = gen.feasibleSet(hist);
            tc.verifyEqual(lb(1), ub(1), "the first variable must be pinned");

            [best, diag] = rid.opt.GridSearch("nLevels", 5).solve( ...
                gen, hist, @(t) prob.withGradient(t));

            tc.verifyEqual(diag.startInfo.nGridPoints, 5, ...
                "one collapsed axis times five levels, not 25");
            tc.verifyEqual(best.theta, [0; 0], "AbsTol", 1e-12);
        end

        function itRefusesAGridItShouldNotBeAskedToEnumerate(tc)
            % The default problem is continuous and enormous; the oracle exists
            % for a DELIBERATELY tiny one. An exhaustive search that quietly took
            % an hour would be indistinguishable from a hung run.
            %
            % Fails if: maxGridPoints is advisory, or the size is computed after
            % the grid is built.
            gen  = OptTestProblem.generator(8, [0 4], [-4 4]);
            hist = wellFormedHistory(2);
            tc.verifyError(@() rid.opt.GridSearch("nLevels", 5, ...
                "maxGridPoints", 1000).solve(gen, hist, @(t) sum(t.^2)), ...
                "rid:opt:gridTooLarge");
        end

        function itReportsNoStartsRatherThanNoDeadStarts(tc)
            % An exhaustive search has no starting points, so the dead-on-arrival
            % count does not apply. Zero would be the claim that no start was
            % dead, which is a different statement.
            %
            % Fails if: deadOnArrival is initialized to 0 and left there.
            gen  = OptTestProblem.generator(2, [0 4], [-4 4]);
            hist = wellFormedHistory(2);
            [~, diag] = rid.opt.GridSearch("nLevels", 3).solve( ...
                gen, hist, @(t) sum(t.^2));

            tc.verifyTrue(isnan(diag.deadOnArrival));
            tc.verifyEqual(diag.nStarts, 0);
            tc.verifyEqual(diag.solver, "rid.opt.GridSearch");
        end

        function aStartsProviderIsAcceptedAndIgnored(tc)
            % The loop hands every optimizer the same four arguments. An
            % exhaustive search has no use for starting points, and taking them
            % would suggest the result depended on them.
            %
            % Fails if: solve rejects the fourth argument, or seeds its grid from it.
            gen  = OptTestProblem.generator(2, [0 4], [-4 4]);
            hist = wellFormedHistory(2);
            sol  = rid.opt.GridSearch("nLevels", 3);

            withoutStarts = sol.solve(gen, hist, @(t) sum(t.^2));
            withStarts    = sol.solve(gen, hist, @(t) sum(t.^2), ...
                                OptTestProblem.fixedStarts([0.1 0.2; 3.9 3.8]));

            tc.verifyEqual(withStarts.theta, withoutStarts.theta);
        end

    end

end
