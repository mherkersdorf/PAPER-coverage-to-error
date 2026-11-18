classdef OptimizerContractTest < matlab.unittest.TestCase
%OPTIMIZERCONTRACTTEST  The four-point contract of the method description.
%
%   Each test names the input that would make it fail, because a contract test
%   that no wrong implementation fails is decoration. In particular:
%
%     - the ORACLE problem is built so that the optimum is at NO starting point.
%       Without that, an optimizer that never moves passes item 4.
%     - DETERMINISM is asserted in both directions. Same seed identical is half
%       the statement; without "different seed distinct", an optimizer that
%       ignored its starts entirely would pass.
%     - DEAD ON ARRIVAL is asserted against an EXACT expected count. A test that
%       checks only `>= 0` is vacuous, and this diagnostic exists precisely
%       because the count is informative. It is READ from
%       `info.nShell`, so the objective reports it and the optimizer must not
%       second-guess it: the test therefore also pins a start whose OBJECTIVE
%       gradient is nonzero while its criterion gradient is zero, which is the
%       case that separates "read" from "estimate by differencing".
%     - the CONSTRAINT case has its unconstrained optimum outside the feasible
%       set, so a solver that ignored lb/ub/A/b must fail rather than merely
%       return something admissible by luck.
%
%   The oracle problem is deliberately SMOOTH and unimodal. Item 4 tests the
%   SEARCH against exhaustive enumeration, not the criterion's landscape; the flat
%   regions that make the real landscape hard are tested where they belong, in the
%   dead-on-arrival case.

    properties (Constant)
        % A positive-definite quadratic, so the minimum is unique and exact.
        M      = [2 0.5; 0.5 1]
        Offset = 12345.678       % the secret of OptTestProblem
        Tol    = 1e-6            % the solver's own OptimalityTolerance scale
    end

    methods (Test)

        % --- item 4: agreement with exhaustive search --------------------------

        function multiStartFindsTheGridSearchOptimumOnATinyProblem(tc)
            % The oracle of section 4 item 4. Two decision variables on a coarse
            % five-level grid -- the only regime where exhaustive search is
            % possible at all, since the default problem is continuous.
            %
            % The minimum sits ON a grid point, at [1; 3], and at NO start: a
            % solver that returned its own starting point, or that never moved,
            % fails here.
            %
            % Fails if: the local search ignores the gradient sign, or the grid is
            % not built from the declared feasible set.
            c      = [1; 3];
            gen    = OptTestProblem.generator(2, [0 4], [-4 4]);
            hist   = wellFormedHistory(2);
            starts = [0.2 0.7; 3.5 0.4; 0.8 3.9; 2.5 2.5];
            tc.verifyGreaterThan(min(vecnorm(starts.' - c, 2, 1)), 0.5, ...
                "the optimum must not be at a start, or a solver that never " + ...
                "moves would pass");

            [bestMs, diagMs] = rid.opt.MultiStart("nStarts", size(starts, 1)).solve( ...
                gen, hist, tc.scoreHandle(c), OptTestProblem.fixedStarts(starts));
            [bestGs, diagGs] = rid.opt.GridSearch("nLevels", 5).solve( ...
                gen, hist, tc.scoreHandle(c));

            tc.verifyEqual(bestGs.theta, c, "AbsTol", 1e-12, ...
                "the grid contains the optimum exactly");
            tc.verifyEqual(bestMs.theta, bestGs.theta, "AbsTol", 1e-5);
            tc.verifyEqual(bestMs.score, bestGs.score, "RelTol", 1e-9);
            tc.verifyEqual(diagMs.gradientPath, "analytic");
            tc.verifyEqual(diagGs.gradientPath, "none");
        end

        function bothAgreeWhenTheOptimumIsPushedOntoTheBoundary(tc)
            % The same oracle with the unconstrained minimum outside the box, so
            % the answer is a CORNER of the feasible set rather than an interior
            % stationary point. Both must still agree.
            %
            % Fails if: either optimizer treats the bounds as advisory, or the
            % grid is built over a box that is not the declared one.
            c      = [9; 9];
            gen    = OptTestProblem.generator(2, [0 4], [-4 4]);
            hist   = wellFormedHistory(2);
            starts = [0.2 0.7; 3.5 0.4; 0.8 3.9; 2.5 2.5];

            bestMs = rid.opt.MultiStart("nStarts", size(starts, 1)).solve( ...
                gen, hist, tc.scoreHandle(c), OptTestProblem.fixedStarts(starts));
            [bestGs, diagGs] = rid.opt.GridSearch("nLevels", 5).solve( ...
                gen, hist, tc.scoreHandle(c));

            tc.verifyEqual(bestGs.theta, [4; 4], "AbsTol", 1e-12);
            tc.verifyEqual(bestMs.theta, bestGs.theta, "AbsTol", 1e-5);
            tc.verifyTrue(diagGs.hitBound);
        end

        % --- item 1: scoreFcn is a black box -----------------------------------

        function everyEvaluationGoesThroughTheSuppliedScoreFcn(tc)
            % The objective counts its own calls and adds a secret offset. An
            % optimizer that reimplemented the criterion could not produce that
            % offset, and one that evaluated anything outside the handle -- or
            % cached calls it claimed to make -- could not match the count.
            %
            % Fails if: the optimizer inspects the criterion, or reports the
            % solver's own funcCount instead of the calls it actually made.
            c      = [1; 3];
            gen    = OptTestProblem.generator(2, [0 4], [-4 4]);
            hist   = wellFormedHistory(2);
            starts = [0.2 0.7; 3.5 0.4; 0.8 3.9; 2.5 2.5];
            prob   = tc.problem(c);

            [best, diag] = rid.opt.MultiStart("nStarts", size(starts, 1)).solve( ...
                gen, hist, @(t) prob.withGradient(t), OptTestProblem.fixedStarts(starts));

            tc.verifyEqual(diag.nEvaluations, prob.count, ...
                "nEvaluations must equal the calls the objective actually saw");
            tc.verifyGreaterThan(diag.nEvaluations, size(starts, 1));
            tc.verifyEqual(best.score, OptimizerContractTest.Offset, ...
                "AbsTol", 1e-8, ...
                "the best score must carry the offset no reimplementation knows");
        end

        function theCriterionAndThePenaltyAreReportedSeparatelyAtTheOptimum(tc)
            % Contract item 3 requires the penalty contribution at the optimum,
            % and the criterion/objective split requires J and the penalty to be logged apart: neither is
            % recoverable from the objective value alone.
            %
            % Fails if: the optimizer reports the objective as J, or never asks
            % for the third output.
            c    = [1; 3];
            gen  = OptTestProblem.generator(2, [0 4], [-4 4]);
            hist = wellFormedHistory(2);
            prob = tc.problem(c);

            [best, diag] = rid.opt.MultiStart("nStarts", 4).solve( ...
                gen, hist, @(t) prob.withGradient(t), ...
                OptTestProblem.fixedStarts([0.2 0.7; 3.5 0.4; 0.8 3.9; 2.5 2.5]));

            tc.verifyTrue(diag.penaltyReported);
            tc.verifyEqual(diag.penalty, OptimizerContractTest.Offset, "AbsTol", 1e-12);
            tc.verifyEqual(diag.J, best.score - OptimizerContractTest.Offset, ...
                "AbsTol", 1e-8);
        end

        function anObjectiveWithoutTheThirdOutputLeavesTheAbsenceVisible(tc)
            % The optimizer can only report what the objective tells it. When the
            % third output does not exist the fields stay NaN and the flag says
            % so, rather than a zero penalty being invented.
            %
            % Fails if: a missing info struct is reported as penalty = 0.
            gen  = OptTestProblem.generator(2, [0 4], [-4 4]);
            hist = wellFormedHistory(2);
            prob = tc.problem([1; 3]);

            [~, diag] = rid.opt.MultiStart("nStarts", 4).solve( ...
                gen, hist, @(t) prob.valueOnly(t), ...
                OptTestProblem.fixedStarts([0.2 0.7; 3.5 0.4; 0.8 3.9; 2.5 2.5]));

            tc.verifyFalse(diag.penaltyReported);
            tc.verifyTrue(isnan(diag.penalty));
            tc.verifyTrue(isnan(diag.J));
        end

        % --- item 2: deterministic given seed ----------------------------------

        function theSameSeedGivesTheSameSearchAndADifferentSeedDoesNot(tc)
            % BOTH halves are load-bearing. Same-seed-identical alone is passed by
            % an optimizer that ignores its starts entirely and always runs the
            % same single search; the per-start scores are what expose that.
            %
            % The objective is non-convex on purpose, so different starts land in
            % different basins and the spread is real rather than numerical.
            %
            % Fails if: the optimizer reads the global RNG stream, or discards the
            % start set.
            gen  = OptTestProblem.generator(3, [-2 3], [-1 1]);
            hist = wellFormedHistory(0);

            dA = tc.wigglyRun(gen, hist, 11);
            dB = tc.wigglyRun(gen, hist, 11);
            dC = tc.wigglyRun(gen, hist, 12);

            tc.verifyEqual(dA.scores, dB.scores, ...
                "the same seed must reproduce every start exactly");
            tc.verifyEqual(dA.bestScore, dB.bestScore);
            tc.verifyNotEqual(dA.scores, dC.scores, ...
                "a different seed must move the starts; an optimizer that " + ...
                "ignored them would return the same scores");
        end

        % --- item 3: the dead-on-arrival diagnostic ----------------------------

        function deadOnArrivalCountsExactlyTheStartsWithAnEmptyShell(tc)
            % a start whose trajectory
            % leaves every supporting point outside the shell cannot move, and the
            % search reports convergence at iteration zero having done nothing.
            % `grad J` is a sum over shell points, so `nShell == 0` IS a zero
            % criterion gradient -- exactly, by construction.
            %
            % Two of the four starts below have an empty shell and two do not, so
            % the expected count is 2 and nothing else. A test asserting only
            % `>= 0` would be vacuous.
            %
            % Fails if: the count is derived from the objective gradient rather
            % than read from info.nShell -- see the next test, which is what makes
            % that failure unavoidable rather than incidental.
            gen    = OptTestProblem.generator(2, [-2 2], [-4 4]);
            hist   = wellFormedHistory(0);
            starts = [ 0.0  0.0;      % empty shell
                       0.5 -0.5;      % empty shell
                       2.0  2.0;      % in the shell
                      -1.5  1.8];     % in the shell
            prob   = tc.flatProblem();

            [~, diag] = rid.opt.MultiStart("nStarts", 4).solve( ...
                gen, hist, @(t) prob.withGradient(t), ...
                OptTestProblem.fixedStarts(starts));

            tc.verifyTrue(diag.deadOnArrivalReported);
            tc.verifyEqual(diag.deadOnArrival, 2);
        end

        function aStartDeadInTheCriterionCountsEvenWhereThePenaltyStillMoves(tc)
            % THIS is the test that forces "read, not estimate". The objective here
            % has a NONZERO gradient at both flat-region starts -- a constant slope
            % standing in for a penalty pushing toward feasibility -- while
            % info.nShell says the criterion is flat there. Dead-on-arrival is a
            % statement about the CRITERION, which is what psiMax governs,
            %   so the count must still be 2.
            %
            % Fails if: the count is `all(g == 0)` on the objective gradient, which
            % returns 0 here -- plausible, and answering the wrong question.
            gen    = OptTestProblem.generator(2, [-2 2], [-4 4]);
            hist   = wellFormedHistory(0);
            starts = [0.0 0.0; 0.5 -0.5; 2.0 2.0; -1.5 1.8];
            inShell = @(t) double(any(abs(t) > 1));
            prob   = OptTestProblem( ...
                @(t) sum(max(0, abs(t) - 1).^2) + 0.05 * sum(t), ...
                @(t) 2 * max(0, abs(t) - 1) .* sign(t) + 0.05, 0, inShell);

            [~, diag] = rid.opt.MultiStart("nStarts", 4).solve( ...
                gen, hist, @(t) prob.withGradient(t), ...
                OptTestProblem.fixedStarts(starts));

            tc.verifyEqual(diag.deadOnArrival, 2);
        end

        function anObjectiveThatCannotReportTheShellLeavesTheCountUnstated(tc)
            % The optimizer can only report what the objective tells it, and
            % differencing would answer a different question (the objective's
            % gradient, not the criterion's). NaN with the flag false says
            % "unstated"; 0 would claim that no start was dead.
            %
            % Fails if: a missing info struct is reported as deadOnArrival = 0, or
            % if the optimizer falls back to a difference estimate.
            gen    = OptTestProblem.generator(2, [-2 2], [-4 4]);
            hist   = wellFormedHistory(0);
            starts = [0.0 0.0; 0.5 -0.5; 2.0 2.0; -1.5 1.8];
            prob   = tc.flatProblem();

            [~, diag] = rid.opt.MultiStart("nStarts", 4).solve( ...
                gen, hist, @(t) prob.valueOnly(t), ...
                OptTestProblem.fixedStarts(starts));

            tc.verifyFalse(diag.deadOnArrivalReported);
            tc.verifyTrue(isnan(diag.deadOnArrival));
            tc.verifyEqual(diag.gradientPath, "finiteDifference");
        end

        % --- the feasible set ---------------------------------------------------

        function theOptimumSatisfiesEveryDeclaredConstraint(tc)
            % The unconstrained minimum is at [10; 10], far outside the box AND
            % outside the rate band, so a solver that ignored either must fail.
            % The constrained optimum is computable by hand: the rate limit caps
            % the first sample at u(k-1) + duMax = 0.4 and the second at 0.8.
            %
            % Fails if: A and b are dropped, or only the bounds are passed on.
            gen  = OptTestProblem.generator(2, [-2 3], [-0.4 0.4]);
            hist = wellFormedHistory(0);
            [lb, ub, A, b] = gen.feasibleSet(hist);

            best = rid.opt.MultiStart("nStarts", 4, "seed", 7).solve( ...
                gen, hist, tc.scoreHandle([10; 10]));

            tol = 1e-8;
            tc.verifyTrue(all(best.theta >= lb - tol), "lower bound violated");
            tc.verifyTrue(all(best.theta <= ub + tol), "upper bound violated");
            tc.verifyTrue(all(A * best.theta <= b + tol), "rate band violated");
            tc.verifyEqual(best.theta, [0.4; 0.8], "AbsTol", 1e-5, ...
                "the rate-limited corner is the constrained optimum");
        end

        function anInfeasibleStartIsRefusedAndNeverRepaired(tc)
            % starts are feasible BY
            % CONSTRUCTION and there is no repair path, so a violation means the
            % construction is wrong. Projecting it would hide that.
            %
            % Fails if: the optimizer clamps theta0 into the box before searching,
            % which is what the underlying solver would do on its own.
            gen  = OptTestProblem.generator(2, [0 4], [-4 4]);
            hist = wellFormedHistory(2);
            bad  = OptTestProblem.fixedStarts([0.5 0.5; 99 0.5]);

            tc.verifyError(@() rid.opt.MultiStart("nStarts", 2).solve( ...
                gen, hist, tc.scoreHandle([1; 3]), bad), "rid:opt:infeasiblePoint");
        end

        % --- the gradient-free path --------------------------------------------

        function theSolverConvergesAndSaysSoWithoutAnAnalyticGradient(tc)
            % Tier 3 of the method description: with no gradient
            % available the solver differences the objective itself. It must still
            % converge, and the diagnostics must SAY which path was taken so a
            % fall-back cannot masquerade as the fast path.
            %
            % Fails if: the optimizer assumes a second output exists (an error
            % rather than a fall-back), or reports "analytic" regardless.
            c      = [1; 3];
            gen    = OptTestProblem.generator(2, [0 4], [-4 4]);
            hist   = wellFormedHistory(2);
            starts = [0.2 0.7; 3.5 0.4; 0.8 3.9; 2.5 2.5];
            prob   = tc.problem(c);

            [best, diag] = rid.opt.MultiStart("nStarts", 4).solve( ...
                gen, hist, @(t) prob.valueOnly(t), OptTestProblem.fixedStarts(starts));

            tc.verifyEqual(diag.gradientPath, "finiteDifference");
            tc.verifyEqual(best.theta, c, "AbsTol", 1e-4);
            tc.verifyEqual(diag.nEvaluations, prob.count);
        end

    end

    methods (Access = private)

        function prob = problem(tc, c)
            %PROBLEM  The quadratic (t-c)'M(t-c) plus the secret offset.
            M = OptimizerContractTest.M;
            prob = OptTestProblem(@(t) (t - c).' * M * (t - c), ...
                                  @(t) 2 * M * (t - c), ...
                                  OptimizerContractTest.Offset);
        end

        function h = scoreHandle(tc, c)
            %SCOREHANDLE  A fresh quadratic objective, [f, g, info].
            prob = tc.problem(c);
            h    = @(t) prob.withGradient(t);
        end

        function prob = flatProblem(tc) %#ok<MANU>
            %FLATPROBLEM  Exactly zero, with exactly zero gradient, inside |t| <= 1.
            %
            %   The analogue of the criterion's resolved and saturated regions:
            %   `max(0, |t| - 1)^2` is not merely small there, it is the constant
            %   zero, so both the analytic gradient and any difference quotient
            %   are exact zeros.
            excess = @(t) max(0, abs(t) - 1);
            % nShell mirrors the flat region exactly: no coordinate outside the
            % box means no supporting point in the shell.
            prob   = OptTestProblem( ...
                @(t) sum(excess(t).^2), ...
                @(t) 2 * excess(t) .* sign(t), 0, ...
                @(t) double(any(abs(t) > 1)));
        end

        function diag = wigglyRun(tc, gen, hist, seed)
            %WIGGLYRUN  A non-convex objective through the DEFAULT starts provider.
            %
            %   The default provider is where the seed acts: it is
            %   gen.initialGuesses on its own stream, so this is the path
            %   contract item 2 is a statement about.
            c    = 0.5 * ones(gen.nParams(), 1);
            prob = OptTestProblem( ...
                @(t) sum((t - c).^2) + 0.5 * sum(sin(3 * t)), ...
                @(t) 2 * (t - c) + 1.5 * cos(3 * t), 0);
            [~, diag] = rid.opt.MultiStart("nStarts", 4, "seed", seed).solve( ...
                gen, hist, @(t) prob.withGradient(t));
        end

    end

end
