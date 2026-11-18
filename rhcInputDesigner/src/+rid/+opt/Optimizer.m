classdef Optimizer
%OPTIMIZER  The contract every optimizer over the candidate space must satisfy.
%
%   Choosing the best candidate is a separate concern from generating candidates
%   and from scoring them, so the loop talks to this interface and never to a
%   particular search:
%
%       sol          = rid.opt.MultiStart(...)   % the DEFAULT
%       sol          = rid.opt.GridSearch(...)   % the correctness oracle
%       [best, diag] = sol.solve(gen, hist, scoreFcn, startsFcn)
%
%   THE OPTIMIZER IS HANDED THE GENERATOR, NOT A CANDIDATE SET. For the default
%   parametrization no enumerable set exists -- the decision space is
%   `nHorizon * nInputs` dimensional and a coarse five-level grid over a
%   twenty-sample horizon has 5^20 ~ 1e14 members -- so `solve` obtains
%   `feasibleSet` and `toTrajectory` from `gen`, and its starting points from
%   `startsFcn`.
%
%   `startsFcn` IS A CLOSURE, AND THAT IS WHAT KEEPS BOTH HALVES OF SECTION 4
%   TRUE. `gen.initialGuesses` needs `thetaPrev`, `targetLevels` and
%   `roiInputRange`, which the LOOP owns and `+opt` must not learn about; but
%   "how many starts" is the optimizer's decision, not the loop's. A closure
%   settles both:
%
%       startsFcn = @(nStarts) gen.initialGuesses(hist, thetaPrev, ...
%                       "nStarts", nStarts, "targetLevels", L, "roiInputRange", R);
%
%   -- the optimizer decides the count and asks; the loop supplied its own data
%   by closure; the optimizer never learns what a target level is. This is the
%   device `rid.gen.affineMap` already uses for `znormOf`, for the same reason:
%   a handle carries data across a dependency boundary that a package must not
%   cross. Omitting it means "first step, unguided".
%
%   THE FOUR-POINT CONTRACT, and where each is enforced:
%
%     1. `scoreFcn` is a BLACK BOX. No optimizer may inspect or reimplement the
%        criterion. Every call goes through `countedScore`, so the count is
%        checkable against a scoreFcn that counts itself.
%     2. DETERMINISTIC given `seed`. No implementation reads the global stream;
%        randomness enters only through `startsFcn`, which carries its own.
%     3. RETURNS DIAGNOSTICS -- the fields of `emptyDiagnostics` below.
%     4. CHECKABLE AGAINST EXHAUSTIVE SEARCH on a small discrete problem. This is
%        a correctness requirement, not a benchmark, and it is the only exhaustive
%        check available once the default problem is continuous.
%
%   `scoreFcn` RETURNS THE OBJECTIVE, NEVER `J` ALONE: the criterion plus
%   the soft-constraint penalties. Its three forms, in the order of preference an
%   implementation should ask for them:
%
%       f            = scoreFcn(theta)   % gradient-free; the solver differences it
%       [f, g]       = scoreFcn(theta)   % g (nTheta x 1), the analytic path
%       [f, g, info] = scoreFcn(theta)   % .J, .penalty and .nShell, SEPARATELY
%
%   The third form is asked for at the optimum, because contract item 3 requires
%   the penalty contribution there and the criterion/objective split requires `J` and the penalty to be
%   logged apart; and once per START, because item 3 also requires the
%   dead-on-arrival count and `info.nShell` is where it is READ from. `grad J` is a
%   sum over shell points, so `nShell == 0` IS a zero criterion gradient, exactly
%   and by construction -- no differencing, no tolerance. The HOT PATH still never
%   asks for the third output: once per start is not once per evaluation.
%
%   A scoreFcn that cannot supply it leaves `.penalty`, `.J` and `.deadOnArrival`
%   as NaN with `.penaltyReported` and `.deadOnArrivalReported` false, so the
%   absence is visible rather than silent. Estimating the count by differencing the
%   objective instead would answer a different question: a start dead in the
%   CRITERION may still move under the penalty gradient, toward feasibility rather
%   than toward coverage.
%
%   OPTIMIZERS ARE VALUE OBJECTS. They hold settings and no run state, so the two
%   intended homes of mutable state stay `rid.crit.CoverageState` and
%   the update schedulers.

    methods (Abstract)
        % [best, diag] = solve(obj, gen, hist, scoreFcn, startsFcn)
        %   best  (1x1) struct  .theta (nParams x 1) .uSeq (nHorizon x nInputs)
        %                       .score (1x1)
        %   diag  (1x1) struct  see emptyDiagnostics
        [best, diag] = solve(obj, gen, hist, scoreFcn, startsFcn)
    end

    methods (Static, Access = protected)

        function diag = emptyDiagnostics()
            %EMPTYDIAGNOSTICS  The diagnostic shape of contract item 3, in one place.
            %
            %   Written once so that a new optimizer cannot quietly report a
            %   different set of fields, which would make two runs incomparable
            %   for a reason that has nothing to do with the design.
            %
            %   `deadOnArrival` is NaN in two different situations and
            %   `deadOnArrivalReported` tells them apart: the concept does not
            %   apply (an exhaustive search has no starts), or the objective did
            %   not report `nShell`. NaN and not 0 in both: zero would claim that
            %   no start was dead, which is a different statement from "there were
            %   no starts" and from "the objective could not say".
            diag = struct( ...
                "nEvaluations",    0, ...
                "bestScore",       NaN, ...
                "scores",          zeros(0, 1), ...
                "scoreSpread",     NaN, ...
                "hitBound",        false, ...
                "nAtBound",        0, ...
                "J",               NaN, ...
                "penalty",         NaN, ...
                "penaltyReported", false, ...
                "deadOnArrival",   NaN, ...
                "deadOnArrivalReported", false, ...
                "nStarts",         0, ...
                "gradientPath",    "none", ...
                "exitFlags",       zeros(0, 1), ...
                "iterations",      zeros(0, 1), ...
                "startDisplacement", zeros(0, 1), ...
                "startInfo",       struct(), ...
                "solver",          "");
        end

        function [hitBound, nAtBound] = boundsActive(theta, lb, ub, tol)
            %BOUNDSACTIVE  Did the search end against a bound?
            %
            %   The tolerance is the solver's own documented `ConstraintTolerance`
            %   (rid.defaultOptions), not a number invented here: "active" has to
            %   mean the same thing to this report as it did to the search.
            atBound  = (theta <= lb + tol) | (theta >= ub - tol);
            hitBound = any(atBound);
            nAtBound = sum(atBound);
        end

        function assertFeasible(theta, lb, ub, A, b, tol, what)
            %ASSERTFEASIBLE  A point is inside the declared feasible set, or it fails.
            %
            %   Never repaired, never projected: the method description makes every start feasible BY CONSTRUCTION, so a
            %   violation means the construction is wrong and silently fixing it
            %   would hide that. `rid.gen.Direct.initialGuesses` carries the same
            %   guard; this one also covers a starts-provider that is not it.
            viol = max([0; lb(:) - theta(:); theta(:) - ub(:)]);
            if ~isempty(b)
                viol = max(viol, max([0; A * theta(:) - b(:)]));
            end
            if viol > tol
                error("rid:opt:infeasiblePoint", ...
                    "%s violates the feasible set by %g (tolerance %g). Starts " + ...
                    "are feasible by construction and there is no repair path.", ...
                    what, viol, tol);
            end
        end

        function [J, penalty, reported] = splitAtOptimum(scoreFcn, theta)
            %SPLITATOPTIMUM  Ask scoreFcn once, at the optimum, for J and the penalty.
            %
            %   Contract item 3 requires the penalty contribution at the optimum,
            %   and the criterion/objective split requires `J` and the penalty to be reported separately --
            %   they are different objects and only one of them is comparable
            %   across runs. Neither is recoverable from the objective value alone.
            %
            %   ONE call, after the search: the hot path never asks for the third
            %   output, so an objective that builds an info struct pays for it once
            %   per step rather than once per evaluation.
            J        = NaN;
            penalty  = NaN;
            reported = false;
            [values, nOut] = probeArity(scoreFcn, theta, 3);
            if nOut < 3 || ~isstruct(values{3})
                return
            end
            info = values{3};
            if isfield(info, "J") && isfield(info, "penalty")
                J        = info.J;
                penalty  = info.penalty;
                reported = true;
            end
        end

    end

end
