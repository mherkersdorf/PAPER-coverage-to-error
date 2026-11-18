classdef GridSearch < rid.opt.Optimizer
%GRIDSEARCH  Exhaustive search over a level grid: the correctness ORACLE.
%
%   NOT THE DEFAULT SOLVER, and no longer a candidate for it. The default
%   parametrization makes every future input sample a continuous decision variable,
%   so there is no set to enumerate: a coarse five-level grid over a
%   twenty-sample horizon has 5^20 ~ 1e14 members. This class survives for two
%   reasons:
%
%     1. it is the exhaustive check of contract item 4 -- on a small DISCRETE
%        problem any optimizer must find the grid-search optimum. That is a
%        correctness requirement, not a benchmark, and it is now the only
%        exhaustive check available anywhere in the design;
%     2. it is the solver for the discrete parametrizations.
%
%   BECAUSE THE DEFAULT PROBLEM IS CONTINUOUS, THE ORACLE PROBLEM MUST BE MADE
%   TINY ON PURPOSE -- a two-sample horizon on a coarse level grid -- where
%   exhaustive search is actually possible. `maxGridPoints` refuses anything
%   larger rather than appearing to work: an exhaustive search that quietly took
%   an hour would be indistinguishable from a hung run.
%
%   IT BUILDS ITS OWN CANDIDATE SET. `startsFcn` is accepted and ignored, and the
%   argument exists only so that the loop can hand any optimizer the same four
%   things. Starting points are meaningless to an exhaustive search, which is also
%   why `diag.deadOnArrival` is NaN here: there are no starts to be dead, and
%   reporting 0 would claim that none were.
%
%   THE GRID IS BUILT FROM THE DECLARED FEASIBLE SET, per variable, and then
%   FILTERED by the linear inequalities. Bounds define the grid because they are a
%   box; the rate constraints cannot, because they couple variables -- so they are
%   applied as a filter and the infeasible combinations are never scored, which is
%   loop invariant 4 (no infeasible candidate is ever scored).
%
%   A variable whose bounds coincide contributes ONE level, not `nLevels`
%   identical ones. `feasibleSet` produces that case whenever the rate limit pins
%   the first sample exactly, and scoring the same point five times would inflate
%   `nEvaluations` without changing anything.
%
%   Construction:
%       sol = rid.opt.GridSearch()
%       sol = rid.opt.GridSearch("nLevels", 3)
%
%   Inputs:
%       opts.nLevels             (1x1)  grid levels per decision variable
%       opts.maxGridPoints       (1x1)  refuse a grid larger than this
%       opts.constraintTolerance (1x1)

    properties (SetAccess = immutable)
        nLevels             (1,1) double
        maxGridPoints       (1,1) double
        constraintTolerance (1,1) double
    end

    methods

        function obj = GridSearch(opts)
            arguments
                opts.nLevels (1,1) double {mustBeInteger, mustBePositive} ...
                    = rid.defaultOptions().opt.gridLevels
                opts.maxGridPoints (1,1) double {mustBePositive} ...
                    = rid.defaultOptions().opt.maxGridPoints
                opts.constraintTolerance (1,1) double {mustBePositive} ...
                    = rid.defaultOptions().opt.constraintTolerance
            end
            obj.nLevels             = opts.nLevels;
            obj.maxGridPoints       = opts.maxGridPoints;
            obj.constraintTolerance = opts.constraintTolerance;
        end

        function [best, diag] = solve(obj, gen, hist, scoreFcn, startsFcn)
            %SOLVE  Score every feasible grid point and return the best.
            %
            %   Inputs:
            %       gen        (1x1)            the generator
            %       hist       (1x1) struct     the signal history
            %       scoreFcn   function_handle  the OBJECTIVE, J + penalties
            %       startsFcn                   accepted and IGNORED; see the header
            %   Outputs:
            %       best  (1x1) struct  .theta .uSeq .score
            %       diag  (1x1) struct  rid.opt.Optimizer.emptyDiagnostics
            arguments
                obj       (1,1) rid.opt.GridSearch
                gen       (1,1)
                hist      (1,1) struct
                scoreFcn  (1,1) function_handle
                startsFcn                       = []      %#ok<INUSA> exhaustive
            end

            [lb, ub, A, b] = gen.feasibleSet(hist);
            nTheta         = gen.nParams();

            axes  = cell(1, nTheta);
            total = 1;
            for iTheta = 1:nTheta
                if lb(iTheta) > ub(iTheta)
                    error("rid:opt:emptyFeasibleSet", ...
                        "Variable %d has lb = %g above ub = %g.", ...
                        iTheta, lb(iTheta), ub(iTheta));
                end
                if lb(iTheta) == ub(iTheta)
                    axes{iTheta} = lb(iTheta);
                else
                    axes{iTheta} = linspace(lb(iTheta), ub(iTheta), obj.nLevels);
                end
                total = total * numel(axes{iTheta});
                if total > obj.maxGridPoints
                    error("rid:opt:gridTooLarge", ...
                        "The grid has more than %g points (%d variable(s) at " + ...
                        "up to %d level(s) each). Exhaustive search is the " + ...
                        "correctness oracle of a DELIBERATELY TINY problem, " + ...
                        "not a solver for the default parametrization.", ...
                        obj.maxGridPoints, nTheta, obj.nLevels);
                end
            end

            mesh = cell(1, nTheta);
            [mesh{:}] = ndgrid(axes{:});
            Theta = zeros(total, nTheta);
            for iTheta = 1:nTheta
                Theta(:, iTheta) = mesh{iTheta}(:);
            end

            % Loop invariant 4: an infeasible candidate is never scored. The
            % bounds already hold by construction of the axes; the coupled rate
            % constraints are what this removes.
            keep = true(total, 1);
            if ~isempty(b)
                keep = all(A * Theta.' <= b(:) + obj.constraintTolerance, 1).';
            end
            Theta = Theta(keep, :);
            if isempty(Theta)
                error("rid:opt:noFeasibleGridPoint", ...
                    "None of the %d grid point(s) satisfies the linear " + ...
                    "constraints. The grid is too coarse to contain a feasible " + ...
                    "combination, not the feasible set empty.", total);
            end

            [counted, readCount] = countedScore(scoreFcn);

            nKept  = size(Theta, 1);
            scores = zeros(nKept, 1);
            for iCand = 1:nKept
                scores(iCand) = counted(Theta(iCand, :).');
            end

            [bestScore, iBest] = min(scores);
            thetaBest = Theta(iBest, :).';

            [J, penalty, reported] = rid.opt.Optimizer.splitAtOptimum(counted, thetaBest);
            [hitBound, nAtBound]   = rid.opt.Optimizer.boundsActive( ...
                thetaBest, lb, ub, obj.constraintTolerance);

            best = struct("theta", thetaBest, ...
                          "uSeq",  gen.toTrajectory(thetaBest), ...
                          "score", bestScore);

            diag = rid.opt.Optimizer.emptyDiagnostics();
            diag.nEvaluations    = readCount();
            diag.bestScore       = bestScore;
            diag.scores          = scores;
            diag.scoreSpread     = max(scores) - min(scores);
            diag.hitBound        = hitBound;
            diag.nAtBound        = nAtBound;
            diag.J               = J;
            diag.penalty         = penalty;
            diag.penaltyReported = reported;
            diag.deadOnArrival         = NaN;    % no starts; see the header
            diag.deadOnArrivalReported = false;
            diag.nStarts         = 0;
            diag.gradientPath    = "none";
            diag.exitFlags       = zeros(0, 1);
            diag.startInfo       = struct("nGridPoints", total, "nFeasible", nKept, ...
                                          "nLevels", obj.nLevels);
            diag.solver          = "rid.opt.GridSearch";
        end

    end

end
