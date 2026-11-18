classdef MultiStart < rid.opt.Optimizer
%MULTISTART  The default optimizer: one local search per starting point.
%
%   `Znorm` is affine in `theta` for the default configuration, but `J` is a
%   MINIMUM over samples of a Euclidean norm, and a minimum of convex functions is
%   not convex. The non-convexity is real, so
%   multi-start is not optional: one local search is run from each row of
%   the start set, and the best result wins.
%
%       sol          = rid.opt.MultiStart()
%       [best, diag] = sol.solve(gen, hist, scoreFcn, startsFcn)
%
%   THE STARTS MATTER MORE HERE THAN IN ORDINARY MULTI-START. Beyond escaping
%   local optima, large parts of this landscape have EXACTLY ZERO gradient: a
%   supporting point contributes nothing when it is resolved (`dTilde <= 1`) or
%   saturated (`dTilde >= 1 + psiMax`). A start whose trajectory leaves every
%   supporting point outside the shell `1 < dTilde < 1 + psiMax` cannot move at
%   all, and the search reports convergence at iteration zero having done nothing.
%
%   That is why `diag.deadOnArrival` is a REQUIRED diagnostic rather than a
%   curiosity. A high count is the run saying that `psiMax` is too small for the
%   optimizer -- the tension in the psiMax = 1 default made observable -- and without the count
%   it presents as a design that merely searches badly.
%
%   DEAD ON ARRIVAL IS READ, NOT ESTIMATED. `grad J` is a sum over shell points, so
%   `nShell == 0` IS a zero criterion gradient -- exactly, and by construction, with
%   no tolerance and no difference quotient anywhere near it. The count comes from
%   the objective's `info.nShell`, which is free on the analytic path
%   because the evaluation at `theta0` happens anyway, and costs one evaluation per
%   start otherwise -- against `nTheta` per start for a difference estimate.
%
%   IT IS DEAD WITH RESPECT TO THE CRITERION, which is what
%   the method description means and what `psiMax` governs. A start
%   with `nShell == 0` may still move under the PENALTY gradient, toward
%   feasibility rather than toward coverage; conflating the two would blunt the
%   diagnostic into "did the solver move", which is a different and much less
%   useful question.
%
%   An objective that does not report `nShell` leaves the count NaN, with
%   `diag.deadOnArrivalReported` false. Estimating it from the objective gradient
%   instead would silently answer the wrong question.
%
%   THE GRADIENT PATH IS DETECTED, NOT ASSUMED, AND ALWAYS REPORTED. The objective
%   is asked for three outputs once, at the first start; what it supplies settles
%   both questions at once -- two or more outputs means the analytic path, three
%   means `info.nShell` is available as well. Otherwise the solver differences the
%   objective itself. `diag.gradientPath` records which happened,
%   so a fall-back to finite differences -- the tier-3 case of
%   the method description, where the affine check did not pass --
%   cannot masquerade as the fast path.
%
%   Detection is by CALLING, never by `nargout(scoreFcn)`: an anonymous handle
%   reports -1 there and passes its own nargout through to whatever it wraps, so
%   the declared arity classifies nothing. See `probeArity`.
%
%   NOTHING HERE IS RANDOM. Contract item 2 -- deterministic given `seed` -- holds
%   because every search is deterministic and the only randomness in the step
%   enters through `startsFcn`, which carries its own stream. The seed
%   held here is used only by the default starts provider.
%
%   Construction:
%       sol = rid.opt.MultiStart()
%       sol = rid.opt.MultiStart("nStarts", 8, "seed", 3)
%
%   Inputs:
%       opts.nStarts                       (1x1)  local searches per step
%       opts.seed                          (1x1)  for the DEFAULT starts provider
%                                                 (the +gen stream)
%       opts.algorithm                     (1x1) string
%       opts.maxIterations                 (1x1)
%       opts.maxEvaluationsPerVariable     (1x1)
%       opts.optimalityTolerance           (1x1)
%       opts.stepTolerance                 (1x1)
%       opts.constraintTolerance           (1x1)

    properties (SetAccess = immutable)
        nStarts                    (1,1) double
        seed                       (1,1) double
        algorithm                  (1,1) string
        maxIterations              (1,1) double
        maxEvaluationsPerVariable  (1,1) double
        optimalityTolerance        (1,1) double
        stepTolerance              (1,1) double
        constraintTolerance        (1,1) double
    end

    methods

        function obj = MultiStart(opts)
            arguments
                % The count is a SEARCH-STRATEGY setting, so it is declared where
                % the decision lives. rid.gen.Direct.initialGuesses takes
                % it as an argument with no default of its own: +opt decides and
                % +gen answers, so there is one declaration and one decider.
                opts.nStarts (1,1) double {mustBeInteger, mustBePositive} ...
                    = rid.defaultOptions().opt.nStarts
                opts.seed    (1,1) double {mustBeInteger, mustBeNonnegative} ...
                    = rid.defaultOptions().gen.seed
                opts.algorithm (1,1) string = rid.defaultOptions().opt.algorithm
                opts.maxIterations (1,1) double {mustBePositive} ...
                    = rid.defaultOptions().opt.maxIterations
                opts.maxEvaluationsPerVariable (1,1) double {mustBePositive} ...
                    = rid.defaultOptions().opt.maxEvaluationsPerVariable
                opts.optimalityTolerance (1,1) double {mustBePositive} ...
                    = rid.defaultOptions().opt.optimalityTolerance
                opts.stepTolerance (1,1) double {mustBePositive} ...
                    = rid.defaultOptions().opt.stepTolerance
                opts.constraintTolerance (1,1) double {mustBePositive} ...
                    = rid.defaultOptions().opt.constraintTolerance
            end
            obj.nStarts                   = opts.nStarts;
            obj.seed                      = opts.seed;
            obj.algorithm                 = opts.algorithm;
            obj.maxIterations             = opts.maxIterations;
            obj.maxEvaluationsPerVariable = opts.maxEvaluationsPerVariable;
            obj.optimalityTolerance       = opts.optimalityTolerance;
            obj.stepTolerance             = opts.stepTolerance;
            obj.constraintTolerance       = opts.constraintTolerance;
        end

        function [best, diag] = solve(obj, gen, hist, scoreFcn, startsFcn)
            %SOLVE  One local search per start; the best result wins.
            %
            %   Inputs:
            %       gen        (1x1)            the generator, e.g. rid.gen.Direct
            %       hist       (1x1) struct     the signal history
            %       scoreFcn   function_handle  the OBJECTIVE, J + penalties
            %       startsFcn  function_handle  nStarts -> [Theta0, info]; omitted
            %                                   or empty means first step, unguided
            %   Outputs:
            %       best  (1x1) struct  .theta .uSeq .score
            %       diag  (1x1) struct  rid.opt.Optimizer.emptyDiagnostics
            arguments
                obj       (1,1) rid.opt.MultiStart
                gen       (1,1)
                hist      (1,1) struct
                scoreFcn  (1,1) function_handle
                startsFcn                       = []
            end

            [lb, ub, A, b] = gen.feasibleSet(hist);
            nTheta         = gen.nParams();

            % The default provider is the first-step, unguided case: no previous
            % solution to shift and no target levels, which is exactly what
            % initialGuesses does with an empty thetaPrev.
            if isempty(startsFcn)
                startsFcn = @(n) gen.initialGuesses(hist, zeros(0, 1), ...
                                     "nStarts", n, "seed", obj.seed);
            elseif ~isa(startsFcn, "function_handle")
                error("rid:opt:invalidStartsProvider", ...
                    "startsFcn must be a function handle nStarts -> " + ...
                    "[Theta0, info], or empty; got %s.", class(startsFcn));
            end

            [Theta0, startInfo] = startsFcn(obj.nStarts);
            if ~isequal(size(Theta0), [obj.nStarts, nTheta])
                error("rid:opt:invalidStarts", ...
                    "startsFcn returned %s; %d start(s) of %d parameter(s) " + ...
                    "were requested.", mat2str(size(Theta0)), obj.nStarts, nTheta);
            end

            [counted, readCount] = countedScore(scoreFcn);

            thetaOpt    = zeros(nTheta, obj.nStarts);
            scores      = inf(obj.nStarts, 1);
            exitFlags   = zeros(obj.nStarts, 1);
            iterations  = zeros(obj.nStarts, 1);
            % max |theta - theta0| per start: zero is a start the solver
            % never moved, which a flat objective produces.
            displaced   = zeros(obj.nStarts, 1);
            dead        = false(obj.nStarts, 1);
            nOut        = 0;

            for iStart = 1:obj.nStarts
                theta0 = Theta0(iStart, :).';
                rid.opt.Optimizer.assertFeasible(theta0, lb, ub, A, b, ...
                    obj.constraintTolerance, sprintf("Start %d", iStart));

                % One call at theta0 settles both questions: how many outputs the
                % objective has (hence the gradient path) and what nShell is at the
                % start (hence dead-on-arrival). The arity is probed once, on the
                % first start, and reused -- an objective does not change shape
                % between starts.
                if iStart == 1
                    [values, nOut] = probeArity(counted, theta0, 3);
                else
                    values = cell(1, nOut);
                    [values{:}] = counted(theta0);
                end

                if nOut >= 2 && numel(values{2}) ~= nTheta
                    error("rid:opt:invalidGradient", ...
                        "scoreFcn returned a gradient with %d entries; " + ...
                        "theta has %d.", numel(values{2}), nTheta);
                end
                if nOut >= 3
                    info = values{3};
                    if ~isstruct(info) || ~isfield(info, "nShell")
                        error("rid:opt:invalidScoreInfo", ...
                            "The third output of scoreFcn must be a struct " + ...
                            "carrying nShell.");
                    end
                    % grad J is a sum over shell points, so nShell == 0 IS a zero
                    % criterion gradient. Exact, by construction.
                    dead(iStart) = info.nShell == 0;
                end

                [thetaOpt(:, iStart), scores(iStart), solveInfo] = ...
                    rid.ext.constrainedMinimize(counted, theta0, lb, ub, A, b, ...
                        "useGradient",                nOut >= 2, ...
                        "algorithm",                  obj.algorithm, ...
                        "maxIterations",              obj.maxIterations, ...
                        "maxEvaluationsPerVariable",  obj.maxEvaluationsPerVariable, ...
                        "optimalityTolerance",        obj.optimalityTolerance, ...
                        "stepTolerance",              obj.stepTolerance);
                exitFlags(iStart)  = solveInfo.exitFlag;
                iterations(iStart) = solveInfo.iterations;
                displaced(iStart)  = max(abs(thetaOpt(:, iStart) - theta0));
            end

            % min() takes the lowest index among ties, which keeps the choice
            % deterministic when two starts converge to the same value.
            [bestScore, iBest] = min(scores);
            thetaBest = thetaOpt(:, iBest);

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
            % NaN, not 0, when the objective does not report nShell: zero would
            % claim that no start was dead, which is a different statement from
            % "the objective could not say".
            diag.deadOnArrival         = NaN;
            diag.deadOnArrivalReported = nOut >= 3;
            if nOut >= 3
                diag.deadOnArrival = sum(dead);
            end
            diag.nStarts      = obj.nStarts;
            diag.gradientPath = "finiteDifference";
            if nOut >= 2
                diag.gradientPath = "analytic";
            end
            diag.exitFlags = exitFlags;
            diag.iterations        = iterations;
            diag.startDisplacement = displaced;
            diag.startInfo = startInfo;
            diag.solver    = "rid.opt.MultiStart";
        end

    end

end
