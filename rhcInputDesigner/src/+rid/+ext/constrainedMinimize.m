function [theta, fVal, info] = constrainedMinimize(objFcn, theta0, lb, ub, A, b, opts)
%CONSTRAINEDMINIMIZE  Local minimization under bounds and linear inequalities.
%
%   TOOLBOX: Optimization Toolbox (`fmincon`, `optimoptions`). This is the ONLY
%   place in the project that calls it, so the project's dependence on a local
%   solver can be enumerated by listing +ext and replaced by editing one file.
%   +opt never names fmincon.
%
%       min_theta  objFcn(theta)   s.t.   lb <= theta <= ub,  A*theta <= b
%
%   The functionality isolated here is "a local search over a polytope", not "a
%   call to fmincon": a replacement needs to honour this signature and nothing
%   else. In particular the gradient enters through `objFcn` returning a second
%   output, which is how every reasonable local solver takes it.
%
%   THE ALGORITHM IS SQP, NOT THE SOLVER'S OWN DEFAULT. The soft-constraint
%   penalty is deliberately C1 and NOT C2 -- the second derivative of
%   `max(0, v)^2` jumps from 0 to 2 -- and
%   the method description records that SQP tolerates exactly
%   that. Every other setting below is the solver's documented default, taken
%   rather than invented (see rid.defaultOptions).
%
%   `MaxFunctionEvaluations` is expressed as a count PER VARIABLE because that is
%   how the solver documents its own default (`100*numberOfVariables`). The
%   decision vector is `nHorizon * nInputs` long and the horizon is a run-level
%   choice, so a fixed absolute cap would silently tighten as the horizon grew.
%
%   NOTHING IS PRINTED. Display is forced off: +ext is on the no-output list,
%   and a solver banner inside the candidate loop would be printed
%   nStarts times per step.
%
%   Inputs:
%       objFcn      function_handle  theta (nTheta x 1) -> f, or -> [f, g] when
%                                    opts.useGradient is true
%       theta0      (nTheta x 1)     starting point; must already be feasible
%       lb, ub      (nTheta x 1)     bounds
%       A           (nRows x nTheta) linear inequality matrix; may be empty
%       b           (nRows x 1)      right-hand side; may be empty
%       opts.useGradient              (1x1) logical
%       opts.algorithm                (1x1) string
%       opts.maxIterations            (1x1) double
%       opts.maxEvaluationsPerVariable(1x1) double
%       opts.optimalityTolerance      (1x1) double
%       opts.stepTolerance            (1x1) double
%   Outputs:
%       theta   (nTheta x 1)  the local minimizer
%       fVal    (1x1)         objFcn(theta)
%       info    (1x1) struct  .exitFlag .iterations .funcCount .algorithm
%                             .message .firstOrderOpt

arguments
    objFcn                            (1,1) function_handle
    theta0                            (:,1) double {mustBeReal, mustBeFinite}
    lb                                (:,1) double {mustBeReal}
    ub                                (:,1) double {mustBeReal}
    A                                 (:,:) double {mustBeReal}
    b                                 (:,1) double {mustBeReal}
    opts.useGradient                  (1,1) logical = false
    opts.algorithm                    (1,1) string = rid.defaultOptions().opt.algorithm
    opts.maxIterations                (1,1) double {mustBePositive} ...
        = rid.defaultOptions().opt.maxIterations
    opts.maxEvaluationsPerVariable    (1,1) double {mustBePositive} ...
        = rid.defaultOptions().opt.maxEvaluationsPerVariable
    opts.optimalityTolerance          (1,1) double {mustBePositive} ...
        = rid.defaultOptions().opt.optimalityTolerance
    opts.stepTolerance                (1,1) double {mustBePositive} ...
        = rid.defaultOptions().opt.stepTolerance
end

nTheta = numel(theta0);
if numel(lb) ~= nTheta || numel(ub) ~= nTheta
    error("rid:ext:dimensionMismatch", ...
        "lb has %d and ub has %d entries; theta0 has %d.", ...
        numel(lb), numel(ub), nTheta);
end
if ~isempty(A) && size(A, 2) ~= nTheta
    error("rid:ext:dimensionMismatch", ...
        "A has %d column(s); theta0 has %d entries.", size(A, 2), nTheta);
end
if size(A, 1) ~= numel(b)
    error("rid:ext:dimensionMismatch", ...
        "A has %d row(s) and b has %d entries.", size(A, 1), numel(b));
end

solverOpts = optimoptions("fmincon", ...
    "Algorithm",                 opts.algorithm, ...
    "SpecifyObjectiveGradient",  opts.useGradient, ...
    "MaxIterations",             opts.maxIterations, ...
    "MaxFunctionEvaluations",    opts.maxEvaluationsPerVariable * nTheta, ...
    "OptimalityTolerance",       opts.optimalityTolerance, ...
    "StepTolerance",             opts.stepTolerance, ...
    "Display",                   "off");

% The generator builds its rate band with kron and spdiags, so A arrives
% sparse; the sqp algorithm densifies it and warns while doing so, and +ext is
% on the no-output list. Converting here says what the solver is
% going to do anyway, in the one place that knows it is fmincon.
Ause = full(A);
buse = b;
if isempty(A)
    Ause = [];
    buse = [];
end

[theta, fVal, exitFlag, output] = fmincon(objFcn, theta0, Ause, buse, ...
    [], [], lb, ub, [], solverOpts);

theta = theta(:);
info  = struct( ...
    "exitFlag",      exitFlag, ...
    "iterations",    output.iterations, ...
    "funcCount",     output.funcCount, ...
    "algorithm",     string(output.algorithm), ...
    "message",       string(output.message), ...
    "firstOrderOpt", output.firstorderopt);

end
