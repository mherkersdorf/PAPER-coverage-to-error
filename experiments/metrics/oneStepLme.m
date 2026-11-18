function lme = oneStepLme(stepFcn, p, grid)
%ONESTEPLME  The one-step-increment LME used to choose and to report alpha.
%
%   For each row of `grid`, Delta x = stepFcn(x, u, p) - x is the one-step state
%   increment. The best AFFINE model of Delta x as a function of [u, x] -- least
%   squares WITH AN INTERCEPT, fitted once over the whole grid -- is the linear
%   reference; calculate_lme(deltaX, deltaXLinear) (default, no offset option:
%   the intercept already sits in the affine fit) returns the LME PER COMPONENT.
%
%   `grid` is [u, x] in that column order, matching rid.rep.StateSpace's Zphys
%   (the joint-space representation contract) -- one row per supporting point, and this
%   function does not itself decide what the grid is; the caller states its rule
%   (e.g. a 21-per-axis full-factorial grid over the ROI).
%
%   Normalized by Delta x's own variance, this LME tends to the LME of the
%   CONTINUOUS dynamics as tS -> 0; reported here at whatever
%   tS is in `p`. It always refers to `grid`'s region -- a grid over the ROI
%   reports the ROI's LME, and nowhere else.
%
%   Inputs:
%       stepFcn  (1x1) function_handle  xNext = stepFcn(x, u, p),
%                                       (1 x nStates), (1x1), struct -> (1 x nStates)
%       p        (1x1) struct           passed to stepFcn unchanged
%       grid     (nPoints x (1+nStates)) double  columns [u, x1, ..., xn]
%   Outputs:
%       lme  (1 x nStates) double  LME [%] per state component, grid column order
arguments
    stepFcn (1,1) function_handle
    p       (1,1) struct
    grid    (:,:) double {mustBeReal, mustBeFinite}
end

if size(grid, 2) < 2
    error("experiments:oneStepLme:dimensionMismatch", ...
        "grid must have at least 2 columns [u, x1, ...]; got %s.", ...
        mat2str(size(grid)));
end

nPoints = size(grid, 1);
nStates = size(grid, 2) - 1;

u = grid(:, 1);
x = grid(:, 2:end);

deltaX = zeros(nPoints, nStates);
for iPoint = 1:nPoints
    xNext = stepFcn(x(iPoint, :), u(iPoint), p);
    deltaX(iPoint, :) = xNext - x(iPoint, :);
end

Phi = [grid, ones(nPoints, 1)];
theta = Phi \ deltaX;
deltaXLinear = Phi * theta;

lme = calculate_lme(deltaX, deltaXLinear);
end
