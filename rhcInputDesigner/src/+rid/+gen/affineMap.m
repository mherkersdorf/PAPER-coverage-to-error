function [A, c, info] = affineMap(znormOf, theta0, opts)
%AFFINEMAP  The map theta -> Znorm, MEASURED by probe, never read off the model.
%
%       c   = Znorm(theta0)
%       A_i = [ Znorm(theta0 + h*e_i) - Znorm(theta0) ] / h      i = 1 .. nTheta
%
%   so that  Znorm(theta)(:) = A*theta + c(:), and
%   dJ/dtheta = (dJ/dZnorm)(:).' * A -- assembled by +opt, never here.
%
%   FOR AN AFFINE MAP THE PROBE IS EXACT FOR ANY h ~= 0. It evaluates the linear
%   part rather than approximating it, so there is no truncation error and no step
%   size to tune; `h` is chosen in the natural scale of the input for
%   floating-point conditioning only. This is NOT finite differencing and the
%   section 8 rules about step size and tolerance do not apply to it.
%
%   MEASURED, NOT EXTRACTED. Reading the Markov parameters
%   `g_m = C_d A_d^(m-1) B_d` would require every surrogate to expose its
%   state-space matrices, which contradicts the minimal interface of the minimal surrogate interface -- a
%   vendor FMU can supply none of them. Probing through the caller's own closure
%   keeps the surrogate and the representation black boxes, so lags, differences and
%   the normalization are captured without `+gen` knowing anything about them.
%
%   `znormOf` IS A CLOSURE, NOT A MODEL. It is supplied by the loop, which owns the
%   surrogate, the representation and the normalization. That is what keeps `+gen`
%   free of dependencies on `+rep` and `+roi`, and it is what lets this function
%   probe any composition at all.
%
%   A CARRIES THE HORIZON ROWS ONLY. Rows of Znorm built from the realized history
%   do not move with theta (section 2), so `znormOf` must return the horizon block
%   alone. Carrying history rows as zeros is waste; forgetting they are constant is
%   worse, and is the error section 3 warns cannot be caught by a finite-difference
%   reference.
%
%   AFFINENESS IS FALSIFIED, NEVER PROVED. Finitely many checks cannot establish
%   affineness -- a nonlinear map may coincide with an affine one at the points
%   tried -- so `info.tier` reports what the checks found and the caller must take
%   the affine path ONLY when it reads "affine". The check points must be GENERIC:
%   the probe's own points are axis-aligned, so an axis-aligned validation point
%   tests nothing new. Multi-start rows are free generic points already inside the
%   feasible set.
%
%   Inputs:
%       znormOf            function_handle  theta (nTheta x 1) -> (nRows x nDim)
%                                           normalized HORIZON coordinates
%       theta0             (nTheta x 1)     the probe centre
%       opts.h             (1x1)            probe step; conditioning only.
%                                           `opts.gen.affineProbeStep`
%       opts.checkPoints   (nChecks x nTheta) GENERIC validation points
%       opts.checkTol      (1x1)            relative tolerance for the check, and
%                                           therefore the affine/notAffine VERDICT.
%                                           `opts.gen.affineCheckTol`
%   Outputs:
%       A     (nRows*nDim x nTheta)  Znorm(:) = A*theta + c(:)
%       c     (nRows x nDim)         the constant part, Znorm(theta0) - A*theta0
%       info  (1x1) struct  .tier ("affine" | "notAffine" | "unchecked")
%                           .maxRelError .nChecks .h .nProbes .nRows .nDim
%
%   BOTH SETTINGS COME FROM `rid.defaultOptions` and neither is a literal here, so
%   a record's effective options say what the probe actually did. They are not the
%   same kind of number: `h` is conditioning and cannot change an answer, while
%   `checkTol` decides the verdict and therefore which gradient path the run took.

arguments
    znormOf          (1,1) function_handle
    theta0           (:,1) double {mustBeReal, mustBeFinite}
    opts.h           (1,1) double {mustBeReal, mustBeNonzero} ...
        = rid.defaultOptions().gen.affineProbeStep
    opts.checkPoints (:,:) double {mustBeReal} = zeros(0, 0)
    opts.checkTol    (1,1) double {mustBeReal, mustBePositive} ...
        = rid.defaultOptions().gen.affineCheckTol
end

nTheta = numel(theta0);
Z0     = znormOf(theta0);
[nRows, nDim] = size(Z0);

A = zeros(nRows * nDim, nTheta);
for iTheta = 1:nTheta
    thetaProbe          = theta0;
    thetaProbe(iTheta)  = thetaProbe(iTheta) + opts.h;
    Zi                  = znormOf(thetaProbe);
    if ~isequal(size(Zi), [nRows, nDim])
        error("rid:gen:inconsistentMap", ...
            "znormOf returned %s at probe %d but %s at theta0.", ...
            mat2str(size(Zi)), iTheta, mat2str([nRows, nDim]));
    end
    A(:, iTheta) = (Zi(:) - Z0(:)) / opts.h;
end

% c is the constant part of the affine map, not Znorm(theta0).
c = reshape(Z0(:) - A * theta0, nRows, nDim);

info = struct("tier", "unchecked", "maxRelError", NaN, ...
              "nChecks", size(opts.checkPoints, 1), "h", opts.h, ...
              "nProbes", nTheta, "nRows", nRows, "nDim", nDim);

if isempty(opts.checkPoints)
    return
end
if size(opts.checkPoints, 2) ~= nTheta
    error("rid:gen:invalidCheckPoints", ...
        "checkPoints has %d column(s); theta has %d entries.", ...
        size(opts.checkPoints, 2), nTheta);
end

maxRel = 0;
for iCheck = 1:size(opts.checkPoints, 1)
    thetaCheck = opts.checkPoints(iCheck, :).';
    actual     = znormOf(thetaCheck);
    predicted  = reshape(A * thetaCheck, nRows, nDim) + c;
    denom      = max(1, max(abs(actual(:))));
    maxRel     = max(maxRel, max(abs(actual(:) - predicted(:))) / denom);
end
info.maxRelError = maxRel;
if maxRel <= opts.checkTol
    info.tier = "affine";
else
    info.tier = "notAffine";
end

end
