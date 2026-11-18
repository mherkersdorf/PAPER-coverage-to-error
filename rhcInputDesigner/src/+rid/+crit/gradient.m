function [gradZnorm, info] = gradient(Znorm, Q, epsReq, opts)
%GRADIENT  dJ/dZnorm, the criterion's derivative with respect to its own input.
%
%   Returns one row per SAMPLE: `gradZnorm(i, :)` is dJ/d(Znorm(i, :)). It knows
%   nothing about the parametrization, the surrogate or the representation -- it
%   needs only the shell membership and the nearest-neighbour data, which is why it
%   lives in `+crit` and not in `+gen`. The chain to the decision vector,
%
%       dJ/dtheta = (dJ/dZnorm)(:).' * A
%
%   is assembled by `+opt`; nothing here forms an objective.
%
%   For supporting point q_j with nearest sample i*(j):
%
%       uvec_j          = (z_{i*} - q_j) / dRaw_j
%       dJ/dz_{i*}     += dJ/dpsi_j * (1/epsReq_j) * uvec_j     for j in the shell
%       dJ/dpsi_j       = (1/nSup) * (psi_j / J)^(pNorm - 1)
%
%   THE FORM OF dJ/dpsi MATTERS. The algebraically identical
%   `psi_j^(p-1) * S^(1/p-1) / nSup` forms `0 * Inf` at `J = 0` and is poorly
%   conditioned near it. The form used here also degenerates correctly at
%   `pNorm = 1`, returning `(psi/J)^0 = 1`, so one expression covers `p >= 1`.
%
%   pNorm = Inf IS REFUSED. The worst-case limit makes J a max over supporting
%   points, whose derivative is that of the argmax term alone. The p-norm form
%   below would evaluate there to 1/nSup at the max point, where 1 is correct --
%   wrong by a factor of nSup, finite, and entirely plausible-looking. Failing
%   loudly is the project's rule for a case that is documented but not yet built.
%
%   J == 0 => GRADIENT EXACTLY ZERO, on exact equality and with NO tolerance.
%   Exact zeros genuinely arise: `max(dTilde - 1, 0)` returns an exact 0, so J is
%   exactly 0 when every point is resolved -- which is the termination condition of
%   certification mode, reached by design. A `J < tol` guard would silently zero
%   real gradients wherever J is small but the design is still moving, which is the
%   endgame of a certification run and exactly where a false zero does most damage.
%   Zero is usually the TRUE gradient there, not a substitute: J is identically
%   zero on a neighbourhood unless some dTilde sits exactly at 1.
%
%   The three other nonsmooth cases of section 5, all returning a finite element of
%   the Clarke generalized gradient and never NaN or Inf:
%
%       argmin tie          lowest index, deterministically. `min` already does
%                           this, and the value is a GENUINE element of the Clarke
%                           gradient -- section 5 derives the hull equality -- not
%                           a pragmatic approximation.
%       dRaw_j = 0          contributes zero: uvec is 0/0, and psi_j = 0 there
%                           anyway, so the term is doubly inert.
%       dTilde at a clip    one-sided, taken consistently as the OUTER value: a
%       boundary            point exactly at 1 or at 1 + psiMax is treated as
%                           outside the shell and contributes zero.
%
%   A SUPPORTING POINT WHOSE NEAREST SAMPLE IS NOT A DECISION VARIABLE contributes
%   nothing to dJ/dtheta. That restriction is not applied here -- this function
%   reports the derivative for every sample row it was given -- it falls out of the
%   composition, because `A` carries the horizon rows only (section 2). Callers must
%   therefore compose with an `A` built the same way; `info.nearest` is returned so
%   a test can assert which samples are active.
%
%   Inputs:
%       Znorm         (nSamples x nDim)  normalized coordinates, ALL samples
%       Q             (nSup x nDim)      supporting points, normalized
%       epsReq        (nSup x 1)         demanded resolution per supporting point
%       opts.psiMax   (1x1)
%       opts.pNorm    (1x1)
%   Outputs:
%       gradZnorm     (nSamples x nDim)  dJ/dZnorm
%       info          (1x1) struct  .J .nearest (nSup x 1) .isShell (nSup x 1)
%                                   .nShell .dRaw .dTilde .psi

arguments
    Znorm       (:,:) double {mustBeReal, mustBeFinite}
    Q           (:,:) double {mustBeReal, mustBeFinite}
    epsReq      (:,1) double {mustBeReal}
    opts.psiMax (1,1) double {mustBeReal} = rid.defaultOptions().crit.psiMax
    opts.pNorm  (1,1) double {mustBeReal} = rid.defaultOptions().crit.pNorm
end

nSup     = size(Q, 1);
nSamples = size(Znorm, 1);
nDim     = size(Q, 2);

if nSup < 1
    error("rid:crit:invalidResolution", ...
        "Q must have at least one supporting point.");
end
if numel(epsReq) ~= nSup
    error("rid:crit:invalidResolution", ...
        "epsReq has %d entries; Q has %d supporting points.", ...
        numel(epsReq), nSup);
end
if nSamples > 0 && size(Znorm, 2) ~= nDim
    error("rid:crit:dimensionMismatch", ...
        "Znorm has %d column(s); Q has %d.", size(Znorm, 2), nDim);
end

% REFUSED, NOT APPROXIMATED. At the worst-case limit J is a MAX over supporting
% points, so its derivative is that of the argmax term alone -- a different
% expression, with the Clarke hull and the lowest-index tie-break of section 5
% applying to ties. Evaluating the p-norm form below at pNorm = Inf silently
% returns 1/nSup at the max point where 1 is correct, and 0 elsewhere: a gradient
% wrong by a factor of nSup that looks entirely plausible. The method description documents
% p -> Inf as the worst-case mode, so this has to be implemented eventually; until
% it is, it fails loudly.
if isinf(opts.pNorm)
    error("rid:crit:gradientNormUnsupported", ...
        "The gradient is not implemented for pNorm = Inf. J is a max there, " + ...
        "so the derivative is that of the argmax term alone; the p-norm form " + ...
        "used here would return 1/nSup at the max point where 1 is correct.");
end

gradZnorm = zeros(nSamples, nDim);

% The criterion itself comes from the permanent oracle, so the gradient cannot
% drift from the value it is the gradient OF.
[J, evalInfo] = rid.crit.evaluate(Znorm, Q, epsReq, ...
    "psiMax", opts.psiMax, "pNorm", opts.pNorm);

info = struct("J", J, "nearest", zeros(nSup, 1), ...
              "isShell", false(nSup, 1), "nShell", 0, ...
              "dRaw", evalInfo.dRaw, "dTilde", evalInfo.dTilde, "psi", evalInfo.psi);

if nSamples < 1
    return
end

% min() returns the LOWEST index among ties, which is the deterministic tie-break
% section 5 requires; the selected branch derivative is a genuine element of the
% Clarke generalized gradient.
[~, nearest] = min(rid.ext.pairwiseDistance(Q, Znorm), [], 2);
info.nearest = nearest;

% Strict inequalities: a point exactly at a clip boundary is outside the shell and
% contributes zero -- one of the two one-sided values, taken consistently.
isShell      = evalInfo.dTilde > 1 & evalInfo.dTilde < 1 + opts.psiMax;
info.isShell = isShell;
info.nShell  = sum(isShell);

% J == 0 on EXACT equality, no tolerance. See the header.
if J == 0
    return
end

for jSup = 1:nSup
    if ~isShell(jSup)
        continue
    end
    iNear = nearest(jSup);
    dRaw  = evalInfo.dRaw(jSup);
    if dRaw == 0
        % uvec is 0/0; psi is 0 here anyway, so the term is doubly inert.
        continue
    end
    uvec   = (Znorm(iNear, :) - Q(jSup, :)) / dRaw;
    dJdPsi = (1 / nSup) * (evalInfo.psi(jSup) / J)^(opts.pNorm - 1);
    gradZnorm(iNear, :) = gradZnorm(iNear, :) + ...
        dJdPsi * (1 / epsReq(jSup)) * uvec;
end

end
