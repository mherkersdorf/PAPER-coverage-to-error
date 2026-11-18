function gradTheta = gradientToTheta(gradPhys, A, scale, dimWeights)
%GRADIENTTOTHETA  A physical-space gradient, converted and composed.
%
%   The penalty works in PHYSICAL units -- gMid, gHalf and the
%   admissible width are physical -- while `A` maps into NORMALIZED
%   coordinates. The two do not compose as written. Normalization is
%   affine and diagonal, `z_n = (z_p - centre) ./ scale .* dimWeights`,
%   so one diagonal factor converts the gradient:
%
%       dPen/dZnorm = (dPen/dZphys) * diag(scale ./ dimWeights)
%
%   and only then composes with `A`. The correction is arithmetically
%   trivial and structurally easy to omit, which is why the two arrays
%   no longer share a name. With `scale` and
%   `dimWeights` both 1 the factor is the identity and a test that uses
%   those values proves nothing.
%
%   This composes a SINGLE term, not an objective: assembling
%   `J + penalty` is +opt's.
%
%   IT IS A `+gen` PACKAGE FUNCTION, not a static of `rid.gen.Direct`.
%   Selecting rows and composing with an affine map is generic composition maths
%   and belongs to no parametrization; `+rhc` calling it was the symptom. It is
%   not a static of `rid.gen.Generator` either -- a base class holding utilities
%   no subclass overrides is a namespace wearing a class's clothes -- and `+gen`
%   already has the right pattern, in `affineMap` and `reachabilityReport`.
%
%   Inputs:
%       gradPhys    (nRows x nDim)       d(term)/d(Zphys)
%       A           (nRows*nDim x nTheta) from rid.gen.affineMap
%       scale       (1 x nDim)           the normalization span
%       dimWeights  (1 x nDim)
%   Outputs:
%       gradTheta   (nTheta x 1)

arguments
    gradPhys   (:,:) double {mustBeReal}
    A          (:,:) double {mustBeReal}
    scale      (1,:) double {mustBeReal, mustBePositive}
    dimWeights (1,:) double {mustBeReal, mustBePositive}
end
if size(gradPhys, 2) ~= numel(scale) || numel(scale) ~= numel(dimWeights)
    error("rid:gen:dimensionMismatch", ...
        "gradPhys has %d column(s), scale %d, dimWeights %d.", ...
        size(gradPhys, 2), numel(scale), numel(dimWeights));
end
gradNorm = gradPhys .* (scale ./ dimWeights);
if size(A, 1) ~= numel(gradNorm)
    error("rid:gen:dimensionMismatch", ...
        "A has %d row(s); the gradient has %d entries.", ...
        size(A, 1), numel(gradNorm));
end
gradTheta = (gradNorm(:).' * A).';

end
