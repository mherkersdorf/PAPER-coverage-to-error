function [gradTheta, info] = horizonGradient(gradZnorm, kIndex, kHorizon, A)
%HORIZONGRADIENT  Select the horizon rows, then compose with A.
%
%   `rid.crit.gradient` returns a row for EVERY sample it was given,
%   history included. `A` carries the HORIZON rows only. The two do not
%   line up, so the caller must select -- and a wrong selection yields a
%   gradient of plausible magnitude and wrong content, which no
%   tolerance will catch.
%
%   Two rules, both:
%
%   - the rows are identified from `kIndex`, NOT from an assumed
%     position. That horizon rows come last is a convention of
%     `rep.build`; relying on it silently is what makes the failure
%     undetectable.
%   - the count is ASSERTED against `size(A,1)/nDim`. The convention is
%     safe to rely on only because this check exists.
%
%   IT IS A `+gen` PACKAGE FUNCTION, not a static of `rid.gen.Direct`.
%   Selecting rows and composing with an affine map is generic composition maths
%   and belongs to no parametrization; `+rhc` calling it was the symptom. It is
%   not a static of `rid.gen.Generator` either -- a base class holding utilities
%   no subclass overrides is a namespace wearing a class's clothes -- and `+gen`
%   already has the right pattern, in `affineMap` and `reachabilityReport`.
%
%   Inputs:
%       gradZnorm  (nSamples x nDim)        from rid.crit.gradient
%       kIndex     (nSamples x 1)           time index of each row
%       kHorizon   (nHorizon x 1)           time indices of the horizon
%       A          (nHorizon*nDim x nTheta) from rid.gen.affineMap
%   Outputs:
%       gradTheta  (nTheta x 1)
%       info       (1x1) struct  .nSelected .nExpected .rows

arguments
    gradZnorm (:,:) double {mustBeReal}
    kIndex    (:,1) double {mustBeReal}
    kHorizon  (:,1) double {mustBeReal}
    A         (:,:) double {mustBeReal}
end
nDim = size(gradZnorm, 2);
if numel(kIndex) ~= size(gradZnorm, 1)
    error("rid:gen:dimensionMismatch", ...
        "kIndex has %d entries; gradZnorm has %d row(s).", ...
        numel(kIndex), size(gradZnorm, 1));
end
rows      = find(ismember(kIndex, kHorizon));
nExpected = size(A, 1) / nDim;
if numel(rows) ~= nExpected
    error("rid:gen:horizonRowMismatch", ...
        "Selected %d horizon row(s) from kIndex, but A carries %d " + ...
        "(size(A,1) = %d over nDim = %d). A wrong selection gives a " + ...
        "gradient of plausible magnitude and wrong content.", ...
        numel(rows), nExpected, size(A, 1), nDim);
end
gradSel   = gradZnorm(rows, :);
gradTheta = (gradSel(:).' * A).';
info = struct("nSelected", numel(rows), "nExpected", nExpected, ...
              "rows", rows(:).');

end
