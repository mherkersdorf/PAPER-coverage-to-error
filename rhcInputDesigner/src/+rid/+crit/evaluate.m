function [J, info] = evaluate(Znorm, Q, epsReq, opts)
%EVALUATE  The resolution-relative space-filling criterion, computed from scratch.
%
%   For a candidate distribution Znorm and supporting points Q:
%
%       dTilde = min_i || (z(i) - q_j) / epsReq(j) ||_2
%       psi    = min( (dTilde - 1)_+ , psiMax )
%       J      = ( mean( psi.^pNorm ) )^(1/pNorm)
%
%   J is MINIMIZED. Only positive deviations count: resolving a region finer than
%   demanded yields no further gain, which is what makes a non-uniform target
%   expressible at all. J = 0 certifies the demanded resolution everywhere in the
%   ROI, up to fillDistanceQ (section 5).
%
%   THIS IMPLEMENTATION IS PERMANENT, NOT A STEPPING STONE. It is the
%   naive, obviously-correct version and it stays as the oracle that
%   rid.crit.CoverageState is tested against. The failure mode of the cached path
%   is a SILENT stale cache -- no error, just slightly worse decisions,
%   indistinguishable from "the criterion is weak" -- so an independent
%   recomputation is the only thing that catches it.
%
%   Independence is the whole value of that oracle, so this function shares no
%   code with the cached path, one level further down than is comfortable:
%
%     - the distances are computed HERE, not through rid.ext.pairwiseDistance. A
%       bug in the shared wrapper would otherwise be invisible to the equivalence
%       test. +ext carries its own tests; this is not a statement about it.
%     - the aggregation dTilde -> psi -> mean is written out again in
%       CoverageState rather than factored into a helper both call, for the same
%       reason.
%
%   Speed is deliberately not a concern here: one loop over the supporting points,
%   vectorized over samples. The vectorization rule of the project's test plan governs the hot
%   path, which is CoverageState.
%
%   Limits, made explicit rather than left to the arithmetic:
%     - an EMPTY Znorm gives dRaw = Inf at every supporting point, hence
%       psi = psiMax and J = psiMax: nothing is covered before the first sample;
%     - pNorm = Inf takes its own branch, J = max(psi). The formula would give
%       Inf^(1/Inf) = NaN;
%     - psiMax = Inf removes the UPPER clip only: psi = max(dTilde - 1, 0), so
%       J = mean((dTilde - 1)_+). That equals mean(dTilde) - 1 -- the plain mean
%       nearest-neighbour form of the earlier criteria -- only while dTilde >= 1 at
%       EVERY supporting point (section 10); with an empty Znorm it is legitimately Inf.
%
%   Inputs:
%       Znorm          (nSamples x nDim)  candidate or realized distribution in
%                                        normalized coordinates; may be empty
%       Q             (nSup x nDim)      supporting points, nSup >= 1
%       epsReq        (nSup x 1)         demanded resolution, finite and > 0
%       opts.psiMax   (1x1)              > 0; Inf allowed
%       opts.pNorm    (1x1)              >= 1; Inf allowed
%
%   Outputs:
%       J     (1x1)
%       info  (1x1) struct with fields
%           .dRaw          (nSup x 1)  raw, unnormalized nearest-neighbour distance
%           .dTilde        (nSup x 1)  dRaw ./ epsReq
%           .psi           (nSup x 1)  the clipped penalty
%           .fracResolved  (1x1)       fraction with dTilde <= 1

arguments
    Znorm          (:,:) double {mustBeReal, mustBeFinite}
    Q             (:,:) double {mustBeReal, mustBeFinite}
    epsReq        (:,1) double {mustBeReal}
    opts.psiMax   (1,1) double {mustBeReal} = rid.defaultOptions().crit.psiMax
    opts.pNorm    (1,1) double {mustBeReal} = rid.defaultOptions().crit.pNorm
end

nSup     = size(Q, 1);
nSamples = size(Znorm, 1);

if nSup < 1
    error("rid:crit:invalidResolution", ...
        "Q must have at least one supporting point.");
end
if numel(epsReq) ~= nSup
    error("rid:crit:invalidResolution", ...
        "epsReq must have one entry per supporting point (got %d for nSup = %d).", ...
        numel(epsReq), nSup);
end
if ~all(isfinite(epsReq)) || any(epsReq <= 0)
    error("rid:crit:invalidResolution", ...
        "epsReq must be finite and strictly positive.");
end
if nSamples > 0 && size(Znorm, 2) ~= size(Q, 2)
    error("rid:crit:dimensionMismatch", ...
        "Znorm and Q must have the same number of columns (got %d and %d).", ...
        size(Znorm, 2), size(Q, 2));
end
if isnan(opts.psiMax) || opts.psiMax <= 0
    error("rid:crit:invalidPenaltyClip", ...
        "psiMax must be strictly positive (Inf is allowed); got %g.", opts.psiMax);
end
if isnan(opts.pNorm) || opts.pNorm < 1
    error("rid:crit:invalidNorm", ...
        "pNorm must be >= 1 (Inf is allowed); got %g.", opts.pNorm);
end

% Raw, unnormalized nearest-neighbour distance of every supporting point. Inf
% where there is no sample at all: before the first sample nothing is covered.
dRaw = inf(nSup, 1);
if nSamples > 0
    for jSup = 1:nSup
        dRaw(jSup) = min(vecnorm(Znorm - Q(jSup, :), 2, 2));
    end
end

dTilde = dRaw ./ epsReq;
psi    = min(max(dTilde - 1, 0), opts.psiMax);

if isinf(opts.pNorm)
    J = max(psi);                       % the worst-case limit p -> Inf
else
    J = mean(psi.^opts.pNorm)^(1 / opts.pNorm);
end

if nargout > 1
    info = struct( ...
        "dRaw",         dRaw, ...
        "dTilde",       dTilde, ...
        "psi",          psi, ...
        "fracResolved", mean(dTilde <= 1));
end

end
