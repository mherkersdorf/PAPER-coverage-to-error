function path = saveSupportingSet(Q, nrm, fillDistanceQ, info, name)
%SAVESUPPORTINGSET  Persist a supporting-point set for reuse.
%
%   Amortizes the cost of an expensively well-distributed `Q` -- typically
%   "lhsImprovedFull" in a high dimension -- across every run that reuses it:
%   `opts.roi.method = "stored"` loads what this function writes. It does NOT
%   reduce the cost of EVALUATING the criterion, which `nSup` sets regardless
%   of where `Q` came from.
%
%   THE PROVENANCE TRAVELS INSIDE THE SAME FILE, so it cannot drift from the
%   points it describes: `info`, exactly the third output of
%   rid.roi.supportingPoints, is stored alongside `Q`, and `fillDistanceQ`'s own
%   MEASUREMENT provenance -- the bracket of the probe bracket -- is what loading later
%   re-checks against, never the cached value taken on faith.
%
%   ONLY A MONTE-CARLO fillDistanceQ CAN BE STORED (`info.fillDistance.method`
%   must be `"monteCarlo"`, i.e. `info.fillDistance.bracket` must exist): the
%   analytic grid value has no bracket to re-check a later measurement against,
%   and reuse is bought for the expensive designs, never the cheap grid.
%
%   RECORDS ARE IMMUTABLE, on the `rid.io.saveRun` precedent: an existing file at the same name is never overwritten,
%   because a stored set silently replaced out from under a frozen
%   `fillDistanceQ` elsewhere is exactly the anisotropy-by-stealth the Latin-hypercube method option warns
%   against, one level up.
%
%   Format and location: `.mat` under a git-ignored `supportingSets/`, resolved
%   through `rid.io.projectRoot` -- never a second path definition, and
%   never JSON: `%.17g` at `1e5 x 6` exceeds a hundred megabytes.
%
%   Inputs:
%       Q             (nSup x nDim) double  supporting points, NORMALIZED,
%                                           inside [0, nrm.dimWeights]
%       nrm           (1x1) rid.roi.Normalization  the space Q discretizes
%       fillDistanceQ (1x1) double  the MEASURED fill distance of Q
%       info          (1x1) struct  exactly the third output of
%                                   rid.roi.supportingPoints: .method, .seed,
%                                   .origin, .fillDistance (with .bracket), ...
%       name          (1x1) string  the file's basename under supportingSets/
%   Outputs:
%       path          (1x1) string  where it was written

arguments
    Q             (:,:) double {mustBeReal, mustBeFinite}
    nrm           (1,1) rid.roi.Normalization
    fillDistanceQ (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
    info          (1,1) struct
    name          (1,1) string
end

if size(Q, 2) ~= nrm.nDim
    error("rid:roi:dimensionMismatch", ...
        "Q has %d column(s); the normalization declares %d. Supporting " + ...
        "points are NORMALIZED coordinates, samples in rows.", ...
        size(Q, 2), nrm.nDim);
end
if isempty(Q)
    error("rid:roi:emptyPointSet", "Q is empty; there is nothing to store.");
end
extent  = nrm.dimWeights;
outside = any(Q < 0 | Q > extent, 2);
if any(outside)
    error("rid:roi:supportingSetOutOfBounds", ...
        "%d of %d point(s) lie outside [0, dimWeights] = %s. A stored set " + ...
        "is refused rather than saved out of its own declared box.", ...
        sum(outside), size(Q, 1), mat2str(extent));
end
if ~isfield(info, "fillDistance") || ~isfield(info.fillDistance, "method") || ...
        info.fillDistance.method ~= "monteCarlo" || ...
        ~isfield(info.fillDistance, "bracket")
    error("rid:roi:unsupportedFillDistanceProvenance", ...
        "info.fillDistance.method must be ""monteCarlo"" with a .bracket " + ...
        ": only a Monte-Carlo fillDistanceQ has the bracket a later " + ...
        "load re-checks against. An analytic ""grid"" design is cheap to " + ...
        "regenerate and is not what reuse is bought for.");
end

path = supportingSetPath(name);
if isfile(path)
    error("rid:roi:supportingSetExists", ...
        "%s already exists. Stored sets are immutable, on the run-record " + ...
        "precedent: a set silently replaced out from under a frozen " + ...
        "fillDistanceQ elsewhere is exactly the failure the stored-Q reuse exists to " + ...
        "prevent. Choose a different name.", path);
end
parent = fileparts(path);
if ~isfolder(parent)
    mkdir(parent);
end

stored = struct( ...
    "formatVersion",   1, ...
    "Q",               Q, ...
    "nDim",            nrm.nDim, ...
    "dimWeights",      extent, ...
    "fillDistanceQ",   fillDistanceQ, ...
    "fillDistanceInfo", info.fillDistance, ...
    "generation",      info, ...
    "savedAt",         string(datetime("now", "Format", "uuuu-MM-dd'T'HHmmss")));

save(path, "stored", "-v7.3");

end
