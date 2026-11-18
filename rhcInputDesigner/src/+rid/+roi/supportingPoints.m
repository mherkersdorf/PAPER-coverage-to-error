function [Q, fillDistanceQ, info] = supportingPoints(nrm, nSup, opts)
%SUPPORTINGPOINTS  Generate the supporting points Q discretizing the ROI.
%
%   Q discretizes the ROI for the criterion: it is fixed for the duration of a run
%   -- changing Q changes the meaning of J, so a resampled Q is a new run -- and it
%   is returned in NORMALIZED coordinates, as is everything the criterion sees.
%   `toPhysical` is used only where physical units are required: signal
%   generation, constraints on u, and plotting.
%
%   `nSup` is an INPUT. Q is generated from it and `fillDistanceQ` is
%   MEASURED, never predicted: there is no calibrated constant and no inverse
%   "size nSup to reach a target fill distance" mode. `fillDistanceQ` is the
%   certifiable floor -- the finest resolution this discretization can certify,
%   whatever the data -- and the headroom `epsReq0 / fillDistanceQ` is reported at
%   run start by the caller.
%
%   This function uses only `nrm.nDim` and `nrm.normalizedBox()`; it never touches
%   the physical bounds.
%
%   Methods:
%
%   "lhs" (default) -- a maximin Latin hypercube, scaled onto the normalized ROI.
%       For fillDistanceQ -> 0 the coverage bound of the method description
%       becomes tight, which makes a dense uniform sample the rational default.
%       `fillDistanceQ` is a Monte-Carlo estimate and is biased low; its
%       provenance, including the convergence trace, is in `info.fillDistance`.
%
%   "lhsImproved" / "lhsImprovedFull" -- the SAME design as "lhs", then
%       post-processed by the working group's rid.ext.improveLatinHypercube,
%       which swaps coordinates to raise the minimum inter-point distance:
%       moderately for "lhsImproved" (one distance level per pass) and over all
%       distance levels for "lhsImprovedFull". They cost more computation.
%
%       WHAT THEY GUARANTEE, AND WHAT THEY ONLY TEND TO DO (the Latin-hypercube method option). At equal nSup and seed they raise the MINIMUM INTER-POINT
%       DISTANCE -- the criterion the improvement actually maximizes -- without
%       exception, and that is the test contract:
%
%           minDist("lhsImprovedFull") >= minDist("lhsImproved") > minDist("lhs")
%
%       They USUALLY also lower fillDistanceQ, which is the headroom under
%       epsReq0 they are bought for: over 60 measured cases "lhsImprovedFull"
%       beat "lhs" in 58, median 22.7% smaller. That is NOT guaranteed per case,
%       and the chain through "lhsImproved" is not monotone -- the improvement
%       optimizes a PACKING radius (maximin) while fillDistanceQ is a COVERING
%       radius (minimax), and the two are different functionals. fillDistanceQ is
%       measured on every run anyway, so read the measurement, not the method
%       name.
%
%       THE DEFAULT STAYS "lhs": a default that changed every Q would
%       cost a regeneration of everything frozen for a gain nobody has
%       measured here yet.
%
%       Determinism: the improvement contains no randomness, so it is a
%       deterministic function of the "lhs" design and the seed discipline of
%       the Latin-hypercube seed handling is untouched -- UNLESS `opts.improveMaxTime` bites, which returns
%       the best design found so far and makes the result machine-dependent (a
%       faster machine performs more swaps). That is reported in
%       `info.timeCapHit` and never silent.
%
%       THE METRIC (the anisotropy question, settled empirically). The
%       improvement is applied to the design ALREADY SCALED by dimWeights, with
%       the external function's own normalization switched OFF, so the distance
%       it maximizes is the one +crit measures in. This is legitimate because a
%       swap exchanges two points' coordinates WITHIN one dimension: every
%       column keeps its multiset of values, so the Latin structure and the box
%       survive the scaling untouched. Improving the unit cube and scaling
%       afterwards would maximize the wrong metric under non-unit weights; it
%       was measured to be worse there and identical at dimWeights = 1.
%       `info.improveMetric` records which metric was used.
%
%   "grid" -- a regular CELL-CENTRED grid, available for analytic test cases. The
%       coordinates in dimension k are dimWeights(k) * (0.5 : nPerDim-0.5)/nPerDim,
%       so the grid does NOT include the ROI boundary and the analytic fill
%       distance is the half-diagonal of one cell,
%
%           fillDistanceQ = 0.5 * norm(dimWeights ./ nPerDim)
%
%       which reduces to sqrt(nDim) * h / 2 at dimWeights = 1. A grid including the
%       endpoints would cover the corners exactly and give a different number; the
%       convention is fixed here and the tests assert this one.
%
%       The grid uses an EQUAL COUNT PER DIMENSION, hence anisotropic spacing
%       under non-unit dimWeights: nPerDim is the largest integer with
%       nPerDim^nDim <= nSup, so nSup is treated as a BUDGET and never overspent.
%       The achieved nPerDim^nDim therefore generally differs from the requested
%       nSup -- 2000 in 3-D gives 12 per dimension, hence 1728 -- and BOTH are
%       returned in `info`. The count is never silently adjusted: a smaller Q
%       changes the measured fillDistanceQ, and an unreported discrepancy would
%       surface much later as an unexplained difference between the grid and LHS
%       paths.
%
%       nPerDim is computed by exact integer comparison rather than as
%       floor(nSup^(1/nDim)), which loses a layer on perfect powers:
%       1000^(1/3) evaluates to 9.999999999999998.
%
%   Reproducibility of Q comes from STORING it, not from regenerating it: the LHS
%   generator wraps a toolbox algorithm that may change between releases.
%   The generator name, its seeds and the
%   measured fillDistanceQ are recorded as provenance regardless.
%
%   "stored" -- LOADS a set rid.roi.saveSupportingSet wrote, named by
%       opts.supportingSet, rather than generating one. This is what makes an
%       expensively well-distributed Q -- typically "lhsImprovedFull" in a high
%       dimension -- affordable to reuse: the design cost is paid once. It does
%       NOT reduce the cost of evaluating the criterion, which nSup sets
%       regardless of where Q came from.
%
%       THREE THINGS ARE CHECKED, NONE CORRECTED (the pattern). nSup is a
%       DECLARATION and is compared against the file's point count, which is
%       never silently adopted. nDim and the normalized box extent
%       (nrm.dimWeights) must match the file's; loading NEVER RESCALES, because
%       a maximin design stretched per dimension is no longer maximin and
%       nothing about the result looks wrong (the anisotropy argument).
%       Every loaded point must lie inside [0, dimWeights].
%
%       THE STORED fillDistanceQ IS A CHECKSUM, never the live value: it is
%       RE-MEASURED here with this call's own probe settings, and the two
%       measurements' brackets must overlap. Disjoint brackets mean a
%       different box or damaged points -- not a probe seed the caller happened
%       to change, which is why the check is bracket overlap and not a bare
%       tolerance on the two point estimates.
%
%       IDENTITY IS THE CONTENT HASH, not the path: `info.hash` is computed
%       from the loaded Q, `info.origin` is "stored", and `info.seed` is NaN --
%       no generator ran this call.
%
%   Generation is kept behind this one function so that adaptive placement can
%   replace it without touching callers.
%
%   Inputs:
%       nrm                 (1x1)  rid.roi.Normalization
%       nSup                (1x1)  requested number of supporting points, positive
%                                  integer
%       opts.method         (1x1)  string, "lhs" (default), "lhsImproved",
%                                  "lhsImprovedFull", "grid" or "stored"
%       opts.supportingSet  (1x1)  string, the basename of a file under
%                                  supportingSets/; REQUIRED by "stored" and
%                                  refused by every other method
%       opts.improveMaxTime (1x1)  positive, seconds allowed to the improvement
%       opts.progress       (1x1)  string, "auto" (default) prints one liveness
%                                  line for the improving methods only, "off" is
%                                  silent
%       opts.seed           (1x1)  nonnegative integer, supporting-point generator
%       opts.probeSeed      (1x1)  nonnegative integer, fill-distance probes
%       opts.nProbes        (1x1)  positive integer, fill-distance probes
%       opts.probeMethod    (1x1)  string, "sobol" (default) or "uniform"; the
%                                  probe set of the fill-distance estimate
%       opts.blockSize      (1x1)  positive integer, probes per block
%       opts.traceTolerance (1x1)  nonnegative, the relative climb across the
%                                  fill-distance trace above which the estimate is
%                                  reported unconverged
%       opts.lhsCriterion   (1x1)  string, passed to rid.ext.latinHypercube
%       opts.lhsIterations  (1x1)  nonnegative integer, passed through
%       opts.lhsSmooth      (1x1)  logical, passed through
%
%   Outputs:
%       Q              (nSup x nDim)  supporting points, normalized coordinates,
%                                     inside [0, dimWeights]
%       fillDistanceQ  (1x1)          measured fill distance of Q over the ROI
%       info           (1x1) struct with fields
%           .method          (1x1) string, the method used
%           .nDim            (1x1)
%           .nSupRequested   (1x1)  as asked for
%           .nSup            (1x1)  as achieved; equals size(Q,1)
%           .seed            (1x1)  generator seed ("grid" ignores it)
%           .improveElapsed  (1x1)  seconds spent in the improvement; 0 when the
%                                   method does not improve
%           .timeCapHit      (1x1)  logical, true when opts.improveMaxTime cut the
%                                   improvement short -- the design is then
%                                   machine-dependent; false otherwise
%           .improveMetric   (1x1)  string, the metric the improvement maximized:
%                                   "dimWeighted" or "none"
%           .nPerDim         (1x1)  "grid" only, points per dimension
%           .fillDistance    (1x1) struct, provenance of fillDistanceQ:
%                                  the rid.roi.estimateFillDistance info for "lhs",
%                                  or struct("method","gridAnalytic","spacing",...,
%                                  "converged",true). Both carry `.converged`, the
%                                  convergence-trace safeguard verdict on the measurement.
%           .origin          (1x1) string, "generated" for every method above,
%                                  "stored" for "stored". `info` is what
%                                  rid.rhc.Setup's `qInfo` argument is
%                                  typed against; a Q that never went through
%                                  this function is tagged "supplied" by the
%                                  caller instead.
%           .hash            (1x1) string, "" except for "stored", where it is
%                                  the content hash of the loaded Q

arguments
    nrm  (1,1) rid.roi.Normalization
    nSup (1,1) double {mustBeReal}
    opts.method          (1,1) string  = rid.defaultOptions().roi.method
    opts.supportingSet   (1,1) string  = rid.defaultOptions().roi.supportingSet
    opts.improveMaxTime  (1,1) double {mustBeReal, mustBePositive} ...
                                     = rid.defaultOptions().roi.improveMaxTime
    opts.progress        (1,1) string  = rid.defaultOptions().roi.progress
    opts.seed          (1,1) double  = rid.defaultOptions().roi.seed
    opts.probeSeed     (1,1) double  = rid.defaultOptions().roi.probeSeed
    opts.nProbes       (1,1) double  = rid.defaultOptions().roi.nProbes
    opts.probeMethod   (1,1) string  = rid.defaultOptions().roi.probeMethod
    opts.blockSize     (1,1) double  = rid.defaultOptions().roi.blockSize
    opts.traceTolerance (1,1) double = rid.defaultOptions().roi.traceTolerance
    opts.lhsCriterion  (1,1) string  = rid.defaultOptions().roi.lhsCriterion
    opts.lhsIterations (1,1) double  = rid.defaultOptions().roi.lhsIterations
    opts.lhsSmooth     (1,1) logical = rid.defaultOptions().roi.lhsSmooth
end

if ~isfinite(nSup) || mod(nSup, 1) ~= 0 || nSup < 1
    error("rid:roi:invalidSupportingPointCount", ...
        "nSup must be a positive integer (got %g).", nSup);
end

if ~ismember(opts.progress, ["auto", "off"])
    error("rid:roi:invalidProgress", ...
        "progress must be ""auto"" or ""off"" (got ""%s"").", opts.progress);
end

% supportingSet is REQUIRED by "stored" and REFUSED by every other method
%: the illegal combination is unrepresentable rather than sayable and
% then caught, the same discipline the Latin-hypercube method option applies to method itself.
if opts.method == "stored"
    if opts.supportingSet == ""
        error("rid:roi:missingSupportingSet", ...
            "method = ""stored"" requires opts.supportingSet, the basename " + ...
            "of a file rid.roi.saveSupportingSet wrote under " + ...
            "supportingSets/.");
    end
elseif opts.supportingSet ~= ""
    error("rid:roi:supportingSetNotApplicable", ...
        "opts.supportingSet is set but method is ""%s"", not ""stored"". A " + ...
        "supporting set that would not be used is refused rather than " + ...
        "silently ignored.", opts.method);
end

nDim   = nrm.nDim;
extent = nrm.dimWeights;

info = struct("method", opts.method, "nDim", nDim, ...
    "nSupRequested", nSup, "nSup", nSup, "seed", opts.seed, ...
    "improveElapsed", 0, "timeCapHit", false, "improveMetric", "none", ...
    "origin", "generated", "hash", "");

switch opts.method

    case {"lhs", "lhsImproved", "lhsImprovedFull"}
        Q = rid.ext.latinHypercube(nSup, nDim, opts.seed, ...
            "criterion",  opts.lhsCriterion, ...
            "iterations", opts.lhsIterations, ...
            "smooth",     opts.lhsSmooth) .* extent;

        isImproving = opts.method ~= "lhs";
        announce    = isImproving && opts.progress == "auto";

        % The liveness signal, and the ONE place this project prints from a
        % library package: a struct handed back after half an hour
        % cannot tell anybody that the half hour was expected.
        if announce
            announceStart(nSup, nDim, opts.method, opts.improveMaxTime);
        end

        if isImproving
            [Q, info.improveElapsed, info.timeCapHit] = ...
                improveDesign(Q, opts.method, opts.improveMaxTime);
            info.improveMetric = "dimWeighted";
        end

        [fillDistanceQ, fillInfo] = rid.roi.estimateFillDistance(Q, nrm, ...
            "nProbes",        opts.nProbes, ...
            "probeMethod",    opts.probeMethod, ...
            "seed",           opts.probeSeed, ...
            "blockSize",      opts.blockSize, ...
            "traceTolerance", opts.traceTolerance);

        if announce
            announceFinish(info.improveElapsed, fillDistanceQ, ...
                info.timeCapHit, opts.improveMaxTime);
        end

    case "grid"
        % Largest nPerDim with nPerDim^nDim <= nSup. The floating-point root is
        % only a starting guess; the acceptance test is exact integer arithmetic.
        nPerDim = max(1, floor(nSup^(1 / nDim)));
        while (nPerDim + 1)^nDim <= nSup
            nPerDim = nPerDim + 1;
        end
        while nPerDim > 1 && nPerDim^nDim > nSup
            nPerDim = nPerDim - 1;
        end

        % Cell-centred: the grid excludes the ROI boundary, and one cell has
        % spacing dimWeights(k)/nPerDim in dimension k.
        axisCoords = cell(1, nDim);
        for kDim = 1:nDim
            axisCoords{kDim} = extent(kDim) * ((0.5:nPerDim - 0.5).' / nPerDim);
        end
        gridCols = cell(1, nDim);
        [gridCols{:}] = ndgrid(axisCoords{:});
        Q = zeros(nPerDim^nDim, nDim);
        for kDim = 1:nDim
            Q(:, kDim) = gridCols{kDim}(:);
        end

        spacing       = extent ./ nPerDim;
        fillDistanceQ = 0.5 * norm(spacing);

        info.nPerDim = nPerDim;
        info.nSup    = size(Q, 1);          % achieved; may differ from requested
        % `converged` is true by construction: the grid value is analytic, so
        % there is no probe set and nothing to converge. The field is present so
        % that every fillDistance provenance struct answers the question,
        % whichever path produced it.
        fillInfo     = struct("method", "gridAnalytic", "spacing", spacing, ...
                              "converged", true);

    case "stored"
        setPath = supportingSetPath(opts.supportingSet);
        if ~isfile(setPath)
            error("rid:roi:supportingSetNotFound", ...
                "No stored supporting set at %s.", setPath);
        end
        loaded = load(setPath, "stored");
        stored = loaded.stored;

        % THREE CHECKS, NONE CORRECTED (the pattern) -- nSup, nDim and the
        % normalized box extent are all DECLARATIONS compared against the
        % file, never adopted from it.
        if stored.nDim ~= nDim
            error("rid:roi:supportingSetDimensionMismatch", ...
                "%s was saved for nDim = %d; this normalization declares " + ...
                "%d.", opts.supportingSet, stored.nDim, nDim);
        end
        if size(stored.Q, 1) ~= nSup
            error("rid:roi:supportingSetSizeMismatch", ...
                "%s holds %d point(s); nSup declares %d. Two declarations, " + ...
                "compared and neither corrected.", ...
                opts.supportingSet, size(stored.Q, 1), nSup);
        end
        if ~isequal(stored.dimWeights, extent)
            error("rid:roi:supportingSetBoxMismatch", ...
                "%s was saved over a normalized extent of %s; this run " + ...
                "declares %s. Loading never rescales: a maximin design " + ...
                "stretched per dimension is no longer maximin, and nothing " + ...
                "about the result looks wrong.", ...
                opts.supportingSet, mat2str(stored.dimWeights), mat2str(extent));
        end
        outside = any(stored.Q < 0 | stored.Q > extent, 2);
        if any(outside)
            error("rid:roi:supportingSetOutOfBounds", ...
                "%s carries %d of %d point(s) outside [0, dimWeights]. " + ...
                "Loading never rescales or clips.", ...
                opts.supportingSet, sum(outside), size(stored.Q, 1));
        end

        Q = stored.Q;

        % THE STORED VALUE IS A CHECKSUM, re-measured and never trusted on its
        % own. The two the probe bracket BRACKETS must overlap: each is
        % [estimate, estimate + probeFillDistance], biased low in one
        % direction, so a disjoint pair cannot both describe the same true
        % fill distance -- a different box or damaged points, not a probe
        % seed the caller happened to change between save and load.
        [fillDistanceQ, fillInfo] = rid.roi.estimateFillDistance(Q, nrm, ...
            "nProbes",        opts.nProbes, ...
            "probeMethod",    opts.probeMethod, ...
            "seed",           opts.probeSeed, ...
            "blockSize",      opts.blockSize, ...
            "traceTolerance", opts.traceTolerance);

        storedBracket = stored.fillDistanceInfo.bracket;
        liveBracket   = fillInfo.bracket;
        disjoint = liveBracket(1) > storedBracket(2) || ...
                   storedBracket(1) > liveBracket(2);
        if disjoint
            error("rid:roi:supportingSetFillDistanceMismatch", ...
                "%s declares fillDistanceQ = %.6g (bracket [%.6g, %.6g]); " + ...
                "re-measured now at %.6g (bracket [%.6g, %.6g]) -- the two " + ...
                "brackets do not overlap. The stored value is a CHECKSUM, " + ...
                "never the live one: a material disagreement means a " + ...
                "different box or damaged points.", ...
                opts.supportingSet, stored.fillDistanceQ, storedBracket(1), ...
                storedBracket(2), fillDistanceQ, liveBracket(1), liveBracket(2));
        end

        % IDENTITY IS THE CONTENT HASH, not the path; nothing about a
        % generator is invented -- no generator ran THIS call.
        info.origin = "stored";
        info.seed   = NaN;
        info.hash   = contentHash(Q);
        info.nSup   = size(Q, 1);

    otherwise
        error("rid:roi:unknownMethod", ...
            "method must be ""lhs"", ""lhsImproved"", ""lhsImprovedFull"", " + ...
            """grid"" or ""stored"" (got ""%s"").", opts.method);
end

info.fillDistance = fillInfo;

end

% ---------------------------------------------------------------------------

function [Q, elapsed, timeCapHit] = improveDesign(Q, method, maxTime)
%IMPROVEDESIGN  Raise the minimum inter-point distance of Q in place.
%
%   Q (nSup x nDim), already scaled by dimWeights. `normalize` is switched OFF so
%   the criterion maximized is the dimWeights-scaled one; see the header.
%
%   The vendored rid.ext.improveLatinHypercube is NOT edited -- an updated copy
%   from the working group must drop in unchanged. Its four print sites
%   are gated by its own `display` option and its plots by `plot`, both passed
%   off here; its ungated warnings are captured and discarded by evalc.

improveArgs = {'plot', 0, 'display', false, 'normalize', false, ...
               'maxCalculationTime', maxTime};
if method == "lhsImprovedFull"
    improveArgs = [improveArgs, {'worstDistance2improve', Inf}];
end

% This timer strictly encloses the external function's own, so a cap that bit
% inside always shows up as elapsed > maxTime here: the test has no false
% negatives, and its false positives (a natural finish in the last instant
% before the cap) only over-report machine dependence.
improveTimer = tic;
discardedOutput = evalc('Q = rid.ext.improveLatinHypercube(Q, improveArgs{:});'); %#ok<NASGU>
elapsed    = toc(improveTimer);
timeCapHit = elapsed > maxTime;

end

% ---------------------------------------------------------------------------

function announceStart(nSup, nDim, method, maxTime)
%ANNOUNCESTART  One liveness line before a design that can take minutes.

if method == "lhsImprovedFull"
    cheaper = """lhsImproved"" or ""lhs""";
else
    cheaper = """lhs""";
end

fprintf("[rid.roi] Designing supporting points -- nSup = %d, nDim = %d, " + ...
    "method = ""%s"", time cap %g s. This runs BEFORE the first input step " + ...
    "and is not a stall; for a faster design use method = %s, and " + ...
    "opts.roi.progress = ""off"" silences this line.\n", ...
    nSup, nDim, method, maxTime, cheaper);

end

% ---------------------------------------------------------------------------

function announceFinish(elapsed, fillDistanceQ, timeCapHit, maxTime)
%ANNOUNCEFINISH  The closing half of the liveness signal.

if timeCapHit
    capText = sprintf("the %g s TIME CAP WAS HIT, so this design is " + ...
        "machine-dependent", maxTime);
else
    capText = "the time cap was not hit";
end

fprintf("[rid.roi] Supporting points designed in %.1f s -- " + ...
    "fillDistanceQ = %.6g; %s.\n", elapsed, fillDistanceQ, capText);

end
