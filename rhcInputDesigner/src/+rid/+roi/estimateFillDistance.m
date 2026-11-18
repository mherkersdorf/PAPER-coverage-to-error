function [fillDistance, info] = estimateFillDistance(P, nrm, opts)
%ESTIMATEFILLDISTANCE  Monte-Carlo estimate of the fill distance of a point set over the ROI.
%
%   The fill distance is the worst-case distance from any point of the ROI to its
%   nearest point of P,
%
%       fillDistance = sup_{z in ROI} min_{p in P} || z - p ||_2 ,
%
%   so the region -- not P alone -- defines it. `nrm` supplies that region: this
%   function uses only `nrm.nDim` and `nrm.normalizedBox()`, never the physical
%   bounds. P is in normalized coordinates, as is the returned value; all
%   fill-distance quantities live in the normalized coordinates of the active
%   representation.
%
%   THE ESTIMATE IS BIASED LOW, systematically and in one direction. It is
%   a maximum over finitely many probes of a supremum, so it can only ever
%   under-report; no probe set can over-report. No corner augmentation is applied,
%   so the corners of the ROI -- typically the worst-covered points -- are reached
%   only in the limit. The consequence is not a wrong claim but the criterion's
%   blind spot: a demanded resolution below what Q can resolve means new samples
%   fall between supporting points and register no gain. In the normal regime
%   `fillDistanceQ << epsReq0` this cannot arise.
%
%   THE PROBES ARE LOW-DISCREPANCY, AND THE SHORTFALL IS BRACKETED. This
%   estimator takes a MAXIMUM, and a maximum is decided entirely by coverage: a
%   gap in the probe set is precisely where the deepest hole of P hides. An i.i.d.
%   uniform sample leaves gaps of its own, so the default probe set is a digitally
%   shifted Sobol set (rid.ext.lowDiscrepancySet); "uniform" remains selectable
%   for comparison. Because the distance-to-nearest function is 1-Lipschitz, a
%   probe set whose own fill distance is hProbe cannot miss by more than hProbe,
%
%       hEst  >=  h - hProbe  ,   hence   h  in  [hEst, hEst + hProbe] ,
%
%   so the reported quantity is a BRACKET, not a point estimate with an unknown
%   bias. `info.probeFillDistance` and `info.bracket` carry it.
%
%   THE PROBE SET IS SEEDED WITH THE ROI'S CORNERS (corner seeding).
%   "sobolCorners" is the DEFAULT: the 2^nDim vertices and the 2*nDim face centres
%   come first, then the Sobol set. A design that fills the bulk of a box leaves
%   its deepest holes towards the vertices, and the volume within reach of any
%   vertex collapses with dimension, so a probe set that does not visit them
%   cannot report them. Measured on an lhs design at nSup = 2000, seeding recovers
%   +8.7 %% (nDim 2), +15.5 %% (3), +11.8 %% (4), +30.6 %% (5) and +22.1 %% (6) over
%   the better of "uniform" and "sobol". There is no dimension at which it begins
%   to matter; it matters at all of them.
%
%   THE CEILING IS REPORTED, NOT REFUSED. The seeding costs 2^nDim probes, free to
%   about nDim = 12 and unaffordable beyond. Where 2^nDim + 2*nDim > nProbes the
%   method FALLS BACK to "sobol" and says so, in `info.probeMethod` (the effective
%   method), `info.probeMethodRequested` and `info.probeSeedingFellBack`, all three
%   of which reach the run record. This is the precedent -- a count is
%   adjusted, both values are returned, nothing is silent -- and it is why the
%   fallback is not an error: refusing would make the default unusable in high
%   dimension for a reason the caller cannot act on. The corner set is never
%   TRUNCATED to fit: an arbitrary subset of the vertices is a different probe
%   set, not a smaller one, so the fallback is exactly "sobol".
%
%   HOW HONEST THAT BRACKET IS. hProbe is itself measured by this same estimator,
%   with an independent probe set of the same size drawn from the same stream, so
%   it is a Monte-Carlo estimate and is biased low exactly as the main value is:
%   the upper end of the bracket is the best available statement, NOT a rigorous
%   upper bound. A rigorous hProbe needs a construction whose covering radius is
%   analytic -- a deterministic grid, at a cost of (1/s)^nDim. The bracket
%   is still worth its cost: it converts "biased low by an unknown amount" into
%   "biased low by about this much", which is what a reader needs in order to
%   judge a headroom.
%
%   CONVERGENCE IS DECIDED HERE, NOT BY A READER. The probe maximum is
%   recorded at a quarter, a half and all of the STOCHASTIC probes, at no extra
%   cost, and the relative climb across that trace is compared against
%   `traceTolerance`. A running maximum over a growing probe set can only rise, so
%   the increments measure how far the estimate still has to go; the verdict
%   `info.converged` is read by rid.crit.resolutionReport, because a headroom
%   epsReq0/fillDistanceQ computed from an unconverged estimate is a headroom
%   computed from the wrong number.
%
%   THE TRACE EXCLUDES THE DETERMINISTIC SEEDS, AND THE ESTIMATE DOES NOT (the convergence-trace safeguard). The seeds are probed first and are placed exactly where the
%   deepest holes are, so a trace over all probes would be flat from the first
%   quarter onwards and `converged` would be true by construction -- a safeguard
%   retired by a feature added a week after it. Restricting the trace to the
%   random probes restores the question it was built to ask: HAVE WE PROBED ENOUGH
%   BEYOND THE SEEDS? The reported `fillDistance` remains the maximum over ALL
%   probes, seeds included. For "sobol" and "uniform" nothing changes: they have
%   no seeded prefix, so the trace is over the whole probe set exactly as before.
%
%   The decision is taken on the LARGEST of the increments, not the last one. Both
%   steps of the trace are a doubling of the probe count, so they are comparable
%   and neither is privileged; and a running maximum is a STEP function that sits
%   flat between rare jumps, so a last increment of zero is not evidence of
%   convergence. Measured on a 3-D cell-centred grid whose supremum is known
%   exactly, the last-increment rule called n = 40 converged at 12 % below the
%   truth and n = 200 converged at 9 % below it, while the maximum rule rejected
%   both. A trace with fewer than two entries -- nProbes so small that the
%   quarter, half and full counts collapse -- is NOT converged: nothing was
%   compared.
%
%   Nothing is raised. An under-resolved probe set makes fillDistanceQ too small,
%   which makes the certifiable floor look finer than it is -- worth saying
%   loudly, not worth refusing to run over. The probe count and seed are
%   returned with the value and belong in the run record.
%
%   A rigorous upper bound -- a deterministic probe grid of spacing s, giving
%   max + sqrt(nDim)*s/2 -- is deliberately not implemented.
%
%   Randomness. Probes are drawn from a private RandStream built from `seed`. The
%   global stream is never read and never written, so the value is reproducible in
%   isolation and inside parfor.
%
%   THE CENTRE. `info.centre` is the probe at which the maximum is attained
%   -- the first such probe in probe order -- in normalized coordinates: the
%   estimated centre of the largest empty ball, inside the normalized box by
%   construction. It is where the probe maximum sits, not the exact maximizer:
%   with several holes of nearly equal depth it may lie in any of them.
%
%   PREFIXES IN ONE PASS. With `opts.prefixSizes` = [n_1 < ... < n_m], the
%   estimate is returned for each leading block P(1:n_i, :) -- the fill distance of
%   a design as it grows -- from ONE probe set: each probe's distance to its nearest
%   point is carried from one prefix to the next and lowered by the new points only,
%   and the probe set's own fill distance is measured once, since it does not depend
%   on P. Every per-prefix quantity is bit-identical to a separate call on
%   P(1:n_i, :) with the same options: a minimum and a maximum are exact in any
%   grouping, and each distance is computed pair by pair. Left empty (the default),
%   the call is the single estimate of all of P, exactly as before the prefix estimation.
%
%   Memory. The probe set itself is materialized -- (nProbes x nDim), and measuring
%   its own fill distance requires it as a point set -- but the DISTANCE matrix is
%   built in blocks of `blockSize`, so it is at most (blockSize x nPoints)
%   regardless of `nProbes`. Probe order is preserved, so `blockSize` changes the
%   memory footprint and nothing else: the result is bit-identical across block
%   sizes. Across prefixes only per-prefix maxima are kept, so the memory does not
%   grow with m beyond O(m): every prefix of an 8000-point design costs no more
%   memory than one estimate.
%
%   Cost. The main estimate is nProbes x nPoints distances; the probe set's own
%   fill distance adds nProbes x nProbes, which dominates whenever nProbes >
%   nPoints. Both are paid once, at run start -- and once for all m prefixes, which
%   cost no more distances than the largest of them alone.
%
%   Inputs:
%       P                (nPoints x nDim)  point set in normalized coordinates,
%                                          samples in rows; must be non-empty
%       nrm              (1x1)             rid.roi.Normalization
%       opts.nProbes     (1x1)             positive integer, probe count
%       opts.probeMethod (1x1)             string, "sobolCorners" (default),
%                                          "sobol" or "uniform"
%       opts.seed        (1x1)             nonnegative integer probe seed
%       opts.blockSize   (1x1)             positive integer, probes per block
%       opts.traceTolerance (1x1)          nonnegative, the relative climb across
%                                          the trace above which the estimate is
%                                          reported unconverged
%       opts.prefixSizes (1 x m)           strictly increasing positive integers,
%                                          at most nPoints: estimate P(1:n, :) for
%                                          each n. Default zeros(1, 0): one
%                                          estimate of all of P, m = 1
%
%   Outputs (m = numel(prefixSizes), 1 by default):
%       fillDistance     (1 x m)  per prefix, the probe maximum over ALL probes,
%                                 seeds included
%       info             (1x1) struct with fields
%           .method            (1x1) string, "monteCarlo"
%           .nProbes           (1x1)
%           .probeMethod       (1x1) string, the EFFECTIVE method
%           .probeMethodRequested (1x1) string, as asked for; differs from the
%                                    above only when the ceiling below fired
%           .probeSeedingFellBack (1x1) logical, true when "sobolCorners" could
%                                    not afford its 2^nDim + 2*nDim seeds
%           .nSeededProbes     (1x1)  deterministic probes, 0 unless seeded
%           .nStochasticProbes (1x1)  nProbes - nSeededProbes; the trace's domain
%           .probeFillDistance (1x1) hProbe, the probe set's own fill distance,
%                                    itself estimated and itself biased low
%           .bracket           (m x 2) [fillDistance, fillDistance + hProbe]
%           .centre            (m x nDim) the probe attaining fillDistance
%           .prefixSizes       (1 x m) the prefix lengths; nPoints by default
%           .seed              (1x1)
%           .blockSize         (1x1)
%           .traceTolerance    (1x1)
%           .converged         (1 x m) logical, the verdict on each trace
%           .trace.nProbes     (1x3)  [ceil(s/4), ceil(s/2), s] counted in
%                                     STOCHASTIC probes, s = .nStochasticProbes
%           .trace.value       (m x 3)  running maximum over the stochastic probes
%                                     alone at those counts, non-decreasing per row
%           .trace.relativeClimb (m x 2)  consecutive relative increments of .value;
%                                       the verdict is max(.) <= traceTolerance

arguments
    P   (:,:) double {mustBeReal, mustBeFinite}
    nrm (1,1) rid.roi.Normalization
    opts.nProbes   (1,1) double {mustBeReal} = rid.defaultOptions().roi.nProbes
    opts.probeMethod (1,1) string = rid.defaultOptions().roi.probeMethod
    opts.seed      (1,1) double {mustBeReal} = rid.defaultOptions().roi.probeSeed
    opts.blockSize (1,1) double {mustBeReal} = rid.defaultOptions().roi.blockSize
    opts.traceTolerance (1,1) double {mustBeReal} ...
                                   = rid.defaultOptions().roi.traceTolerance
    opts.prefixSizes (1,:) double {mustBeReal} = zeros(1, 0)
end

if size(P, 2) ~= nrm.nDim
    error("rid:roi:dimensionMismatch", ...
        "P must have nDim = %d columns (got %d).", nrm.nDim, size(P, 2));
end
if isempty(P)
    error("rid:roi:emptyPointSet", ...
        "The fill distance of an empty point set is unbounded; P must be non-empty.");
end
if ~isfinite(opts.nProbes) || mod(opts.nProbes, 1) ~= 0 || opts.nProbes < 1
    error("rid:roi:invalidProbeCount", ...
        "nProbes must be a positive integer (got %g).", opts.nProbes);
end
if ~isfinite(opts.blockSize) || mod(opts.blockSize, 1) ~= 0 || opts.blockSize < 1
    error("rid:roi:invalidProbeCount", ...
        "blockSize must be a positive integer (got %g).", opts.blockSize);
end
if ~ismember(opts.probeMethod, ["sobol", "uniform", "sobolCorners"])
    error("rid:roi:unknownProbeMethod", ...
        "probeMethod must be ""sobol"", ""uniform"" or ""sobolCorners"" " + ...
        "(got ""%s"").", opts.probeMethod);
end
if ~isfinite(opts.traceTolerance) || opts.traceTolerance < 0
    error("rid:roi:invalidTraceTolerance", ...
        "traceTolerance must be finite and nonnegative (got %g).", ...
        opts.traceTolerance);
end
prefixSizes = opts.prefixSizes;
if isempty(prefixSizes)
    prefixSizes = size(P, 1);
elseif any(~isfinite(prefixSizes)) || any(mod(prefixSizes, 1) ~= 0) || ...
        prefixSizes(1) < 1 || any(diff(prefixSizes) <= 0) || prefixSizes(end) > size(P, 1)
    error("rid:roi:invalidPrefixSizes", ...
        "prefixSizes must be strictly increasing integers in [1, %d] (got %s).", ...
        size(P, 1), mat2str(prefixSizes));
end

nProbes = opts.nProbes;
extent  = nrm.dimWeights;                      % normalizedBox is [0, extent]
stream  = RandStream("twister", "Seed", opts.seed);

[probes, seeding] = probeSet(opts.probeMethod, nProbes, nrm.nDim, stream);
probes  = probes .* extent;
nSeeded = seeding.nSeeded;
nStoch  = nProbes - nSeeded;

% Probe counts at which the running maximum is recorded, so that a still-climbing
% estimate is visible in the report rather than assumed away. They count
% STOCHASTIC probes: the deterministic seeds sit at the deepest holes by design,
% so including them would make the trace flat and the verdict free.
traceCounts = unique([ceil(nStoch / 4), ceil(nStoch / 2), nStoch]);
traceCounts = traceCounts(traceCounts >= 1);

% Per prefix, the maximum over ALL probes and where it sits, and the running
% maximum over the stochastic probes at the trace counts -- from one pass.
[fillDistance, iCentre, traceValues] = prefixMaxima(probes, P, opts.blockSize, ...
    prefixSizes, nSeeded, traceCounts);
centre = probes(iCentre, :);

% The probe set's own fill distance, measured the same way with an INDEPENDENT
% set of the same size drawn from the same stream. This is what turns the
% estimate into a bracket: the distance-to-nearest function is 1-Lipschitz, so no
% probe set of fill distance hProbe can miss the true value by more than hProbe.
% It does not depend on P, so it is measured once for every prefix.
verification      = probeSet(opts.probeMethod, nProbes, nrm.nDim, stream) .* extent;
probeFillDistance = prefixMaxima(verification, probes, opts.blockSize, nProbes, ...
    0, zeros(1, 0));

nPrefix       = numel(prefixSizes);
converged     = false(1, nPrefix);
relativeClimb = zeros(nPrefix, max(0, numel(traceCounts) - 1));
for iPrefix = 1:nPrefix
    [converged(iPrefix), relativeClimb(iPrefix, :)] = ...
        traceVerdict(traceValues(iPrefix, :), opts.traceTolerance);
end

info = struct( ...
    "method",                "monteCarlo", ...
    "nProbes",               nProbes, ...
    "probeMethod",           seeding.effective, ...
    "probeMethodRequested",  seeding.requested, ...
    "probeSeedingFellBack",  seeding.fellBack, ...
    "nSeededProbes",         nSeeded, ...
    "nStochasticProbes",     nStoch, ...
    "probeFillDistance",     probeFillDistance, ...
    "bracket",               [fillDistance.', fillDistance.' + probeFillDistance], ...
    "centre",                centre, ...
    "prefixSizes",           prefixSizes, ...
    "seed",                  opts.seed, ...
    "blockSize",             opts.blockSize, ...
    "traceTolerance",        opts.traceTolerance, ...
    "converged",             converged, ...
    "trace",                 struct("nProbes", traceCounts, "value", traceValues, ...
                                    "relativeClimb", relativeClimb));

end

% ---------------------------------------------------------------------------

function [X, seeding] = probeSet(method, nProbes, nDim, stream)
%PROBESET  The probe points, in the unit hypercube, in probe order.
%
%   X (nProbes x nDim). Consecutive calls on one stream give independent sets, so
%   the caller can measure the probe set with a second one.
%
%   seeding (1x1) struct: .requested, .effective, .nSeeded, .fellBack. The two
%   method names differ only when "sobolCorners" could not afford its seeds and
%   fell back -- reported, never silent (corner seeding, the precedent).

seeding = struct("requested", method, "effective", method, ...
                 "nSeeded", 0, "fellBack", false);

switch method
    case "sobol"
        X = rid.ext.lowDiscrepancySet(nProbes, nDim, stream);
    case "sobolCorners"
        nSeeded = 2^nDim + 2 * nDim;
        if nSeeded > nProbes
            % The ceiling. Falling back to the whole "sobol" set rather than to a
            % truncated corner set: an arbitrary subset of the vertices is a
            % different probe set, not a smaller one.
            seeding.effective = "sobol";
            seeding.fellBack  = true;
            X = rid.ext.lowDiscrepancySet(nProbes, nDim, stream);
        else
            seeding.nSeeded = nSeeded;
            X = cornerSeededSet(nProbes, nDim, nSeeded, stream);
        end
    case "uniform"
        % Drawn transposed: the stream is consumed in probe order, one probe's
        % coordinates at a time, so the probe sequence of a shorter run is a
        % PREFIX of a longer one's at the same seed -- which is what makes the
        % convergence trace comparable across probe counts.
        X = rand(stream, nDim, nProbes).';
    otherwise
        error("rid:roi:unknownProbeMethod", ...
            "probeMethod must be ""sobol"", ""uniform"" or ""sobolCorners"" " + ...
            "(got ""%s"").", method);
end

end

% ---------------------------------------------------------------------------

function X = cornerSeededSet(nProbes, nDim, nSeeded, stream)
%CORNERSEEDEDSET  Low-discrepancy probes preceded by the ROI corners.
%
%   X (nProbes x nDim), in the unit hypercube, with nSeeded <= nProbes. The 2^nDim
%   vertices come first, then the 2*nDim face centres, then Sobol points for the
%   remainder. The order is deliberate twice over: it puts the geometrically
%   deepest candidates first, and it makes the seeded prefix a contiguous block
%   that the trace can exclude by index alone.

% Vertices, as the binary expansions of 0 : 2^nDim-1.
corners = zeros(2^nDim, nDim);
for kVertex = 0:2^nDim - 1
    corners(kVertex + 1, :) = bitget(kVertex, 1:nDim);
end

% Face centres: the box centre with one coordinate pushed onto each face.
faces = repmat(0.5, 2 * nDim, nDim);
for kDim = 1:nDim
    faces(2 * kDim - 1, kDim) = 0;
    faces(2 * kDim,     kDim) = 1;
end

X = [corners; faces];
if nProbes > nSeeded
    X = [X; rid.ext.lowDiscrepancySet(nProbes - nSeeded, nDim, stream)];
end

end

% ---------------------------------------------------------------------------

function [runningMax, iMax, traceValues] = prefixMaxima(probes, P, blockSize, ...
                                                      prefixSizes, nSeeded, traceCounts)
%PREFIXMAXIMA  Per prefix of P, the probe maximum of the distance to the nearest point.
%
%   probes (nProbes x nDim), P (nPoints x nDim), prefixSizes (1 x m) increasing,
%   traceCounts (1 x nTrace) counted in the probes AFTER the first nSeeded.
%   runningMax (1 x m)       max over all probes of the distance to P(1:n, :)
%   iMax       (1 x m)       the first probe, in probe order, attaining it
%   traceValues (m x nTrace) max over probes nSeeded + (1:c), for each count c
%
%   The probes are swept in blocks of blockSize; within a block, each probe's
%   distance is carried from one prefix to the next and lowered only by the points
%   the prefix adds, so the distance matrix is at most (blockSize x
%   nPoints) and only O(m) maxima are kept. Maxima are exact and probe order is
%   preserved, so the result does not depend on blockSize or on the prefixes asked.

nProbes     = size(probes, 1);
nPrefix     = numel(prefixSizes);
runningMax  = zeros(1, nPrefix);
iMax        = ones(1, nPrefix);
traceValues = zeros(nPrefix, numel(traceCounts));
for iStart = 1:blockSize:nProbes
    rows   = iStart:min(iStart + blockSize - 1, nProbes);
    dBlock = inf(numel(rows), 1);
    nDone  = 0;
    % The block's probes inside each trace window, fixed across prefixes.
    inTrace = cell(1, numel(traceCounts));
    for iTrace = 1:numel(traceCounts)
        inTrace{iTrace} = find(rows > nSeeded & rows <= nSeeded + traceCounts(iTrace));
    end
    for iPrefix = 1:nPrefix
        added  = nDone + 1:prefixSizes(iPrefix);
        dBlock = min(dBlock, min(rid.ext.pairwiseDistance(probes(rows, :), P(added, :)), [], 2));
        nDone  = prefixSizes(iPrefix);

        % Strictly greater, so a tie keeps the earlier probe.
        [blockMax, iBlock] = max(dBlock);
        if blockMax > runningMax(iPrefix)
            runningMax(iPrefix) = blockMax;
            iMax(iPrefix)       = rows(iBlock);
        end
        for iTrace = 1:numel(traceCounts)
            if ~isempty(inTrace{iTrace})
                traceValues(iPrefix, iTrace) = max(traceValues(iPrefix, iTrace), ...
                    max(dBlock(inTrace{iTrace})));
            end
        end
    end
end

end

% ---------------------------------------------------------------------------

function [converged, relativeClimb] = traceVerdict(traceValues, tol)
%TRACEVERDICT  Draw the conclusion from the convergence trace.
%
%   traceValues (1 x nTrace), non-decreasing. Returns the consecutive relative
%   increments, each taken against the LATER (larger) value so the quantity is a
%   fraction of what is currently reported, and the verdict on their maximum.

relativeClimb = zeros(1, max(0, numel(traceValues) - 1));
for kStep = 1:numel(relativeClimb)
    later = traceValues(kStep + 1);
    if later > 0
        relativeClimb(kStep) = (later - traceValues(kStep)) / later;
    end
end

% An empty trace of increments means the quarter, half and full counts collapsed
% onto one another -- a probe count so small that nothing was compared. That is
% not a converged estimate.
converged = ~isempty(relativeClimb) && max(relativeClimb) <= tol;

end
