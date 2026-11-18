function [levels, info] = targetLevelsFor(srg, rep, roi, nInputs, nLevels, nProbes)
%TARGETLEVELSFOR  Input levels whose resting OUTPUTS -- or resting STATES -- tile the ROI.
%
%   A start that spans the input range says nothing about the coordinates the
%   criterion actually scores. The process low-passes, and where the static map
%   is nonlinear a uniform draw in `u` is badly clustered in the predicted
%   coordinates -- the regressor cloud comes out a flat sheet and the start is
%   worthless however well it spans the input. Those coordinates are half the
%   regressor space, so aiming at them is the point.
%
%   TWO BASES, ONE CONSTRUCTION. Which coordinates are aimed at follows from the
%   representation's KINDS, never from its class:
%
%       "output"      an output coordinate exists: tile the ROI's output range,
%                     the original construction
%       "stateLocus"  state coordinates and NO output coordinate: tile the
%                     STEADY-STATE LOCUS inside the declared state box
%                    
%
%   The locus case is what `rid.rep.StateSpace` needs. It has no output
%   coordinate at all, so the output branch would find nothing to tile and every
%   state-space run would start unguided -- the fallback the surrogate-guided initialization exists to replace.
%   The study compares the two spaces, so an initialization that is guided in one and
%   blind in the other would confound exactly the comparison the milestone is
%   for. A representation with BOTH kinds keeps the output basis: that is the
%   case the surrogate-guided initialization was written for, and the stored canary runs must not move (the first-canary rule,
%   the path canary).
%
%   THE LOCUS IS BUILT THROUGH `rep.build`, one probe at a time, rather than by
%   reading `srg.steadyState` into state columns directly. The k-th state
%   coordinate of a representation need not be the k-th state of the surrogate,
%   and inventing that rule here is exactly the kind of second mapping that
%   drifts. Arc length is measured in the DECLARED box, so it does not
%   depend on the physical units of the individual states, and only probes inside
%   the box count -- a level whose rest state the ROI excludes is not a level.
%   Consecutive probes only: a gap left by an excursion outside the box is a
%   piece of locus that was not probed, so its chord is not added.
%
%   IT LIVES IN `+rhc`, AND IS PASSED TO `+gen` AS A VECTOR OF NUMBERS. Handing
%   `initialGuesses` a surrogate would couple `+gen` to `+surrogate` and `+roi`;
%   a vector keeps it pure and deterministically testable, puts the
%   model-dependent step where the model lives, and makes the guided and blind
%   cases the same code path with a different vector.
%
%   IT USES ONLY THE MINIMAL SURROGATE INTERFACE -- `steadyState` and one
%   `simulate` step -- so it works for a black-box model, which has no time
%   constant to read but always has a rest state and a step response that can be
%   MEASURED. A non-monotone static map makes this a lookup rather
%   than an inverse, which is what it already is.
%
%   CALLED AT SETUP AND THEN ONCE PER SURROGATE UPDATE, never per step: the map it
%   inverts changes only when the model does, and `Setup` holds a current model from
%   the first step. Called only in the update branch it was never built at all in a
%   default run, since `updateInterval` is Inf.
%
%   WHAT RETURNS NOTHING IS REPORTED, never worked around. With several input
%   channels a target is an input COMBINATION and only a whole row has a defined
%   resting value, so `rid.gen.Direct` refuses `targetLevels` there and MISO runs
%   initialize unguided. A representation with neither an output
%   coordinate nor state coordinates has nothing to aim at. And a locus with
%   fewer than two probes inside the declared box has no arc length to spread
%   levels along -- the ROI and the model's static map disagree, which
%   `reachabilityReport` is the place to point at. Nor has a locus whose
%   inside probes span ZERO arc length, which is the same refusal for a different
%   reason: nothing disagrees, the locus is simply a point.
%
%   THE MAPPED LEVELS ARE NOT CLIPPED HERE. The ROI declares its `u` and its
%   predicted bounds independently, so a level's input can fall outside
%   the `u` range; `rid.gen.Direct.initialGuesses` clips into the ROI range and
%   the amplitude box, and `reachabilityReport` reports the inconsistency.
%   Clipping twice would hide which declaration was the loose one.
%
%   Inputs:
%       srg      (1x1) object   the current surrogate
%       rep      (1x1) object   the representation; coordinateInfo is read, never
%                               a coordinate NAME
%       roi      (1x1) struct   the ROI declaration
%       nInputs  (1x1) double   input channels
%       nLevels  (1x1) double   levels to tile the range or the locus with
%       nProbes  (1x1) double   probes of the static map per level search
%   Outputs:
%       levels   (nLevels x nInputs) double, or (0 x 0) when unavailable
%       info     (1x1) struct  .available .reason .basis .yRange .uRange
%                              .yAchieved .xAchieved

levels = zeros(0, 0);
info   = struct("available", false, "reason", "", "basis", "", ...
                "yRange", [NaN NaN], "uRange", [NaN NaN], ...
                "yAchieved", zeros(0, 1), "xAchieved", zeros(0, 0));

if nInputs > 1
    info.reason = "targetLevels are input COMBINATIONS and only a whole row " + ...
        "has a defined resting output, so multi-input runs initialize " + ...
        "unguided.";
    return
end

coords = rep.coordinateInfo();
kinds  = [coords.kind];

% Classify by KIND, never by parsing a name.
outDim   = find(kinds == "output", 1);
stateDim = find(kinds == "state");
if isempty(outDim) && isempty(stateDim)
    info.reason = "the representation has no output coordinate and no state " + ...
        "coordinate, so there is no predicted range to tile.";
    return
end

inputDims = find(kinds == "inputLevel" & [coords.channel] == 1);
if isempty(inputDims)
    info.reason = "the representation has no input-level coordinate for " + ...
        "channel 1, so there is no u range to search.";
    return
end
% The lowest lag is the coordinate closest to the value being decided.
[~, iLowest] = min([coords(inputDims).lag]);
inDim        = inputDims(iLowest);

box         = rep.bounds(roi);
info.uRange = box(inDim, :);

% The probe grid over the declared u range, shared by both bases: the static map
% is PROBED rather than inverted, which is all the minimal surrogate interface
% guarantees.
uProbe = linspace(info.uRange(1), info.uRange(2), nProbes).';

if isempty(outDim)
    info.basis     = "stateLocus";
    [levels, info] = tileStateLocus(srg, rep, box, stateDim, uProbe, nLevels, info);
else
    info.basis     = "output";
    info.yRange    = box(outDim, :);
    [levels, info] = tileOutputRange(srg, uProbe, nLevels, info);
end

end

% ---------------------------------------------------------------------------

function [levels, info] = tileOutputRange(srg, uProbe, nLevels, info)
%TILEOUTPUTRANGE  Levels spread uniformly over the ROI's output range.
%
%   One rest state and one `simulate` step per probe, so a black-box model needs
%   nothing but the minimal interface to be probed this way.
%
%   Inputs:
%       uProbe  (nProbes x 1)  the probe grid over the ROI's u range
%   Outputs:
%       levels  (nLevels x 1)
%       info    (1x1) struct   .available .yAchieved filled in

nProbes = numel(uProbe);
yProbe  = zeros(nProbes, 1);
for iProbe = 1:nProbes
    xRest          = srg.steadyState(uProbe(iProbe));
    yProbe(iProbe) = srg.simulate(uProbe(iProbe), xRest);
end

yTarget = linspace(info.yRange(1), info.yRange(2), nLevels).';
levels  = zeros(nLevels, 1);
for iLevel = 1:nLevels
    [~, iNear]     = min(abs(yProbe - yTarget(iLevel)));
    levels(iLevel) = uProbe(iNear);
end

info.available = true;
info.yAchieved = interp1(uProbe, yProbe, levels, "linear");

end

% ---------------------------------------------------------------------------

function [levels, info] = tileStateLocus(srg, rep, box, stateDim, uProbe, nLevels, info)
%TILESTATELOCUS  Levels spread equally in arc length along the rest locus.
%
%   Inputs:
%       box       (nDim x 2)     the declared bounds, in coordinate order
%       stateDim  (1 x nState)   which coordinates are states
%       uProbe    (nProbes x 1)  the probe grid over the ROI's u range
%   Outputs:
%       levels    (nLevels x 1)  or (0 x 0) when the locus cannot be tiled
%       info      (1x1) struct   .available .reason .xAchieved filled in

levels  = zeros(0, 0);
nProbes = numel(uProbe);
nState  = numel(stateDim);
xProbe  = zeros(nProbes, nState);

for iProbe = 1:nProbes
    xRest = srg.steadyState(uProbe(iProbe));
    % THROUGH THE REPRESENTATION. A one-sample history at rest, whose measured
    % and simulated states are the same state, so the row comes out the same in
    % either pastSource -- and no rule mapping the k-th state coordinate to
    % surrogate state k is invented here.
    histProbe = struct("u",  uProbe(iProbe), "x0", xRest, ...
                       "x",  xRest,          "tS", srg.tS);
    Zprobe    = rep.build(histProbe, srg);
    if isempty(Zprobe)
        info.reason = "the representation yields no row for a one-sample " + ...
            "history, so its steady-state locus cannot be probed through it.";
        return
    end
    xProbe(iProbe, :) = Zprobe(end, stateDim);
end

% Scaled by the DECLARED box, so arc length does not depend on the physical units
% of the individual states, and the ROI is what decides "inside".
lo     = box(stateDim, 1).';
hi     = box(stateDim, 2).';
inside = all(xProbe >= lo & xProbe <= hi, 2);
if sum(inside) < 2
    info.reason = sprintf("%d of %d probed rest states lie inside the " + ...
        "declared state box; two are needed to measure arc length along the " + ...
        "locus. The ROI and the model's static map disagree.", ...
        sum(inside), nProbes);
    return
end

xNorm = (xProbe - lo) ./ (hi - lo);
iIn   = find(inside);
s     = zeros(numel(iIn), 1);
for ii = 2:numel(iIn)
    chord = 0;
    if iIn(ii) == iIn(ii - 1) + 1
        % CONSECUTIVE PROBES ONLY: across a gap the locus left the box, and what
        % happened in between was not probed, so its chord is not arc length.
        chord = norm(xNorm(iIn(ii), :) - xNorm(iIn(ii - 1), :));
    end
    s(ii) = s(ii - 1) + chord;
end

% A LOCUS OF ZERO ARC LENGTH IS THE SAME REFUSAL AS FEWER THAN TWO PROBES
% Two probes inside the box are necessary but not
% sufficient: if the input does not move the rest state, or if no two inside
% probes are consecutive, every target lands on the first probe and the tiling
% returns nLevels copies of one input while reporting them available. That is a
% guided initialization in name only -- worse than the unguided fallback, which
% at least spreads its starts. What cannot be delivered is reported.
if s(end) <= 0
    levels      = zeros(0, 0);
    info.reason = sprintf("the %d probed rest states inside the declared " + ...
        "state box span zero arc length, so there is nothing to spread " + ...
        "levels along: the input does not move the rest state, or no two " + ...
        "probes inside the box are consecutive.", sum(inside));
    return
end

sTarget   = linspace(s(1), s(end), nLevels).';
levels    = zeros(nLevels, 1);
xAchieved = zeros(nLevels, nState);
for iLevel = 1:nLevels
    [~, iNear] = min(abs(s - sTarget(iLevel)));
    levels(iLevel)       = uProbe(iIn(iNear));
    xAchieved(iLevel, :) = xProbe(iIn(iNear), :);
end

info.available = true;
info.xAchieved = xAchieved;

end
