function [fig, info] = coverage(source, opts)
%COVERAGE  The regressor-space scatter: supporting points in their THREE states.
%
%       [fig, info] = rid.viz.coverage(result)
%       [fig, info] = rid.viz.coverage(rid.io.loadRun(runId), "visible", "off")
%
%   THREE STATES, NOT TWO. "resolved" (dTilde <= 1), "active"
%   (1 < dTilde < 1 + psiMax) and "saturated" (dTilde >= 1 + psiMax) are
%   different things, and the ACTIVE SHELL is the only place design pressure
%   exists: psi is zero on the first and clipped on the third, so neither moves
%   the optimizer. A resolved/unresolved picture merges the last two and hides
%   the distinction that explains why a run stalls -- a run whose remaining
%   points are all saturated has an empty gradient for a completely different
%   reason than one whose points are all resolved. The realized samples
%   are drawn OVER the supporting points, so what covered what is readable.
%
%   THE DIMENSION IS NEVER SILENTLY DROPPED. `nDim = 2` is a scatter, `nDim = 3`
%   a scatter3, and above that a matrix of ALL pairwise projections -- never the
%   first two coordinates of a higher-dimensional space, which is a picture of a
%   different problem. Each projection panel is LABELLED AS A PROJECTION,
%   because resolution is a property of the FULL space: a point resolved in
%   nDim can look far from everything in a 2-D shadow, and a point that looks
%   covered there can be unresolved. A reader who takes a projection for the
%   truth is the failure that label prevents.
%
%   A RECORD WITHOUT `Q` CANNOT PRODUCE THIS FIGURE. `opts.io.storeQ = false`
%   makes a record non-replayable and it makes it
%   un-drawable for the same reason: `Q` is not regenerated from a seed, because
%   the generating algorithm belongs to a toolbox and may change between
%   releases. This reports what is missing and stops, in the manner of
%   `rid.io.replay`, rather than drawing a picture of points the run never used.
%
%   THE REALIZED LAYER IS ALWAYS DRAWN. It used to be withheld for a
%   run whose surrogate had been refitted: `Z` is derived and never stored,
%   a `Setup` is immutable so `setup.srg` is the model the run BEGAN
%   with, and invariant 3 remaps every past coordinate under the new one. The
%   record now carries the FINAL surrogate in `outputs`, which is exactly
%   sufficient, so there is no longer a case to withhold.
%
%   NO CRITERION MATHS AT ALL. The three states come from
%   `rid.crit.coverageReport`, which owns the classification; this file draws
%   the answer and forms no part of the definition. Nothing in `+viz` computes
%   dTilde, psi or a state boundary.
%
%   Inputs:
%       source        (1x1) struct  a rid.rhc.run result or a run record
%       opts.visible  (1x1) matlab.lang.OnOffSwitchState = "on"
%       opts.showRealized (1x1) logical = true   draw the realized samples
%       opts.maxPanels    (1x1) double  = 15     projection panels above which
%                                                the matrix is refused rather
%                                                than drawn illegibly
%   Outputs:
%       fig   (1x1) matlab.ui.Figure
%       info  (1x1) struct
%           .runId        (1x1) string
%           .nSup         (1x1) double
%           .count        (1 x 3) double   resolved, active, saturated
%           .label        (1 x 3) string
%           .nRealized    (1x1) double
%           .isProjection (1x1) logical    true when nDim > 3
%           .nAxes        (1x1) double     panels drawn
%           .pairs        (nAxes x 2) double  the coordinate pair per panel

arguments
    source                (1,1) struct
    opts.visible          (1,1) matlab.lang.OnOffSwitchState = "on"
    opts.showRealized     (1,1) logical = true
    opts.maxPanels        (1,1) double {mustBeInteger, mustBePositive} = 15
end

v = viewOf(source);

if ~v.hasQ
    error("rid:viz:recordWithoutQ", ...
        "This record does not store Q (opts.io.storeQ was false), so the " + ...
        "regressor-space scatter cannot be drawn. Q is not regenerated from " + ...
        "its seed: the generating algorithm belongs to a toolbox and may " + ...
        "change between releases, so the points would be a picture of a " + ...
        "design this run never used. " + ...
        "Re-record with opts.io.storeQ.");
end

% THE CLASSIFICATION IS NOT THIS PACKAGE'S TO MAKE. `+viz` asks the
% criterion what the states are and draws the answer; it forms no part of the
% definition itself.
st    = rid.crit.coverageReport(v.state);
pairs = panelPairs(v.nDim);

Zreal = zeros(0, v.nDim);
if opts.showRealized
    Zreal = realizedPoints(v);
end

if v.nDim > 3 && size(pairs, 1) > opts.maxPanels
    error("rid:viz:tooManyPanels", ...
        "nDim = %d gives %d pairwise projections, above maxPanels = %d. " + ...
        "Raise maxPanels deliberately, or select coordinates first: a " + ...
        "matrix nobody can read is not a diagnostic.", ...
        v.nDim, size(pairs, 1), opts.maxPanels);
end

% THE ROI IS THE FRAME, not the data's own extent. Coverage is a statement
% about a REGION, and axes cropped to the points drawn make a design that
% reached a third of the space look like one that filled it.
roiBox = v.nrm.normalizedBox();

fig = newFigure("coverage - " + v.runId, opts.visible);

if v.nDim <= 3
    ax = axes(fig);
    drawPanel(ax, v, st, Zreal, pairs(1, :), false, roiBox);
    axList = ax;
else
    % ALL pairs, never the first two. The layout is square-ish for legibility
    % only; the panel order is the pair order, which is coordinate order.
    nCols  = ceil(sqrt(size(pairs, 1)));
    nRows  = ceil(size(pairs, 1) / nCols);
    tl     = tiledlayout(fig, nRows, nCols, "TileSpacing", "compact", ...
                         "Padding", "compact");
    axList = gobjects(size(pairs, 1), 1);
    for iPanel = 1:size(pairs, 1)
        axList(iPanel) = nexttile(tl);
        drawPanel(axList(iPanel), v, st, Zreal, pairs(iPanel, :), true, roiBox);
    end
    title(tl, sprintf("coverage - %s", v.runId), "Interpreter", "none");
    % THE LABEL IS THE POINT, not decoration: resolution is a property of the
    % full space and every panel below is a shadow of it.
    subtitle(tl, sprintf("PAIRWISE PROJECTIONS of a %d-D space. Resolution " + ...
        "is a property of the FULL space: a point resolved in %d-D can look " + ...
        "far from everything here, and the reverse.", v.nDim, v.nDim));
end

legendOn(axList(1));

info = struct( ...
    "runId",        v.runId, ...
    "nSup",         v.nSup, ...
    "count",        st.count, ...
    "label",        st.stateName, ...
    "nRealized",    size(Zreal, 1), ...
    "isProjection", v.nDim > 3, ...
    "nAxes",        numel(axList), ...
    "pairs",        pairs);

end

% ---------------------------------------------------------------------------

function drawPanel(ax, v, st, Zreal, pair, isProjection, roiBox)
%DRAWPANEL  One panel: the three states, then the realized samples over them.
%
%   Every supporting point is drawn exactly once, in exactly one series, so the
%   three series together hold nSup points in every panel.
style = panelStyle();
hold(ax, "on");

is3D = numel(pair) == 3;
for iState = 1:3
    sel = st.(style(iState).mask);
    P   = v.Q(sel, pair);
    if is3D
        scatter3(ax, P(:, 1), P(:, 2), P(:, 3), style(iState).size, ...
            "filled", "MarkerFaceColor", style(iState).color, ...
            "MarkerFaceAlpha", style(iState).alpha, ...
            "DisplayName", style(iState).name, "Tag", style(iState).tag);
    else
        scatter(ax, P(:, 1), P(:, 2), style(iState).size, ...
            "filled", "MarkerFaceColor", style(iState).color, ...
            "MarkerFaceAlpha", style(iState).alpha, ...
            "DisplayName", style(iState).name, "Tag", style(iState).tag);
    end
end

if ~isempty(Zreal)
    R = Zreal(:, pair);
    if is3D
        scatter3(ax, R(:, 1), R(:, 2), R(:, 3), 14, "filled", ...
            "MarkerFaceColor", [0.10 0.10 0.10], "MarkerFaceAlpha", 0.85, ...
            "DisplayName", "realized", "Tag", "realized");
    else
        scatter(ax, R(:, 1), R(:, 2), 14, "filled", ...
            "MarkerFaceColor", [0.10 0.10 0.10], "MarkerFaceAlpha", 0.85, ...
            "DisplayName", "realized", "Tag", "realized");
    end
end

xlabel(ax, v.dimName(pair(1)), "Interpreter", "none");
ylabel(ax, v.dimName(pair(2)), "Interpreter", "none");
if is3D
    zlabel(ax, v.dimName(pair(3)), "Interpreter", "none");
    view(ax, 3);
end
grid(ax, "on");
box(ax, "on");
axis(ax, "equal");

% The declared extent, never the data's. A margin only, so a point ON the
% boundary is still drawn whole.
margin = 0.02 * max(roiBox(pair, 2) - roiBox(pair, 1));
xlim(ax, roiBox(pair(1), :) + [-margin margin]);
ylim(ax, roiBox(pair(2), :) + [-margin margin]);
if is3D
    zlim(ax, roiBox(pair(3), :) + [-margin margin]);
end

if isProjection
    title(ax, sprintf("projection %s-%s", v.dimName(pair(1)), ...
                      v.dimName(pair(2))), "Interpreter", "none", ...
          "FontWeight", "normal");
else
    title(ax, sprintf("coverage - %s", v.runId), "Interpreter", "none");
end
hold(ax, "off");

end

% ---------------------------------------------------------------------------

function pairs = panelPairs(nDim)
%PANELPAIRS  Which coordinates each panel draws.
%
%   Two and three dimensions are drawn whole; above that every pair is drawn,
%   in coordinate order, so no coordinate is silently omitted.
switch nDim
    case 2
        pairs = [1 2];
    case 3
        pairs = [1 2 3];
    otherwise
        pairs = nchoosek(1:nDim, 2);
end

end

% ---------------------------------------------------------------------------

function style = panelStyle()
%PANELSTYLE  One entry per state, in the order resolved, active, saturated.
%
%   The active shell is the loud one on purpose: it is the only state that can
%   still move the design.
style = struct( ...
    "tag",   {"resolved", "active", "saturated"}, ...
    "mask",  {"isResolved", "isActive", "isSaturated"}, ...
    "name",  {"resolved (dTilde <= 1)", ...
              "active (1 < dTilde < 1+psiMax)", ...
              "saturated (dTilde >= 1+psiMax)"}, ...
    "color", {[0.62 0.72 0.62], [0.85 0.33 0.10], [0.20 0.29 0.55]}, ...
    "size",  {36, 64, 46}, ...
    "alpha", {0.55, 0.95, 0.80});

end

% ---------------------------------------------------------------------------

function legendOn(ax)
%LEGENDON  One legend, on the first panel only.
lg = legend(ax, "Location", "bestoutside");
lg.Interpreter = "none";

end
