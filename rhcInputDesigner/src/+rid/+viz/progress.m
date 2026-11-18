function [fig, info] = progress(source, opts)
%PROGRESS  J over time, drawn as the segments within which it is comparable.
%
%       [fig, info] = rid.viz.progress(result)
%       [fig, info] = rid.viz.progress(rid.io.loadRun(runId), "visible", "off")
%
%   J IS MONOTONE ONLY WITHIN A SEGMENT. The loop resets its monotonicity
%   baseline at every applied surrogate update, every applied resolution change
%   and every the coverage reset coverage reset, and its assert never spans one -- because
%   the values on either side measure different things: a rebuilt cache, a
%   different epsReq, or a cache deliberately emptied. So the trace is drawn as
%   SEPARATE segments with the boundaries visible, never as one continuous
%   line, which would imply a comparability that does not hold and would show a
%   reader a rise the loop never saw.
%
%   THE MARKS ARE GENERIC. Every entry in the event log is marked, whichever of
%   the seven kinds it is, not only the kinds
%   that break a segment and not only `"surrogateUpdate"`. A jump in J with no
%   mark beside it is the puzzling plot months later that the canaries exist to
%   prevent, and a mark that appears only for the kinds this file happens to
%   know would make a newly-emitted kind invisible at exactly the moment
%   somebody needs to see it. Warning-severity events are drawn distinctly:
%   consecutive coverage resets are a statement about the CONFIGURATION
%   and are the one thing on this axes a reader must not skim past.
%
%   EVERY DRAWN OBJECT STAYS HANDLE-VISIBLE. The legend is assembled from an
%   explicit list of representative handles rather than by hiding the repeats,
%   because `HandleVisibility = "off"` also hides an object from `findobj` --
%   so a test counting what was drawn would see one segment and one mark per
%   kind, and would pass while the plot silently drew nothing else.
%
%   Inputs:
%       source        (1x1) struct  a rid.rhc.run result or a run record
%       opts.visible  (1x1) matlab.lang.OnOffSwitchState = "on"
%   Outputs:
%       fig   (1x1) matlab.ui.Figure
%       info  (1x1) struct
%           .runId      (1x1) string
%           .nSteps     (1x1) double
%           .nSegments  (1x1) double
%           .segments   (nSeg x 1) struct  .first .last
%           .nMarks     (1x1) double       one per event, all kinds
%           .marks      (nMarks x 1) struct .step .kind .status .severity .breaks
%           .kinds      (1 x 7) string     the closed set, in section 1's order
%           .status     (1x1) string       the run's termination status

arguments
    source       (1,1) struct
    opts.visible (1,1) matlab.lang.OnOffSwitchState = "on"
end

v            = viewOf(source);
[seg, marks] = jSegments(v.Jtrace, v.eventLog);

fig = newFigure("progress - " + v.runId, opts.visible);
ax  = axes(fig);
hold(ax, "on");

key = struct("handle", {}, "name", {});   % one legend entry per distinct thing

for iSeg = 1:numel(seg)
    kSteps = (seg(iSeg).first:seg(iSeg).last).';
    h = plot(ax, kSteps, v.Jtrace(kSteps), "-o", ...
        "Color", [0.00 0.35 0.62], "MarkerFaceColor", [0.00 0.35 0.62], ...
        "MarkerSize", 4, "LineWidth", 1.4, "Tag", "segment");
    key = remember(key, h, "J (comparable within a segment)");
end

key = drawMarks(ax, key, marks, v.Jtrace);

xlabel(ax, "step");
ylabel(ax, "J");
title(ax, sprintf("progress - %s", v.runId), "Interpreter", "none");
subtitle(ax, sprintf("%d segment(s): J is comparable only WITHIN one. " + ...
    "Terminated %s.", numel(seg), v.status));
grid(ax, "on");
box(ax, "on");
if ~isempty(v.Jtrace)
    xlim(ax, [0.5, numel(v.Jtrace) + 0.5]);
end
if ~isempty(key)
    lg = legend(ax, [key.handle], {key.name}, "Location", "bestoutside");
    lg.Interpreter = "none";
end
hold(ax, "off");

info = struct( ...
    "runId",     v.runId, ...
    "nSteps",    numel(v.Jtrace), ...
    "nSegments", numel(seg), ...
    "segments",  seg, ...
    "nMarks",    numel(marks), ...
    "marks",     marks, ...
    "kinds",     rid.eventKinds(), ...
    "status",    v.status);

end

% ---------------------------------------------------------------------------

function key = drawMarks(ax, key, marks, Jtrace)
%DRAWMARKS  One vertical mark per event, for every kind in the closed set.
%
%   The mark is drawn at the step the event took effect, spanning the axes, so
%   it sits between the value it ends and the value it begins.
if isempty(marks)
    return
end

yTop = 1;
finiteJ = Jtrace(isfinite(Jtrace));
if ~isempty(finiteJ) && max(finiteJ) > 0
    yTop = max(finiteJ);
end

for iMark = 1:numel(marks)
    m = marks(iMark);

    lineStyle = "--";
    width     = 1.0;
    if m.severity == "warning"
        % A statement about the configuration, not a routine event.
        lineStyle = "-";
        width     = 2.0;
    end

    h = plot(ax, [m.step m.step], [0 yTop * 1.05], lineStyle, ...
        "Color", colorFor(m.kind), "LineWidth", width, "Tag", char(m.kind));
    key = remember(key, h, m.kind + labelSuffix(m));
end

end

% ---------------------------------------------------------------------------

function key = remember(key, handle, name)
%REMEMBER  Keep the FIRST handle of each distinct legend entry, and only that.
if ~isempty(key) && any(string({key.name}) == name)
    return
end
key(end + 1) = struct("handle", handle, "name", name);

end

% ---------------------------------------------------------------------------

function suffix = labelSuffix(m)
%LABELSUFFIX  What distinguishes one occurrence of a kind from another.
%
%   `status` and `breaks` are on the label because they are the difference
%   between a logged event and a discontinuity: a rejected fit left the model
%   in place, so the values on either side are still comparable.
suffix = "";
if m.status ~= ""
    suffix = suffix + " (" + m.status + ")";
end
if m.severity == "warning"
    suffix = suffix + " [warning]";
end
if ~m.breaks
    suffix = suffix + " - no baseline break";
end

end

% ---------------------------------------------------------------------------

function color = colorFor(kind)
%COLORFOR  A colour per kind, and a fallback so a NEW kind is still visible.
%
%   An unrecognized kind is drawn in grey, never dropped. `rid.eventKinds` is
%   the set's single owner and is used here only to assign colours, so
%   a kind emitted before that list catches up is still visible.
palette = [0.85 0.33 0.10      % surrogateUpdate
           0.47 0.67 0.19      % resolutionChange
           0.64 0.08 0.18      % coverageReset
           0.30 0.75 0.93      % supportChange
           0.49 0.18 0.56      % normalizationChange
           0.93 0.69 0.13      % modeChange
           0.25 0.25 0.25];    % seed

idx = find(rid.eventKinds() == kind, 1);
if isempty(idx)
    color = [0.40 0.40 0.40];
else
    color = palette(idx, :);
end

end
