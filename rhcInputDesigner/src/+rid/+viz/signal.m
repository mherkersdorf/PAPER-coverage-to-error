function [fig, info] = signal(source, opts)
%SIGNAL  The designed input trajectory, its measured output, and its bounds.
%
%       [fig, info] = rid.viz.signal(result)
%       [fig, info] = rid.viz.signal(rid.io.loadRun(runId), "visible", "off")
%
%   THE BOUNDS ARE DRAWN, so a saturated run is visible AS SATURATION rather
%   than as a flat line. An input pinned at its amplitude limit and an input
%   the design simply stopped moving produce the same picture without them, and
%   they call for opposite responses: the first says the box is the binding
%   constraint, the second says the criterion has no gradient left. The rate
%   box is drawn for the same reason -- a rate-limited signal is flat between
%   samples for a reason nothing on the amplitude axes reveals. Samples sitting
%   ON a bound are marked, and counted in `info`.
%
%   THE INITIAL CONDITION IS SEPARATED FROM THE DESIGN. `hist.u` opens with the
%   pre-run actuator position and whatever priming samples the experiment
%   started with; only the rows the loop committed are the designed
%   signal. The split is read from the per-step diagnostics -- `steps.uApplied`
%   is what was committed, step by step -- and never guessed from a
%   horizon length, because only the APPLIED UNIT is committed and the rest of
%   each horizon was discarded.
%
%   THE BOX IS ALWAYS THERE. `constraints` is on the
%   `rid.gen.Generator` contract, so a generator that cannot state its actuator
%   limits fails at construction and never reaches a figure. This file asks no
%   question about whether bounds exist; it draws them.
%
%   THE MEASURED OUTPUT IS DRAWN WHEN THERE IS ONE. An absent `y` is how a
%   history says it has none (an offline run has nothing to measure), so its
%   absence is stated on the figure, not filled in from the surrogate.
%
%   Inputs:
%       source         (1x1) struct  a rid.rhc.run result or a run record
%       opts.visible   (1x1) matlab.lang.OnOffSwitchState = "on"
%       opts.boundTol  (1x1) double = 1e-6  a sample counts as ON a bound when
%                                   it is within this FRACTION of the channel's
%                                   box width of it. A display threshold, not an
%                                   algorithm parameter, which is why it lives
%                                   here and not in rid.defaultOptions.
%   Outputs:
%       fig   (1x1) matlab.ui.Figure
%       info  (1x1) struct
%           .runId        (1x1) string
%           .nSamples     (1x1) double   rows of hist.u
%           .nDesigned    (1x1) double   rows the loop committed
%           .nPriming     (1x1) double   nSamples - nDesigned
%           .nInputs      (1x1) double
%           .hasOutput    (1x1) logical
%           .nAtAmplitude (1 x nInputs) double  designed samples on an amplitude bound
%           .nAtRate      (1 x nInputs) double  designed steps on a rate bound
%           .nAxes        (1x1) double
%           .missing      (1 x n) string

arguments
    source        (1,1) struct
    opts.visible  (1,1) matlab.lang.OnOffSwitchState = "on"
    opts.boundTol (1,1) double {mustBeReal, mustBePositive} = 1e-6
end

v = viewOf(source);

u        = v.hist.u;
nSamples = size(u, 1);
nInputs  = size(u, 2);
tS       = v.hist.tS;
t        = (0:nSamples - 1).' * tS;

nDesigned = designedCount(v.steps);
firstDesigned = nSamples - nDesigned + 1;

b       = appliedBounds(v.gen, nInputs);
missing = strings(1, 0);

hasOutput = isfield(v.hist, "y") && ~isempty(v.hist.y);
if ~hasOutput
    missing(end + 1) = "no measured output: this run declared no process, so " + ...
                       "there is nothing measured to draw.";
end

du = diff(u, 1, 1);   % (nSamples-1 x nInputs), aligned with t(2:end)

fig    = newFigure("signal - " + v.runId, opts.visible);
nAxes  = 2 + double(hasOutput);
tl     = tiledlayout(fig, nAxes, 1, "TileSpacing", "compact", "Padding", "compact");

axU = nexttile(tl);
nAtAmplitude = drawChannels(axU, t, u, firstDesigned, b, "amplitude", ...
                            opts.boundTol, "u");
ylabel(axU, "u");

axD = nexttile(tl);
nAtRate = drawChannels(axD, t(2:end), du, max(firstDesigned - 1, 1), b, ...
                       "rate", opts.boundTol, "du");
ylabel(axD, "du");

if hasOutput
    axY = nexttile(tl);
    plot(axY, t, v.hist.y, "-", "LineWidth", 1.3, ...
        "Color", [0.20 0.29 0.55], "Tag", "y");
    ylabel(axY, "y");
    grid(axY, "on");
    box(axY, "on");
    xlabel(axY, "time [s]");
else
    xlabel(axD, "time [s]");
end

% ONE TIME AXIS FOR ALL PANELS. Stacked panels whose x-limits differ invite
% exactly the misreading they exist to prevent: the same horizontal position
% meaning two different instants, so an input excursion appears to line up
% with the wrong measurement.
linkaxes(findobj(fig, "Type", "axes"), "x");

title(tl, sprintf("signal - %s", v.runId), "Interpreter", "none");
subtitle(tl, sprintf("%d designed sample(s) after %d of initial condition; " + ...
    "%d on an amplitude bound, %d on a rate bound.", ...
    nDesigned, nSamples - nDesigned, sum(nAtAmplitude), sum(nAtRate)));

if ~isempty(missing)
    % On the figure as well as in the struct: a figure saved to disk outlives
    % the struct and must still say what it is not showing.
    subtitle(axU, strjoin(missing, " "));
end

info = struct( ...
    "runId",        v.runId, ...
    "nSamples",     nSamples, ...
    "nDesigned",    nDesigned, ...
    "nPriming",     nSamples - nDesigned, ...
    "nInputs",      nInputs, ...
    "hasOutput",    hasOutput, ...
    "nAtAmplitude", nAtAmplitude, ...
    "nAtRate",      nAtRate, ...
    "nAxes",        nAxes, ...
    "missing",      missing);

end

% ---------------------------------------------------------------------------

function nAtBound = drawChannels(ax, t, X, firstDesigned, b, which, tol, tag)
%DRAWCHANNELS  One stairs per channel, its box, and the samples sitting on it.
%
%   STAIRS, NOT PLOT: the input is held between samples, and a straight line
%   between two commits draws a signal the actuator never produced.
nChannels = size(X, 2);
nAtBound  = zeros(1, nChannels);
hold(ax, "on");

for iCh = 1:nChannels
    stairs(ax, t, X(:, iCh), "-", "LineWidth", 1.3, ...
        "Color", channelColor(iCh), "Tag", char(tag));

    box2 = b.(which)(iCh, :);
    for iSide = 1:2
        yline(ax, box2(iSide), "--", "Color", [0.64 0.08 0.18], ...
            "LineWidth", 1.0, "Tag", char(tag + "Bound"), ...
            "Label", sprintf("%s %s", b.channelName(iCh), ...
                             sideName(which, iSide)), ...
            "LabelHorizontalAlignment", "left", "Interpreter", "none");
    end

    % ON a bound, within a fraction of the channel's own box width: an absolute
    % threshold would mean something different on every channel.
    width  = box2(2) - box2(1);
    atLow  = abs(X(:, iCh) - box2(1)) <= tol * width;
    atHigh = abs(X(:, iCh) - box2(2)) <= tol * width;
    at     = atLow | atHigh;
    at(1:min(firstDesigned - 1, numel(at))) = false;   % the design's doing only

    nAtBound(iCh) = sum(at);
    if any(at)
        plot(ax, t(at), X(at, iCh), "o", "MarkerSize", 7, ...
            "MarkerEdgeColor", [0.64 0.08 0.18], "LineWidth", 1.2, ...
            "Tag", char(tag + "Saturated"));
    end
end

if firstDesigned > 1 && firstDesigned <= numel(t)
    % Where the initial condition ends and the design begins.
    xline(ax, t(firstDesigned), "-", "Color", [0.45 0.45 0.45], ...
        "LineWidth", 1.2, "Label", "design begins", ...
        "LabelVerticalAlignment", "bottom", "Tag", "designStart");
end

grid(ax, "on");
box(ax, "on");
hold(ax, "off");

end

% ---------------------------------------------------------------------------

function n = designedCount(steps)
%DESIGNEDCOUNT  How many samples the loop committed, from the per-step log.
%
%   Only the APPLIED UNIT of each winning candidate is committed, so
%   this counts what `steps.uApplied` actually holds rather than multiplying a
%   horizon by a step count.
n = 0;
for kStep = 1:numel(steps)
    n = n + size(steps(kStep).uApplied, 1);
end

end

% ---------------------------------------------------------------------------

function name = sideName(which, iSide)
%SIDENAME  The label a bound carries.
if which == "amplitude"
    names = ["min", "max"];
else
    names = ["rate min", "rate max"];
end
name = names(iSide);

end

% ---------------------------------------------------------------------------

function c = channelColor(iCh)
%CHANNELCOLOR  A colour per input channel, cycling for many-input problems.
palette = [0.00 0.35 0.62
           0.85 0.33 0.10
           0.47 0.67 0.19
           0.49 0.18 0.56];
c = palette(mod(iCh - 1, size(palette, 1)) + 1, :);

end
