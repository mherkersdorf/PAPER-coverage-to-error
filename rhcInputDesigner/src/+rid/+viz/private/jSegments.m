function [seg, marks] = jSegments(Jtrace, eventLog)
%JSEGMENTS  Where the J trace is comparable with itself, and where it is not.
%
%   J IS MONOTONE ONLY WITHIN A SEGMENT. `rid.rhc.run` sets its monotonicity
%   baseline to Inf at three places -- an APPLIED surrogate update, an APPLIED
%   resolution change, and a coverage reset -- and the assert never spans
%   one, because the values on either side measure different things: a rebuilt
%   cache, a different epsReq, or a cache deliberately emptied. Drawing the
%   trace as one continuous line would imply a comparability that does not
%   hold, and a reader would see a rise where the loop saw none.
%
%   THE BOUNDARY LIES AFTER THE STEP THAT LOGGED THE EVENT. Within step k the
%   loop stores JAfter into `Jtrace(k)` and only then polls the schedulers and
%   applies the reset, so an event at step k separates `Jtrace(k)` from
%   `Jtrace(k+1)`.
%
%   `.status` IS WHAT DECIDES, NOT THE KIND. A `"surrogateUpdate"` is logged
%   whenever the scheduler had something to say -- `"applied"`, `"rejected"` or
%   `"unsupported"` (rid.surrogate.Scheduler) -- and only `"applied"` replaced
%   the model. Breaking the trace on a rejected fit would invent a
%   discontinuity the run does not have. `"coverageReset"` carries no status
%   and always breaks.
%
%   EVERY EVENT IS MARKED, whatever its kind. The closed set has seven members
%   and only three are emitted today; a mark
%   that appears only for the kinds this file happens to know would make a new
%   kind invisible, and a jump in J with no mark beside it is exactly the
%   puzzling plot months later that the canaries exist to prevent. Marking is
%   therefore generic and segmenting is the narrow, named exception.
%
%   Inputs:
%       Jtrace    (nSteps x 1) double
%       eventLog  (nEvents x 1) struct  .step .kind .detail
%   Outputs:
%       seg    (nSeg x 1) struct   .first .last (1x1) step indices, inclusive
%       marks  (nEvents x 1) struct
%           .step      (1x1) double
%           .kind      (1x1) string
%           .status    (1x1) string   "" when the detail carries none
%           .severity  (1x1) string   "info" | "warning" | ""
%           .breaks    (1x1) logical  did this event break the J baseline?

nSteps = numel(Jtrace);

marks = struct("step", {}, "kind", {}, "status", {}, "severity", {}, ...
               "breaks", {});
marks = reshape(marks, 0, 1);
for iEvent = 1:numel(eventLog)
    e      = eventLog(iEvent);
    status = detailField(e.detail, "status");
    marks(end + 1, 1) = struct( ...
        "step",     e.step, ...
        "kind",     string(e.kind), ...
        "status",   status, ...
        "severity", detailField(e.detail, "severity"), ...
        "breaks",   breaksBaseline(string(e.kind), status)); %#ok<AGROW>
end

if nSteps == 0
    seg = reshape(struct("first", {}, "last", {}), 0, 1);
    return
end

isBreakAfter = false(nSteps, 1);
for iMark = 1:numel(marks)
    step = marks(iMark).step;
    if marks(iMark).breaks && step >= 1 && step <= nSteps
        isBreakAfter(step) = true;
    end
end

seg   = reshape(struct("first", {}, "last", {}), 0, 1);
first = 1;
for kStep = 1:nSteps
    if isBreakAfter(kStep) || kStep == nSteps
        seg(end + 1, 1) = struct("first", first, "last", kStep); %#ok<AGROW>
        first = kStep + 1;
    end
end

end

% ---------------------------------------------------------------------------

function tf = breaksBaseline(kind, status)
%BREAKSBASELINE  Does this event make the values on either side incomparable?
%
%   The three places rid.rhc.run resets JBaseline to Inf, and only those.
switch kind
    case "coverageReset"
        tf = true;
    case {"surrogateUpdate", "resolutionChange"}
        tf = status == "applied";
    otherwise
        % supportChange, normalizationChange, modeChange, seed: marked, but
        % not treated as a discontinuity, because the loop does not reset its
        % baseline on them. When one of them starts being emitted, this is the
        % line that must be revisited -- deliberately, not by default.
        tf = false;
end

end

% ---------------------------------------------------------------------------

function value = detailField(detail, name)
%DETAILFIELD  One optional field of a kind-specific detail, as a string.
value = "";
if isstruct(detail) && isscalar(detail) && isfield(detail, name)
    value = string(detail.(name));
end

end
