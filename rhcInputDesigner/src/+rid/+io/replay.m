function [ok, difference] = replay(rec)
%REPLAY  Re-run a record's algorithm and compare against its stored history.
%
%   A record is only worth keeping if it can be acted on, and this is the act:
%   re-derive the run from `(seed, opts, hist, eventLog)` and check that the
%   reproduced history matches bit-exactly.
%
%   IT NEVER TOUCHES HARDWARE. For an offline run the process is `[]`, as it was.
%   For an ONLINE run the measurements come from the record's own `y` -- and its
%   `x`, when the representation read measured states -- through a
%   playback source, never from a device: `process` is the one thing a replay must
%   not have, and what is being reproduced is the DESIGN, not the plant.
%   The playback is faithful exactly while the reproduced input matches -- and the
%   step at which it stops matching is what `difference` names.
%
%   THE MUTABLE STATE IS REBUILT, NOT REUSED. `setup.state` is a handle the
%   original run advanced, so its cache and its `epsReq` are what the run FINISHED
%   with; the schedulers likewise carry an interval clock that would change WHEN
%   the replay updates. A fresh `CoverageState` is built from the stored
%   `epsReqInitial` -- or, for a criterion with no resolution, the stored one's
%   `freshCopy()` -- and the schedulers are `reset`. Reusing either would produce a
%   replay that diverges for a reason having nothing to do with the code.
%
%   A RECORD WITHOUT `Q` IS DECLARED NON-REPLAYABLE and this stops, rather than
%   regenerating `Q` from its seed (section 7). Regeneration is not a
%   reproducibility mechanism: the generator belongs to a toolbox and may change
%   between releases, so a regenerated `Q` would produce a mismatch that looks
%   exactly like a code regression.
%
%   THE EVENT LOG IS COMPARED, NOT MERELY THE HISTORY. With the blocking scheduler
%   of v1 the update policy is deterministic, so re-running it reproduces the swap
%   instants; comparing the reproduced log against the stored one is what makes
%   that checkable rather than assumed. When the asynchronous scheduler arrives
%   the instants depend on wall-clock fit time and must be RE-APPLIED from the
%   log rather than re-timed -- which needs a scheduler that plays the log
%   back, and this comparison is what will fail the day one is needed and absent.
%
%   THE WALL-CLOCK GUARD IS NOT REPLAYED. `maxWallTime` is ignored and
%   the replay repeats the recorded number of steps; a run stopped by its guard
%   is replayed to that step and keeps the status "guardTripped".
%
%   Inputs:
%       rec         (1x1) struct  a record from rid.rhc.run or rid.io.loadRun
%   Outputs:
%       ok          (1x1) logical  true when the history reproduced exactly
%       difference  (1x1) struct
%           .replayable  (1x1) logical  false when the record lacks Q
%           .reason      (1x1) string   why, when not replayable or not ok
%           .quantity    (1x1) string   "u" | "y" | "x" | "status" | "eventLog" | ""
%           .step        (1x1) double   the FIRST divergence, NaN when none
%           .expected .actual           the two values at that step
%           .result      (1x1) struct   the replayed result, for inspection

arguments
    rec (1,1) struct
end

difference = struct("replayable", false, "reason", "", "quantity", "", ...
                    "step", NaN, "expected", [], "actual", [], ...
                    "result", struct());
ok = false;

if ~isfield(rec, "setup") || ~isfield(rec.setup, "setup")
    error("rid:io:invalidRecord", ...
        "The record carries no setup group, so there is nothing to re-run. " + ...
        "Records written by rid.rhc.run store the Setup object itself, " + ...
        "because a description string cannot be consumed.");
end

% ON THE CONTENTS, NOT ON THE FLAG. The two used to be able to disagree
% -- the record carried `Q` inside the stored `Setup` while `storedQ` said it did
% not -- and refusing on the announcement was therefore stricter than the record
% required, for a reason that was not true. What decides is whether the points
% are there.
if isempty(rec.problem.Q)
    difference.reason = "The record does not store Q, so it is NOT exactly " + ...
        "replayable. Q is not regenerated from " + ...
        "its seed here: the generator belongs to a toolbox and may change " + ...
        "between releases, so a regenerated Q would produce a mismatch that " + ...
        "looks exactly like a code regression. Re-record with opts.io.storeQ.";
    return
end
difference.replayable = true;

old   = rec.setup.setup;
Q     = rec.problem.Q;
% A criterion WITH a resolution is rebuilt from the snapshot the run started
% with, exactly as before the criterion interface, so older records and the canaries replay
% unchanged. One WITHOUT a resolution has no configuration the loop can move
% (`epsReqInitial` is NaN), so a fresh copy of the stored object is the run's
% starting state.
if ~isempty(old.state.resolution())
    state = rid.crit.CoverageState(Q, rec.setup.epsReqInitial, ...
                "psiMax", rec.configuration.effective.crit.psiMax, ...
                "pNorm",  rec.configuration.effective.crit.pNorm);
else
    state = old.state.freshCopy();
end

old.sch.reset();
old.fieldSch.reset();

% The SAME solver, from the record: a replay with a different search does not
% reproduce, which is why `sol` is inside `setup` at all.
% `qInfo` is `old.qInfo` unchanged: a replay reproduces the run, it does not
% re-derive Q's provenance.
setup = rid.rhc.Setup(old.srg, old.rep, old.nrm, old.roi, Q, old.fields, ...
                      old.gen, old.sol, state, old.sch, old.fieldSch, old.qInfo);

% The stored EFFECTIVE options, never defaults plus overrides: a record means what
% it meant when it was written, and merging defaults would let a later change of
% default alter an old run.
overrides    = asOverrides(rec.configuration.effective);
overrides.io.saveRun = false;   % a replay is a check, not a new experiment

% THE CLOCK IS NOT PART OF THE DESIGN. `maxWallTime` stops a run after a
% wall time, so the step it stops at depends on the machine and its load, and a
% replay holding the stored value stops wherever the replay's own clock says.
% The replay repeats the RECORDED number of steps instead: a run stopped by its
% guard is replayed to exactly its `nSteps` through `maxSteps`, the step guard,
% so it ends with the same status. Every other run stops for its own,
% deterministic reason at its `nSteps`; `maxSteps` is capped one step later only
% so that the replay stays finite without the clock, and a replay that reaches
% the cap has diverged and is reported as such.
overrides.rhc.maxWallTime = Inf;
if rec.outputs.status == "guardTripped"
    overrides.rhc.maxSteps = rec.outputs.nSteps;
else
    overrides.rhc.maxSteps = min(overrides.rhc.maxSteps, rec.outputs.nSteps + 1);
end

% THE INITIAL HISTORY IS EVERY ROW BEFORE THE FIRST DESIGNED STEP, not row 1
% alone. Row 1 is the pre-run actuator position and its measurement; a
% run may also start from more -- the NARX warm-up, or a priming signal -- and
% whatever `rid.rhc.run` was handed is the initial condition. The loop
% appends one row per designed step, so the initial condition is the first
% size(u,1) - nSteps rows, and playback starts after them.
nInit = size(rec.outputs.u, 1) - rec.outputs.nSteps;

process = [];
if size(rec.outputs.y, 1) > nInit
    % THE MEASURED STATES TRAVEL WITH IT, split at the same row:
    % `x` is a second measurement of the same applied samples, so it is offset
    % exactly as `y` is. A measured-state run replayed without them would be
    % anchored somewhere else at every step.
    xStored = rec.outputs.x;
    if ~isempty(xStored)
        xStored = xStored(nInit + 1:end, :);
    end
    process = rid.io.PlaybackProcess(rec.outputs.y(nInit + 1:end, :), xStored);
end

% The initial condition is reconstructed from the record's own first rows --
% the pre-run actuator position and its measurement, and any priming
% after it -- rather than from options, because that is where it lives now
%.
hist0 = struct("u", rec.outputs.u(1:nInit, :), "x0", rec.outputs.x0, ...
               "tS", rec.outputs.tS);
if ~isempty(rec.outputs.y)
    hist0.y = rec.outputs.y(1:nInit, :);
end
% The states measured BEFORE the run began. Present only when the
% run measured states: an absent `x` is how a history says it measured none, and
% seeding an empty one would turn a replay of an offline run into a differently
% specified run.
if ~isempty(rec.outputs.x)
    hist0.x = rec.outputs.x(1:nInit, :);
end

result            = rid.rhc.run(process, setup, hist0, overrides);
difference.result = result;

[ok, difference] = compare(rec, result, difference);

end

% ---------------------------------------------------------------------------

function [ok, difference] = compare(rec, result, difference)
%COMPARE  Name the FIRST divergence, not merely report inequality.
%
%   "The histories differ" is useless to someone deciding whether a change was
%   intended; "u diverged at step 12" says where to look.
ok = false;

[same, step] = firstDivergence(rec.outputs.u, result.hist.u);
if ~same
    difference.quantity = "u";
    difference.step     = step;
    difference.expected = rowAt(rec.outputs.u, step);
    difference.actual   = rowAt(result.hist.u, step);
    difference.reason   = sprintf("The designed input diverged at step %d.", step);
    return
end

replayY = zeros(0, 0);
if isfield(result.hist, "y")
    replayY = result.hist.y;
end
[same, step] = firstDivergence(rec.outputs.y, replayY);
if ~same
    difference.quantity = "y";
    difference.step     = step;
    difference.expected = rowAt(rec.outputs.y, step);
    difference.actual   = rowAt(replayY, step);
    difference.reason   = sprintf("The output diverged at step %d.", step);
    return
end

% AFTER `y`, AND FOR THE SAME REASON. The measured states are the
% second half of what an online run observed, and they are what each horizon was
% anchored on -- so a replay that reproduced the input and the output but not the
% states reproduced the signal by coincidence rather than by mechanism.
replayX = zeros(0, 0);
if isfield(result.hist, "x")
    replayX = result.hist.x;
end
[same, step] = firstDivergence(rec.outputs.x, replayX);
if ~same
    difference.quantity = "x";
    difference.step     = step;
    difference.expected = rowAt(rec.outputs.x, step);
    difference.actual   = rowAt(replayX, step);
    difference.reason   = sprintf("The measured state diverged at step %d.", step);
    return
end

if rec.outputs.status ~= result.status
    difference.quantity = "status";
    difference.reason   = "The run stopped for a different reason: """ + ...
        rec.outputs.status + """ then, """ + result.status + """ now.";
    difference.expected = rec.outputs.status;
    difference.actual   = result.status;
    return
end

% The instants an update took effect at, which invariant 5 is about.
[same, step] = compareEvents(rec.outputs.eventLog, result.eventLog);
if ~same
    difference.quantity = "eventLog";
    difference.step     = step;
    difference.reason   = sprintf("The event log diverged at entry %d: an " + ...
        "update took effect at a different step, or a different kind of " + ...
        "update happened.", step);
    return
end

ok = true;

end

% ---------------------------------------------------------------------------

function [same, step] = firstDivergence(expected, actual)
%FIRSTDIVERGENCE  Bit-for-bit, with the first differing row named.
%
%   `isequal` on doubles is EXACT, which is what section 4 requires: a replay that
%   agreed only to a tolerance would hide precisely the reorderings a canary
%   exists to catch.
same = true;
step = NaN;

if size(expected, 2) ~= size(actual, 2)
    same = false;
    step = 1;
    return
end
nCommon = min(size(expected, 1), size(actual, 1));
for iRow = 1:nCommon
    if ~isequal(expected(iRow, :), actual(iRow, :))
        same = false;
        step = iRow;
        return
    end
end
if size(expected, 1) ~= size(actual, 1)
    same = false;
    step = nCommon + 1;
end

end

% ---------------------------------------------------------------------------

function [same, index] = compareEvents(expected, actual)
%COMPAREEVENTS  Same kinds, at the same steps, in the same order.
%
%   The DETAIL is not compared: a scheduler's message may legitimately carry a
%   file path or a timing. What must reproduce is WHEN each change took effect and
%   WHAT it was.
same  = true;
index = NaN;

if numel(expected) ~= numel(actual)
    same  = false;
    index = min(numel(expected), numel(actual)) + 1;
    return
end
for iEvent = 1:numel(expected)
    if expected(iEvent).step ~= actual(iEvent).step || ...
            expected(iEvent).kind ~= actual(iEvent).kind
        same  = false;
        index = iEvent;
        return
    end
end

end

% ---------------------------------------------------------------------------

function row = rowAt(matrix, step)
%ROWAT  One row, or an explicit absence when the matrix is shorter.
row = [];
if ~isnan(step) && step <= size(matrix, 1)
    row = matrix(step, :);
end

end

% ---------------------------------------------------------------------------

function overrides = asOverrides(effective)
%ASOVERRIDES  The effective options, in the shape rid.rhc.run accepts.
%
%   `run` resolves overrides onto the current defaults, so handing it the whole
%   effective struct pins every value the original run used -- including the ones
%   that were defaults then and may not be now, which is the entire reason the
%   record stores the EFFECTIVE struct rather than the overrides alone.
overrides = effective;

end
