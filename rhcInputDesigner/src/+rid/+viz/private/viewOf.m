function v = viewOf(source)
%VIEWOF  One plain view of a run, from a result OR from a run record.
%
%   Every +viz entry point takes either. The difference between them is one
%   field -- a result CARRIES its record (rid.rhc.run, `result.record`) -- so
%   this reduces both to the record and reads everything from there. Two code
%   paths reading two field layouts would be two chances to draw one of them
%   wrongly.
%
%   IT READS ONLY DECLARED ACCESSORS. `Q` comes from `problem.Q`, which is what
%   `rid.io.replay` refuses on. Under
%   `storeQ = false` there is nowhere else to read it from: the storeQ rule removed `Q`
%   from the stored `Setup` too, `setup.Q` and `setup.state` alike, so the flag
%   and the contents can no longer disagree. This reads the declared place
%   regardless, because that is the one the specification names.
%
%   THE COVERAGE STATE IS THE RUN'S FINAL, HONEST ONE. `setup.setup.state` is
%   the handle the loop advanced, and rid.rhc.run rebuilds it over the full
%   history at run end when a reset occurred, precisely so that the
%   stored object is not the deliberate fiction. Its `dRaw`, `epsReq`,
%   `psiMax` and `isActive` are public reads.
%
%   Inputs:
%       source  (1x1) struct  a rid.rhc.run result or a run record
%   Outputs:
%       v  (1x1) struct
%           .runId       (1x1) string
%           .nDim .nSup  (1x1) double
%           .dimName     (1 x nDim) string   coordinate names, coordinate order
%           .hasQ        (1x1) logical
%           .Q           (nSup x nDim) normalized, empty when ~hasQ
%           .state       (1x1) rid.crit.CoverageState
%           .Jtrace      (nSteps x 1)
%           .eventLog    (nEvents x 1) struct .step .kind .detail
%           .steps       (nSteps x 1) struct
%           .hist        (1x1) struct .u .x0 .tS, plus .y and .x when the run
%                                              measured them
%           .rep .nrm .gen                     the stored objects
%           .srg         (1x1) Surrogate       the run's FINAL surrogate, from
%                                              `outputs` -- NOT `setup.srg`,
%                                              which is the initial one
%           .status .mode (1x1) string

arguments
    source (1,1) struct
end

rec = recordOf(source);

setup = rec.setup.setup;

v.runId   = string(rec.identity.runId);
v.state   = setup.state;
v.nSup    = setup.state.nSup;
v.nDim    = setup.state.nDim;
v.dimName = setup.nrm.dimName;

% Section 7, through the declared accessor. `storedQ` false is a refusal, not a
% cue to reconstruct: the generating algorithm belongs to a toolbox and may
% change between releases, so regenerated points would be a picture of a design
% the run never used.
v.hasQ = ~isempty(rec.problem.Q);
if v.hasQ
    v.Q = rec.problem.Q;
else
    v.Q = zeros(0, v.nDim);
end

v.Jtrace   = rec.outputs.Jtrace;
v.eventLog = rec.outputs.eventLog;
v.steps    = rec.outputs.steps;
v.status   = string(rec.outputs.status);
v.mode     = string(rec.outputs.certification.mode);

% THE HISTORY IS REASSEMBLED, NOT SUMMARIZED, because `realizedPoints` calls
% `rep.build` on it: a "measured" representation reads its past coordinates from
% the history itself, so a field left out here is not a plot drawn slightly
% wrongly -- it is a run that cannot be plotted at all.
v.hist = struct("u",  rec.outputs.u,  "y", rec.outputs.y, ...
                "x",  rec.outputs.x, ...
                "x0", rec.outputs.x0, "tS", rec.outputs.tS);
if isempty(v.hist.y)
    v.hist = rmfield(v.hist, "y");   % an absent y is how a history says it has none
end
if isempty(v.hist.x)
    v.hist = rmfield(v.hist, "x");   % and an absent x, that it measured no states
end

v.rep = setup.rep;
% THE FINAL SURROGATE, NOT THE ONE THE RUN STARTED WITH. `setup.srg` is
% immutable and holds the initial model; rebuilding coordinates with it after a
% refit would show points the run never visited.
v.srg = rec.outputs.srg;
v.nrm = setup.nrm;
v.gen = setup.gen;

end

% ---------------------------------------------------------------------------

function rec = recordOf(source)
%RECORDOF  Reduce a result or a record to the record, or say which it is not.
if isfield(source, "record") && isstruct(source.record) && ...
        isfield(source.record, "outputs")
    rec = source.record;      % a result carries its record
    return
end
if isfield(source, "outputs") && isfield(source, "setup") && ...
        isfield(source, "problem") && isfield(source, "identity")
    rec = source;
    return
end
error("rid:viz:unknownSource", ...
    "The source is neither a rid.rhc.run result (it has no .record with " + ...
    ".outputs) nor a run record (it lacks the identity/problem/setup/outputs " + ...
    "groups of a run record).");

end
