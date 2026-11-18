function generateCanary(name)
%GENERATECANARY  Produce a canary file: a declarative spec plus expected outputs.
%
%   RE-RUNNING THIS IS A RE-RECORDING. The method description and canary replay
%   permit it only in a commit that touches the canary files alone, with a
%   one-line reason, and only after deciding that the behaviour change was
%   INTENDED. A canary regenerated reflexively when it goes red is a mirror of the
%   code and can never catch anything again.
%
%   WHY IT EXISTS AT ALL. The golden-file rule was written after a golden layer sat inert for its
%   whole existence because nothing about how it had been made was recorded. The
%   generator beside the artefact is what makes its provenance
%   checkable rather than asserted.
%
%   IT USES `lhsdesign` ONCE, HERE, AT FREEZE TIME, and writes the resulting `Q` as
%   explicit numbers. Nothing in the replay path regenerates it: the canary specification forbids a
%   canary that is in part a test of the Statistics Toolbox's generator.
%
%   Every double is written at `%.17g`, which round-trips exactly, and the
%   generator verifies that by decoding what it wrote and rebuilding the run
%   before it returns.
%
%   TWO CANARIES, TWO ACCEPTANCE CRITERIA. `sisoNarxBudget` is the VALUE
%   canary: `rid.opt.GridSearch`, no toolbox anywhere in its replay path, every
%   number compared at `1e-12`. `sisoNarxMultiStart` is the PATH canary: the same
%   configuration with `rid.opt.MultiStart`, so `fmincon` IS in the replay path,
%   and what is frozen is the affine tier, the gradient path, monotonicity within
%   a segment, the reset count, evaluation counts to an order of magnitude, and
%   the final `J` within a declared band -- deliberately not `u`, `theta` or `J`.
%
%   Usage (from the project root, after `setupPath`):
%       rid.io.projectRoot();  % assert we are in the tree
%       generateCanary("sisoNarxBudget")
%       generateCanary("sisoNarxMultiStart")
%
%   Inputs:
%       name  (1x1) string  the canary's basename in tests/data/canaries/
%   Outputs: none; writes <name>.json beside this file.

arguments
    name (1,1) string {mustBeMember(name, ["sisoNarxBudget", "sisoNarxMultiStart"])} ...
        = "sisoNarxBudget"
end

here = string(fileparts(mfilename("fullpath")));

% TWO CANARIES, TWO ACCEPTANCE CRITERIA. The value canary is built by
% `valueCanaryNode`, the path canary by `pathCanaryNode`; the serializer below is
% shared, so both files are written by the same %.17g round-tripping encoder and a
% change to the format cannot reach one artefact without reaching the other.
switch name
    case "sisoNarxBudget"
        node = valueCanaryNode(name);
    case "sisoNarxMultiStart"
        node = pathCanaryNode(name);
end

path = fullfile(here, name + ".json");
fid  = fopen(path, "w");
if fid < 0
    error("rid:io:canaryWriteFailed", "Cannot write %s.", path);
end
closer = onCleanup(@() fclose(fid));
fprintf(fid, "%s\n", encode(node, 0));
clear closer

fprintf("wrote %s\n", path);

end

% --- the VALUE canary: rid.opt.GridSearch, every number at 1e-12 -----------

function node = valueCanaryNode(name)
%VALUECANARYNODE  The spec and expected outputs of `sisoNarxBudget`.

% --- the configuration, and why it is this one ------------------------------
% Twelve supporting points, epsReq0 = 0.30, a four-sample horizon and a six-step
% budget. Checked against the coverage reset: the run ends at J = 0.24 with an EMPTY event log,
% so it never certifies, never resets, and is not the degenerate regime in which a
% canary is a constant signal that passes forever.
nSup     = 12;
epsReq0  = 0.30;
nHorizon = 4;
nSteps   = 6;
seedQ    = 0;      % used HERE only; Q is frozen as numbers below

roi     = struct("u1_k", [-1 1], "du1_k", [-0.5 0.5], "y_km1", [-2 2]);
rep     = rid.rep.Narx(2, 1, "pastSource", "simulated");
[box, dimName] = rep.bounds(roi);
nrm     = rid.roi.Normalization(box, dimName);
srg     = rid.surrogate.LinearFirstOrder(0.5, 0.1, 1);
Q       = rid.roi.supportingPoints(nrm, nSup, "seed", seedQ);
fields  = rid.fields.nominalFields(nSup);
epsReq  = rid.crit.demandedResolution(fields, "epsReq0", epsReq0, "mode", "budget", ...
              "fillDistanceQ", rid.roi.estimateFillDistance(Q, nrm));

% THE SOLVER IS THE EXHAUSTIVE ORACLE, NOT THE DEFAULT SEARCH (the single-setup rule amendment). rid.opt.MultiStart would put `fmincon` in the
% replay path, and its iterates are the one thing in this project a MATLAB release
% may legitimately move -- which is the failure the canary specification names about `lhsdesign`, one
% size larger. rid.opt.GridSearch decides by exhaustive comparison over a fixed
% grid: no seed and no toolbox chooses the committed input.
sol = rid.opt.GridSearch();

% qInfo IS "supplied": this Q is frozen as explicit numbers below and
% no generator runs at replay, so nothing about its provenance may be invented.
qInfo = struct("origin", "supplied", "method", "", "seed", NaN, "hash", "");

setup = rid.rhc.Setup(srg, rep, nrm, roi, Q, fields, ...
                      rid.gen.Direct(rid.gen.Constraints(rep, roi), nHorizon), ...
                      sol, rid.crit.CoverageState(Q, epsReq), ...
                      rid.surrogate.ImmediateScheduler(), ...
                      rid.fields.ImmediateScheduler(), qInfo);

hist0     = struct("u", 0, "x0", srg.steadyState(0), "tS", srg.tS);
overrides = struct("crit", struct("epsReq0", epsReq0, "mode", "budget"), ...
                   "rhc",  struct("nSteps", nSteps), ...
                   "io",   struct("saveRun", false));

r = rid.rhc.run([], setup, hist0, overrides);

% --- the file ---------------------------------------------------------------
node = struct();
node.canary      = name;
node.description = "SISO first-order surrogate in a three-dimensional NARX " + ...
    "regressor space; exhaustive grid search; six-step budget run.";
node.tolerance   = jStruct(struct("relative", 1e-12));

node.data = jStruct(struct( ...
    "Q",             jMatrix(Q), ...
    "fields",        jStruct(struct("relevance",  jColumn(fields.relevance), ...
                                    "complexity", jColumn(fields.complexity), ...
                                    "noiseLevel", jColumn(fields.noiseLevel))), ...
    "epsReqInitial", jColumn(epsReq), ...
    "hist0",         jStruct(struct("u",  jMatrix(hist0.u), ...
                                    "x0", jMatrix(hist0.x0), ...
                                    "tS", hist0.tS))));

% The twelve members, each with its CLASS and its
% CONSTRUCTOR ARGUMENTS. They are listed in the order the harness builds them,
% so that a {"ref": ...} never points forward.
node.setup = jStruct(struct( ...
    "roi",      jStruct(struct("u1_k",  jRow(roi.u1_k), ...
                               "du1_k", jRow(roi.du1_k), ...
                               "y_km1", jRow(roi.y_km1))), ...
    "rep",      jClass("rid.rep.Narx", {2, 1, "pastSource", "simulated"}), ...
    "nrm",      jClass("rid.roi.Normalization", {jMatrix(box), jStrings(dimName)}), ...
    "srg",      jClass("rid.surrogate.LinearFirstOrder", {0.5, 0.1, 1}), ...
    "fields",   jRef("fields"), ...
    "Q",        jRef("Q"), ...
    "gen",      jClass("rid.gen.Direct", ...
                    {jClass("rid.gen.Constraints", {jRef("rep"), jRef("roi")}), nHorizon}), ...
    "sol",      jClass("rid.opt.GridSearch", {}), ...
    "state",    jClass("rid.crit.CoverageState", {jRef("Q"), jRef("epsReqInitial")}), ...
    "sch",      jClass("rid.surrogate.ImmediateScheduler", {}), ...
    "fieldSch", jClass("rid.fields.ImmediateScheduler", {}), ...
    "qInfo",    jStruct(struct("origin", "supplied", "method", "", ...
                               "seed", NaN, "hash", ""))));

node.options = jStruct(struct( ...
    "crit", jStruct(struct("epsReq0", epsReq0, "mode", "budget")), ...
    "rhc",  jStruct(struct("nSteps", nSteps)), ...
    "io",   jStruct(struct("saveRun", false))));

node.expected = jStruct(struct( ...
    "status",  r.status, ...
    "nSteps",  r.nSteps, ...
    "J",       r.J, ...
    "Jtrace",  jColumn(r.Jtrace), ...
    "u",       jMatrix(r.hist.u), ...
    "steps",   jStruct(struct( ...
        "JBefore",       jColumn([r.steps.JBefore].'), ...
        "JAfter",        jColumn([r.steps.JAfter].'), ...
        "nActive",       jColumn([r.steps.nActive].'), ...
        "uApplied",      jMatrix(vertcat(r.steps.uApplied)), ...
        "theta",         jMatrix(horzcat(r.steps.theta).'), ...
        "coverageReset", jColumn(double([r.steps.coverageReset].')))), ...
    "certification", jStruct(struct( ...
        "isSatisfied",   double(r.certification.isSatisfied), ...
        "epsReqMax",     r.certification.epsReqMax, ...
        "fillDistanceQ", r.certification.fillDistanceQ, ...
        "bound",         r.certification.bound, ...
        "mode",          r.certification.mode)), ...
    "nEvents",         numel(r.eventLog), ...
    "nCoverageResets", sum(string({r.eventLog.kind}) == "coverageReset"), ...
    "nWarnings",       numel(r.warnings)));


end


% --- the PATH canary: rid.opt.MultiStart, discrete facts about the path ----

function node = pathCanaryNode(name)
%PATHCANARYNODE  The spec and expected outputs of `sisoNarxMultiStart`.
%
%   THE CONFIGURATION IS THE VALUE CANARY'S, WITH ONE MEMBER CHANGED. Everything
%   above -- ROI, representation, surrogate, Q, fields, horizon, budget, initial
%   history -- is identical; `sol` is `rid.opt.MultiStart` instead of
%   `rid.opt.GridSearch`. The two canaries therefore differ in exactly one member
%   of the single-setup rule's eleven, so what this one adds is exactly the optimizer and the
%   analytic-gradient path the other deliberately leaves out.
%
%   WHAT IS FROZEN IS THE PATH THE RUN TOOK, NOT THE VALUES IT REACHED.
%   `J` is a minimum over samples and is NOT convex, which is why
%   multi-start exists at all; a changed `fmincon` iterate path across a MATLAB
%   release therefore does not perturb the answer slightly, it can land in a
%   DIFFERENT LOCAL OPTIMUM, and the committed input jumps discontinuously. No
%   tolerance covers a basin switch, so `u`, `theta` and `J` are not stored here
%   and adding them "for strength" is the failure the path canary predicts, not an
%   improvement.
%
%   Measured, not assumed: over ten legitimate perturbations of the iterate path
%   (start seeds 0-9) the committed input changed completely at every step while
%   the affine tier, the gradient path, monotonicity and the reset count did not
%   move at all.

% --- the configuration -----------------------------------------------------
% Checked against the coverage reset exactly as the value canary is: the run ends at
% J = 0.2555 with six of twelve supporting points still unresolved and an EMPTY
% event log, so it never certifies, never resets, and is not the degenerate
% regime in which a canary is a constant signal that passes forever.
nSup     = 12;
epsReq0  = 0.30;
nHorizon = 4;
nSteps   = 6;
seedQ    = 0;      % used HERE only; Q is frozen as numbers below

roi     = struct("u1_k", [-1 1], "du1_k", [-0.5 0.5], "y_km1", [-2 2]);
rep     = rid.rep.Narx(2, 1, "pastSource", "simulated");
[box, dimName] = rep.bounds(roi);
nrm     = rid.roi.Normalization(box, dimName);
srg     = rid.surrogate.LinearFirstOrder(0.5, 0.1, 1);
Q       = rid.roi.supportingPoints(nrm, nSup, "seed", seedQ);
fields  = rid.fields.nominalFields(nSup);
epsReq  = rid.crit.demandedResolution(fields, "epsReq0", epsReq0, "mode", "budget", ...
              "fillDistanceQ", rid.roi.estimateFillDistance(Q, nrm));

% THE SOLVER IS THE DEFAULT SEARCH, KNOWINGLY PUTTING `fmincon` IN THE REPLAY
% PATH. That is precisely why the acceptance criterion below is a set of discrete
% facts rather than values: the toolbox is admitted, and the criterion is chosen
% so that what the toolbox may legitimately move is not what the canary reads.
sol = rid.opt.MultiStart();

% qInfo IS "supplied": this Q is frozen as explicit numbers below and
% no generator runs at replay, so nothing about its provenance may be invented.
qInfo = struct("origin", "supplied", "method", "", "seed", NaN, "hash", "");

setup = rid.rhc.Setup(srg, rep, nrm, roi, Q, fields, ...
                      rid.gen.Direct(rid.gen.Constraints(rep, roi), nHorizon), ...
                      sol, rid.crit.CoverageState(Q, epsReq), ...
                      rid.surrogate.ImmediateScheduler(), ...
                      rid.fields.ImmediateScheduler(), qInfo);

hist0     = struct("u", 0, "x0", srg.steadyState(0), "tS", srg.tS);
overrides = struct("crit", struct("epsReq0", epsReq0, "mode", "budget"), ...
                   "rhc",  struct("nSteps", nSteps), ...
                   "io",   struct("saveRun", false));

r = rid.rhc.run([], setup, hist0, overrides);

% --- the file ---------------------------------------------------------------
node = struct();
node.canary      = name;
node.criterion   = "path";
node.description = "SISO first-order surrogate in a three-dimensional NARX " + ...
    "regressor space; multi-start SQP on the analytic gradient; six-step " + ...
    "budget run. Freezes the PATH, not the values.";

% THE ACCEPTANCE CRITERION, DECLARED IN THE ARTEFACT. A path canary has no single
% tolerance, so the file carries the three numbers its criterion needs instead of
% the value canary's `tolerance.relative`.
%
% `JFinalBand` IS CHOSEN, NOT DERIVED.
% `evaluationDecadeSlack` is the order of magnitude of the path canary read as one decade of
% slack either side of the stored count, because the measured per-step spread
% under a changed iterate path is a factor of ~29 -- a ratio band of ten would
% flap, and a bare decade bucket would sit a factor of two from its own boundary.
% `evaluationTotalFactor` is the sharper aggregate check the per-step slack cannot
% be: the measured spread of the TOTAL is only a factor of ~2.2.
node.acceptance = jStruct(struct( ...
    "JFinalBand",            jRow([0.10, 0.60]), ...
    "evaluationDecadeSlack", 1, ...
    "evaluationTotalFactor", 10));

node.data = jStruct(struct( ...
    "Q",             jMatrix(Q), ...
    "fields",        jStruct(struct("relevance",  jColumn(fields.relevance), ...
                                    "complexity", jColumn(fields.complexity), ...
                                    "noiseLevel", jColumn(fields.noiseLevel))), ...
    "epsReqInitial", jColumn(epsReq), ...
    "hist0",         jStruct(struct("u",  jMatrix(hist0.u), ...
                                    "x0", jMatrix(hist0.x0), ...
                                    "tS", hist0.tS))));

% The twelve members, each with its CLASS and its
% CONSTRUCTOR ARGUMENTS, in the order the harness builds them so that a
% {"ref": ...} never points forward.
node.setup = jStruct(struct( ...
    "roi",      jStruct(struct("u1_k",  jRow(roi.u1_k), ...
                               "du1_k", jRow(roi.du1_k), ...
                               "y_km1", jRow(roi.y_km1))), ...
    "rep",      jClass("rid.rep.Narx", {2, 1, "pastSource", "simulated"}), ...
    "nrm",      jClass("rid.roi.Normalization", {jMatrix(box), jStrings(dimName)}), ...
    "srg",      jClass("rid.surrogate.LinearFirstOrder", {0.5, 0.1, 1}), ...
    "fields",   jRef("fields"), ...
    "Q",        jRef("Q"), ...
    "gen",      jClass("rid.gen.Direct", ...
                    {jClass("rid.gen.Constraints", {jRef("rep"), jRef("roi")}), nHorizon}), ...
    "sol",      jClass("rid.opt.MultiStart", {}), ...
    "state",    jClass("rid.crit.CoverageState", {jRef("Q"), jRef("epsReqInitial")}), ...
    "sch",      jClass("rid.surrogate.ImmediateScheduler", {}), ...
    "fieldSch", jClass("rid.fields.ImmediateScheduler", {}), ...
    "qInfo",    jStruct(struct("origin", "supplied", "method", "", ...
                               "seed", NaN, "hash", ""))));

node.options = jStruct(struct( ...
    "crit", jStruct(struct("epsReq0", epsReq0, "mode", "budget")), ...
    "rhc",  jStruct(struct("nSteps", nSteps)), ...
    "io",   jStruct(struct("saveRun", false))));

% EVERY ENTRY BELOW IS A DISCRETE FACT ABOUT THE PATH, or a count read to an order
% of magnitude. `affineTier` is per STEP rather than per surrogate epoch because
% the step vector states both halves at once: which step probed, and that no other
% step re-probed. With no surrogate update in this run there is exactly one epoch.
%
% `nActive`, `bestStart` and `deadOnArrival` are deliberately absent. The first two
% were measured to move with the iterate path (nActive ran from [10 9 9 9 7 6] to
% [9 8 6 6 6 6] over the ten start seeds), and the third is a count whose starts
% depend on the previous solution from step 2 on. The path criterion asserts the
% PROPERTIES of those instead -- non-increasing, never certified, always reported.
node.expected = jStruct(struct( ...
    "status",            r.status, ...
    "nSteps",            r.nSteps, ...
    "solver",            r.steps(1).optDiag.solver, ...
    "nStarts",           r.steps(1).optDiag.nStarts, ...
    "affineTier",        jStrings(string({r.steps.affineTier})), ...
    "gradientPath",      jStrings(string({r.steps.gradientPath})), ...
    "nEvaluations",      jColumn([r.steps.nEvaluations].'), ...
    "nEvaluationsTotal", sum([r.steps.nEvaluations]), ...
    "certification",     jStruct(struct( ...
        "mode",        r.certification.mode, ...
        "isSatisfied", double(r.certification.isSatisfied))), ...
    "nEvents",           numel(r.eventLog), ...
    "nCoverageResets",   sum(string({r.eventLog.kind}) == "coverageReset"), ...
    "nWarnings",         numel(r.warnings)));

end

% --- the JSON tree, tagged so the harness knows what to build ---------------

function n = jClass(cls, args),  n = struct("class", string(cls), "args", {args}); end
function n = jRef(name),         n = struct("ref", string(name));                  end
function n = jRow(v),            n = struct("row", reshape(v, 1, []));             end
function n = jColumn(v),         n = struct("column", reshape(v, [], 1));          end
function n = jMatrix(M),         n = struct("matrix", M);                          end
function n = jStrings(s),        n = struct("strings", {reshape(string(s), 1, [])}); end
function n = jStruct(s),         n = struct("plain", s);                           end

% --- the serializer: every double at %.17g, which round-trips exactly -------

function txt = encode(node, depth)
pad  = string(repmat('  ', 1, depth + 1));
padC = string(repmat('  ', 1, depth));

if isstruct(node) && isscalar(node) && isfield(node, "plain")
    txt = encodeObject(node.plain, depth);
elseif isstruct(node) && isscalar(node) && isfield(node, "class")
    inner = ["""class"": """ + node.class + """", ...
             """args"": " + encodeList(node.args, depth + 1)];
    txt = "{" + newline + pad + strjoin(inner, "," + newline + pad) + newline + padC + "}";
elseif isstruct(node) && isscalar(node) && isfield(node, "ref")
    txt = "{""ref"": """ + node.ref + """}";
elseif isstruct(node) && isscalar(node) && isfield(node, "row")
    txt = "{""row"": " + numberList(node.row) + "}";
elseif isstruct(node) && isscalar(node) && isfield(node, "column")
    txt = "{""column"": " + numberList(node.column) + "}";
elseif isstruct(node) && isscalar(node) && isfield(node, "strings")
    txt = "{""strings"": [" + strjoin("""" + node.strings + """", ", ") + "]}";
elseif isstruct(node) && isscalar(node) && isfield(node, "matrix")
    M    = node.matrix;
    rows = strings(size(M, 1), 1);
    for iRow = 1:size(M, 1)
        rows(iRow) = numberList(M(iRow, :));
    end
    txt = "{""matrix"": [" + newline + pad + ...
          strjoin(rows, "," + newline + pad) + newline + padC + "]}";
elseif isstruct(node) && isscalar(node)
    txt = encodeObject(node, depth);
elseif isstring(node) || ischar(node)
    txt = """" + string(node) + """";
elseif islogical(node)
    if node
        txt = "true";
    else
        txt = "false";
    end
elseif isnumeric(node) && isscalar(node)
    txt = num(node);
else
    error("rid:io:canaryEncodeFailed", "Unencodable node of class %s.", class(node));
end
end

function txt = encodeObject(s, depth)
pad   = string(repmat('  ', 1, depth + 1));
padC  = string(repmat('  ', 1, depth));
names = string(fieldnames(s));
parts = strings(numel(names), 1);
for iName = 1:numel(names)
    parts(iName) = """" + names(iName) + """: " + encode(s.(names(iName)), depth + 1);
end
txt = "{" + newline + pad + strjoin(parts, "," + newline + pad) + newline + padC + "}";
end

function txt = encodeList(items, depth)
if isempty(items)
    txt = "[]";
    return
end
pad   = string(repmat('  ', 1, depth + 1));
padC  = string(repmat('  ', 1, depth));
parts = strings(numel(items), 1);
for iItem = 1:numel(items)
    parts(iItem) = encode(items{iItem}, depth + 1);
end
txt = "[" + newline + pad + strjoin(parts, "," + newline + pad) + newline + padC + "]";
end

function txt = numberList(v)
parts = arrayfun(@num, reshape(v, 1, []));
txt   = "[" + strjoin(parts, ", ") + "]";
end

function s = num(x)
s = string(sprintf("%.17g", x));
end
