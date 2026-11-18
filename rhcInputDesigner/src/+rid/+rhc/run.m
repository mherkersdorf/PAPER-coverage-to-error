function result = run(process, setup, hist0, overrides)
%RUN  The receding-horizon loop: design, commit, measure, update.
%
%       setup  = rid.rhc.Setup(srg, rep, nrm, roi, Q, fields, gen, sol, state, ...
%                              sch, fieldSch, qInfo)
%       result = rid.rhc.run(process, setup, hist0, overrides)
%
%   At each step the whole future input trajectory over the horizon is optimized --
%   every future sample a decision variable -- by simulating it through the
%   surrogate and scoring the resulting points with `J`; then ONLY THE APPLIED UNIT
%   of the winning candidate is committed. The rest of the horizon is
%   look-ahead, never a plan.
%
%   TWO THINGS THIS FUNCTION IS NOT ALLOWED TO KNOW (invariant 6). It contains no
%   criterion maths -- certification is `state.isSatisfied()`, not a comparison
%   this file writes -- and no model maths: the surrogate is simulated
%   through its own interface and never inspected. A formula here would be in the
%   wrong package.
%
%   WHAT IS COMMITTED IS NOT WHAT WAS SCORED. The candidate's points are a
%   PREDICTION; what enters the coverage state are the REALIZED points, built from
%   the grown history according to `pastSource`. Committing the simulated
%   candidate points instead would let the design believe its own forecast.
%
%   THE INITIAL HISTORY COUNTS TOO. Whatever `rep.build` yields from
%   `hist0` is realized data and is added BEFORE step 1, so that the incremental
%   count and a full rebuild agree for ANY initial condition -- not merely for the
%   one-row history, which yields no regressor rows and hid the gap.
%
%   ONLINE AND OFFLINE ARE THE SAME LOOP. `isempty(process)` is the one permitted
%   branch on mode and the whole budget of it; everything else is carried by
%   `pastSource` and the surrogate (section 8). Declaring no process together with
%   `pastSource = "measured"` is incoherent -- an offline run has nothing to
%   measure -- and is refused at assembly, comparing two declarations the user made
%   and correcting neither.
%
%   THE STATE IS RE-DERIVED EVERY STEP, never carried. `xEnd` from the
%   previous step is not merely stale after a surrogate update, it is
%   type-incorrect: the new model's state lives in a different space. Re-deriving
%   from `hist.x0` with the current surrogate is automatically correct across an
%   update, and it keeps the distance cache the only persistent derived state.
%
%   WITH MEASURED STATES THE HORIZON IS ANCHORED ON A MEASUREMENT.
%   A representation that reads measured states looks back zero samples, so no
%   measured value would reach the horizon through the splice: the
%   measurements would change the realized points and never the design. The loop
%   then asks the process for its states as a second output, appends them to
%   `hist.x`, and each horizon starts one surrogate step from the last measured
%   state -- still re-derived every step, never carried. Which state each
%   horizon started from is in `steps(k).xNow`, the one place the anchor is
%   observable.
%
%   CANDIDATE POINTS COME FROM A WINDOW, the last `w` applied samples
%   followed by the horizon, through the ordinary `rep.build`. `w` is read from
%   `kIndex` -- the warm-up depth `build` itself discards -- and never from
%   a representation-specific property, so this works for any representation. The
%   alternative, passing the whole history, costs `O(N)` per OBJECTIVE EVALUATION
%   rather than per step.
%
%   FOUR WAYS TO STOP, and it never stops silently (section 5):
%
%       "certified"      certification mode, `state.isSatisfied()`
%       "budgetReached"  the prescribed experiment length
%       "stalled"        no improvement in `J` for `stallSteps` steps
%       "guardTripped"   maximum steps or maximum wall time
%
%   In BUDGET mode `J` reaching zero does not stop the run -- the length is the
%   commitment -- but it does trigger the coverage RESET, because a run
%   with every point resolved has an empty active set and a gradient that is zero
%   everywhere: every start would report dead on arrival and the remaining budget
%   would produce arbitrary input. The reset is a deliberate fiction, so the cache
%   understates coverage from then on, and the honest final coverage is recovered
%   by one `rebuild` over the full history at the end.
%
%   BOTH SCHEDULERS, IN ONE ORDER: surrogate first, then fields, because
%   estimators take `srg` and must see the new one; and the cache follows the same
%   order, `rebuild` before `setResolution`, each recomputing the active set.
%   Every event that makes earlier values incomparable is logged, and the
%   monotonicity assert never spans one. THE RESET COMES AFTER BOTH: a
%   rebuild in the same step would otherwise put the whole history back and annul
%   it inside the step that performed it.
%
%   THE CRITERION IS ANY rid.crit.Criterion, and the loop never asks
%   which class it holds. What is specific to the coverage criterion is reached
%   through the interface's defaults: certification and the reset through
%   `isSatisfied()`, invariant 3 through `isMonotone()`, the resolution check,
%   the `epsReqInitial` snapshot, field updates and the certification bound
%   through `resolution()` -- recorded as NaN where it is empty -- and
%   `nActive` through `nActive()`.
%
%   THE GRADIENT PATH IS DETECTED, NEVER ASSUMED. The affine map is MEASURED by
%   probe through the same `rep` and `srg` the loop uses, and the affine path is taken only when the check points pass; `pNorm = Inf`
%   is discovered by CALLING `probeGradient` and catching exactly
%   `rid:crit:gradientNormUnsupported`. Which path ran is in the
%   diagnostics, so a fall-back cannot masquerade as the fast path.
%
%   Inputs:
%       process    (1x1) object or []  online: `yMeasured = process.apply(uApplied)`,
%                                      (nApplied x nInputs) in, (nApplied x nOutputs)
%                                      out; and
%                                      `[yMeasured, xMeasured] = process.apply(uApplied)`
%                                      when the representation reads measured
%                                      states, xMeasured (nApplied x nStates)
%                                     . `[]` is offline.
%       setup      (1x1) rid.rhc.Setup the problem
%       hist0      (1x1) struct        the experiment's INITIAL CONDITION, in the
%                                      plain history struct: `.u` with
%                                      `u(1,:)` the pre-run actuator position
%                                     , `.x0`, `.tS`, `.y` when a process
%                                      will supply measurements, and `.x` when the
%                                      representation reads measured states
%       overrides  (1x1) struct        nested like rid.defaultOptions, partial;
%                                      must carry crit.epsReq0, plus rhc.nSteps
%                                      in budget mode
%   Outputs:
%       result  (1x1) struct
%           .hist      (1x1) struct        the designed signal: .u .y .x0 .tS,
%                                          and .x for a measured-state run
%           .srg       (1x1) Surrogate      the FINAL surrogate, which is not
%                                           `setup.srg` once a refit has been
%                                           applied
%           .status    (1x1) string        one of the four above
%           .J         (1x1) double        the final, honest criterion value
%           .nSteps    (1x1) double        steps designed
%           .steps     (nSteps x 1) struct per-step diagnostics (section 6)
%           .eventLog  (nEvents x 1) struct .step .kind .detail
%           .Jtrace    (nSteps x 1) double
%           .certification (1x1) struct    .isSatisfied .epsReqMax .fillDistanceQ
%                                          .bound .mode
%           .warnings  (n x 1) string      configuration diagnostics the run
%                                          produced; never printed
%           .opts      (1x1) struct        the EFFECTIVE options
%           .record    (1x1) struct        the run record
%           .recordPath (1x1) string       where it was written, "" if not
%
%   THE HORIZON HAS A FLOOR. `rep.minHorizon()` is the shortest horizon
%   that scores the applied input, and assembly refuses anything shorter with
%   `rid:rhc:horizonTooShort`. The representation states it, because only the
%   representation knows where its input block starts and this file must not name
%   `nk`.

arguments
    process
    setup     (1,1) rid.rhc.Setup
    hist0     (1,1) struct
    overrides (1,1) struct = struct()
end

[opts, explicitOverrides] = resolveOptions(overrides);

% `storeQ`'s CONSTRUCTION default is false, but section 7 promises "always for
% a user-supplied Q" and the stored-Q reuse extends that to a reused one -- a Q that is not
% this run's own generation has no seed to regenerate it from, so a record
% that dropped it would name points nobody can recover. That promise is
% enforced HERE, once, rather than left to whoever assembles overrides to
% remember, and only when the caller left `io.storeQ` unset: an EXPLICIT
% override, true or false, is never touched.
if setup.qInfo.origin ~= "generated" && ...
        ~(isfield(explicitOverrides, "io") && isfield(explicitOverrides.io, "storeQ"))
    opts.io.storeQ = true;
end

srg      = setup.srg;
rep      = setup.rep;
nrm      = setup.nrm;
Q        = setup.Q;
fields   = setup.fields;
gen      = setup.gen;
sol      = setup.sol;
state    = setup.state;
sch      = setup.sch;
fieldSch = setup.fieldSch;

nSup     = size(Q, 1);
nDim     = size(Q, 2);
nHorizon = gen.nHorizon;
nInputs  = gen.nInputs();

% --- mode, and the one branch on it (section 8) ----------------------
isOffline = isempty(process);
if isOffline && rep.pastSource == "measured"
    error("rid:rhc:incoherentDeclaration", ...
        "No process was supplied and the representation declares " + ...
        "pastSource = ""measured"". An offline run has nothing to measure. " + ...
        "Both are declarations you made and neither is corrected here: either " + ...
        "supply a process, or declare pastSource = ""simulated"".");
end
if ~isOffline
    assertProcess(process);
end

% --- the horizon must be long enough to SCORE the applied input -----
% Two declarations the user made, compared and neither corrected. A
% representation whose input block starts at u(k-nk) contains the applied sample
% first in the horizon's (nk+1)-th row, so a shorter horizon leaves `J` constant
% in the very sample about to be committed: every start scores alike, the
% optimizer returns whatever its start held, and an arbitrary input is applied
% with nothing raised. A warning was rejected -- the committed input would still
% be arbitrary. The floor is asked of the REPRESENTATION, never read as `nk`,
% which is a property of one representation class that +rhc must not name.
if nHorizon < rep.minHorizon()
    error("rid:rhc:horizonTooShort", ...
        "nHorizon is %d and this representation needs at least %d: its " + ...
        "input block does not start at u(k), so no scored point of a shorter " + ...
        "horizon depends on the sample about to be committed. J would be " + ...
        "constant in it and the loop would apply whatever the winning start " + ...
        "held. The floor is necessary, not sufficient -- a surrogate " + ...
        "with an internal dead time delays the predicted coordinates further.", ...
        nHorizon, rep.minHorizon());
end

isBudget = opts.crit.mode == "budget";

% --- pre-flight capture, before anything is computed ----------------
codeState   = rid.io.gitState(rid.io.projectRoot());
environment = rid.io.environment();
if opts.io.requireCleanTree && codeState.isDirty
    error("rid:io:dirtyWorkingTree", ...
        "requireCleanTree is set and the working tree has uncommitted " + ...
        "changes at %s. A run destined for a publication must name a commit " + ...
        "a reader can check out.", codeState.shortCommit);
end

% --- the experiment's initial condition --------------
% It ARRIVES as a history rather than being assembled here from two options. An
% initial condition is data, the object for it already exists, and
% `rid.assertHistory` is its one structural validator -- so `u(1,:)` being
% the pre-run actuator position and `x0` coming from `srg.steadyState` are
% the caller's declarations, made where the experiment is set up.
requireOption(opts, "crit", "epsReq0");
rid.assertHistory(hist0);
hist = hist0;

if size(hist.u, 2) ~= nInputs
    error("rid:rhc:invalidInitialHistory", ...
        "hist0.u has %d column(s); the parametrization declares %d input " + ...
        "channel(s). Samples are in rows.", size(hist.u, 2), nInputs);
end
if size(hist.u, 1) < 1
    error("rid:rhc:invalidInitialHistory", ...
        "hist0.u has no rows. A run starts from the actuator's position " + ...
        "before it began -- unrecorded, not undefined -- so there is " + ...
        "always at least one.");
end
% ONLINE MEANS THE HISTORY CARRIES MEASUREMENTS FROM ITS FIRST ROW. The loop
% appends one per step, and `y` with fewer rows than `u` is a malformed history,
% so the absence has to be caught here rather than at the first append.
if ~isOffline && ~isfield(hist, "y")
    error("rid:rhc:invalidInitialHistory", ...
        "A process was supplied and hist0 has no ""y"". An online run's " + ...
        "history carries a measurement for every applied sample, including " + ...
        "the pre-run one; an absent y is how an OFFLINE history says it has " + ...
        "none.");
end
% AND THE SAME ARGUMENT FOR THE STATES. A representation that reads
% measured states needs `hist.x(1,:)` -- the state measured before the run began --
% because the loop appends one row per step from the first and an `x` with fewer
% rows than `u` is a malformed history. It is asked of the REPRESENTATION, never
% inferred from the presence of the field: an absent `x` is how a history says no
% states were measured, and inferring the mode from it would silently turn a
% mis-specified run into a different one.
if ~isOffline && rep.needsMeasuredStates() && ~isfield(hist, "x")
    error("rid:rhc:invalidInitialHistory", ...
        "The representation reads measured states and hist0 has no ""x"". " + ...
        "The run's history carries the state measured at every applied sample, " + ...
        "including the pre-run one, so hist0.x(1,:) is where it starts " + ...
        ". An absent x is how a history says NO states were " + ...
        "measured; declare pastSource = ""simulated"" for such a run.");
end
% How many states the process reports is the HISTORY's declaration, read once.
% Zero means the one-output the process contract form, and is the only thing the loop branches on
% when it applies.
nStatesMeasured = 0;
if ~isOffline && rep.needsMeasuredStates()
    nStatesMeasured = size(hist.x, 2);
end
% The sample time is the HISTORY's. A surrogate discretized at a
% different one is a real inconsistency and silent otherwise: every simulated
% horizon would advance by a different amount of physical time than the data.
if isfinite(srg.tS) && srg.tS ~= hist.tS
    error("rid:rhc:inconsistentSampleTime", ...
        "hist0.tS is %g and the surrogate is discretized at %g. The two " + ...
        "describe the same experiment.", hist.tS, srg.tS);
end

% --- at least one stopping rule must be finite (section 5) ------------------
% Before any work: a run that cannot end, or that omits its budget, is refused
% here rather than after the fill-distance measurement and the resolution report.
if isBudget
    requireOption(opts, "rhc", "nSteps");
end
assertStoppable(opts, isBudget);

% --- what this discretization can resolve ---------------------------
[fillDistanceQ, fillInfo] = rid.roi.estimateFillDistance(Q, nrm, ...
    "nProbes", opts.roi.nProbes, "seed", opts.roi.probeSeed, ...
    "probeMethod", opts.roi.probeMethod, "blockSize", opts.roi.blockSize, ...
    "traceTolerance", opts.roi.traceTolerance);
% The estimate's own verdict on its probe trace travels with the headroom it was
% used to compute: an unconverged fillDistanceQ makes the certifiable
% floor look finer than it is, and the run-start statement says so.
resolution = rid.crit.resolutionReport(opts.crit.epsReq0, fillDistanceQ, ...
                                       opts.crit.mode, ...
                                       "estimateConverged", fillInfo.converged, ...
                                       "probeFillDistance", fillInfo.probeFillDistance, ...
                                       "relativeClimb", fillInfo.trace.relativeClimb);

% --- the fields the run's epsReq was derived from -------------------
% They arrive in `setup`, beside `Q`, because they are data indexed by the
% supporting points rather than settings. The check below is
% what makes that safe: the field scheduler recomputes `epsReq` from these, so a
% state built from different ones would have its resolution silently replaced by
% one derived from fields the run never saw. A criterion with no demanded
% resolution has nothing to compare; a field update for it raises at
% `setResolution` instead.
epsReq = demandedFrom(fields, opts, fillDistanceQ);
hasResolution = ~isempty(state.resolution());
if hasResolution && ~isequal(epsReq, state.resolution())
    error("rid:rhc:inconsistentResolution", ...
        "The coverage state's epsReq is not what these fields produce " + ...
        "(largest difference %g). The field scheduler overrides fields and " + ...
        "recomputes epsReq from them, so a run whose state was built from " + ...
        "different fields would have its resolution silently replaced by one " + ...
        "derived from fields it never saw. Pass the same fields to rid.rhc.Setup, " + ...
        "or build the state from rid.crit.demandedResolution with the same " + ...
        "options.", max(abs(epsReq - state.resolution())));
end

% --- run state (section 6) --------------------------------------------------
kStep        = 0;
thetaPrev    = zeros(0, 1);
events       = emptyEventLog();
steps        = emptyStepLog();
Jtrace       = zeros(0, 1);
JBaseline    = Inf;     % the monotonicity reference; Inf until the first step
Jbest        = Inf;     % the stall reference
stallCount   = 0;
didReset     = false;
resetAt      = 0;      % samples before this index were discarded by a reset
nConsecutive = 0;      % resets in a row; the degenerate-freeze diagnostic
warnings     = strings(0, 1);
needAffine   = true;
affineMatrix = zeros(0, 0);
% BUILT HERE, NOT ONLY AFTER AN UPDATE. The levels invert a STATIC map,
% and `setup` holds a surrogate from the first step -- nothing about the
% construction requires the map to have just changed, only that it be current,
% which the update branch below keeps true. Called only there, the vector was
% never built at all in a default run, because `updateInterval` is Inf:
% every run initialized unguided, which is the fallback the surrogate-guided initialization exists to replace.
[targetLevels, levelInfo] = targetLevelsFor(srg, rep, setup.roi, nInputs, ...
    opts.rhc.nTargetLevels, opts.rhc.nStaticProbes);
roiRange     = roiInputRanges(rep, setup.roi, nInputs);

% `pNorm = Inf` is a supported mode of the criterion whose GRADIENT is not
% implemented. Found by calling, not by reading opts.crit.pNorm: the
% loop asks the criterion and believes the answer.
useAnalytic = gradientAvailable(state, nDim);
penaltyOpts = struct("violationScale", opts.gen.violationScale, ...
                     "tailWeight",     opts.gen.tailWeight, ...
                     "psiMax",         opts.crit.psiMax);

% What replay must rebuild the state from; NaN where the criterion has no
% resolution, which then has no configuration the run can move.
epsReqInitial = resolutionOrNaN(state);

% --- the initial history is DATA, and it counts ---------------------
% The step loop admits only rows created by the step just committed, so rows
% `rep.build` already yields from `hist0` would be added by nothing: the run
% would overstate the work remaining and spend budget re-covering ground the
% supplied data already covered. The trivial one-row history yields ZERO rows
% and this is then a no-op, which is exactly why the defect was
% invisible. Re-adding a point the state already holds cannot corrupt it: `add`
% takes a minimum against the cached nearest-neighbour distances.
[Zprimed, ~] = rep.build(hist, srg);
if ~isempty(Zprimed)
    state.add(nrm.toNormalized(Zprimed));
end

runTimer      = tic;
status        = "";

while true

    % --- stop? ------------------------------------------------------------
    if ~isBudget && state.isSatisfied()
        status = "certified";
        break
    end
    if isBudget && kStep >= opts.rhc.nSteps
        status = "budgetReached";
        break
    end
    if kStep >= opts.rhc.maxSteps || toc(runTimer) >= opts.rhc.maxWallTime
        status = "guardTripped";
        break
    end
    if stallCount >= opts.rhc.stallSteps
        status = "stalled";
        break
    end

    stepTimer = tic;
    JBefore   = state.value();

    % --- the geometry of this step -------------------------
    [~, kIndexNow] = rep.build(hist, srg);
    w   = size(hist.u, 1) - numel(kIndexNow);
    ctx = struct("gen", gen, "rep", rep, "nrm", nrm, "srg", srg, ...
                 "state", state, "win", stepWindow(hist, srg, rep, w), ...
                 "nHorizon", nHorizon, "A", affineMatrix, ...
                 "penaltyOpts", penaltyOpts);

    % --- the starts provider: +opt decides the count, +gen answers -
    stepSeed  = opts.gen.seed + kStep;
    startsFcn = @(n) gen.initialGuesses(hist, thetaPrev, "nStarts", n, ...
                    "seed", stepSeed, "targetLevels", targetLevels, ...
                    "roiInputRange", roiRange);

    % --- the affine map, once per surrogate update
    affineInfo = struct("tier", "notProbed", "maxRelError", NaN);
    if useAnalytic && needAffine
        % The count comes from opts.opt.nStarts, where the starts-provider rule declares it --
        % not from the solver object, which need not have starts at all:
        % rid.opt.GridSearch reports nStarts = 0 and is now selectable (the single-setup rule
        % amendment). One declaration, two readers, no second source to drift.
        Theta0 = startsFcn(opts.opt.nStarts);
        % The probe's own points are axis-aligned, so a validation point must be
        % GENERIC. The remaining start rows are free generic points already
        % inside the feasible set.
        % The EFFECTIVE probe settings, not affineMap's own defaults: checkTol
        % decides the affine/notAffine verdict and hence which gradient path this
        % run takes, so a record that did not carry the value it ran under could
        % not explain its own `gradientPath`.
        [affineMatrix, ~, affineInfo] = rid.gen.affineMap( ...
            @(theta) znormOnly(ctx, theta), Theta0(1, :).', ...
            "checkPoints", Theta0(2:end, :), ...
            "h",           opts.gen.affineProbeStep, ...
            "checkTol",    opts.gen.affineCheckTol);
        useAnalytic = affineInfo.tier == "affine";
        needAffine  = false;
        ctx.A       = affineMatrix;
    end

    % --- the objective, in the arity the detected tier supports ----
    if useAnalytic
        scoreFcn = @(theta) stepObjective(ctx, theta);
    else
        scoreFcn = @(theta) valueOnly(ctx, theta);
    end

    [best, optDiag] = sol.solve(gen, hist, scoreFcn, startsFcn);

    % --- commit the APPLIED UNIT, whatever the parametrization says
    uApplied  = gen.applied(best.uSeq);
    nApplied  = size(uApplied, 1);
    yMeasured = zeros(0, 0);
    xMeasured = zeros(0, 0);
    if ~isOffline
        [yMeasured, xMeasured] = applyToProcess(process, uApplied, nApplied, ...
            size(hist.y, 2), nStatesMeasured);
    end

    hist.u = [hist.u; uApplied];
    if ~isOffline
        hist.y = [hist.y; yMeasured];
    end
    if nStatesMeasured > 0
        hist.x = [hist.x; xMeasured];
    end
    kStep     = kStep + 1;
    thetaPrev = best.theta;

    % --- the REALIZED points, identified from kIndex --------------
    [Zall, kIndex] = rep.build(hist, srg);
    kNow    = size(hist.u, 1);
    newRows = kIndex > kNow - nApplied;
    state.add(nrm.toNormalized(Zall(newRows, :)));

    JAfter = state.value();
    % Invariant 3 is a property of the criterion, asserted only where it holds:
    % a criterion that is not monotone in the samples may rise when a sample
    % lands where the data are already dense.
    if state.isMonotone()
        assertMonotone(opts, JBaseline, JAfter, kStep);
    end
    JBaseline = JAfter;
    Jtrace(kStep, 1) = JAfter; %#ok<AGROW>

    if JAfter < Jbest - opts.rhc.stallTolerance
        Jbest      = JAfter;
        stallCount = 0;
    else
        stallCount = stallCount + 1;
    end

    % --- the schedulers: surrogate FIRST, then fields -------------
    [srgNew, schInfo] = sch.poll(hist, srg);
    if schInfo.status ~= "none"
        events(end + 1, 1) = event(kStep, "surrogateUpdate", schInfo); %#ok<AGROW>
    end
    if ~isempty(srgNew)
        srg = srgNew;
        % THE NEW MODEL MUST STILL FIT THE SPACE. `Setup` asked
        % this at assembly, but an update may return a model of a different state
        % dimension, and the coordinates ARE the surrogate's states. Asked again
        % here so the refusal names the swap rather than surfacing one `build`
        % later as a shape complaint about a trajectory.
        rep.assertCompatible(srg);
        % A surrogate update moves historical coordinates when the coordinates or
        % the past outputs depend on the model (invariant 3). The representation
        % is asked; the loop never infers it from the mode.
        %
        % IT REBUILDS FROM THE POST-RESET HISTORY, NOT THE WHOLE OF IT.
        % A reset is the deliberate fiction that no data existed, and the coverage reset says
        % the cache understates coverage FROM THEN ON -- which a rebuild over the
        % full history would silently revoke, putting every discarded sample back
        % and returning the run to the empty active set the reset exists to
        % escape. `resetAt` is that boundary; the honest final coverage is still
        % recovered from the full history at the end, which is the one place the
        % fiction is meant to be dropped.
        if rep.requiresRebuildOnSurrogateUpdate()
            [Zrebuilt, kIdxRebuilt] = rep.build(hist, srg);
            kept = kIdxRebuilt > resetAt;
            state.rebuild(nrm.toNormalized(Zrebuilt(kept, :)));
        end
        % `A` is constant only BETWEEN surrogate updates, and the guided target
        % levels invert a static map that has just changed. This REFRESHES
        % what was built at setup; it is not the only construction.
        needAffine   = true;
        useAnalytic  = gradientAvailable(state, nDim);
        [targetLevels, levelInfo] = targetLevelsFor(srg, rep, setup.roi, ...
            nInputs, opts.rhc.nTargetLevels, opts.rhc.nStaticProbes);
        JBaseline  = Inf;
        Jbest      = Inf;
        stallCount = 0;
    end

    [fieldsNew, fieldInfo] = fieldSch.poll(hist, srg, rep, nrm, Q, fields);
    if fieldInfo.status ~= "none"
        events(end + 1, 1) = event(kStep, "resolutionChange", fieldInfo); %#ok<AGROW>
    end
    if ~isempty(fieldsNew)
        fields = fieldsNew;
        epsReq = demandedFrom(fields, opts, fillDistanceQ);
        state.setResolution(epsReq);
        JBaseline  = Inf;
        Jbest      = Inf;
        stallCount = 0;
    end

    % --- the budget-mode coverage reset, LAST in the step -
    % It follows both schedulers so that nothing in the step can re-add the
    % history behind it. Reset first and rebuilt second, a step doing both put
    % the whole history back and restored full coverage inside the step that had
    % just cleared it: the run then sat permanently in the state the coverage reset exists to
    % prevent -- empty active set, zero gradient, every start dead on arrival --
    % resetting and being annulled at every step. `requiresRebuildOnSurrogateUpdate`
    % is true for any simulated-past run, so that is ordinary offline design with
    % updates on, not a corner.
    %
    % Placing it here also asks the right QUESTION: whether the state is certified
    % under the model and the resolution now in force, rather than under the ones
    % the step began with. The active-set rule's order is untouched -- `rebuild` still precedes
    % `setResolution`, and both precede this.
    didResetHere = false;
    if isBudget && opts.rhc.resetOnCertification && state.isSatisfied()
        state.rebuild(zeros(0, nDim));
        nConsecutive = nConsecutive + 1;
        % CONSECUTIVE RESETS ARE A STATEMENT ABOUT THE CONFIGURATION, not about
        % the run. If a single horizon can cover the whole ROI the cache
        % is empty again at every step, the problem posed to the optimizer is
        % identical each time, and the input freezes into a repeat. That is
        % information the user needs -- nSup too small, epsReq0 too loose, or
        % nHorizon too long against the ROI.
        %
        % It is carried in the log and in `result.warnings`, NOT raised with
        % `warning`: +rhc may not write to the console, and
        % diagnostics travel in returned structs. No fifth termination status
        % either -- the stall guard genuinely will not fire, since J jumps at
        % every reset, and a status for a mis-specification is machinery.
        severity = "info";
        if nConsecutive > 1
            severity = "warning";
            warnings(end + 1, 1) = sprintf("Step %d: %d consecutive coverage " + ...
                "resets. A single horizon is covering the whole ROI, so the " + ...
                "design problem is identical at every step and the input " + ...
                "freezes. nSup may be too small, epsReq0 too loose, or " + ...
                "nHorizon too long against the ROI.", ...
                kStep, nConsecutive); %#ok<AGROW>
        end
        events(end + 1, 1) = event(kStep, "coverageReset", ...
            struct("nSup", nSup, "nSamplesDiscarded", kNow, ...
                   "consecutive", nConsecutive, "severity", severity, ...
                   "reason", "J reached zero in budget mode; every point is " + ...
                   "marked unresolved so the design keeps a gradient")); %#ok<AGROW>
        didReset     = true;
        didResetHere = true;
        resetAt      = kNow;   % everything up to here is discarded
        % J before and after a reset is not comparable -- a discontinuity of the
        % same species as a resolution change.
        JBaseline  = Inf;
        Jbest      = Inf;
        stallCount = 0;
    else
        nConsecutive = 0;
    end

    steps(kStep, 1) = stepRecord(kStep, JBefore, JAfter, state, best, ...
        optDiag, uApplied, yMeasured, ctx.win.xNow, w, affineInfo, ...
        useAnalytic, didResetHere, toc(stepTimer)); %#ok<AGROW>
end

% --- the honest final coverage --------------------------------------
if didReset
    % The reset is a deliberate fiction and the cache understates coverage from
    % then on. One rebuild over the full history, O(nSup * nAll) once, at a point
    % where that is nothing.
    [Zfinal, ~] = rep.build(hist, srg);
    state.rebuild(nrm.toNormalized(Zfinal));
end

certification = struct( ...
    "mode",          opts.crit.mode, ...
    "isSatisfied",   state.isSatisfied(), ...
    "epsReqMax",     max(resolutionOrNaN(state)), ...
    "fillDistanceQ", fillDistanceQ, ...
    "bound",         max(resolutionOrNaN(state)) + fillDistanceQ, ...
    "report",        resolution);

result = struct( ...
    "hist",          hist, ...
    "srg",           srg, ...
    "status",        status, ...
    "J",             state.value(), ...
    "nSteps",        kStep, ...
    "steps",         steps, ...
    "eventLog",      events, ...
    "Jtrace",        Jtrace, ...
    "certification", certification, ...
    "warnings",      warnings, ...
    "opts",          opts, ...
    "record",        struct(), ...
    "recordPath",    "");

result.record = buildRecord(result, setup, opts, explicitOverrides, ...
    codeState, environment, fillDistanceQ, fillInfo, fields, levelInfo, ...
    epsReqInitial);

if opts.io.saveRun
    % The EFFECTIVE runsDir, not the default: an option recorded in the effective
    % struct and then ignored is worse than one that does not exist.
    % The EFFECTIVE runsDir, through the one resolution `loadRun` also uses: an
    % option recorded in the effective struct and then ignored is worse than one
    % that does not exist, and a write the reader cannot find is worse still.
    result.recordPath = rid.io.saveRun(result.record, opts.io.note, ...
        "root", rid.io.runsRoot(opts.io.runsDir));
end

end

% ---------------------------------------------------------------------------

function assertProcess(process)
%ASSERTPROCESS  The contract, checked before a run rather than mid-run.
%
%   Duck-typed on purpose: the process is the caller's adapter to their own
%   hardware or simulator, and the project declares no base class for it. What is
%   checked is the one method the contract names.
if ~(isobject(process) || isstruct(process)) || ...
        ~(ismethod(process, "apply") || (isstruct(process) && isfield(process, "apply")))
    error("rid:rhc:invalidProcess", ...
        "process must be [] for an offline run, or an object with an " + ...
        "`apply` method: yMeasured = process.apply(uApplied), " + ...
        "(nApplied x nInputs) in and (nApplied x nOutputs) out. " + ...
        "Got %s.", class(process));
end

end

% ---------------------------------------------------------------------------

function [yMeasured, xMeasured] = applyToProcess(process, uApplied, nApplied, ...
                                                 nOutputs, nStatesMeasured)
%APPLYTOPROCESS  The applied unit in, its measurement out.
%
%   ERRORS RAISE AND ARE NOT CAUGHT. What to DO about a failing process -- retry,
%   abort gracefully, park the actuator -- is a policy, and policies stay out until
%   there is real hardware to have opinions about, exactly as updateInterval = Inf
%   keeps update policy out.
%
%   THE SECOND OUTPUT IS ASKED FOR ONLY WHEN THE REPRESENTATION READS MEASURED
%   STATES. A process that measures none keeps the process contract's one-output form untouched,
%   which is what lets every existing adapter go on working; requesting two from
%   all of them would break each one at its first step.
%
%   A PROCESS THAT CANNOT SUPPLY THEM IS A CONTRACT VIOLATION, NOT A MATLAB
%   ACCIDENT. `[y, x] = process.apply(u)` against a one-output process raises
%   MATLAB's own "Too many output arguments" -- `MATLAB:TooManyOutputs` for a class
%   method, `MATLAB:maxlhs` for a struct carrying a function handle -- which says
%   nothing about the declaration that asked for them and points the reader at the
%   loop rather than at the pairing they chose. Exactly those two identifiers are
%   translated and everything else is rethrown, so a genuine failure INSIDE the
%   process is never disguised as a shape complaint. That is the discipline the scoreFcn arity rule
%   applies to the objective's arity: classify by calling, and catch precisely.
xMeasured = zeros(0, 0);

if nStatesMeasured == 0
    yMeasured = process.apply(uApplied);
else
    try
        [yMeasured, xMeasured] = process.apply(uApplied);
    catch err
        if ~ismember(string(err.identifier), ...
                ["MATLAB:TooManyOutputs", "MATLAB:maxlhs"])
            rethrow(err)
        end
        error("rid:rhc:invalidMeasurement", ...
            "The representation reads measured states, so the loop asks this " + ...
            "process for them -- [yMeasured, xMeasured] = " + ...
            "process.apply(uApplied), with xMeasured (%d x %d) (measured states, as the representation " + ...
            "requires) -- and it returns one output only. Either measure the states, " + ...
            "or declare pastSource = ""simulated"".", nApplied, nStatesMeasured);
    end
end

if ~isnumeric(yMeasured) || ~isreal(yMeasured) || ~ismatrix(yMeasured) || ...
        ~all(isfinite(yMeasured), "all")
    error("rid:rhc:invalidMeasurement", ...
        "process.apply returned a %s; the contract is a real, finite " + ...
        "(nApplied x nOutputs) numeric matrix.", class(yMeasured));
end
if ~isequal(size(yMeasured), [nApplied, nOutputs])
    error("rid:rhc:invalidMeasurement", ...
        "process.apply returned %s for an applied unit of %d sample(s) and " + ...
        "%d output channel(s). It measures exactly what was applied.", ...
        mat2str(size(yMeasured)), nApplied, nOutputs);
end

% THE STATES ARE CHECKED LIKE THE OUTPUT, and separately, so the message names
% which of the two is wrong. Row k is the state at the instant u(k) was applied
%, so the row count is the applied unit's, exactly as for y.
if nStatesMeasured > 0
    if ~isnumeric(xMeasured) || ~isreal(xMeasured) || ~ismatrix(xMeasured) || ...
            ~all(isfinite(xMeasured), "all")
        error("rid:rhc:invalidMeasurement", ...
            "process.apply returned a %s as its measured states; the " + ...
            "contract is a real, finite (nApplied x nStates) numeric matrix.", ...
            class(xMeasured));
    end
    if ~isequal(size(xMeasured), [nApplied, nStatesMeasured])
        error("rid:rhc:invalidMeasurement", ...
            "process.apply returned %s as its measured states for an applied " + ...
            "unit of %d sample(s), where the history declares %d state(s). " + ...
            "Row k is the state at the instant u(k) was applied.", ...
            mat2str(size(xMeasured)), nApplied, nStatesMeasured);
    end
end

end

% ---------------------------------------------------------------------------

function tf = gradientAvailable(state, nDim)
%GRADIENTAVAILABLE  Is the cached gradient implemented for this criterion?
%
%   ASKED BY CALLING, and only `rid:crit:gradientNormUnsupported` is caught --
%   the same discipline the scoreFcn arity rule applies to the objective's arity. Rethrowing
%   everything else is the load-bearing half: without it a genuine bug inside the
%   criterion would be indistinguishable from an unimplemented gradient, and the
%   run would silently drop onto the finite-difference path.
try
    state.probeGradient(zeros(0, nDim));
    tf = true;
catch err
    if err.identifier ~= "rid:crit:gradientNormUnsupported"
        rethrow(err)
    end
    tf = false;
end

end

% ---------------------------------------------------------------------------

function Znorm = znormOnly(ctx, theta)
%ZNORMONLY  The composition the affine probe measures: theta -> normalized rows.
[~, Znorm] = candidatePoints(ctx, theta);

end

% ---------------------------------------------------------------------------

function f = valueOnly(ctx, theta)
%VALUEONLY  A genuinely single-output objective, for the tier-3 path.
%
%   It must have no second output to give: `probeArity` classifies by CALLING, and
%   an anonymous handle wrapping the three-output objective would pass its own
%   nargout through and be misread as the analytic path.
f = stepObjective(ctx, theta);

end

% ---------------------------------------------------------------------------

function epsReq = resolutionOrNaN(state)
%RESOLUTIONORNAN  The criterion's demanded resolution, or NaN where it has none.
%
%   NaN encodes the ABSENCE of the concept, never a default: a criterion
%   without a resolution has no bound and no initial resolution to snapshot.
epsReq = state.resolution();
if isempty(epsReq)
    epsReq = NaN;
end

end

% ---------------------------------------------------------------------------

function epsReq = demandedFrom(fields, opts, fillDistanceQ)
%DEMANDEDFROM  epsReq from the fields, with this run's options. One call site each.
epsReq = rid.crit.demandedResolution(fields, ...
    "epsReq0",       opts.crit.epsReq0, ...
    "beta",          opts.crit.beta, ...
    "betaNoise",     opts.crit.betaNoise, ...
    "mode",          opts.crit.mode, ...
    "fillDistanceQ", fillDistanceQ);

end

% ---------------------------------------------------------------------------

function assertMonotone(opts, JBaseline, JAfter, kStep)
%ASSERTMONOTONE  Invariant 3, and only where it holds.
%
%   `J` is monotone non-increasing BETWEEN surrogate and resolution updates: a
%   committed sample can only shrink a nearest-neighbour distance, so every psi_j
%   is non-increasing and so is their mean. Across an update, a resolution change
%   or a coverage reset it is not comparable at all, which is why the baseline is
%   set to Inf at each of those and this assert then passes vacuously.
%
%   THE ALLOWANCE IS ROUNDING, NOT TOLERANCE. Summing smaller terms can round up
%   by an ulp; a real violation -- a stale cache, points added that were never
%   committed, the wrong rows selected -- is a jump, not an ulp.
if ~opts.rhc.assertMonotonicity || ~isfinite(JBaseline)
    return
end
allowance = 8 * eps(max(JBaseline, 1));
if JAfter > JBaseline + allowance
    error("rid:rhc:monotonicityViolated", ...
        "J rose from %.17g to %.17g at step %d with no surrogate update, " + ...
        "resolution change or coverage reset in between. Adding samples can " + ...
        "only shrink nearest-neighbour distances, so this is a stale cache or " + ...
        "the wrong rows committed, not a numerical artefact (algorithmic " + ...
        "invariant 3).", JBaseline, JAfter, kStep);
end

end

% ---------------------------------------------------------------------------

function assertStoppable(opts, isBudget)
%ASSERTSTOPPABLE  Refuse a run that has no finite way to end.
%
%   Certification may never be reached -- the remaining supporting points can be
%   dynamically unreachable, which is what "stalled" exists to report -- so a
%   certification run with every guard at Inf would loop forever. Refusing is the
%   honest answer; inventing a default cap would be a number nobody chose.
finite = (isBudget && isfinite(opts.rhc.nSteps)) || ...
         isfinite(opts.rhc.maxSteps) || isfinite(opts.rhc.maxWallTime) || ...
         isfinite(opts.rhc.stallSteps);
if ~finite
    error("rid:rhc:noStoppingRule", ...
        "Every stopping rule is infinite: nSteps, maxSteps, maxWallTime and " + ...
        "stallSteps. Certification may never be reached -- the remaining " + ...
        "supporting points can be dynamically unreachable -- so set at least " + ...
        "one finite guard. No default cap is invented here.");
end

end

% ---------------------------------------------------------------------------

function requireOption(opts, section, name)
%REQUIREOPTION  A commitment with no defensible default, on the precedent.
if ~isfield(opts.(section), name)
    reasons = struct( ...
        "epsReq0", "the experiment's baseline commitment -- the coarsest " + ...
                   "resolution demanded anywhere", ...
        "nSteps",  "the prescribed experiment length, which IS the commitment " + ...
                   "in budget mode");
    error("rid:rhc:missingOption", ...
        "opts.%s.%s is required and has no default: it is %s.", ...
        section, name, reasons.(name));
end

end

% ---------------------------------------------------------------------------

function ranges = roiInputRanges(rep, roi, nInputs)
%ROIINPUTRANGES  The ROI's own u range per channel, (nInputs x 2).
%
%   The unguided start row draws its targets over this rather than over the
%   amplitude box: the ROI is what the design is trying to cover. Coordinates are
%   classified by `.kind` and `.channel` and the lowest `.lag` is taken -- never by
%   parsing a name.
box    = rep.bounds(roi);
coords = rep.coordinateInfo();
kinds  = [coords.kind];
ranges = zeros(nInputs, 2);

for cIn = 1:nInputs
    dims = find(kinds == "inputLevel" & [coords.channel] == cIn);
    if isempty(dims)
        % No input level in the space: let +gen fall back to the amplitude box,
        % which it records in info.targetRange rather than substituting silently.
        ranges = zeros(0, 2);
        return
    end
    [~, iLowest]   = min([coords(dims).lag]);
    ranges(cIn, :) = box(dims(iLowest), :);
end

end

% ---------------------------------------------------------------------------

function e = event(step, kind, detail)
%EVENT  One entry of the event log.
%
%   `kind` is drawn from a CLOSED set and code branches on it, never on a parsed
%   string. The set has ONE owner, `rid.eventKinds`, read here and by
%   `+viz`: it used to be written out in both, and two copies of one normative
%   list are two things that can drift.
%
%   `"coverageReset"` IS ITS OWN KIND and not a mislabelled `"resolutionChange"`:
%   a budget-mode reset clears the coverage cache without `epsReq` changing
%   at all, and calling it a resolution change would put a false statement in the
%   one record a later reader trusts.
kinds = rid.eventKinds();
if ~ismember(kind, kinds)
    error("rid:rhc:unknownEventKind", ...
        "``%s'' is not one of the closed set %s.", kind, ...
        strjoin("""" + kinds + """", ", "));
end
e = struct("step", step, "kind", string(kind), "detail", detail);

end

% ---------------------------------------------------------------------------

function log = emptyEventLog()
%EMPTYEVENTLOG  A zero-row event log with the right fields.
log = struct("step", {}, "kind", {}, "detail", {});
log = reshape(log, 0, 1);

end

% ---------------------------------------------------------------------------

function log = emptyStepLog()
%EMPTYSTEPLOG  A zero-row step log with the right fields.
log = struct("step", {}, "JBefore", {}, "JAfter", {}, "nActive", {}, ...
             "theta", {}, "uApplied", {}, "yMeasured", {}, "xNow", {}, ...
             "bestStart", {}, ...
             "scoreSpread", {}, "nEvaluations", {}, "penalty", {}, ...
             "deadOnArrival", {}, "gradientPath", {}, "affineTier", {}, ...
             "windowDepth", {}, "coverageReset", {}, "wallTime", {}, ...
             "optDiag", {});
log = reshape(log, 0, 1);

end

% ---------------------------------------------------------------------------

function s = stepRecord(kStep, JBefore, JAfter, state, best, optDiag, ...
                        uApplied, yMeasured, xNow, w, affineInfo, ...
                        useAnalytic, didResetHere, wallTime)
%STEPRECORD  The per-step diagnostics of the method description.
%
%   `nActive` is free -- CoverageState maintains it anyway, and a criterion
%   without an active set reports NaN -- and is the most
%   direct read there is on how much of the ROI is still in play: it falls
%   monotonically between updates and reaches zero exactly at certification.
%   LOGGED, NEVER TERMINATED ON: `nActive == 0` and `J == 0` are equivalent, and
%   the run stops on `J` so that one fact has one source.
%
%   `bestStart` is the winning multi-start row, which is what "chosen candidate
%   index" means for a continuous parametrization: there is no candidate set to
%   index into.
%
%   `xNow` IS THE STATE THIS STEP'S HORIZON STARTED FROM (section 6): re-derived
%   from `hist.x0`, or one surrogate step from the last measured state when the
%   representation reads measured states. It is the ONLY place the anchor
%   part 6 is observable -- the anchor changes the design and nothing else in the
%   record would show whether it acted.
[~, bestStart] = min(optDiag.scores);
s = struct( ...
    "step",          kStep, ...
    "JBefore",       JBefore, ...
    "JAfter",        JAfter, ...
    "nActive",       state.nActive(), ...
    "theta",         best.theta, ...
    "uApplied",      uApplied, ...
    "yMeasured",     yMeasured, ...
    "xNow",          xNow, ...
    "bestStart",     bestStart, ...
    "scoreSpread",   optDiag.scoreSpread, ...
    "nEvaluations",  optDiag.nEvaluations, ...
    "penalty",       optDiag.penalty, ...
    "deadOnArrival", optDiag.deadOnArrival, ...
    "gradientPath",  optDiag.gradientPath, ...
    "affineTier",    string(affineInfo.tier) + ...
                     ternary(useAnalytic, "", " (analytic path off)"), ...
    "windowDepth",   w, ...
    "coverageReset", didResetHere, ...
    "wallTime",      wallTime, ...
    "optDiag",       optDiag);

end

% ---------------------------------------------------------------------------

function out = ternary(condition, ifTrue, ifFalse)
%TERNARY  A one-line conditional, for building a diagnostic string.
if condition
    out = ifTrue;
else
    out = ifFalse;
end

end

% ---------------------------------------------------------------------------

function rec = buildRecord(result, setup, opts, overrides, codeState, ...
                           environment, fillDistanceQ, fillInfo, fields, ...
                           levelInfo, epsReqInitial)
%BUILDRECORD  The run record of the method description.
%
%   ALWAYS BUILT, as part of `result`; writing it to disk is a separate decision
%   (section 2 rule 1). The EFFECTIVE options struct is the load-bearing item: a
%   record storing only "defaults plus epsReq0 = 0.05" changes meaning the day a
%   default changes, and defaults will change while the psiMax = 1 default and the betaNoise default are open. The
%   overrides are kept alongside because they preserve the INTENT.
%
%   THE SETUP IS STORED IN FULL, in its own group, not as a summary. The single-setup rule makes
%   `setup` exactly what `replay` consumes, and a `describe()` string cannot be
%   consumed -- rebuilding objects by parsing one would make a display convention
%   load-bearing, which is the defect the coordinate-meaning rule exists to prevent. The human-readable
%   strings are kept in `problem` alongside, for the JSON summary and for a reader.
%
%   `epsReqInitial` IS STORED BESIDE IT, and the distinction is the whole of it:
%   `setup.state` is a handle the loop ADVANCED, so at the end of a run `setup` no
%   longer holds the resolution the run started from. A replay must start from what
%   the run started with. That is bounded -- `rid.crit.CoverageState` is the
%   project's only intentional handle class -- but the rule generalises, and is
%   written down in the single-setup rule: any handle member of `Setup` needs an explicit initial
%   snapshot here, and a second one arriving requires the same treatment.
%
%   `qOrigin` AND `qHash` ARE READ FROM `setup.qInfo`, NEVER FROM `opts.roi.*`.
%   `opts.roi.method` and `opts.roi.seed` describe a call to
%   rid.roi.supportingPoints that may never have happened -- Q can arrive
%   supplied by the caller or loaded from a stored set -- and qInfo is
%   the one place that actually observed which. `qHash` is "" unless
%   `qOrigin == "stored"`: identity is the CONTENT HASH, not the path a set
%   happened to be loaded from, for the same reason the record stores
%   a git SHA rather than a working-tree location.
%
%   `Q` is stored only when asked for (section 7) -- and when it is not, it is
%   removed from the stored `Setup` too: the flag exists for size, and
%   while `Q` sat inside `setup` the file never got smaller while the record
%   announced that it had. Every record keeps the generator, its seed and the
%   measured fillDistanceQ as provenance whether or not the points are stored. Regeneration from the seed is NOT a reproducibility
%   mechanism -- the generating algorithm belongs to a toolbox and may change
%   between releases -- so a record without `Q` is declared non-replayable rather
%   than silently reconstructed.
rec.identity = struct( ...
    "runId",     runId(codeState, opts.io.slug), ...
    "sweepId",   "", ...
    "createdAt", environment.capturedAt, ...
    "host",      environment.host, ...
    "platform",  environment.platform, ...
    "matlabRelease", environment.matlabRelease, ...
    "toolboxes", environment.toolboxes);

rec.code = codeState;

rec.configuration = struct( ...
    "effective",     opts, ...
    "overrides",     overrides, ...
    "defaultsHash",  optionsHash(rid.defaultOptions()));

% READ FROM WHAT ACTED, NOT FROM WHAT WAS CONFIGURED. `optimizer` was
% `opts.gen.seed` -- the generator's seed under a second name -- so the record
% asserted an equality that is not a fact and stayed plausible for as long as the
% two happened to coincide. The solver carries its own.
%
% `supportingPoints` comes from `setup.qInfo`, NEVER from `opts.roi.seed`
%: that option describes a call to rid.roi.supportingPoints that may
% never have happened -- Q can arrive supplied or stored -- and qInfo is the
% one place that actually observed which. NaN when no generator ran.
rec.seeds = struct( ...
    "supportingPoints", setup.qInfo.seed, ...
    "fillDistanceProbe", opts.roi.probeSeed, ...
    "generator",        opts.gen.seed, ...
    "optimizer",        seedOf(setup.sol));

rec.problem = struct( ...
    "roi",            setup.roi, ...
    "representation", setup.rep.describe(), ...
    "normalization",  setup.nrm.describe(), ...
    "surrogate",      setup.srg.describe(), ...
    "generator",      setup.gen.describe(), ...
    "optimizer",      string(class(setup.sol)), ...
    "schedulers",     setup.sch.describe() + " | " + setup.fieldSch.describe(), ...
    "updateInterval", struct( ...
        "surrogate", intervalOf(setup.sch), ...
        "fields",    intervalOf(setup.fieldSch)), ...
    "nSup",           size(setup.Q, 1), ...
    "nDim",           size(setup.Q, 2), ...
    "fillDistanceQ",  fillDistanceQ, ...
    "fillDistanceInfo", fillInfo, ...
    "fields",         fields, ...
    "criterion",      setup.state.describe(), ...
    "epsReq",         resolutionOrNaN(setup.state), ...
    "targetLevels",   levelInfo, ...
    "storedQ",        opts.io.storeQ, ...
    "qOrigin",        setup.qInfo.origin, ...
    "qHash",          setup.qInfo.hash, ...
    "Q",              []);
if opts.io.storeQ
    rec.problem.Q = setup.Q;
end

% Its own group, because it is not provenance about the
% problem: it IS the problem, in the form `replay` consumes.
%
% "IN FULL" HAS EXACTLY ONE EXCEPTION. With
% `storeQ = false` the stored Setup has `Q` removed from it -- from `setup.Q` and
% from `setup.state` alike -- so that `rec.problem.storedQ` describes the record
% truthfully and the file actually gets smaller, which is the only thing the flag
% is for. It is stripped HERE rather than on the way to disk so that
% `result.record` makes the same statement as the record read back.
storedSetup = setup;
if ~opts.io.storeQ
    storedSetup = setup.withoutQ();
end
rec.setup = struct("setup", storedSetup, "epsReqInitial", epsReqInitial);

% THE FINAL SURROGATE, beside the history. `setup.srg` is the model the
% run BEGAN with and a Setup is immutable, so after a refit the record held no
% model corresponding to the coordinates that existed at the end: `Z` is derived
% and never stored, and invariant 3 remaps every past coordinate under
% the new model. One object closes that, with nothing left over -- which is why
% it is preferred to storing `Znorm`, a derived quantity the history-as-single-source rule keeps out of
% storage and a second source for something the record can compute.
%
% The rule, both halves: a member the loop ADVANCES in place needs its starting
% value recorded (`epsReqInitial`); a member the loop REPLACES needs its
% ending value recorded (this).
rec.outputs = struct( ...
    "u",             result.hist.u, ...
    "y",             historyOutput(result.hist), ...
    "x",             historyStates(result.hist), ...
    "tS",            result.hist.tS, ...
    "x0",            result.hist.x0, ...
    "srg",           result.srg, ...
    "eventLog",      result.eventLog, ...
    "steps",         result.steps, ...
    "Jtrace",        result.Jtrace, ...
    "J",             result.J, ...
    "status",        result.status, ...
    "nSteps",        result.nSteps, ...
    "certification", result.certification);

rec.note = opts.io.note;

end

% ---------------------------------------------------------------------------

function seed = seedOf(sol)
%SEEDOF  The seed the OPTIMIZER used, or an explicit absence.
%
%   Read from the object that acted, never copied from a setting that may not
%   have reached it. `rid.opt.GridSearch` is exhaustive and has no starts at all,
%   so it has no seed; NaN says so. Zero would name a seed the optimizer
%   never used, which is the plausible-but-false entry this exists to remove --
%   the same reason `emptyDiagnostics` leaves `deadOnArrival` NaN when the concept
%   does not apply.
seed = NaN;
if isprop(sol, "seed")
    seed = sol.seed;
end

end

% ---------------------------------------------------------------------------

function interval = intervalOf(sch)
%INTERVALOF  The update interval the SCHEDULER holds, or an absence.
%
%   Schedulers read their interval from `rid.defaultOptions` at CONSTRUCTION, so
%   an already-built scheduler keeps its own value while a run's override sits in
%   the effective options unread. The record must say which one polled, so it asks
%   the object. A policy with no interval at all reports NaN rather than a number
%   nobody chose.
interval = NaN;
if isprop(sch, "updateInterval")
    interval = sch.updateInterval;
end

end

% ---------------------------------------------------------------------------

function y = historyOutput(hist)
%HISTORYOUTPUT  The measured output, or an explicit absence.
%
%   `hist.y` is ABSENT offline and that is the only way to say it. The
%   record is a serialized snapshot rather than a history, so it stores an empty
%   matrix and `isempty` there means "no measurements were taken" -- a distinction
%   the record's own reader must not confuse with "not yet".
y = zeros(0, 0);
if isfield(hist, "y")
    y = hist.y;
end

end

% ---------------------------------------------------------------------------

function x = historyStates(hist)
%HISTORYSTATES  The measured states, or an explicit absence.
%
%   `hist.x` is absent unless the representation reads measured states, and the
%   record says so exactly as it says it for `y`: an empty matrix, whose
%   `isempty` means "no states were measured". Without this a measured-state run
%   is not replayable at all -- every horizon is anchored on `hist.x(end,:)`, so
%   a record that dropped the states would reproduce a different design from its
%   first step, with nothing left to compare it against.
x = zeros(0, 0);
if isfield(hist, "x")
    x = hist.x;
end

end

% ---------------------------------------------------------------------------

function id = runId(codeState, slug)
%RUNID  <ISO-8601 timestamp>_<short git hash>_<slug>.
stamp = string(datetime("now", "Format", "uuuuMMdd'T'HHmmss"));
id    = stamp + "_" + codeState.shortCommit + "_" + slug;

end

% ---------------------------------------------------------------------------

function h = optionsHash(opts)
%OPTIONSHASH  A stable hash of the defaults, so a record says which they were.
try
    h = string(lower(reshape(dec2hex(typecast( ...
        java.security.MessageDigest.getInstance("SHA-256").digest( ...
            uint8(char(jsonencode(opts)))), "uint8")).', 1, [])));
catch
    h = "unknown";
end

end
