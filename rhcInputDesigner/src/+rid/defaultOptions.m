function opts = defaultOptions()
%DEFAULTOPTIONS  The project's option defaults, in one place.
%
%   Every default value used anywhere in rid lives here and nowhere else, so that
%   the effective configuration of a run is one struct and a change of default is
%   one edit (the project's coding standard). Functions read their
%   name-value defaults from this struct in their `arguments` block; no function
%   carries a literal default of its own.
%
%   Sections are named after the package they configure. Only the sections of
%   implemented packages are present: an option is added here in the same change
%   that first uses it, never in advance.
%
%   Inputs:  none
%   Outputs:
%       opts  (1x1) struct  one field per configured package

% --- +roi: supporting points and the fill-distance measurement ---------------
%
% The seeds are defaults, not a reproducibility mechanism: they are always
% returned in the `info` struct of the function that used them and are stored as
% run-record provenance.

% The default method stays "lhs" deliberately: a default that changes
% every Q costs a regeneration of everything frozen, and the gain of the
% improving variants is not yet measured in this project.

opts.roi.method         = "lhs";      % "lhs", "lhsImproved", "lhsImprovedFull", "grid",
                                      % "stored"
opts.roi.supportingSet  = "";         % basename of a file under supportingSets/; required by
                                      % "stored" and refused by every other method
opts.roi.improveMaxTime = 1800;       % s allowed to the improvement; a cap that BITES makes the
                                      % design machine-dependent, so info.timeCapHit reports it
opts.roi.progress       = "auto";     % "auto" prints one liveness line for the improving methods
                                      % only; "off" is silent
opts.roi.seed           = 0;          % supporting-point generator seed
opts.roi.probeSeed      = 0;          % Monte-Carlo probe seed, independent of the above
opts.roi.probeMethod    = "sobolCorners";  % "sobolCorners", "sobol" or "uniform". Coverage decides
                                      % a maximum, so the probe set is not i.i.d. uniform,
                                      % and it is seeded with the ROI's corners. Above the
                                      % 2^nDim ceiling this falls back to "sobol" and reports it.
opts.roi.nProbes        = 10000;      % probes for the fill-distance estimate
opts.roi.traceTolerance = 0.02;       % relative climb across the quarter/half/full trace above
                                      % which the estimate is reported UNCONVERGED
opts.roi.blockSize      = 2048;       % probes per distance-matrix block; memory, not maths
opts.roi.lhsCriterion   = "maximin";  % passed through to rid.ext.latinHypercube
opts.roi.lhsIterations  = 5;
opts.roi.lhsSmooth      = true;

% --- +crit: the criterion and the demanded resolution ------------------------
%
% `epsReq0` deliberately has NO default: it is the experiment's baseline
% commitment -- the coarsest resolution demanded anywhere -- and no defensible
% default exists. Inventing one would hide the single most consequential choice in
% the design, so rid.crit.demandedResolution errors when it is absent.
%
% `beta` and `betaNoise` are DERIVED HERE from `convOrder`. Overriding
% `convOrder` at a call site therefore does not change `beta`: a caller wanting a
% different convergence order must override `beta` as well. The alternative --
% deriving beta inside demandedResolution when it is not supplied -- needs
% nargin-style defaulting, which the project bans, and would split one number's
% definition across two files.

opts.crit.pNorm     = 1;       % aggregation exponent; 1 = mean, Inf = worst case
opts.crit.psiMax    = 1;       % penalty clip; OPEN, to be revisited
opts.crit.convOrder = 2;       % convergence order r of the local error; OPEN
opts.crit.beta      = 1 / opts.crit.convOrder;
opts.crit.betaNoise = opts.crit.beta;   % OPEN, pragmatic not derived
opts.crit.mode      = "budget";% "budget" or "certification"; budget certifies nothing

% --- +rep: the coordinates coverage is measured in ---------------------------
%
% The default coordinates are [u(k), du(k), yHat(k-1)]: nk = 0 so the
% current input is represented, and the difference basis so that a rate limit
% makes the reachable set an axis-aligned box rather than a diagonal band. The
% lag orders are NOT here -- they come from the surrogate's order through
% `defaultRepresentation`, or from the caller.
%
% `pastSource` is never inferred from online/offline: an offline run must set
% "simulated" explicitly.

opts.rep.pastSource      = "measured";  % "measured" or "simulated"
opts.rep.nk              = 0;           % input delay; 0 represents u(k)
opts.rep.inputDifference = true;        % difference basis for the input block
opts.rep.nInputs         = 1;           % input channels p

% --- +surrogate: the model and its update policy -----------------------------
%
% `updateInterval` is Inf on purpose: the update MECHANISM exists from the
% start, the POLICY does not, and a default that silently refits would presuppose
% one. A run that wants adaptation constructs its scheduler with a finite one.
%
% IT IS THE SCHEDULER'S CONSTRUCTION DEFAULT, NOT A RUN-TIME OVERRIDE.
% rid.surrogate.ImmediateScheduler reads it in its `arguments` block and then holds
% its own copy, so setting it in a run's overrides reaches an already-built
% scheduler not at all. Naming it here as the construction default leaves ONE
% declaration: the run record reports the interval read from the scheduler that
% actually polled, never this number. A refusal on mismatch (the pattern) was
% rejected -- it would oblige the user to keep two things in step where deleting
% the second costs nothing.
%
% The steadyState settings serve only the BASE-CLASS fallback, which simulates a
% constant input to convergence for a model that cannot state its rest state
% analytically. rid.surrogate.LinearFirstOrder overrides it and never
% reads them.

opts.surrogate.updateInterval        = Inf;    % scheduler CONSTRUCTION default
opts.surrogate.minUpdateSamples      = 20;     % rows required before a fit is worth it
opts.surrogate.steadyStateTol        = 1e-10;  % relative state change at settling
opts.surrogate.steadyStateMaxSamples = 1e6;    % give up rather than loop forever
opts.surrogate.steadyStateBlock      = 1000;   % samples per convergence check

% --- +fields: when the refinement fields are re-estimated --------------------
%
% `updateInterval` is Inf for exactly the reason the surrogate's is: the
% MECHANISM exists from the start, the POLICY does not. With the default nothing
% is re-estimated, `epsReq` stays at `epsReq0` everywhere, and the run is the pure
% space-filling design of rid.fields.Nominal. It is likewise the scheduler's
% CONSTRUCTION default rather than a run-time override, and the record reports what
% rid.fields.ImmediateScheduler held.
%
% There is deliberately no default ESTIMATOR here. The estimators are objects, not
% values, and they travel in rid.fields.ImmediateScheduler's constructor, which
% names rid.fields.Nominal as its own default. An options field holding a class
% name would make "which estimator" a string the code branches on, which the coordinate-meaning rule
% exists to prevent.

opts.fields.updateInterval = Inf;   % scheduler CONSTRUCTION default

% --- +gen: the soft-constraint penalty and multi-start initialization --------
%
% The penalty WEIGHT is deliberately NOT here. It is derived where it is used, as
% `w = psiMax / violationScale^2`, reading `psiMax` from `opts.crit` and
% never duplicating it in this section. Storing a raw weight would decalibrate the
% penalty silently the moment `psiMax` changed -- which is the exact failure the
% `violationScale` parametrization exists to prevent.
%
% The RATE LIMIT is deliberately NOT here either. It is part of the constraint
% DECLARATION and travels as a name-value argument to rid.gen.Constraints:
% every bound is declared and none is derived, and the default comes from the
% declared admissible region, not from an options struct.
%
% `nHorizon` deliberately has NO default, on the same footing as `epsReq0`:
%   the method description makes it user-defined and fixed for the run,
% and the time constant is guidance rather than a formula, so no defensible default
% exists. It is a REQUIRED constructor argument of rid.gen.Direct, which is its
% single holder -- the loop reads `gen.nHorizon` rather than keeping a copy.
%
% `nStarts` is NOT here either, and moved out of this section: the number of local
% searches is a search-strategy setting, so it is declared where the decision
% lives, at `opts.opt.nStarts`. `rid.gen.Direct.initialGuesses` takes it as
% an argument with no default of its own -- +opt decides the count and +gen
% answers. One declaration, one decider, no second source to drift.
%
% Both values below are OPEN: interpretable defaults with no study behind
% them, and `violationScale` couples to `psiMax` through the psiMax = 1 default.

% The two AFFINE-PROBE values move here from rid.gen.affineMap's arguments block,
% and they move for DIFFERENT reasons.
%
% `affineProbeStep` is CONDITIONING and cannot change an answer: for an affine map
% the probe is exact for any h ~= 0 (that function's header), so it evaluates the
% linear part rather than approximating it -- there is no truncation error and no
% step size to tune. It is here because a record's effective options should say
% what the run used, not because the value is consequential.
%
% `affineCheckTol` decides the affine/notAffine VERDICT, and that is behaviour:
% the caller takes the analytic path only when the check passes, so this number
% chooses between the analytic gradient and finite differences. As a literal
% inside the function it was unreconstructable from a record -- a replay could not
% know which tolerance produced the tier it was comparing against.

opts.gen.affineProbeStep = 1e-3;  % probe step; conditioning only, never a verdict
opts.gen.affineCheckTol  = 1e-9;  % relative tolerance of the affineness CHECK

opts.gen.violationScale = 0.05;  % fraction of the admissible width whose SUSTAINED
                                 % violation costs psiMax; OPEN
opts.gen.tailWeight     = 0.1;   % gamma(nHorizon)/gamma(1): a violation at the far
                                 % end of the horizon costs this share of one now
opts.gen.seed           = 0;     % random-walk start seed, own RandStream

% --- +opt: the local search and the exhaustive oracle -------------------------
%
% THE SOLVER SETTINGS ARE THE SOLVER'S OWN DOCUMENTED DEFAULTS, taken rather than
% invented: `fmincon` at R2025b documents MaxIterations = 400,
% MaxFunctionEvaluations = 100*numberOfVariables, OptimalityTolerance = 1e-6,
% StepTolerance = 1e-6 (sqp), ConstraintTolerance = 1e-6 and
% FiniteDifferenceStepSize = sqrt(eps). Citing a source is not inventing a number,
% and reproducing them here rather than leaving them unset means the effective
% options struct of a run says what the search actually did -- which is the whole
% purpose of storing the EFFECTIVE options in the run record.
%
% `maxEvaluationsPerVariable` is per variable because that is how the solver
% expresses its own default. The decision vector is nHorizon*nInputs long and the
% horizon is a run-level choice, so an absolute cap would silently tighten as the
% horizon grew.
%
% THE ALGORITHM IS THE ONE EXCEPTION and it is a project choice, not the solver's
% default (which is interior-point): the soft-constraint penalty is deliberately
% C1 and not C2, and the method description records that SQP
% tolerates exactly that.
%
% `nStarts` lives here rather than under +gen because it is a decision about the
% SEARCH, not about the parametrization: +opt chooses how many local
% searches to run and asks the generator for that many starting points. It is the
% one value in this section that is neither a solver default nor a grid setting.
% The method description explains what the four rows are; it no
% longer declares how many there are.
%
% `constraintTolerance` is used for three things that must agree: the solver's own
% feasibility tolerance, the guard that rejects an infeasible starting point, and
% the test for whether the optimum ended against a bound. "Active" has to mean the
% same to the report as it did to the search.
%
% The two grid values are CHOSEN, not sourced. They configure the correctness
% oracle of the method description and never the default solver, so
% they set how tiny the deliberately-tiny oracle problem may be. `gridLevels = 5`
% matches the coarse five-level grid that section 4 uses in its own 5^20 estimate;
% `maxGridPoints` refuses anything an exhaustive search should not be asked to do.
% They join violationScale and tailWeight as interpretable values with no study
% behind them.

opts.opt.nStarts                   = 4;        % local searches per step
opts.opt.algorithm                 = "sqp";    % project choice; see above
opts.opt.maxIterations             = 400;      % fmincon documented default
opts.opt.maxEvaluationsPerVariable = 100;      % fmincon: 100*numberOfVariables
opts.opt.optimalityTolerance       = 1e-6;     % fmincon documented default
opts.opt.stepTolerance             = 1e-6;     % fmincon documented default (sqp)
opts.opt.constraintTolerance       = 1e-6;     % fmincon documented default
opts.opt.gridLevels                = 5;        % CHOSEN, oracle only
opts.opt.maxGridPoints             = 1e5;      % CHOSEN, oracle only

% --- +rhc: the loop and its stopping rules ----------------------------------
%
% The table of the method description, and nowhere else. Every value there is
% a mechanism switched OFF, or has no default at all, except the last two.
%
% `nSteps` IS ABSENT ON PURPOSE, on the epsReq0 precedent: the prescribed
% experiment length IS the commitment in budget mode, so no default is defensible
% and rid.rhc.run refuses a budget run that omits it.
%
% THE EXPERIMENT'S INITIAL CONDITION IS NOT HERE, and used to be. `u0` and `y0`
% were options; they are data, and they arrive as the `hist0` argument in the history struct's
% plain history struct, where `rid.assertHistory` validates them.
% The three refinement fields left for the same reason -- they are `(nSup x 1)`
% vectors indexed by `Q` and live in `rid.rhc.Setup`.
%
% ASSEMBLY REFUSES A CONFIGURATION IN WHICH NO STOPPING RULE IS FINITE, rather
% than inventing a cap: a certification run may never certify, since the remaining
% supporting points can be dynamically unreachable, and a default cap would be a
% number nobody chose.
%
% `nTargetLevels` and `nStaticProbes` configure the guided initialization of
% the surrogate-guided initialization -- how many output levels the target grid has, and how finely the
% surrogate's static map is probed to find the input reaching each. They are the
% only two numbers in this section nobody can defend: CHOSEN, not sourced, joining
% violationScale, tailWeight and gridLevels.

opts.rhc.maxSteps             = Inf;   % guard: steps
opts.rhc.maxWallTime          = Inf;   % guard: seconds
opts.rhc.stallSteps           = Inf;   % guard: steps without improvement in J
opts.rhc.stallTolerance       = 0;     % what counts as an improvement
opts.rhc.assertMonotonicity   = true;  % invariant 3, between updates only
opts.rhc.resetOnCertification = true;  % budget mode only
opts.rhc.nTargetLevels        = 9;     % CHOSEN, the surrogate-guided initialization initialization
opts.rhc.nStaticProbes        = 201;   % CHOSEN, the surrogate-guided initialization initialization

% --- +io: what is recorded, and where ---------------------------------------
%
% The table of the method description, and nowhere else.
%
% `saveRun` is true because every run produces a record and writing it is on by
% default (section 2 rule 1); false exists for tests and is never used for an
% experiment. `storeQ` is false HERE because Q is only sometimes worth its
% size, and a record without it is declared NON-REPLAYABLE rather than
% silently regenerated (section 7) -- canaries override it and always store Q.
% `rid.rhc.run` overrides this CONSTRUCTION default at run time whenever
% `setup.qInfo.origin` is not "generated": a supplied or
% reused Q has nowhere else to be regenerated from, so section 7's promise
% "always for a user-supplied Q" is enforced there rather than left to
% whoever forgets to set this flag. An EXPLICIT override of this option, in
% either direction, is never touched. `requireCleanTree` is honoured at RUN
% START, not at save time: a refusal is worth having before eight hours of
% experiment and worth nothing after.
%
% `runsDir` is read by `rid.rhc.run`, which passes it to `saveRun`; `saveRun` and
% `loadRun` fall back to this default only when called on their own, with no run
% to take an effective value from. Without that hand-off an overridden `runsDir`
% would be recorded in the effective options and then ignored -- an option that
% exists and does nothing, which is worse than one that does not exist.

opts.io.saveRun          = true;
opts.io.storeQ           = false;
opts.io.requireCleanTree = false;   % true for runs destined for a publication
opts.io.runsDir          = "runs";  % relative to the project root
opts.io.slug             = "run";   % the trailing part of runId
opts.io.note             = "";      % free text stored in the record

end
