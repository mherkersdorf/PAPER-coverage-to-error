function cfg = p3Config()
%P3CONFIG  The ONE configuration of P3 (Bouc-Wen, NARX mode), transcribed.
%
%   Every value marked
%   [pilot] is NaN here until the pilot delivers it and a value is chosen; the ROI bounds are set.
%   Drivers read a [pilot] value through a check that
%   errors by name while it is NaN.
%
%   The process constants (cfg.process) are DEFINED in p3BoucWenSimulate,
%   which cannot read this private file; they are repeated here for the
%   record, and a test checks that both agree.
%
%   Outputs:
%       cfg  (1x1) struct, one field per section; see below.

%% Design spaces: contiguous lags, no feedthrough u(k), input delay 1
cfg.spaces.lagOrders   = [1, 2, 3];   % m
cfg.spaces.nZ          = [2, 4, 6];   % n_z = 2m
cfg.spaces.feedthrough = false;
cfg.spaces.inputDelay  = 1;

%% Surrogate: BLA of order 3, our own estimate
cfg.surrogate.type        = "BLA";
cfg.surrogate.order       = 3;
cfg.surrogate.source      = "own estimate";
cfg.surrogate.band        = [5, 150];   % [Hz], random-phase multisines, within the ROI's u range
cfg.surrogate.file        = "experiments/p3/surrogates/bla.json";
cfg.surrogate.inputPath   = "interp750";   % 750 Hz samples, interp(u, 20) to 15 kHz: every 750 Hz signal's path
cfg.surrogate.nRealizations     = 4;     % M
cfg.surrogate.periodLength      = 8192;  % samples at 750 Hz
cfg.surrogate.nPeriods          = 6;     % 1 + 5; the first is dropped
cfg.surrogate.nTransientPeriods = 1;
cfg.surrogate.phases            = "uniform random, one RandStream(""twister"", seed) per realization";
cfg.surrogate.lineRule          = "every line l with 5 <= l fs/N <= 150 Hz (multisineSignal's band rule)";
cfg.surrogate.scaling           = "per realization: zero mean, largest amplitude inside the ROI's u range";
% The robust method (Pintelon & Schoukens 2012), by the toolbox routines of
% the PNLSS toolbox (not included). maxr, the
% weighting and the feedthrough are the toolbox's own settings (BLA.m and the
% periodic branch of @losmot/bestaffineapproximation.m); the refinement and the
% stability handling are choices of the study.
cfg.surrogate.estimation.nonparametric     = "fCovarFrf";     % G, covGML (total), covGn (noise)
cfg.surrogate.estimation.parametric        = "fLoopSubSpace"; % frequency-domain subspace, McKelvey et al. 1996
cfg.surrogate.estimation.weighting         = "covGML";        % total variance, as bestaffineapproximation
cfg.surrogate.estimation.maxr              = 20;              % BLA.m default
cfg.surrogate.estimation.directFeedthrough = true;            % BLA.m default
% LM refinement of the weighted cost through fLoopSubSpace (fLevMarqFreqSSz).
% Inf selects fLevMarqFreqSSz's own stopping rule: the cost of the last 10
% successful steps differs by less than 0.1 %, or 1000 iterations.
cfg.surrogate.estimation.lmIterations      = Inf;
% No stabilization (no fStabilize): an unstable subspace model stops A3 with an
% error, recorded.
cfg.surrogate.estimation.forceStability    = false;
cfg.surrogate.estimation.keepRule          = "the refined model if it is stable and its weighted cost " + ...
    "is lower than the subspace model's; otherwise the subspace model";
% Conversion check against the toolbox's simulation (simLMSSN).
cfg.surrogate.check.nInputs = 3;
cfg.surrogate.check.N       = 2000;    % samples per input, i.i.d. uniform in the ROI's u range
cfg.surrogate.check.relTol  = 1e-10;   % max|y_rid - y_LMSSN| / max|y_LMSSN|
% Fit on our own data: ROI realization 1 of cfg.roi.record (never test data).
cfg.surrogate.fitCheck.realization = 1;
cfg.surrogate.fitCheck.nrmse       = "rms(y - yHat) / std(y), as experiments/p2";
% Information only: the fit also under the zero-order-hold path,
% G times the zero-order-hold and interpolation responses on |f| <= this limit (interp's passband).
cfg.surrogate.fitCheck.zohConversionMaxHz = 187.5;

%% Excitation signals, sampling, simulation, noise
cfg.process.fs                = 750;        % working sampling frequency [Hz]
cfg.process.fsInt             = 15000;      % integration sampling frequency [Hz] (input upsampled by 20)
cfg.process.upsample          = 20;
cfg.process.decimationFactors = [2, 2, 5];  % decimate, 'fir', in the prime factors
cfg.process.decimationFilter  = "fir";
cfg.process.initialConditions = [0, 0, 0];  % y0, dy0, z0
cfg.process.nExtraBefore      = 20;         % zero input in front, cut
cfg.process.nExtraAfter       = 20;         % hold of u(end), cut
% A 750 Hz signal reaches 15 kHz by band-limited interpolation, interp(u, 20) with
% its default filter (replaces the zero-order hold), so that
% training, validation and BLA data take the test data's path.
cfg.process.interpolation     = struct("function", "interp", "factor", 20, "l", 4, "alpha", 0.5);
cfg.process.code              = "BoucWen_NewmarkIntegration.p (Bouc-Wen benchmark, not included)";
cfg.process.noise.distribution = "Gaussian, band-limited 0-375 Hz";
cfg.process.noise.rmsNominal   = 8e-6;      % [m] = 8e-3 mm, as the benchmark states it
cfg.process.noise.addedBy      = "the integrator, from the global stream; seeded per simulation";
cfg.process.noise.training     = true;      % training and validation data are noisy; noiseless is not allowed
cfg.process.noise.test         = false;     % the test data are noiseless (as provided)
cfg.process.inputNoise         = false;     % the input is noiseless

cfg.signals.band      = [5, 150];   % [Hz] multisine and chirp; the test sweep's 20-50 Hz lies inside
cfg.signals.N         = 5000;       % per realization; shorter lengths are nested prefixes (RHC, APRBS)
cfg.signals.inputRange = "roi.u";   % every signal, validation included, stays in the ROI's u range
cfg.signals.rateLimit = "none";
% Output limit: every signal, validation
% included, keeps the BLA's output from rest, without noise, in the ROI's y range;
% each standard signal gets the largest factor meeting both limits.
%   One factor per signal: per realization for the APRBS, set at
% N = 5000 and kept by its prefixes; per length for multisine and chirp.
cfg.signals.outputRange = "roi.y";

% 40 realizations per seeded method are generated and trained; 1-30 are the
% sample, fixed in advance, and 31-40 enter no reported number.
cfg.signals.realizationsSample = 1:30;
cfg.signals.rhc.realizationsPerLagOrder = 40;   % 120 in all

cfg.signals.aprbs.realizations      = 40;
cfg.signals.aprbs.generator         = "experiments/signals/aprbsSignal.m";
cfg.signals.aprbs.construction      = "expD";   % mean hold 2T
cfg.signals.aprbs.levelConstruction = "rand";   % levels uniform in the input range
% APRBS.m's T = tau/2, tau = -1/ln|z| of the BLA pole of largest |z|, in samples
% at 750 Hz, no rounding; mean hold 2T = tau. Pilot A3:
% p3Bla_20261029T074058, results.json dominantPole.aprbsT (= bla.json aprbsT);
% a check confirms that both agree bit for bit.
cfg.signals.aprbs.T = 18.376993524073175;
% The realization is scaled as a whole about zero, at most by 1: its levels stay uniform over the scaled range.
cfg.signals.aprbs.maxScale = 1;

cfg.signals.multisine.generator    = "experiments/signals/multisineSignal.m";
cfg.signals.multisine.phases       = "Schroeder";
cfg.signals.multisine.perLength    = 1;   % one full period of N samples; deterministic, no seed
cfg.signals.multisine.scaling      = "zero mean, largest amplitude inside the input range";

cfg.signals.chirp.generator = "experiments/signals/chirpSignal.m";
cfg.signals.chirp.sweep     = [5, 150];   % [Hz], linear, over N samples
cfg.signals.chirp.perLength = 1;          % deterministic, no seed
cfg.signals.chirp.scaling   = "zero mean, largest amplitude inside the input range";

% Pilot A5: one signal per standard method at these lengths; the
% APRBS of the shorter length is the prefix of the longer one (nested prefixes).
cfg.amplitude.lengths     = [500, 5000];
cfg.amplitude.lagOrders   = [1, 2, 3];
cfg.amplitude.signalOrder = ["aprbs", "multisine", "chirp"];   % noise seeds in this order, short first

%% ROI and space-filling hyperparameters
cfg.roi.nRealizations = 100;
cfg.roi.class         = "random phase, 5-150 Hz, RMS 50 N, one steady-state period of 8192 samples";
cfg.roi.band          = [5, 150];   % [Hz]
cfg.roi.rms           = 50;         % [N]
cfg.roi.periodLength  = 8192;
cfg.roi.nPeriods      = 5;          % plus the filter period (the benchmark's example integration)
cfg.roi.rule          = "median over realizations of the per-realization minimum and maximum, " + ...
                        "widened by 5 % of the span on each side; every lag coordinate shares its range";
cfg.roi.widening      = 0.05;
% From the record p3Roi_20261028T150220, results.json
% u.bounds.lo/hi and y.bounds.lo/hi, at full precision. Every lag of u (of y)
% shares its range.
cfg.roi.record = "p3Roi_20261028T150220";
cfg.roi.uMin = -204.08499516690307;     % [N]
cfg.roi.uMax = 207.14615989178017;      % [N]
cfg.roi.yMin = -0.0023574812693029994;  % [m]
cfg.roi.yMax = 0.0023417799306230443;   % [m]
% The RHC design's admissible region (rid.gen.Constraints): the same
% limits as every other method. Each lag
% of u gets the u range, enforced hard (lb/ub); each lag of y the y range,
% enforced softly, as a penalty on the BLA-predicted output in the objective.
cfg.admissible.u = [cfg.roi.uMin, cfg.roi.uMax];   % [N]
cfg.admissible.y = [cfg.roi.yMin, cfg.roi.yMax];   % [m]
% The values the multisine test set would give, used ONLY for the reported
% difference; never as bounds.
cfg.roi.testDataReference.u = [-220.98955652282007, 186.00413825851584];        % [N]
cfg.roi.testDataReference.y = [-0.002542998235677952, 0.00225213989637634];     % [m]

cfg.spaceFilling.nZ            = [2, 4, 6];
cfg.spaceFilling.eps0          = [0.010, 0.119, 0.296];   % sqrt(n_z) / (2 N^(1/n_z)), N = 5000
% Supporting points: RULE A, the largest cell-centred full grid
% with at most 22 000 points. Decided from the SMOKE RUNS
% (20 steps per lag order, L = 20,
% laptop, 3 workers), not from the pilot's 1000-step runs: t_step / N_sup was about
% flat in n_z (fitted exponent a ~ 0.15 < 0.5, the threshold used). The pilot
% still reports t_step / N_sup and a, as information. Only n_z = 2 changes against
% rule B (N_sup n_z <= 100 000; 223^2 = 49 729 points, fill distance 0.0032);
% the grid columns below are rule A's.
cfg.spaceFilling.pointsPerAxis = [148, 12, 5];
cfg.spaceFilling.nSup          = [21904, 20736, 15625];
cfg.spaceFilling.epsQ          = [0.0048, 0.083, 0.245];  % sqrt(n_z) / (2 n)
cfg.spaceFilling.r             = [0.42, 0.60, 0.73];      % 0.3 sqrt(n_z)
cfg.spaceFilling.psiMax        = [41.4, 4.05, 1.48];      % r / eps0 - 1
cfg.spaceFilling.supportRule   = "A";                     % N_sup <= 22 000 (decided above)
cfg.spaceFilling.supportBudget = 100000;                  % rule B, no longer used
cfg.spaceFilling.ruleAMaxPoints = 22000;                  % rule A
cfg.spaceFilling.rFactor       = 0.3;
cfg.spaceFilling.mode          = "budget";
cfg.spaceFilling.budgetN       = 5000;
% Horizon L = 40 for every RHC run. Pilot
% B's declared rule gave L = 20 (median epsilon ratio L20/L40 0.962; records
% p3PilotRhc_20261030T075830 and p3PilotRhcDiag_20261030T082354); L = 40 because the
% pilot cannot test whether a longer horizon reaches hard-to-reach regions better
% once the easy ones are covered. The generation job sets read it here.
cfg.spaceFilling.horizon = 40;

%% Training-data lengths
cfg.lengths = [500, 564, 637, 719, 812, 916, 1035, 1168, 1318, 1488, 1680, ...
               1896, 2141, 2416, 2728, 3079, 3476, 3924, 4429, 5000];

%% LMSSN
cfg.lmssn.code                 = "LMSSN (not included)";
cfg.lmssn.localModelCap        = NaN;    % [pilot] the cap method's C, used only if a pilot sends training to it
cfg.lmssn.trainingsPerSignal   = 1;      % one initialization
cfg.lmssn.explicitValidation   = true;   % no internal data splitting: datasets(train, val, empty test)
% The model and the training setting.
% Field names follow the vendor's (infos, optimizationOptions); the training code sets
% exactly these and leaves every other vendor property at its default.
cfg.lmssn.n_x                  = 3;
cfg.lmssn.validity             = "lolimot";
cfg.lmssn.LMN_equations        = {1:3};  % one network over the three state equations; the output equation is affine
cfg.lmssn.optimizer            = "QN";   % fminunc, quasi-Newton
cfg.lmssn.specifyObjectiveGradient = true;
cfg.lmssn.maxIterations        = 2e7;
% Evaluation limit per fminunc call: 5 x the largest
% function-evaluation count of any child optimization of the one-pass trainings
% (121800f, N = 500, 5000, 4429, 3924, 3476; 592 trainings, 14 795 children),
% rounded up to two significant digits: 5 x 150 104 = 750 520 -> 760 000 (record
% p3EvalLimit_20261105T125917; median 447, p99 16 649, second 147 474). The three
% largest were cut by the then 30 h soft time guard (lower bounds); the largest
% that ended on its own is 126 132. The vendor's fminusub hands the line search
% the remaining budget, so the limit also ends a line search.
cfg.lmssn.maxFunctionEvaluations = 760000;
% Training starts like validation and test:
% the training data get the same zero-input pre-roll (cfg.lmssn.preRoll samples,
% output 0, weight 0) and x_0 is FIXED -- not optimized -- at the initial model's
% zero state, so the model settles to its own rest before the scored training
% data. Before: x_0 estimated, which in C1 started the models 0.16-0.55 mm off
% their rest (p3PilotLmssnC1preRoll_20261030T091753). transform_x_0train is the
% vendor's default, true: each state transform carries x_0train
% into the new coordinates. With false it stayed at 0 while the coordinates
% moved, and C1 broke down from split 3 (p3PilotLmssnC1auto_20261031T002002: the
% split-2 model's training NRMSE 0.0188 -> 51.7 under the transform).
cfg.lmssn.optimize_x_0train    = false;
cfg.lmssn.transform_x_0train   = true;
cfg.lmssn.trainingPreRoll      = true;
% Relative tolerance of the training-start checks: the
% state transform must leave each model's training loss unchanged, and split 1's
% candidates must start at the parent's final loss.
cfg.lmssn.continuityRelTol     = 1e-9;
cfg.lmssn.yWeights             = "ones"; % training samples and the scored validation samples weigh 1; the pre-roll 0
cfg.lmssn.SNR                  = 40;     % [dB], sets estimatebestmodel's tolerance sqrt(10^(-SNR/10)) / 4
cfg.lmssn.terminationcriterion = "val";  % automatic stop, selection by validation, no cap
cfg.lmssn.term_deteriorations  = 3;
% Validation pre-roll: zero input of this
% many samples, weight 0, in front of every validation signal, so the model
% settles from its zero state (x_0val = 0 in the toolbox's coordinates, not its
% rest after `transform`) to its own rest before the scored part. 750 samples
% was the starting point; record p3PilotLmssnC1preRoll_20261030T091559 found
% that split 23 of C1 needs 1097 (transient 1.5e-5 m left after 750, noise level
% 8e-6 m), so 1500 (2 s at 750 Hz, 37 % above the longest need); checked again
% on the same models at 1500.
cfg.lmssn.preRoll              = 1500;
cfg.lmssn.preRollTolerance     = cfg.process.noise.rmsNominal;   % [m], 8e-6, the transient left at its end
% Guards on the automatic stop: beyond safetyCap
% local models, or past maxTrainingTime [s], every further optimization ends at
% once, so the 'val' rule ends the training within a few cheap splits and keeps
% the best optimized model (verified on pilot data).
% Values after the pilot rerun (selected 10 local
% models): cap 30. Time guard 30 h: first 10 h, sized to the
% deadline (training Saturday 11:00 to Thursday morning, 5 cluster jobs at a time)
% so that it only stops runaways; raised because the training-start fix (decision
% 3) made a training 5.6x slower (p3PilotLmssnC1auto_20261031T153016: 13 370 s
% on the laptop at N = 5000, about 7.4 h on the cluster at factor 2).
cfg.lmssn.safetyCap            = 30;
% SOFT TIME GUARD OFF: Inf. It acted only at the optimizer's
% calls and let a training overrun it (rhcM1_k24, N = 500); the hard deadline of
% the training client (cfg.training.hardDeadline, 108 000 s = the 30 h above) is
% the only time limit, and every time-limited training ends by the deadline rule.
cfg.lmssn.maxTrainingTime      = Inf;   % [s]
cfg.lmssn.useParallel          = false;  % one computational thread per training
cfg.lmssn.nThreads             = 1;
% One fixed training seed for every training, set before
% each training; 70001 lies in no reserved block (section 4). The training draws
% no random number if the starting model comes from the data alone; the training code
% records whether it did.
cfg.lmssn.trainingSeed         = 70001;

%% Validation: of each method's own kind, 1200 samples
% Per realization: realization k of RHC m and of the APRBS has
% its own validation signal, the same for every length of that realization; an
% RHC validation signal is one run of its lag order's configuration at N = 1200,
% an APRBS one is scaled like its training realization (own factor, cap 1).
% Multisine and chirp keep one validation signal each.
cfg.validation.N     = 1200;
cfg.validation.kinds = ["RHC, per lag order and realization, own seed", "APRBS, per realization, own seed", ...
                        "random-phase multisine 5-150 Hz, own seed", "chirp 150 -> 5 Hz"];

%% Testing: both official test sets, simulation, evaluation only
cfg.test.sets = ["random-phase multisine: one steady-state period of 8192 samples, 5-150 Hz, RMS 50 N", ...
                 "sine sweep: 153 000 samples, 20-50 Hz, amplitude 40 N, from zero initial conditions"];
cfg.test.mode = "simulation";   % the model gets the test input only
cfg.test.use  = "evaluation only; opened by the test evaluation alone";
% The sine sweep is scored after the same zero-input pre-roll as the validation
% data; the multisine test is scored on the
% second of two simulated periods.
cfg.test.sweepPreRoll = cfg.lmssn.preRoll;

% The test sets' lengths as the benchmark states them; the training driver checks
% the loaded test sets against them before any training starts.
cfg.test.nMultisine = 8192;
cfg.test.nSweep     = 153000;

%% Training jobs
% Per length a sample job of 122 trainings: RHC m = 1, 2, 3 and the
% APRBS, realizations 1-30 (first N samples of the 5000-sample data set, the
% realization's own validation data set), multisine and chirp (their data set of
% length N, their one validation data set each); and a spare job of realizations
% 31-40 (40 trainings), submitted later if at all. The data-set records the cluster jobs read.
cfg.training.dataRecords = struct("standard", "p3GenStandard_20261030T104610", ...
                                  "rhc", "p3DataSetsRhc_20261031T094032");
cfg.training.seededMethods = ["rhcM1", "rhcM2", "rhcM3", "aprbs"];
cfg.training.singleMethods = ["multisine", "chirp"];
cfg.training.realizations  = 1:40;
% Panel (d)'s iterations: the iteration at which the training loss has made this
% share of its reduction from the start to the selected model.
cfg.training.convergenceShare = 0.95;
% A test simulation DIVERGES if its peak |yHat| exceeds this factor times the test
% set's peak |y| (a declared setting; it
% only counts failures -- the RMSE enters as it is).
cfg.training.divergenceFactor = 10;
% TWO STAGES: the first five lengths (500, 5000, 4429, 3924,
% 3476) are trained in one pass at this commit; the others in two stages
% (initialization, then training) at a later commit, so that a training
% job holds no Control System Toolbox licence. the collection accepts exactly these
% two commits for the sample, with a passed equivalence record.
cfg.training.onePassCommit = "121800fbe9b7eb57417f287e9284d58bad709e82";
% HARD DEADLINE [s]: the training client stops a training
% that has run this long on its worker (within the client's 60 s poll), no grace,
% and evaluates it from its last-saved lmssn.mat with the vendor's normal
% selection (deadlineHit). The only time limit: the soft
% time guard is off (cfg.lmssn.maxTrainingTime = Inf).
cfg.training.hardDeadline = 108000;

%% Measures
cfg.measures.error  = "RMSE of the simulated output per test set (e_RMSt, eq. 6) and their mean";
cfg.measures.panelD = "local models or iterations to 95 % of the final training-loss reduction; decided after training";
cfg.measures.logged = ["number of local models", "training-loss history", "iterations", ...
                       "termination reason", "final training and validation loss", "training time"];
cfg.measures.maximumError = false;

%% Seeds
cfg.seeds.roiMultisine = 60001:60100;   % A2, realization k -> 60000 + k
cfg.seeds.blaMultisine = 62001:62004;   % A3
cfg.seeds.blaCheck     = 62101:62103;   % A3, the three inputs of the conversion check; not data
cfg.seeds.aprbsA5      = 40001;         % A5
cfg.seeds.noise.A1     = 50001:50020;
cfg.seeds.noise.A2     = 50101:50200;   % realization k -> 50100 + k
cfg.seeds.noise.sweep  = 50201;         % A2's own sine sweep
cfg.seeds.noise.A3     = 50301:50304;
cfg.seeds.noise.A5     = 50401:50406;
cfg.seeds.noise.A5b    = 50501:50506;   % A5b, the signals under both limits
% Pilot RHC runs 61001 up: opts.gen.seed of the pilot RHC runs, by lag order m = 1, 2, 3, and of the
% horizon runs (m = 2; 61011 is also the m = 2 runtime run). The smoke runs before
% the job take 61091-61093.
cfg.seeds.pilotRhc.runtime = [61001, 61011, 61021];
cfg.seeds.pilotRhc.horizon = 61011:61013;
cfg.seeds.pilotRhc.smoke   = 61091:61093;
% One noise seed per simulation of the true process in the pilot evaluation
% (p3PilotRhc), "50601 up"; the record lists which run took which.
cfg.seeds.noise.pilotB = 50601:50611;
% Pilot trainings: APRBS realizations 1 and 2 and the
% validation signals; one noise seed per
% simulation, "50701 up", in the order of cfg.pilotC.noiseOrder.
cfg.seeds.aprbs           = 40001:40040;   % realization k -> 40000 + k
% 49004 is the pilot's APRBS validation seed only: the main
% runs validate per realization (aprbsRealizations below, which
% drops 49001-49004 from the main runs). 49005 is the pilot's and the main runs'.
cfg.seeds.validation.aprbs     = 49004;
cfg.seeds.validation.multisine = 49005;
cfg.seeds.noise.pilotC    = 50701:50707;
% Signal generation and data sets. opts.gen.seed of the RHC realizations, row m, column k =
% realization k; likewise their validation signals (one per realization). Every signal is simulated ONCE (p3BoucWenSimulate) with
% its own noise seed: RHC and APRBS at their full 5000 samples (the training
% takes prefixes), multisine and chirp per length in the order of cfg.lengths,
% the validation signals at 1200.
cfg.seeds.rhc              = [10001:10040; 20001:20040; 30001:30040];
cfg.seeds.validation.rhcRealizations   = [49101:49140; 49201:49240; 49301:49340];
cfg.seeds.validation.aprbsRealizations = 49401:49440;
cfg.seeds.noise.rhc        = [51001:51040; 52001:52040; 53001:53040];
cfg.seeds.noise.aprbs      = 54001:54040;
cfg.seeds.noise.multisine  = 55001:55020;
cfg.seeds.noise.chirp      = 56001:56020;
cfg.seeds.noise.validation.rhc       = [57101:57140; 57201:57240; 57301:57340];
cfg.seeds.noise.validation.aprbs     = 57401:57440;
cfg.seeds.noise.validation.multisine = 57005;
cfg.seeds.noise.validation.chirp     = 57006;

%% Pilot training
cfg.pilotC.N0          = 500;     % C0, the harness
cfg.pilotC.N1          = 5000;    % C1, automatic stop at full length
cfg.pilotC.capC        = 30;      % C1's reference run of the cap method
cfg.pilotC.poolWorkers = 4;       % C0 (b)
cfg.pilotC.otherGlobalRng = struct("seed", 4242, "generator", "philox");   % C0 (c)
cfg.pilotC.maxWallTimeS = 12 * 3600;   % C1 is stopped from outside beyond this; the vendor saves its history after every split
cfg.pilotC.noiseOrder  = ["aprbs 40001, N = 5000", "aprbs 40002, N = 5000", "multisine, N = 500", ...
                          "chirp, N = 500", "validation aprbs", "validation multisine", "validation chirp"];
cfg.validation.chirpSweep = [150, 5];   % [Hz], downwards (section 8)
end
