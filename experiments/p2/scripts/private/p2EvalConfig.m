function cfg = p2EvalConfig()
%P2EVALCONFIG  The settled configuration of the stage (ii) evaluation (one place, every value sourced).
%
%   Nothing here is tuned.
%
%   The field names of the test signal, grid, cubes and LMNtool settings are
%   those of lolimotPilotConfig, so the pilot's helpers read this struct unchanged.
%
%   Outputs:
%       cfg  (1x1) struct
%           .roiLow .roiHigh  (1x3) double  ROI box of [u, x1, x2]
%           .x0               (1x2) double  rest
%           .stateScale       (1x2) double  ROI half-widths of x1, x2: errors are
%                                           divided by them (normalized state units)
%           .families         (1x6) string  "i" "ii" "iii" "aprbs" "multisine" "chirp"
%           .jobFamilies      (1x4) string  one cluster job each: "i" ... "iii", "baselines"
%           .sizes            (1x80) double p2EvaluationSizes
%           .seeds            (1x100) double paired seeds 1-100
%           .nTotal           (1x1) double  the runs' target length, 8000
%           .generation       (1x1) struct  per arm: .commit (full hash) and
%                                           .slugPrefix of the evaluated main runs
%           .baselineRecord .baselineRunId  the baselines' experiment record and the
%                                           runId inside its exports
%           .gp               (1x1) struct  .lengthScale .sensitivityLengthScales
%                                           .sensitivityN .jitterRel .chunkSize
%           .lmn              (1x1) struct  LMNtool settings (lolimotPilotConfig's) and
%                                           .mRule (text) -- the rule for M
%           .models           (1x3) string  "gp" "lmn" "affine": evaluated at every size
%           .ball             (1x1) struct  the largest-empty-ball estimator's options
%           .gap              (1x1) struct  .binEdges (1 x 26) of the pointwise C1 bins,
%                                           .blockSize of the nearest-sample distances
%           .test .grid .cubes .divergenceBound .predictChunk .checkTol
%           .targetScaleExpected (1x2) double  the documented output normalization
%           .pilot            (1x1) struct  the LOLIMOT pilot record of the M = 200 check
cfg.roiLow     = [-1, -0.9, -1.2];                       % [u, x1, x2]
cfg.roiHigh    = [ 1,  0.9,  1.2];
cfg.x0         = [0, 0];
cfg.stateScale = (cfg.roiHigh(2:3) - cfg.roiLow(2:3)) / 2;   % 0.9, 1.2

cfg.families    = ["i", "ii", "iii", "aprbs", "multisine", "chirp"];
cfg.jobFamilies = ["i", "ii", "iii", "baselines"];
cfg.sizes       = p2EvaluationSizes();
cfg.seeds       = 1:100;
cfg.nTotal      = 8000;

% The main runs evaluated: arms (i)-(iii), generated at commit e6f662b. The slug
% prefixes are p2GenConfig's.
e6f = "e6f662bf7b699ce3bea8dcd6452a8e34323a8c28";
cfg.generation = struct( ...
    "i",   struct("commit", e6f, "slugPrefix", "p2ssS0"), ...
    "ii",  struct("commit", e6f, "slugPrefix", "p2ssS1"), ...
    "iii", struct("commit", e6f, "slugPrefix", "p2narxS1"));
cfg.baselineRecord = "p2Baselines_20261007T140042";
cfg.baselineRunId  = "p2Baselines_20261007T140040";

% The reference GP (the coverage gate's
% settings): Matern 5/2, one length scale in [-1, 1]^3, affine mean, jitter 1e-8.
% Sensitivity 0.25 and 1.0 at N = 8000 only.
cfg.gp = struct("lengthScale", 0.5, "sensitivityLengthScales", [0.25, 1.0], ...
    "sensitivityN", 8000, "jitterRel", 1e-8, "chunkSize", 4096);

% LMNtool 1.5.3, the pilot's settings (lolimotPilotConfig): one network with two outputs, premise = consequent =
% [u, x1, x2] without delay, ROI-normalized inputs with scalingMethod 'none',
% targets divided by targetScale, the LAST network of a fixed-M run,
% maxValidationDeterioration = inf, LOOCV off. M = mRule(N) (fixed
% per size, from the LOLIMOT pilot). The single affine model is
% LMNtool's maxNumberOfLM = 1 at every size.
pilot = lolimotPilotConfig();
cfg.lmn       = pilot.lmn;
cfg.lmn.mRule = "M = max(1, round(min(2.25*sqrt(N), N/10)))";
cfg.models    = ["gp", "lmn", "affine"];

% Largest empty ball and measured epsilon: rid.roi.estimateFillDistance
% in [0, 1]^3 of the true [u, x1, x2], 10^5 probes,
% the default probe method and seed, so the laptop and the job agree. blockSize
% lowers the memory of the probe check (blockSize x 10^5 doubles) and nothing
% else: 50 workers start their runs together.
defaults = rid.defaultOptions().roi;
cfg.ball = struct("nProbes", 1e5, "probeMethod", defaults.probeMethod, ...
    "seed", defaults.probeSeed, "blockSize", 256, "traceTolerance", defaults.traceTolerance);

% Gap map and pointwise C1: every
% grid point's distance to its nearest training sample, in the [0, 1]^3 units of
% epsilon and from the same samples. The 25 bins are fixed for every run, size
% and model, so they aggregate: [0, 1e-3), 23 log-spaced bins up to 1, [1, Inf).
% The edges are an implementation choice ("about 25 log-spaced bins"), not a
% sourced value. blockSize: samples per distance block (G x blockSize doubles).
% The grid's covering radius, sqrt(3)/2 * 1/(41 - 1) = 0.0217, is how far the gap
% map's maximum may lie below the true epsilon.
cfg.gap = struct("binEdges", [0, logspace(-3, 0, 24), Inf], "blockSize", 256);

% Test signal, grid and cubes: the pilot's: APRBS (seed 9601) -> Schroeder multisine -> chirp, 20 000
% samples each, one free run from rest; the 41^3 grid with faces; 5^3 cubes, a
% cube with fewer than 30 scored test points not scored. NRMSE and cube NRMSE
% normalized by std(y) over the scored samples.
cfg.test            = pilot.test;
cfg.grid            = pilot.grid;
cfg.cubes           = pilot.cubes;
cfg.divergenceBound = pilot.divergenceBound;
cfg.predictChunk    = pilot.predictChunk;
cfg.checkTol        = pilot.checkTol;
cfg.targetScaleExpected = [0.13355, 0.21508];

% The exact pilot check: history entry M = 200
% of the M = 201 network at N = 8000 against the pilot's c40 network.
cfg.pilot = struct("record", "lolimotPilot_20261009T000250", "variant", "c40", ...
    "M", 200, "N", 8000, "families", ["i", "ii", "aprbs"], "seeds", 1:10);
end
