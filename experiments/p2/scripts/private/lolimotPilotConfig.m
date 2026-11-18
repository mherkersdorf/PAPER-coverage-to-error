function cfg = lolimotPilotConfig()
%LOLIMOTPILOTCONFIG  The settled configuration of the LOLIMOT pilot (one place, every value sourced).
%
%   Nothing here is tuned.
%
%   Outputs:
%       cfg  (1x1) struct
%           .roiLow .roiHigh  (1x3) double  ROI box of [u, x1, x2]
%           .x0               (1x2) double  rest
%           .families         (1x3) string  "i" "ii" "aprbs"
%           .seeds            (1x10) double paired seeds 1-10
%           .sizes            (1x3) double  N = 500, 2000, 8000 designed inputs
%           .c                (1x3) double  c in M = max(1, floor(N/c))
%           .variants         (1x5) string  "c20" "c40" "c80" "aicc" "m1"
%           .aiccCap          (1x1) double  maxNumberOfLM of the AICc variant
%           .lmn              (1x1) struct  LMNtool settings (below)
%           .test             (1x1) struct  test signal (below)
%           .grid             (1x1) struct  .pointsPerAxis = 41 (faces included)
%           .cubes            (1x1) struct  .perAxis = 5, .minPoints = 30
%           .divergenceBound  (1x1) double  diagnostic only (see below)
%           .predictChunk     (1x1) double  rows per LMNtool prediction call
%           .checkTol         (1x1) double  lean predictor vs LMNtool, flagged above
%           .armCommit        (1x1) string  commit of the main runs
%           .baselineRecord   (1x1) string  the baselines' experiment record
%           .baselineRunId    (1x1) string  the runId inside its exports
cfg.roiLow  = [-1, -0.9, -1.2];                 % [u, x1, x2]
cfg.roiHigh = [ 1,  0.9,  1.2];
cfg.x0      = [0, 0];

cfg.families = ["i", "ii", "aprbs"];
cfg.seeds    = 1:10;
cfg.sizes    = [500, 2000, 8000];
cfg.c        = [20, 40, 80];
cfg.variants = ["c20", "c40", "c80", "aicc", "m1"];
cfg.aiccCap  = 400;

% LMNtool 1.5.3. Premise and consequent
% inputs [u, x1, x2] without delay; inputs normalized to the ROI box by us, so
% LMNtool's own scaling is off ('none' is set BEFORE the inputs are assigned:
% LMNtool scales on assignment). Fixed-M variants: maxNumberOfLM = M, early
% termination off, the LAST network taken explicitly. The AICc variant is LMNtool's own selection with its default early stop
% (maxValidationDeterioration = 2) and the cap maxNumberOfLM = 400: its
% suggestedNet is the network; no search for the overall AICc minimum. The
% remaining defaults are asserted, not set. The starting box
% of the partition is LOLIMOT's own, the bounding box of the training data, not
% the ROI (accepted). Targets are normalized by fixed constants
% (targetScale).
cfg.lmn.scalingMethod = 'none';
cfg.lmn.inputDelay    = {0, 0, 0};
cfg.lmn.defaults      = struct("smoothness", 1, "pointsPerLMFactor", 1, "splits", 1, ...
    "xRegressorDegree", 1, "xRegressorType", 'polynomial', "estimationProcedure", 'LS', ...
    "maxValidationDeterioration", 2, "LOOCV", true, ...
    "lossFunctionTermination", 'penaltyLossFunction', "complexityPenalty", 1);
% LOOCV only enters LMNtool's AICc (effective number of parameters), never the
% split choice. Off for the fixed-M variants since a check
% "loocv" showed the identical last network, bit for bit, on arm (i), seed 1,
% N = 2000, M = 50 (record lolimotPilotDryRun_20261008T131735). The AICc variant keeps it on (LMNtool's default).
cfg.lmn.loocvFixedM   = false;

% Test signal: three chunks of 20,000
% samples, generated exactly like the baselines, concatenated
% APRBS -> multisine -> chirp, one free run from rest.
cfg.test.chunkLength = 20000;
cfg.test.order       = ["aprbs", "multisine", "chirp"];
cfg.test.aprbs       = struct("seed", 9601, "range", [-1, 1], "T", 8, ...
    "construction", "expD", "levelConstruction", "linspace");

cfg.grid.pointsPerAxis = 41;
cfg.cubes.perAxis      = 5;
cfg.cubes.minPoints    = 30;

% A free run is counted as diverging at the first step whose predicted state is
% non-finite or lies beyond 10 ROI half-widths (normalized |x~| > 10). A
% diagnostic only: the run itself is not stopped or altered (no special rule
% for diverging runs), except that a non-finite state ends it.
cfg.divergenceBound = 10;

cfg.predictChunk = 10000;
cfg.checkTol     = 1e-10;
% Independent weighted-LS check: largest prediction difference on
% the training points, in normalized target units, counted as agreement.
cfg.wlsTol       = 1e-9;

cfg.armCommit      = "e6f662bf7b699ce3bea8dcd6452a8e34323a8c28";
cfg.baselineRecord = "p2Baselines_20261007T140042";
cfg.baselineRunId  = "p2Baselines_20261007T140040";
end
