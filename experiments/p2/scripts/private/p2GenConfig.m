function cfg = p2GenConfig(arm)
%P2GENCONFIG  The committed configuration of one design arm of A and B (main runs, N = 8000).
%
%   One configuration per arm. Fixed in advance, not tuned.
%
%     "i"    rid.rep.StateSpace (simulated past) + S0, coverage criterion
%     "ii"   rid.rep.StateSpace (simulated past) + S1, coverage criterion
%     "iii"  rid.rep.Narx(2, 2, nk 1, levels, simulated past) + S1, coverage
%            criterion, the 4-D space [u(k-1), u(k-2), y(k-1), y(k-2)]
%
%   RESOLUTION AND RADIUS SCALE WITH THE DIAGONAL: epsReq0 = 0.043 sqrt(d/3) and
%   r = 0.2 sqrt(d/3), so psiMax = r/epsReq0 - 1 is the same in every space. In
%   4-D that is epsReq0 = 0.04965 and r = 0.2309 (rounded in the plan to 0.050
%   and 0.231). The supporting points have the ABSOLUTE fill distance 0.05:
%   18 per axis in 3-D (sqrt(3)/36 = 0.048), 20 per axis in 4-D (0.050).
%
%   Inputs:
%       arm  (1x1) string  "i" | "ii" | "iii"
%   Outputs:
%       cfg  (1x1) struct  .arm .label .rep ("stateSpace" | "narx") .surrogate
%                          ("S0" | "S1") .criterion ("coverage") .roi
%                          .admissible .mode .nSteps .nHorizon .nDim .epsReq0
%                          .outerRadius .psiMax .nPerAxis .nSup .roiSeed
%                          .maxWallTime (s) .slugPrefix .note
arguments
    arm (1,1) string {mustBeMember(arm, ["i", "ii", "iii"])}
end

%% ========================================================================
%  SHARED BY EVERY ARM
%  ========================================================================
cfg.arm         = arm;
cfg.roi         = struct("u", [-1, 1], "x1", [-0.9, 0.9], "x2", [-1.2, 1.2]);
cfg.admissible  = struct("u", [-1, 1], "x1", [-1.42, 1.42], "x2", [-2.00, 2.00]);
cfg.mode        = "budget";
cfg.nSteps      = 8000;                         % N
cfg.nHorizon    = 17;                           % L = T
cfg.roiSeed     = 1;                            % opts.roi.seed, fixed; seeds 1-100 vary gen
cfg.criterion   = "coverage";

epsReq0In3d     = 0.043;                        % kappa = 1 at N = 8000 in 3-D
outerRadiusIn3d = 0.2;
fillDistanceQ   = 0.05;                         % absolute, every space

%% ========================================================================
%  PER ARM
%  ========================================================================
switch arm
    case "i"
        cfg.label     = "stateSpace + S0";
        cfg.rep       = "stateSpace";
        cfg.surrogate = "S0";
        cfg.nDim      = 3;
        cfg.nPerAxis  = 18;
        cfg.maxWallTime = 234 * 3600;           % two runs per worker in a 20-day job
        cfg.slugPrefix  = "p2ssS0";
    case "ii"
        cfg.label     = "stateSpace + S1";
        cfg.rep       = "stateSpace";
        cfg.surrogate = "S1";
        cfg.nDim      = 3;
        cfg.nPerAxis  = 18;
        cfg.maxWallTime = 234 * 3600;
        cfg.slugPrefix  = "p2ssS1";
    case "iii"
        cfg.label     = "NARX(2,2) + S1";
        cfg.rep       = "narx";
        cfg.surrogate = "S1";
        cfg.nDim      = 4;
        cfg.nPerAxis  = 20;
        cfg.maxWallTime = 468 * 3600;           % one run per worker in a 20-day job
        cfg.slugPrefix  = "p2narxS1";
end

%% ========================================================================
%  DERIVED
%  ========================================================================
cfg.epsReq0     = epsReq0In3d * sqrt(cfg.nDim / 3);
cfg.outerRadius = outerRadiusIn3d * sqrt(cfg.nDim / 3);
cfg.psiMax      = cfg.outerRadius / cfg.epsReq0 - 1;
cfg.nSup        = cfg.nPerAxis^cfg.nDim;

% The cell-centred grid's fill distance is sqrt(d)/(2 nPerAxis); nPerAxis is the
% smallest count that reaches the absolute 0.05.
gridFill = sqrt(cfg.nDim) / (2 * cfg.nPerAxis);
if gridFill > fillDistanceQ + 1e-12 || sqrt(cfg.nDim) / (2 * (cfg.nPerAxis - 1)) <= fillDistanceQ
    error("experiments:p2GenConfig:nPerAxis", ...
        "%d points per axis give fill distance %.4f; the smallest count reaching %.2f is wanted.", ...
        cfg.nPerAxis, gridFill, fillDistanceQ);
end

cfg.note = sprintf("Main run, arm (%s): %s, coverage criterion. epsReq0 = %.5g = " + ...
    "0.043 sqrt(%d/3), r = %.5g = 0.2 sqrt(%d/3), psiMax = %.4f; %d supporting " + ...
    "points (%d per axis, fill distance %.4f); N = %d, L = %d, budget mode, start " + ...
    "at rest (u0 = 0).", arm, cfg.label, ...
    cfg.epsReq0, cfg.nDim, cfg.outerRadius, cfg.nDim, cfg.psiMax, cfg.nSup, ...
    cfg.nPerAxis, gridFill, cfg.nSteps, cfg.nHorizon);
end
