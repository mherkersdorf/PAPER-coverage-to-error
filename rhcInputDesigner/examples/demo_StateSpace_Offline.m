close all
clearvars
clc

% Simple introductory demo of the RHC-based excitation design algorithm.
% The example uses the joint input--state-space representation and an
% order-n linear state-space surrogate without an underlying process, i.e., it
% demonstrates offline excitation design. The demanded resolution stays fixed
% in this example, i.e., pure space-filling is pursued.

%% ========================================================================
%  CONFIGURATION
%  ========================================================================
% Reproducibility and simulation settings. Two INDEPENDENT seeds, not one:
% seedROI controls rid.roi.supportingPoints (Q's own randomness, when the
% method uses any -- "grid" below does not) and seedMultiStart controls the
% RHC's per-step multi-start initial guesses (rid.gen.Direct.initialGuesses).
% They may be set equal, or varied one at a time to isolate which source of
% randomness a change in behaviour comes from.
cfg.seedROI        = 1;       % Supporting-point generator seed, redundant for cfg.supMethod = "grid".
cfg.seedMultiStart = 45;      % Per-step multi-start initial-guess seed
cfg.N    = 100;               % Number of samples available for the experiment
cfg.mode = "budget";          % Stop the algorithm when the data budget is reached
% RHC prediction horizon
cfg.nHorizon = 20;            % Number of future input samples optimized at each step
% Required resolution and coverage criterion
cfg.epsReq0 = 0.1;            % Demanded spatial resolution. Larger than the NARX
                              % demo's 0.03: the design space here is
                              % three-dimensional ([u1_k, x1_k, x2_k]) rather than
                              % two-dimensional, and the same fill distance costs
                              % more supporting points in a higher-dimensional ROI.
cfg.psiMax  = 5;              % Maximum penalty for insufficient coverage;
                              % penalty is clipped once the nearest data sample
                              % is (psiMax + 1) times the demanded resolution
                              % away from a supporting point
% Admissible physical region of the regressor
admissibleRegion.u1_k = [-2 2];     % Admissible range of u(k)
admissibleRegion.x1_k = [-3 3];     % Admissible range of x1(k), wider than the ROI
admissibleRegion.x2_k = [-0.6 0.6]; % Admissible range of x2(k), wider than the ROI
% Supporting-point set used to evaluate coverage
cfg.nSup      = 1000;                % Number of supporting points. Larger than the
                                     % NARX demo's 2000 for the same reason as
                                     % epsReq0: three dimensions need more points
                                     % to reach a comparable fill distance.
cfg.supMethod = "grid";       % Method for generating supporting points
                                     % Options: "lhs", "lhsImproved",
                                     %          "lhsImprovedFull", "grid"
% No global rng(...) call: every randomized call below threads its own
% explicit seed and builds its own RandStream -- seedROI and
% seedMultiStart above -- so seeding the global stream would affect nothing.
%% ========================================================================
%  PROCESS AND SURROGATE
%  ========================================================================
process = [];    % No process model or process data is provided
% Order-2 linear state-space surrogate used to predict future system behavior.
% A lightly damped second-order system with unity static gain: x1 tracks u at
% steady state, x2 is its rate.
w0   = 0.2;      % Natural frequency of the second-order surrogate [rad/s]
zeta = 0.5;      % Damping ratio
Ts   = 1;        % Sampling time
Ac = [0 1; -w0^2 -2*zeta*w0];
Bc = [0; w0^2];
Cc = [1 0];
srg = rid.surrogate.LinearStateSpace.fromContinuous(Ac, Bc, Cc, Ts);
%% ========================================================================
%  REGRESSOR REPRESENTATION AND REGION OF INTEREST
%  ========================================================================
% Define the joint input--state-space regressor used as the design space for
% excitation: [u1_k, x1_k, x2_k].
rep = rid.rep.StateSpace(2, "pastSource", "simulated");
% Region of interest where good coverage is desired. Chosen to lie inside the
% region the surrogate can actually reach: probing it with stepped and
% near-resonant square-wave inputs across the admissible range above shows
% x1 reaching about [-2.65, 2.64] and x2 about [-0.51, 0.51]. The ROI below is
% comfortably inside that reachable region; the admissible region above is
% wider still.
roi.u1_k = [-1, 1];      % Target range of u(k)
roi.x1_k = [-1.5, 1.5];  % Target range of x1(k)
roi.x2_k = [-0.3, 0.3];  % Target range of x2(k)
%% ========================================================================
%  NORMALIZATION AND SUPPORTING POINTS
%  ========================================================================
% Determine the regressor-space bounds of the region of interest
[ROIbox, dimName] = rep.bounds(roi);
% Create the normalization object that maps physical to normalized coordinates
nrm = rid.roi.Normalization(ROIbox, dimName);
% Generate supporting points in the normalized region of interest
[Q, ~, qInfo] = rid.roi.supportingPoints(nrm, cfg.nSup, ...
                             "method", cfg.supMethod, ...
                             "seed",   cfg.seedROI);
% Estimate the fill distance of the supporting-point set
epsQ = rid.roi.estimateFillDistance(Q, nrm);
%% ========================================================================
%  DEMANDED RESOLUTION AND COVERAGE STATE
%  ========================================================================
% Define the fields at which the demanded resolution is evaluated
fields = rid.fields.nominalFields(cfg.nSup);
% Determine the required resolution at each supporting point
epsReq = rid.crit.demandedResolution(fields, ...
             "epsReq0",       cfg.epsReq0, ...
             "mode",          cfg.mode, ...
             "fillDistanceQ", epsQ);
% Initialize the state used to track coverage of the supporting points
state = rid.crit.CoverageState(Q, epsReq, ...
                               "psiMax", cfg.psiMax);
%% ========================================================================
%  INPUT GENERATOR AND RHC SETUP
%  ========================================================================
% Define the physical constraints for the generated input trajectory
constraints = rid.gen.Constraints(rep, admissibleRegion);
% Generate future input trajectories over the RHC prediction horizon
gen = rid.gen.Direct(constraints, cfg.nHorizon);
% Assemble all components required for the RHC excitation design
setup = rid.rhc.Setup(srg, rep, nrm, roi, Q, fields, gen, ...
            rid.opt.MultiStart("seed", cfg.seedMultiStart), ...
            state, ...
            rid.surrogate.ImmediateScheduler(), ...
            rid.fields.ImmediateScheduler(), ...
            qInfo);
%% ========================================================================
%  INITIAL CONDITION
%  ========================================================================
u0 = 0;    % Initial input value
% Initialize the past input and surrogate state required by the state-space
% model. This representation has no warm-up (nWarmup = 0): a
% point needs only u(k) and x(k), and x(1) is known. hist0.u therefore carries
% a SINGLE row, the actuator position before the run began -- repmat(u0,
% rep.nWarmup, 1) would be empty and is refused by rid.gen.Direct.
hist0 = struct("u",  u0, ...
               "x0", srg.steadyState(u0), ...
               "tS", srg.tS);
%% ========================================================================
%  RUN OPTIONS
%  ========================================================================
% Optional settings that override the default RHC configuration
defaultOverrides = struct( ...
    "crit", struct( ...
        "epsReq0", cfg.epsReq0, ...
        "psiMax",  cfg.psiMax, ...
        "mode",    cfg.mode), ...
    "rhc", struct( ...
        "nSteps",  cfg.N), ...
    ...  % roi.seed: supporting-point generator seed. cfg.supMethod = "grid" is
    ...  % deterministic and ignores it -- it would matter only for
    ...  % "lhs"/"lhsImproved*". Independent of gen.seed below (see seedROI /
    ...  % seedMultiStart in CONFIGURATION).
    "roi", struct( ...
        "method", cfg.supMethod, ...
        "seed",   cfg.seedROI), ...
    ...  % gen.seed: per-step multi-start initial-guess seed
    ...  % (rid.gen.Direct.initialGuesses, read as opts.gen.seed + kStep by
    ...  % rid.rhc.run). This is the seed that actually varies the designed
    ...  % signal -- rid.opt.MultiStart's own "seed" argument above is unused
    ...  % inside rid.rhc.run, which always supplies its own starts provider.
    "gen", struct( ...
        "seed", cfg.seedMultiStart), ...
    "io", struct( ...
        "saveRun", false, ...    % Do not save the run to disk
        "storeQ",  true));       % Store the supporting points in the result
%% ========================================================================
%  THE RUN
%  ========================================================================
% Run the receding-horizon excitation design
result = rid.rhc.run(process, setup, hist0, defaultOverrides);
% Extract the designed input signal and the corresponding surrogate states
% as well as output.
u    = result.record.outputs.u;
[yhat, ~, xTraj] = srg.simulateStates(u, result.record.outputs.x0);
%% ========================================================================
%  PLOT
%  ========================================================================
fontSize = 15;
colors = [0.000, 0.447, 0.741;     % u-color: blue
         0.850, 0.325, 0.098;      % y-color: red
         0.000, 0.200, 0.450;      % Regressor-space color: dark blue
         0.75, 0.75, 0.75];
figSignalTrajSize = [600 300];
figRSpaceSize     = [600 600];
scatterSize = 50;
showSupportingPoints = true;
suppPointSize = 25;

% Determine plot limits from the generated signals
uextremeval   = [floor(10*min(u))/10,    ceil(10*max(u))/10];
urange    = uextremeval(2) - uextremeval(1);
yhatextremeval = [floor(10*min(yhat))/10, ceil(10*max(yhat))/10];
yhatrange  = yhatextremeval(2) - yhatextremeval(1);
%% ------------------------------------------------------------------------
%  INPUT TRAJECTORY
%  ------------------------------------------------------------------------
figure('Name', sprintf('Input Trajectory'), ...
       'Color', 'w', ...
       'Position', [750 100 figSignalTrajSize]);
plot(1:length(u)-1, u(2:length(u)), ...
     'Color', colors(1,:), ...
     'LineWidth', 1.5);
xlabel('Time Index $k$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);
ylabel('$u(k)$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);
legend('$u(k)$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize, ...
       'Location', 'best');
ylim([uextremeval(1)-0.05*urange, uextremeval(2)+0.05*urange]);
yticks([uextremeval(1), uextremeval(2)]);
set(gca, ...
    'FontSize', fontSize, ...
    'TickLabelInterpreter', 'latex');
box on;
grid on;
%% ------------------------------------------------------------------------
%  SURROGATE OUTPUT TRAJECTORY
%  ------------------------------------------------------------------------
figure('Name', sprintf('Surrogate Output Trajectory'),...
       'Color', 'w', ...
       'Position', [750 100 figSignalTrajSize]);
plot(1:length(yhat)-1, yhat(2:length(yhat)), ...
     'Color', colors(2,:), ...
     'LineWidth', 1.5);
xlabel('Time Index $k$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);
ylabel('$\hat{y}(k)$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);
ylim([yhatextremeval(1)-0.05*yhatrange, ...
      yhatextremeval(2)+0.05*yhatrange]);
yticks([yhatextremeval(1), yhatextremeval(2)]);
legend('$\hat{y}(k)$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize, ...
       'Location', 'best');
set(gca, ...
    'FontSize', fontSize, ...
    'TickLabelInterpreter', 'latex');
box on;
grid on;
%% ------------------------------------------------------------------------
%  REGRESSOR SPACE
%  ------------------------------------------------------------------------
% Generate the regressor-space plot
figure('Name', sprintf('3D Regressor Space'), ...
       'Color', 'w', ...
       'Position', [750 100 figRSpaceSize]);

hold on;
% Plot supporting points in the background
if showSupportingPoints

    Qphys = nrm.toPhysical(Q);
    scatter3(Qphys(:,1), Qphys(:,2), Qphys(:,3), suppPointSize, ...
             'MarkerFaceColor', colors(4,:), ...
             'MarkerEdgeColor', 'none');
end
% Plot the regressor trajectory as an unconnected scatter plot
scatter3(u, xTraj(:,1), xTraj(:,2), scatterSize, ...
         'MarkerFaceColor', colors(3,:), ...
         'MarkerEdgeColor', 'none');
grid on;
box on;
view(3);

% Set axis limits to the global signal bounds
ROIrange = ROIbox(:,2) - ROIbox(:,1);
xlim([ROIbox(1,1)-0.05*ROIrange(1), ...
      ROIbox(1,2)+0.05*ROIrange(1)]);
ylim([ROIbox(2,1)-0.05*ROIrange(2), ...
      ROIbox(2,2)+0.05*ROIrange(2)]);
zlim([ROIbox(3,1)-0.05*ROIrange(3), ...
      ROIbox(3,2)+0.05*ROIrange(3)]);
% Show only the lower and upper boundary values as axis ticks
xticks([ROIbox(1,1), ROIbox(1,2)]);
yticks([ROIbox(2,1), ROIbox(2,2)]);
zticks([ROIbox(3,1), ROIbox(3,2)]);
% Use LaTeX formatting for axis tick labels
set(gca, ...
    'TickLabelInterpreter', 'latex', ...
    'FontSize', fontSize);
% Axis labels
xlabel('$u(k)$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);
ylabel('$x_1(k)$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);
zlabel('$x_2(k)$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);

%% ------------------------------------------------------------------------
%  FREQUENCY SPECTRUM
%  ------------------------------------------------------------------------
% Amplitude spectrum of the excited input
fs = 1/Ts;
df = fs/cfg.N;
fk = 0:df:fs/2;
udetr = detrend(u);
U = fft(udetr);
UamplSp = 2*df*abs(U)/fs;

fTicks = linspace(fk(1), fk(end), 6);
Umax = max(UamplSp) + 0.025*max(UamplSp);

figure('Name', sprintf('Excited Frequency Spectrum'), ...
       'Color', 'w', ...
       'Position', [750 100 figSignalTrajSize]);
plot(fk, UamplSp(1:length(fk)), ...
     'Color', colors(1,:), ...
     'LineWidth', 1.5);
xlabel('$f \ [\mathrm{Hz}]$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);
ylabel('$|U(f)|$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);
xlim([fk(1), fk(end)]);
xticks(fTicks);
ylim([0, Umax]);
yticks([0, Umax]);
set(gca, ...
    'FontSize', fontSize, ...
    'TickLabelInterpreter', 'latex');
box on;
grid on;
