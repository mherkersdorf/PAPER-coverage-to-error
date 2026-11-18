function [y, info] = p3BoucWenSimulate(u, opts)
%P3BOUCWENSIMULATE  The Bouc-Wen benchmark process: force u at 750 Hz -> noisy displacement y at 750 Hz.
%
%   Follows the benchmark's BoucWen_ExampleIntegration.m step by
%   step (cited "l. n" below) and the user guidelines of the benchmark
%   (BoucWenBenchmark.pdf, Noel & Schoukens 2016, section 4); the settings are
%   those of the study:
%
%   1. TO 15 kHz. A signal that exists only at 750 Hz (opts.inputRate = 750,
%      every training, validation and BLA signal) is interpolated
%      band-limited: interp(., 20) with its default filter (Signal Processing
%      Toolbox, replacing the zero-order hold so that these
%      signals take the test data's path). The default is l = 4, alpha = 0.5:
%      a firls lowpass of 2*20*4 + 1 = 161 taps at 15 kHz, gain 20 on
%      0-187.5 Hz and stopbands around the images, delay-compensated. The
%      750 Hz samples pass nearly unchanged -- to about 4e-10 of max|u|
%      (measured; firls does not make it an exact Nyquist filter),
%      i.e. about 1e-7 N at 200 N -- and each 15 kHz sample depends on the
%      750 Hz samples within +-4. At both ends interp
%      extrapolates by point reflection of 2l = 8 samples. The 750 Hz samples
%      stay in the caller's range; the 15 kHz input may overshoot between
%      them (Gibbs ringing at fast changes), which is allowed and RECORDED in
%      info.uFineRangeKept, never clipped. A signal built at 15 kHz
%      (opts.inputRate = 15000, the test-class multisine and sweep) is
%      integrated as built, as in l. 31-39.
%   2. EXTRA SEGMENTS against the edge effects of the filters (l. 28-29 and
%      61-65). decimate(.., 'fir') filters with an order-30 FIR and
%      compensates its delay, so each stage reaches 15 of its own samples
%      ahead and behind: 15 + 30 + 60 = 105 samples at 15 kHz over the stages
%      2-2-5, i.e. 5.25 samples at 750 Hz; within that reach of either end,
%      decimate mirrors the data instead of using it. The example appends one
%      extra period, which is the true continuation of a periodic signal. For
%      any signal, periodic or not, this function adds
%        - IN FRONT, 20 samples at 750 Hz (400 at 15 kHz) of zero input. The
%          process starts at rest (l. 43-45) and stays at rest under zero
%          force;
%        - AT THE END, 20 samples holding u(end).
%      Both are cut after decimation. On the 750 Hz path they are appended at
%      750 Hz, BEFORE the interpolation, so interp's point reflection acts on
%      constant segments (and reproduces them) 20 samples away from the kept
%      part, beyond its reach of 4. The look-aheads add: a kept output sample
%      depends on input up to 4 (interp) + 5.25 (decimate) = 9.25 samples
%      ahead, within the 20 appended. Behind, the last 4 samples of the front
%      segment carry the interpolation's pre-ringing of u(1..4), part of the
%      band-limited signal: the process is exactly at rest until 16 samples
%      into the front segment, and the kept part starts from the state those
%      4 samples leave, not from exact rest. On the 15 kHz path the segments
%      are appended at 15 kHz as before, and the kept output starts from
%      rest. Because the process is causal, the kept output depends on the
%      appended input only through these look-aheads, which exist in the
%      example as well.
%   3. INTEGRATION, l. 43-48: BoucWen_NewmarkIntegration(h, u, 0, 0, 0) with
%      h = 1/15000 and zero initial conditions.
%   4. DECIMATION, l. 53-57: decimate(.., r, 'fir') for r in factor(20) =
%      2, 2, 5.
%
%   THE NOISE IS THE INTEGRATOR'S OWN, NEVER ADDED HERE. The benchmark
%   (section 4) states that Gaussian noise band-limited to 0-375 Hz with RMS
%   8e-3 mm is added automatically to the synthesized y. The integrator does
%   so, drawing from MATLAB's GLOBAL random stream (an unrecorded check). Adding noise here would add it twice. To
%   make the noise depend on opts.noiseSeed alone, the global stream is
%   replaced by a fresh Mersenne twister seeded with opts.noiseSeed for the
%   integrator call, and the caller's stream object and state are restored by
%   onCleanup, on the error path too (the pattern of aprbsSignal.m). After the
%   decimation chain the noise RMS is about 7.36e-6 m, not 8e-6 m: the noise
%   as the benchmark's own chain delivers it, not rescaled.
%   There is no noiseless mode: every simulation carries the benchmark noise.
%   The noise realization depends on the record length, so
%   prefixes of a data set share their noise only when taken from one
%   simulation at full length.
%
%   Requires BoucWen_NewmarkIntegration.p, the benchmark's integrator (not
%   included; THIRD_PARTY.md), on the path.
%
%   Inputs:
%       u               (N x 1) double  force [N] at 750 Hz (interpolated), or
%                                       (20N x 1) at 15 kHz when opts.inputRate = 15000
%       opts.noiseSeed  (1x1) nonnegative integer  seed of the integrator's noise; required
%       opts.inputRate  (1x1) 750 | 15000  rate of u [Hz]; default 750
%   Outputs:
%       y     (N x 1) double  displacement [m] at 750 Hz, with the benchmark noise
%       info  (1x1) struct
%           .fs                 (1x1) 750, working sampling frequency [Hz]
%           .fsInt              (1x1) 15000, integration sampling frequency [Hz]
%           .upsample           (1x1) 20
%           .decimationFactors  (1x3) [2 2 5]
%           .decimationFilter   (1x1) string "fir"
%           .initialConditions  (1x3) [y0, dy0, z0] = [0 0 0]
%           .inputRate          (1x1) opts.inputRate
%           .inputPath          (1x1) string "interp" (750 Hz) | "asBuilt" (15 kHz)
%           .interpolation      (1x1) struct .function "interp", .factor 20, .l 4,
%                                     .alpha 0.5 on the 750 Hz path; [] otherwise
%           .uFineRangeKept     (1x2) [min, max] of the 15 kHz input over the
%                                     kept part: the intersample overshoot
%           .nExtraBefore       (1x1) 20, samples at 750 Hz of zero input in front
%           .nExtraAfter        (1x1) 20, samples at 750 Hz holding the last input
%           .noiseSeed          (1x1) opts.noiseSeed
%           .noiseSource        (1x1) string, where the noise comes from
%           .integratorFile     (1x1) string, which() of the integrator
%           .uFine              (nFine x 1) the input integrated, extra segments included
%           .runtimeS           (1x1) wall time of the integration [s]
arguments
    u              (:,1) double {mustBeReal, mustBeFinite, mustBeNonempty}
    opts.noiseSeed (1,1) double {mustBeInteger, mustBeNonnegative}
    opts.inputRate (1,1) double {mustBeMember(opts.inputRate, [750, 15000])} = 750
end

% The process as the benchmark defines it (l. 12-15, 43-45, 53-57).
fs             = 750;
upsampleFactor = 20;
fsInt          = fs * upsampleFactor;
h              = 1 / fsInt;
nExtra         = 20;                 % per side, at 750 Hz; see step 2
drate          = factor(upsampleFactor); % 2, 2, 5
initialState   = [0, 0, 0];          % y0, dy0, z0

if exist("BoucWen_NewmarkIntegration", "file") ~= 6
    error("experiments:p3BoucWenSimulate:missingIntegrator", ...
        "BoucWen_NewmarkIntegration.p is not on the path. Install it in " + ...
        "a folder on the MATLAB path; it is part of the Bouc-Wen benchmark (THIRD_PARTY.md) " + ...
        "and is not included in this repository.");
end

% Steps 1 and 2: the input at the integration rate, with its extra segments.
nExtraFine = nExtra * upsampleFactor;
if opts.inputRate == fs
    interpolation = struct("function", "interp", "factor", upsampleFactor, "l", 4, "alpha", 0.5);
    uExtended = [zeros(nExtra, 1); u; repmat(u(end), nExtra, 1)];
    uFine = interp(uExtended, upsampleFactor, interpolation.l, interpolation.alpha);
    inputPath = "interp";
    N = numel(u);
else
    if mod(numel(u), upsampleFactor) ~= 0
        error("experiments:p3BoucWenSimulate:notWholeSamples", ...
            "A 15 kHz input must hold a whole number of 750 Hz samples (a multiple of %d); got %d.", ...
            upsampleFactor, numel(u));
    end
    interpolation = [];
    uFine = [zeros(nExtraFine, 1); u; repmat(u(end), nExtraFine, 1)];
    inputPath = "asBuilt";
    N = numel(u) / upsampleFactor;
end
uKept = uFine(nExtraFine + (1 : N*upsampleFactor));

% Step 3: integration, with the noise seeded through the global stream.
callerStream = RandStream.getGlobalStream();
callerState  = callerStream.State;
restorer     = onCleanup(@() restoreGlobalStream(callerStream, callerState));
RandStream.setGlobalStream(RandStream("twister", "Seed", opts.noiseSeed));

tInt  = tic;
yFine = BoucWen_NewmarkIntegration(h, uFine, initialState(1), initialState(2), initialState(3));
runtimeS = toc(tInt);
clear restorer

% Step 4: low-pass filtering and decimation in the prime factors.
yDec = reshape(yFine, [], 1);
for kStage = 1:numel(drate)
    yDec = decimate(yDec, drate(kStage), 'fir');
end
if numel(yDec) ~= N + 2*nExtra
    error("experiments:p3BoucWenSimulate:wrongLength", ...
        "Decimation returned %d samples; expected %d.", numel(yDec), N + 2*nExtra);
end
y = yDec(nExtra + (1:N));

info = struct( ...
    "fs", fs, "fsInt", fsInt, "upsample", upsampleFactor, ...
    "decimationFactors", drate, "decimationFilter", "fir", ...
    "initialConditions", initialState, "inputRate", opts.inputRate, ...
    "inputPath", inputPath, "interpolation", interpolation, ...
    "uFineRangeKept", [min(uKept), max(uKept)], ...
    "nExtraBefore", nExtra, "nExtraAfter", nExtra, ...
    "noiseSeed", opts.noiseSeed, ...
    "noiseSource", "BoucWen_NewmarkIntegration (global stream = twister(noiseSeed) during the call)", ...
    "integratorFile", string(which("BoucWen_NewmarkIntegration")), ...
    "uFine", uFine, "runtimeS", runtimeS);
end

% ---------------------------------------------------------------------------

function restoreGlobalStream(stream, state)
%RESTOREGLOBALSTREAM  Reinstall the caller's global stream object and its state.
RandStream.setGlobalStream(stream);
stream.State = state;
end
