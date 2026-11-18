function [u, info] = multisineSignal(N, fs, opts)
%MULTISINESIGNAL  One complete period of a Schroeder multisine, from multisine.m; full band by default.
%
%   A thin wrapper; the signal is multisine.m's (B. Holmes; not included,
%   see THIRD_PARTY.md). One complete period of N samples per dataset size, never a
%   prefix of a longer period: a Schroeder multisine sweeps its band within the
%   period, and a prefix holds only part of the sweep.
%
%   Settings, all stated: Schroeder phases (InitialPhase 0), flat magnitude
%   (multisine.m's flat default -- its validator rejects an explicit `false`),
%   frequency-domain synthesis (TimeDomain false), StartAtZero true (the period
%   is circularly shifted to its lowest-gradient zero crossing), Normalise true
%   (peak |u| = 1).
%
%   THE BAND IS [fs/N, fs/2 - fs/N], CHOSEN, not [fs/N, fs/2]:
%   multisine.m maps a limit f to the FFT bin 1 + round(f*N/fs) and rejects any
%   bin > N/2, so fs/2 always errors; and a component AT Nyquist is sampled only
%   at its extremes or zero crossings, so its amplitude would depend on its
%   phase. Implemented as lines 1 to floor(N/2) - 1, the highest line
%   multisine.m admits: exactly fs/2 - fs/N for even N; for odd N, where
%   fs/2 - fs/N falls between lines and multisine.m would round it up to a
%   rejected bin, it is fs/2 - 1.5*fs/N.
%
%   AN EXPLICIT BAND (opts.band, P3: 5-150 Hz) excites every line whose frequency lies in [fLow, fHigh]:
%   lines ceil(fLow*N/fs) to floor(fHigh*N/fs), each limit with a tolerance of
%   1e-9 lines so that a limit lying exactly on a line is included despite
%   rounding. The lines must lie within 1 to floor(N/2) - 1, the range above;
%   a band holding no line is an error. The default, band = [], is the full
%   band above and runs the original code path unchanged, so its output is
%   bit-identical to the wrapper before the option existed.
%
%   Inputs:
%       N   (1x1) integer >= 6  period length in samples (= the dataset size)
%       fs  (1x1) positive integer  sampling frequency [Hz]; multisine.m requires
%                                   an integer. It scales frequencies only; the
%                                   samples do not depend on it.
%       opts.band  (1x2) or [] double  [fLow, fHigh] [Hz]; [] = the full band
%   Outputs:
%       u     (N x 1) double  one period, samples in rows, max |u| = 1
%       info  (1x1) struct
%           .band         (1x2) double  [fLow, fHigh] passed to multisine.m [Hz]
%           .lines        (1x2) double  first and last line index (multiples of fs/N)
%           .nLines       (1x1) double  number of excited lines
%           .crestFactor  (1x1) double  max|u| / rms(u)
arguments
    N  (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(N, 6)}
    fs (1,1) double {mustBeInteger, mustBePositive}
    opts.band (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative, mustBeBandOrEmpty} = []
end

if isempty(opts.band)
    lineLow  = 1;
    lineHigh = floor(N/2) - 1;
else
    [lineLow, lineHigh] = bandLines(opts.band, N, fs);
end
band = [lineLow, lineHigh] * fs / N;

if exist("multisine", "file") ~= 2          % not shipped (THIRD_PARTY.md)
    error("experiments:multisineSignal:missingMultisine", ...
        "multisine.m (B. Holmes) is not on the path. It is not included in this repository; see THIRD_PARTY.md.");
end
u = multisine(band, fs, N, ...
    "PhaseResponse", 'Schroeder', ...
    "InitialPhase", 0, ...
    "TimeDomain", false, ...
    "StartAtZero", true, ...
    "Normalise", true);
u = reshape(u, [], 1);

if numel(u) ~= N
    error("experiments:multisineSignal:wrongLength", ...
        "multisine.m returned %d samples for N = %d.", numel(u), N);
end

info = struct("band", band, "lines", [lineLow, lineHigh], ...
    "nLines", lineHigh - lineLow + 1, ...
    "crestFactor", max(abs(u)) / sqrt(mean(u.^2)));
end

% ---------------------------------------------------------------------------

function [lineLow, lineHigh] = bandLines(bandHz, N, fs)
%BANDLINES  First and last line index with frequency inside bandHz, checked against multisine.m's range.
tolLines = 1e-9;
lineLow  = max(1, ceil(bandHz(1) * N / fs - tolLines));
lineHigh = floor(bandHz(2) * N / fs + tolLines);
if lineHigh > floor(N/2) - 1
    error("experiments:multisineSignal:bandAboveRange", ...
        "band %s reaches line %d; the highest admissible line is floor(N/2) - 1 = %d (%g Hz).", ...
        mat2str(bandHz), lineHigh, floor(N/2) - 1, (floor(N/2) - 1) * fs / N);
end
if lineHigh < lineLow
    error("experiments:multisineSignal:emptyBand", ...
        "band %s holds no line of the resolution fs/N = %g Hz.", mat2str(bandHz), fs / N);
end
end

function mustBeBandOrEmpty(band)
%MUSTBEBANDOREMPTY  [] or [fLow, fHigh] with fLow < fHigh.
if ~(isempty(band) || (numel(band) == 2 && band(1) < band(2)))
    error("experiments:multisineSignal:invalidBand", ...
        "band must be [] or [fLow, fHigh] with fLow < fHigh; got %s.", mat2str(band));
end
end
