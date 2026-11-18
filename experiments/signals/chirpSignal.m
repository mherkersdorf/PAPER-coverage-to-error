function u = chirpSignal(N, tS, frequencyRange)
%CHIRPSIGNAL  One complete chirp of N samples, from MATLAB's chirp.
%
%   A thin wrapper around chirp (Signal Processing Toolbox). One complete sweep
%   per dataset size, never a prefix of a longer sweep: the frequency runs from
%   frequencyRange(1) at the first sample to frequencyRange(2) at the last.
%
%   The chirps of the two studies sweep linearly from fs/N to fs/2 - fs/N (the full band up to
%   just below the Nyquist frequency) in the Duffing study and from 5 to 150 Hz in the
%   Bouc-Wen study; the frequency range is a required argument with no default. The
%   sweep method and the phase are chirp's own defaults, stated: 'linear', 0 degrees --
%   a cosine, so u(1) = 1.
%
%   Inputs:
%       N               (1x1) integer >= 2  number of samples
%       tS              (1x1) positive      sample time [s]
%       frequencyRange  (1x2) nonnegative   [fStart, fEnd] [Hz], each below 1/(2*tS)
%   Outputs:
%       u  (N x 1) double  samples in rows, amplitude 1
arguments
    N              (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(N, 2)}
    tS             (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
    frequencyRange (1,2) double {mustBeReal, mustBeFinite, mustBeNonnegative}
end

if any(frequencyRange >= 1 / (2*tS))
    error("experiments:chirpSignal:aboveNyquist", ...
        "frequencyRange %s must lie below the Nyquist frequency %g Hz.", ...
        mat2str(frequencyRange), 1 / (2*tS));
end

t  = (0:N-1).' * tS;
u  = chirp(t, frequencyRange(1), t(end), frequencyRange(2), "linear", 0);
u  = reshape(u, [], 1);
end
