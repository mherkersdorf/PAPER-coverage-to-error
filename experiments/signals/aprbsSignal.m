function u = aprbsSignal(N, range, T, construction, levelConstruction, seed)
%APRBSSIGNAL  APRBS from the vendor APRBS.m, every option stated, global stream restored.
%
%   A thin wrapper; the signal is APRBS.m's (C. Illg; not included,
%   see THIRD_PARTY.md). It exists for two reasons:
%
%   1. EVERY OPTION IS EXPLICIT. APRBS.m's header says the default Construction
%      is 'expD'; its code says 'expD_Tmin'. No default of APRBS.m is relied on.
%      Type is always 'APRBS' -- this is the APRBS wrapper.
%
%   2. THE CALLER'S GLOBAL STREAM IS UNTOUCHED. APRBS.m calls rng(Seed), which
%      replaces MATLAB's global stream. The global stream object and its state
%      are saved here and restored by onCleanup, on the error path too. Before
%      the call the global stream is set to a Mersenne twister, so rng(Seed)
%      inside APRBS.m seeds a twister whatever generator the caller used: the
%      signal depends on `seed` alone (within one MATLAB release).
%
%   LevelConstruction 'linspace' (APRBS.m's default, stated here): for c - 1
%   hold segments, c equally spaced levels over `range`, ends included, are
%   randomly permuted and assigned one per segment, so no level repeats and the
%   last level is unused. 'rand' draws each level i.i.d. uniform in `range`.
%
%   Inputs:
%       N                  (1x1) positive integer  number of samples
%       range              (1x2) double            [min, max] amplitude
%       T                  (1x1) positive double   APRBS.m's T; mean hold 2*T for
%                                                  'expD' and 'expD_Tmin'
%       construction       (1x1) string  "expD" | "expD_Tmin" | "LFSR"
%       levelConstruction  (1x1) string  "linspace" | "rand"
%       seed               (1x1) nonnegative integer
%   Outputs:
%       u  (N x 1) double  samples in rows
arguments
    N                 (1,1) double {mustBeInteger, mustBePositive}
    range             (1,2) double {mustBeReal, mustBeFinite}
    T                 (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
    construction      (1,1) string {mustBeMember(construction, ["expD", "expD_Tmin", "LFSR"])}
    levelConstruction (1,1) string {mustBeMember(levelConstruction, ["linspace", "rand"])}
    seed              (1,1) double {mustBeInteger, mustBeNonnegative}
end

if range(2) <= range(1)
    error("experiments:aprbsSignal:invalidRange", ...
        "range must be [min, max] with min < max; got %s.", mat2str(range));
end

callerStream = RandStream.getGlobalStream();
callerState  = callerStream.State;
restorer     = onCleanup(@() restoreGlobalStream(callerStream, callerState));

RandStream.setGlobalStream(RandStream("twister", "Seed", 0));
u = APRBS(N, range, T, ...
    "Type", 'APRBS', ...
    "Construction", char(construction), ...
    "LevelConstruction", char(levelConstruction), ...
    "Seed", seed);

u = reshape(u, [], 1);
if numel(u) ~= N
    error("experiments:aprbsSignal:wrongLength", ...
        "APRBS.m returned %d samples for N = %d.", numel(u), N);
end
end

% ---------------------------------------------------------------------------

function restoreGlobalStream(stream, state)
%RESTOREGLOBALSTREAM  Reinstall the caller's global stream object and its state.
RandStream.setGlobalStream(stream);
stream.State = state;
end
