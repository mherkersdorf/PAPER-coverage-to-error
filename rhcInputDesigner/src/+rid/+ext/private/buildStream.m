function stream = buildStream(seed)
%BUILDSTREAM  A private RandStream from an integer seed, or the given stream.
%
%   Shared by the +ext samplers, which take "a nonnegative integer seed, or a
%   RandStream object" under one contract: a private stream keeps a call
%   reproducible in isolation and inside parfor, and accepting a stream lets one
%   caller draw several independent point sets without the seeds having to be
%   derived from one another.

if isa(seed, "RandStream")
    stream = seed;
    return
end

isIntegerSeed = isnumeric(seed) && isscalar(seed) && isreal(seed) ...
    && isfinite(seed) && seed >= 0 && mod(seed, 1) == 0;

if ~isIntegerSeed
    error("rid:ext:invalidSeed", ...
        "seed must be a nonnegative integer scalar or a RandStream object.");
end

stream = RandStream("twister", "Seed", seed);

end
