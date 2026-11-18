function X = latinHypercube(nPoints, nDim, seed, opts)
%LATINHYPERCUBE  Latin hypercube sample of the unit hypercube, from an explicit seed.
%
%   Toolbox: Statistics and Machine Learning Toolbox (lhsdesign).
%
%   Isolates the project's Latin-hypercube sampling functionality. The sample is
%   stratified per dimension: each column contains exactly one point in each of
%   the nPoints intervals ((kk-1)/nPoints, kk/nPoints).
%
%   Randomness. lhsdesign takes neither a seed nor a RandStream and draws
%   from the global stream. This function therefore builds its own stream from
%   `seed`, installs it as the global stream for the duration of the lhsdesign
%   call, and restores the caller's stream with onCleanup -- on the error path
%   too. The caller's stream object and its state are bit-identical before and
%   after, and each parfor worker has its own global stream, so the borrow is not
%   observable from outside this call.
%
%   What the seed buys. It reproduces a design *within a MATLAB release only*:
%   lhsdesign is a toolbox algorithm and may change between releases. Seeded
%   regeneration is therefore not a reproducibility mechanism for a run -- run
%   records store Q itself and keep the generator and seed as provenance.
%
%   The three options are a direct pass-through to lhsdesign and carry *its*
%   defaults, not project defaults: the project's values live in
%   rid.defaultOptions and are passed explicitly by rid.roi.supportingPoints.
%   The admissible values of `criterion` are validated by lhsdesign itself, so
%   this wrapper does not duplicate -- and cannot silently diverge from -- the
%   toolbox's own list.
%
%   Inputs:
%       nPoints          (1x1)  number of samples, positive integer
%       nDim             (1x1)  dimension, positive integer
%       seed             (1x1)  nonnegative integer seed, or a RandStream object
%       opts.criterion   (1x1)  string, pass-through: "maximin" (default),
%                               "correlation" or "none"
%       opts.iterations  (1x1)  pass-through, nonnegative integer, default 5
%       opts.smooth      (1x1)  logical, pass-through, default true; false places
%                               points at the stratum midpoints
%
%   Outputs:
%       X  (nPoints x nDim)  sample in the unit hypercube, values in (0,1)

arguments
    nPoints (1,1) double {mustBeInteger, mustBePositive}
    nDim    (1,1) double {mustBeInteger, mustBePositive}
    seed
    opts.criterion  (1,1) string  = "maximin"
    opts.iterations (1,1) double {mustBeInteger, mustBeNonnegative} = 5
    opts.smooth     (1,1) logical = true
end

stream = buildStream(seed);

if opts.smooth
    smoothArg = 'on';
else
    smoothArg = 'off';
end

% Borrow the global stream for the lhsdesign call and restore it afterwards,
% including on the error path -- the Latin-hypercube seed handling.
previousStream = RandStream.getGlobalStream();
restoreStream  = onCleanup(@() RandStream.setGlobalStream(previousStream));
RandStream.setGlobalStream(stream);

X = lhsdesign(nPoints, nDim, ...
    "criterion",  opts.criterion, ...
    "iterations", opts.iterations, ...
    "smooth",     smoothArg);

end
