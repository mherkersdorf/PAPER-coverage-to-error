function X = lowDiscrepancySet(nPoints, nDim, seed)
%LOWDISCREPANCYSET  Randomly digit-shifted Sobol point set of the unit hypercube.
%
%   Toolbox: Statistics and Machine Learning Toolbox (sobolset).
%
%   Isolates the project's low-discrepancy sampling functionality. It exists
%   because rid.roi.estimateFillDistance takes a MAXIMUM, and a maximum is decided
%   entirely by coverage: a gap in the probe set is precisely where the deepest
%   hole of the measured set can hide. An i.i.d. uniform sample leaves
%   gaps of its own; a Sobol sequence is constructed so that every dyadic box of
%   the appropriate volume receives a point.
%
%   Randomness, and why it is a DIGIT shift. A Sobol sequence is
%   deterministic, so the seed must randomize it. The randomization used here is
%   the standard digital shift: each coordinate's binary expansion is XOR-ed with
%   a random bit pattern drawn from this function's own stream. That is a
%   permutation of the dyadic boxes, so it maps the unit cube onto itself and
%   leaves the sequence's equidistribution EXACTLY intact. The two alternatives
%   were rejected:
%
%     - MATLAB's `scramble` draws from the GLOBAL stream, which the project's
%       randomness rule forbids a library function from reading or writing.
%     - A Cranley-Patterson shift, x -> mod(x + u, 1), wraps around the boundary.
%       On the torus that is an isometry, but the ROI is a BOX: a point pushed
%       past the face stops covering what it covered, and the degradation lands
%       exactly at the boundary, which is where the deepest holes of a ROI-filling
%       design sit.
%
%   Sobol coordinates are dyadic rationals, so the 32-bit XOR is exact for every
%   point count this project can reach; anything below 2^-32 is rounded away, a
%   perturbation far under any geometric quantity here.
%
%   The global stream is never read and never written, so a call is reproducible
%   in isolation and inside parfor. Passing a RandStream rather than a seed lets a
%   caller take SEVERAL independent point sets from one stream -- which
%   rid.roi.estimateFillDistance does, for the probe set and for the independent
%   set that measures the probe set's own fill distance.
%
%   What the seed buys is what it buys for rid.ext.latinHypercube: a point set
%   reproducible WITHIN a MATLAB release. sobolset is a toolbox algorithm.
%
%   Inputs:
%       nPoints  (1x1)  number of points, positive integer
%       nDim     (1x1)  dimension, positive integer (sobolset allows up to 1111)
%       seed     (1x1)  nonnegative integer seed, or a RandStream object
%
%   Outputs:
%       X  (nPoints x nDim)  points in the unit hypercube, values in [0,1)

arguments
    nPoints (1,1) double {mustBeInteger, mustBePositive}
    nDim    (1,1) double {mustBeInteger, mustBePositive}
    seed
end

stream = buildStream(seed);

% The digital shift: one random 32-bit word per dimension, applied to every point.
scale     = 2^32;
shiftBits = uint32(randi(stream, [0, scale - 1], 1, nDim));

p = sobolset(nDim);
X = double(bitxor(uint32(p(1:nPoints, :) * scale), shiftBits)) / scale;

end

