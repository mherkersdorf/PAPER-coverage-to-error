function epsReq = spacingScaledResolution(s, nSup, nDim, scale)
%SPACINGSCALEDRESOLUTION  A demanded resolution on the scale of the SAMPLE SPACING.
%
%   Shared by the two randomized oracle suites, CoverageStateOracleTest and
%   CoverageStateGradientTest, which must land in the same regime to be comparable.
%
%   WHY THE DRAW IS SCALED. Which regime a randomized sequence lands in is decided
%   by epsReq relative to the nearest-neighbour distance the samples leave behind,
%   and that distance depends violently on nDim: 25 points sit 0.04 apart in
%   [0,1], 0.2 apart in [0,1]^2 and 0.34 apart in [0,1]^3. A fixed absolute draw
%   therefore resolved almost everything in one dimension and almost nothing in
%   three. Scaling by the nominal spacing M^(-1/nDim), at the M the rebuild branch
%   caps the history at, makes the regime the same in every dimension.
%
%   WHY IT MATTERS THAT IT DOES. A resolved supporting point contributes zero to
%   the cached path AND zero to rid.crit.evaluate, so once a sequence has resolved
%   itself the oracle compares zero with zero and proves nothing
%   (a randomized sample must be asserted non-degenerate).
%   Under the fixed draw this replaced, two of the five oracle seeds fell to a
%   median active fraction of 0.06 and 0.08 and spent 148 of 1000 operations
%   comparing a fully resolved state. Under this one the median is 0.47 to 0.75 on
%   every seed and no operation is fully resolved.
%
%   THIS IS A PROPERTY OF THE TEST'S SAMPLING, NOT OF THE CRITERION. `epsReq` is an
%   experiment's declared commitment and nothing in the design relates it
%   to a sample count; the relation exists here only to keep a randomized sequence
%   in the regime it is meant to exercise.
%
%   Inputs:
%       s       (1x1) RandStream  the sequence's own stream, never the global one
%       nSup    (1x1)  supporting points
%       nDim    (1x1)  dimensions of the regressor space
%       scale   (1x1)  1 is the mixed regime; a small value is the FINE regime,
%                      where almost nothing resolves
%   Outputs:
%       epsReq  (nSup x 1)  finite and strictly positive

arguments
    s     (1,1) RandStream
    nSup  (1,1) double {mustBeInteger, mustBePositive}
    nDim  (1,1) double {mustBeInteger, mustBePositive}
    scale (1,1) double {mustBePositive}
end

spacing = 25^(-1 / nDim);
epsReq  = scale * (0.1 + 0.9 * rand(s, nSup, 1)) * spacing;

end
