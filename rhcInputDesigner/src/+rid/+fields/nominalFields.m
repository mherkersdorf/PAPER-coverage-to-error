function fields = nominalFields(nSup)
%NOMINALFIELDS  The three refinement fields at their nominal value.
%
%   Every field is calibrated so that 1 is nominal: a field at or below
%   nominal leaves the baseline untouched, and only values above nominal refine.
%   This function returns that neutral configuration, which is what
%   `demandedResolution` produces when no field is supplied -- the pure
%   space-filling design and the robust fallback whenever no target information is
%   available.
%
%   It exists so that no call site has to write `ones(nSup, 1)` three times; the
%   literal 1 for "nominal" is defined here and, and
%   nowhere else.
%
%   IT LIVES IN +fields, NOT IN +crit, because it CONSTRUCTS fields and does not
%   evaluate the criterion. Keeping it in
%   +crit would make +fields depend on the criterion for a field constructor,
%   which is the wrong direction. The plural name is deliberate: `nominal.m` and
%   `Nominal.m` cannot coexist in one package folder on a case-insensitive
%   filesystem, and that failure presents as anything but a naming problem.
%
%   Inputs:
%       nSup    (1x1)  positive integer, number of supporting points
%
%   Outputs:
%       fields  (1x1) struct with fields
%           .relevance   (nSup x 1)  ones
%           .complexity  (nSup x 1)  ones
%           .noiseLevel  (nSup x 1)  ones

arguments
    nSup (1,1) double {mustBeReal}
end

if ~isfinite(nSup) || mod(nSup, 1) ~= 0 || nSup < 1
    error("rid:fields:invalidFieldSize", ...
        "nSup must be a positive integer (got %g).", nSup);
end

fields = struct( ...
    "relevance",  ones(nSup, 1), ...
    "complexity", ones(nSup, 1), ...
    "noiseLevel", ones(nSup, 1));

end
