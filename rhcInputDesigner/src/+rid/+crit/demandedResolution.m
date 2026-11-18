function [epsReq, info] = demandedResolution(fields, opts)
%DEMANDEDRESOLUTION  Locally demanded resolution from the three refinement fields.
%
%   The baseline `epsReq0` is the COARSEST resolution demanded anywhere in the ROI.
%   The application relevance a_j, the process complexity c_j and the local noise
%   level v_j may only refine it:
%
%       refinement = min((relevance .* complexity).^(-beta) .* noiseLevel.^(-betaNoise), 1);
%       epsReq     = epsReq0 .* refinement;
%
%   The clip at 1 is load-bearing, not cosmetic. It makes adaptation
%   strictly additive -- a field can add resolution where it is demanded and can
%   never coarsen any region below the space-filling baseline -- and it is what
%   makes `max(epsReq) <= epsReq0` hold, hence the certification bound
%   `eps <= epsReq0 + fillDistanceQ` independent of the fields (section 5).
%   Equality holds whenever at least one supporting point is unrefined; if every
%   point is refined the maximum lies strictly below the baseline and the bound is
%   conservative.
%
%   NO NORMALIZATION IS APPLIED. The fields are not divided by their mean,
%   median or maximum, so `epsReq(j)` depends only on the fields at q_j. Any such
%   reference statistic would destroy locality: enlarging the ROI with a few easy
%   regions would silently refine the hard ones, and two runs on different ROIs
%   would not be comparable. The price is that a miscalibrated estimator is not
%   neutralized automatically, so `info.fracRefined` and `info.minEpsReqRatio`
%   report it -- a field orders of magnitude off shows up immediately in those two
%   numbers.
%
%   THE FLOOR. In certification mode `epsReq` is projected onto
%   [fillDistanceQ, epsReq0]: nothing finer than the fill distance of Q can
%   be certified. The projection is the one clamp in the design and is never
%   silent -- it is counted in `info.nFloored` and raises
%   `rid:crit:resolutionFloorActive`. A large floored fraction is the signal that
%   `nSup` must grow or `epsReq0` must be coarsened. In budget mode no floor is
%   applied: nothing is certified, `epsReq` only shapes the distribution, and
%   clamping would distort that shape for no benefit.
%
%   Inputs:
%       fields  (1x1) struct with any subset of the three members below; a member
%               that is absent defaults to nominal, `ones(nSup, 1)`, and
%               `nSup` is inferred from the members present. An empty struct is an
%               error, since `nSup` is then not inferable, and an unrecognized
%               member is an error (`rid:crit:unknownField`) -- with optional
%               members a misspelt name would otherwise leave a region unrefined
%               with no symptom.
%           .relevance        (nSup x 1)  a_j, strictly positive, 1 = nominal
%           .complexity       (nSup x 1)  c_j, strictly positive, 1 = nominal
%           .noiseLevel       (nSup x 1)  v_j, strictly positive, 1 = nominal
%       opts.epsReq0          (1x1)  REQUIRED, no default; > 0, finite
%       opts.beta             (1x1)  >= 0; 0 disables target shaping entirely
%       opts.betaNoise        (1x1)  >= 0; the noise field's own exponent
%       opts.mode             (1x1) string  "certification" | "budget"
%       opts.fillDistanceQ    (1x1)  required and > 0 in certification mode;
%                                    NaN (the default) means "not supplied"
%
%   Outputs:
%       epsReq  (nSup x 1)  locally demanded resolution
%       info    (1x1) struct with fields
%           .refinement       (nSup x 1)  the factor phi_j, clipped at 1
%           .fracRefined      (1x1)  fraction with refinement < 1
%           .nFloored         (1x1)  count projected onto the floor (0 in budget mode)
%           .fracFloored      (1x1)
%           .minEpsReq        (1x1)
%           .maxEpsReq        (1x1)
%           .minEpsReqRatio   (1x1)  min(epsReq) / epsReq0
%           .headroom         (1x1)  epsReq0 / fillDistanceQ; NaN if not supplied
%           .mode .epsReq0 .fillDistanceQ .beta .betaNoise .nSup   (provenance)

arguments
    fields (1,1) struct
    opts.epsReq0       (1,1) double {mustBeReal}
    opts.beta          (1,1) double {mustBeReal, mustBeNonnegative} = rid.defaultOptions().crit.beta
    opts.betaNoise     (1,1) double {mustBeReal, mustBeNonnegative} = rid.defaultOptions().crit.betaNoise
    opts.mode          (1,1) string = rid.defaultOptions().crit.mode
    opts.fillDistanceQ (1,1) double {mustBeReal} = NaN
end

fieldNames = ["relevance", "complexity", "noiseLevel"];

% --- options ---------------------------------------------------------------
if ~isfield(opts, "epsReq0")
    error("rid:crit:invalidResolution", ...
        "epsReq0 is required and has no default: it is the experiment's " + ...
        "baseline commitment, the coarsest resolution demanded anywhere.");
end
if ~isfinite(opts.epsReq0) || opts.epsReq0 <= 0
    error("rid:crit:invalidResolution", ...
        "epsReq0 must be finite and strictly positive (got %g).", opts.epsReq0);
end
if ~ismember(opts.mode, ["certification", "budget"])
    error("rid:crit:unknownMode", ...
        "mode must be ""certification"" or ""budget"" (got ""%s"").", opts.mode);
end
if ~isnan(opts.fillDistanceQ) && (~isfinite(opts.fillDistanceQ) || opts.fillDistanceQ < 0)
    error("rid:crit:invalidResolution", ...
        "fillDistanceQ must be NaN (not supplied), or finite and nonnegative (got %g).", ...
        opts.fillDistanceQ);
end

isCertification = opts.mode == "certification";
if isCertification && ~(isfinite(opts.fillDistanceQ) && opts.fillDistanceQ > 0)
    error("rid:crit:invalidResolution", ...
        "certification mode requires a measured fillDistanceQ > 0: the floor of " + ...
        "the certifiability floor is defined by it, and nothing finer than it can be certified.");
end

% --- fields ----------------------------------------------------------------
present = string(fieldnames(fields)).';
unknown = setdiff(present, fieldNames);
if ~isempty(unknown)
    error("rid:crit:unknownField", ...
        "unrecognized member(s) of fields: %s. Recognized members are %s. " + ...
        "A misspelt member would otherwise default to nominal and silently " + ...
        "leave a region unrefined.", ...
        strjoin("""" + unknown + """", ", "), strjoin("""" + fieldNames + """", ", "));
end
if isempty(present)
    error("rid:crit:invalidResolution", ...
        "fields must carry at least one of %s: nSup is inferred from the " + ...
        "members present. Use rid.fields.nominalFields(nSup) for the nominal case.", ...
        strjoin("""" + fieldNames + """", ", "));
end

nSup = size(fields.(present(1)), 1);
for name = present
    value = fields.(name);
    if ~isa(value, "double") || ~isreal(value) || ~iscolumn(value) || size(value, 1) ~= nSup
        error("rid:crit:invalidResolution", ...
            "fields.%s must be a real (nSup x 1) double column matching the " + ...
            "other fields (nSup = %d, got %s).", ...
            name, nSup, mat2str(size(value)));
    end
    if ~all(isfinite(value)) || any(value <= 0)
        error("rid:crit:invalidField", ...
            "fields.%s must be finite and strictly positive. A region that is " + ...
            "not to be modelled is excluded through Q, not through a field.", ...
            name);
    end
end

% Absent members are nominal; see the calibration of the fields. Written out rather than folded into the
% product so that the nominal default is visible at the point it is applied.
relevance  = nominalIfAbsent(fields, "relevance",  nSup);
complexity = nominalIfAbsent(fields, "complexity", nSup);
noiseLevel = nominalIfAbsent(fields, "noiseLevel", nSup);

% --- the demanded resolution -----------------------------------------------
refinement = min((relevance .* complexity).^(-opts.beta) .* noiseLevel.^(-opts.betaNoise), 1);
epsReq     = opts.epsReq0 .* refinement;

nFloored = 0;
if isCertification
    isFloored = epsReq < opts.fillDistanceQ;
    nFloored  = sum(isFloored);
    epsReq(isFloored) = opts.fillDistanceQ;
end

% --- diagnostics -----------------------------------------------------------
% The clamp is reported whenever it acts. Warning rather than info alone, because
% an info field a caller may ignore is a silent clamp in practice. It fires on
% every resolution update while the floor is active; warning only on a CHANGE of
% the floored set needs the previous state, which only +rhc holds.
if nFloored > 0
    warning("rid:crit:resolutionFloorActive", ...
        "%d of %d supporting points (%.1f%%) were projected onto the " + ...
        "certifiability floor fillDistanceQ = %g. A large fraction means the " + ...
        "intended target is not representable at this discretization: grow " + ...
        "nSup or coarsen epsReq0.", ...
        nFloored, nSup, 100 * nFloored / nSup, opts.fillDistanceQ);
end

info = struct( ...
    "refinement",     refinement, ...
    "fracRefined",    mean(refinement < 1), ...
    "nFloored",       nFloored, ...
    "fracFloored",    nFloored / nSup, ...
    "minEpsReq",      min(epsReq), ...
    "maxEpsReq",      max(epsReq), ...
    "minEpsReqRatio", min(epsReq) / opts.epsReq0, ...
    "headroom",       opts.epsReq0 / opts.fillDistanceQ, ...
    "mode",           opts.mode, ...
    "epsReq0",        opts.epsReq0, ...
    "fillDistanceQ",  opts.fillDistanceQ, ...
    "beta",           opts.beta, ...
    "betaNoise",      opts.betaNoise, ...
    "nSup",           nSup);

end

% --- internal ---------------------------------------------------------------

function value = nominalIfAbsent(fields, name, nSup)
%NOMINALIFABSENT  The named field, or the nominal column if it was not supplied.
if isfield(fields, name)
    value = fields.(name);
else
    value = ones(nSup, 1);
end
end
