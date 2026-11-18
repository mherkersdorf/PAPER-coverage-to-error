classdef Normalization
%NORMALIZATION  Affine per-dimension map between physical and normalized coordinates.
%
%   The transform is derived from the declared ROI box and is fixed at run start.
%   It is NEVER derived from data statistics: a data-driven scaling would
%   close a feedback loop -- the geometry would become a function of the data the
%   design itself produces -- and would silently invalidate every cached distance,
%   `fillDistanceQ` and `epsReq`.
%
%   With the default `dimWeights = 1` the physical box maps onto the unit cube
%   [0,1]^nDim; non-unit weights stretch it into [0,w_1] x ... x [0,w_nDim].
%
%   The transform NEVER clamps. The ROI defines the normalization; it does not
%   constrain the signal, so the admissible
%   region may be wider than the ROI, normalized values below 0 and above 1 are
%   legitimate, and they round-trip -- up to floating-point rounding, since
%   (x - offset)/span*span + offset is not bit-identical for every x. Tests assert
%   a relative tolerance of 1e-12.
%
%   Scaled by this transform: signal amplitudes, regressor coordinates and every
%   length in that space -- `epsReq`, `epsReq0`, `fillDistanceQ`, all distances.
%   Not scaled: time. Sample time, step and segment lengths, horizons and the
%   surrogate's time constant stay physical.
%
%   Construction:
%       nrm = rid.roi.Normalization(box, dimName)
%       nrm = rid.roi.Normalization(box, dimName, "dimWeights", w)
%
%   `box` and `dimName` come from `rep.bounds(roi)`, which is an ordered lookup over
%   a declaration carrying one range per dimension -- nothing is expanded or derived.
%   This class consumes the box; it is not an ROI object.
%
%   Inputs:
%       box               (nDim x 2)  physical bounds, [lower upper], upper > lower
%       dimName           (1 x nDim)  string, one name per dimension
%       opts.dimWeights   (1 x nDim)  positive weights, default ones(1, nDim)
%
%   This class is the SINGLE OWNER of box validity: ordering and finiteness are
%   checked here and nowhere else, so `rep.bounds` stays a pure lookup and
%   the two cannot drift apart. Its errors name the offending dimension, because a
%   caller reaching them has a bad declaration and a row index does not say which.

    properties (SetAccess = immutable)
        nDim       (1,1) double   % dimension of the regressor space
        box        (:,2) double   % (nDim x 2) physical bounds, [lower upper]
        dimName    (1,:) string   % (1 x nDim) one name per dimension
        dimWeights (1,:) double   % (1 x nDim) normalized extent per dimension
        offset     (1,:) double   % (1 x nDim) box(:,1).'
        span       (1,:) double   % (1 x nDim) (box(:,2) - box(:,1)).'
    end

    methods

        function obj = Normalization(box, dimName, opts)
            arguments
                box     (:,2) double {mustBeReal}
                dimName (1,:) string
                opts.dimWeights (1,:) double {mustBeReal} = ones(1, size(box, 1))
            end

            nDim = size(box, 1);

            if nDim < 1
                error("rid:roi:invalidBox", ...
                    "box must have at least one row (one dimension).");
            end
            % The dimension-name check comes FIRST so that the box errors below can
            % name the offending dimension. This class is the single owner of box
            % validity -- ordering and finiteness are checked here and nowhere else
            %, so its diagnostics have to be readable on their own: a caller
            % that reaches this point has a bad DECLARATION, and "box(3,2)" does not
            % tell them which one.
            if numel(dimName) ~= nDim
                error("rid:roi:dimensionMismatch", ...
                    "dimName must have one entry per dimension (got %d for nDim = %d).", ...
                    numel(dimName), nDim);
            end
            if ~all(isfinite(box), "all")
                bad = find(~all(isfinite(box), 2), 1);
                error("rid:roi:invalidBox", ...
                    "The bounds declared for ""%s"" (dimension %d) are [%g, %g] " + ...
                    "and must be finite; the ROI is a declared specification, " + ...
                    "so an unbounded dimension is not declarable.", ...
                    dimName(bad), bad, box(bad, 1), box(bad, 2));
            end
            if any(box(:, 2) <= box(:, 1))
                bad = find(box(:, 2) <= box(:, 1), 1);
                error("rid:roi:invalidBox", ...
                    "The bounds declared for ""%s"" (dimension %d) are [%g, %g]; " + ...
                    "the upper bound must be strictly greater than the lower.", ...
                    dimName(bad), bad, box(bad, 1), box(bad, 2));
            end
            if numel(opts.dimWeights) ~= nDim
                error("rid:roi:dimensionMismatch", ...
                    "dimWeights must have one entry per dimension (got %d for nDim = %d).", ...
                    numel(opts.dimWeights), nDim);
            end
            if ~all(isfinite(opts.dimWeights)) || any(opts.dimWeights <= 0)
                error("rid:roi:invalidWeights", ...
                    "dimWeights must be finite and strictly positive.");
            end

            obj.nDim       = nDim;
            obj.box        = box;
            obj.dimName    = dimName;
            obj.dimWeights = opts.dimWeights;
            obj.offset     = box(:, 1).';
            obj.span       = (box(:, 2) - box(:, 1)).';
        end

        function Znorm = toNormalized(obj, Zphys)
            %TONORMALIZED  Physical coordinates to normalized coordinates.
            %
            %   Inputs:
            %       Zphys   (nSamples x nDim)  physical coordinates
            %   Outputs:
            %       Znorm  (nSamples x nDim)  normalized coordinates; values outside
            %                              [0, dimWeights] are legal and returned
            %                              unclamped
            arguments
                obj
                Zphys (:,:) double {mustBeReal}
            end
            obj.assertColumns(Zphys, "Zphys");
            Znorm = (Zphys - obj.offset) ./ obj.span .* obj.dimWeights;
        end

        function Zphys = toPhysical(obj, Znorm)
            %TOPHYSICAL  Normalized coordinates to physical coordinates.
            %
            %   Inputs:
            %       Znorm  (nSamples x nDim)  normalized coordinates
            %   Outputs:
            %       Zphys   (nSamples x nDim)  physical coordinates
            arguments
                obj
                Znorm (:,:) double {mustBeReal}
            end
            obj.assertColumns(Znorm, "Znorm");
            Zphys = Znorm ./ obj.dimWeights .* obj.span + obj.offset;
        end

        function bn = normalizedBox(obj)
            %NORMALIZEDBOX  The ROI in normalized coordinates.
            %
            %   This is the extent the supporting points discretize and over which
            %   the fill distance is a supremum. It is the only thing `+roi`'s
            %   geometry functions need from this object, besides `nDim`.
            %
            %   Outputs:
            %       bn  (nDim x 2)  [zeros(nDim,1), dimWeights.']
            bn = [zeros(obj.nDim, 1), obj.dimWeights.'];
        end

        function s = describe(obj)
            %DESCRIBE  Human-readable summary, for logs and run records.
            %
            %   Outputs:
            %       s  (1x1) string
            lines = strings(obj.nDim, 1);
            for kDim = 1:obj.nDim
                lines(kDim) = sprintf("  %-16s [%g, %g]  ->  [0, %g]", ...
                    obj.dimName(kDim), obj.box(kDim, 1), obj.box(kDim, 2), ...
                    obj.dimWeights(kDim));
            end
            s = sprintf("Normalization, nDim = %d\n%s", obj.nDim, ...
                strjoin(lines, newline));
        end

    end

    methods (Access = private)

        function assertColumns(obj, P, name)
            %ASSERTCOLUMNS  Point sets are (nPoints x nDim); empty sets are legal.
            if size(P, 2) ~= obj.nDim
                error("rid:roi:dimensionMismatch", ...
                    "%s must have nDim = %d columns (got %d). " + ...
                    "Point sets are (nPoints x nDim), samples in rows.", ...
                    name, obj.nDim, size(P, 2));
            end
        end

    end

end
