classdef Narx < rid.rep.Representation
%NARX  External-dynamics (NARX) representation: lagged inputs and outputs.
%
%   The coordinates are built from the RAW SIGNALS, so `dependsOnSurrogate` is
%   FALSE -- always, in both pastSource settings. That flag marks a representation
%   whose COORDINATES are a function of the surrogate parameters (an estimated
%   state, an OBF pole, a transformed axis); it does not mark the fact that the
%   surrogate supplies the output trajectory, which every representation uses.
%   The rebuild rule reads BOTH this flag and
%   `pastSource` -- see requiresRebuildOnSurrogateUpdate in the base class.
%
%   COORDINATES. For input delay `nk` and lag orders `nLagU`, `nLagY`:
%
%       z(k) = [ u_1(k-nk) ... u_1(k-nk-nLagU+1),   % channel-major,
%                ...                                 % lag-minor
%                u_p(k-nk) ... u_p(k-nk-nLagU+1),
%                y(k-1) ... y(k-nLagY) ]
%
%   The OUTPUT lags are always strictly past: y(k) is what a downstream model
%   predicts, so it can never be one of its own regressors. The INPUT block starts
%   at `nk`: nk = 0 includes the current input u(k), nk = 1
%   starts at u(k-1). Whether the space contains u(k) is a property of the
%   REPRESENTATION and is independent of whether the surrogate has a direct term --
%   the two deliberately need not coincide, because the input coordinates
%   are decision variables, not predictions. The asymmetric mismatch (surrogate
%   with a direct term, representation without u(k)) is reported, not refused.
%
%   THE DIFFERENCE BASIS (`inputDifference`). With `inputDifference = true` the
%   input block of each channel becomes
%
%       [ u_c(k-nk), du_c(k-nk), du_c(k-nk-1), ... ],  du_c(k) = u_c(k) - u_c(k-1)
%
%   -- the first coordinate is the level, the rest are successive differences. This
%   is a change of BASIS of the same block, not a different representation: nDim,
%   nWarmup and the ROI expansion pattern are identical either way, which is what
%   makes the comparison a comparison of geometry rather than of size.
%
%   Why the difference rather than the second lag: with [u(k), u(k-1)] a
%   rate limit makes the reachable set a diagonal band inside a square box, so much
%   of the discretization can never be visited yet still consumes supporting points
%   and dilutes J. With [u(k), du(k)] the same reachable set is an axis-aligned box
%   and essentially every supporting point is useful.
%
%   The three candidate representations are all expressible here:
%
%       [u(k-1), y(k-1)]           Narx(1, 1, "nk", 1, "inputDifference", false)
%       [u(k), u(k-1), y(k-1)]     Narx(2, 1, "nk", 0, "inputDifference", false)
%       [u(k), du(k), y(k-1)]      Narx(2, 1, "nk", 0, "inputDifference", true)
%
%   COORDINATE NAMES. One name per DIMENSION, in coordinate order:
%   `requiredNames()` returns exactly nDim strings and they are the same strings
%   `bounds` returns as `dimName`. The convention is
%
%       <base>_k        the coordinate at lag 0
%       <base>_km<n>    the coordinate at lag n >= 1
%
%   with `base` one of "u1" ... "up" (input channels, always subscripted because
%   their number is a free parameter), "du1" ... "dup" (difference coordinates), or
%   "y" (unsubscripted: the method description fixes y as (nSteps x 1), so q = 1 project-wide). The default representation therefore declares
%   `u1_k`, `du1_k`, `y_km1`.
%
%   The lag is IN the name rather than implied by position. A name that carried only
%   an ordinal -- `u1`, `u1_1` -- would mean u(k) at nk = 0 and u(k-1) at nk = 1, so
%   changing the input delay would silently rebind an existing declaration to a
%   different lag: an off-by-one that shifts every coordinate with nothing to detect
%   it, the same failure `kIndex` exists to prevent on the other side of this
%   boundary. With the lag in the name, changing nk changes which names are
%   required and a stale declaration fails loudly.
%
%   These are code names, not display labels. Rendering `u1_km1` as "u1(k-1)" is a
%   +viz concern; a second name for the same quantity is not.
%
%   EVERY BOUND IS DECLARED, PER DIMENSION. `bounds` is an ordered lookup
%   over the declaration: no expansion, no inheritance across lags, no derivation.
%   Two lags of the same channel may legitimately carry different ranges, and a lag
%   whose range is absent is an error rather than a copy of its neighbour's.
%
%   THE LAG WARM-UP: rows are DROPPED, never padded. The first
%   `nWarmup = max(nk + nLagU - 1, nLagY)` samples have no complete regressor
%   vector. A padded row would be a fabricated point that enters CoverageState and
%   claims coverage of a region the process never visited -- a data-bias error
%   manufactured by bookkeeping, in a criterion whose whole purpose is measuring
%   where the data actually is. `build` therefore returns nSteps - nWarmup rows
%   together with `kIndex`, the time index of each row, so the correspondence is
%   testable rather than implied. A history shorter than the warm-up returns
%   (0 x nDim) and is NOT an error: at the start of every run there genuinely is no
%   complete regressor yet.
%
%   Construction:
%       rep = rid.rep.Narx(nLagU, nLagY)
%       rep = rid.rep.Narx(nLagU, nLagY, "nk", 0, "nInputs", 1, ...
%                          "inputDifference", true, "pastSource", "measured")
%
%   Inputs:
%       nLagU                 (1x1)  input coordinates per channel, >= 1
%       nLagY                 (1x1)  output lags, >= 0 (0 gives an NFIR space)
%       opts.nk               (1x1)  input delay, >= 0
%       opts.nInputs          (1x1)  input channels p, >= 1
%       opts.inputDifference  (1x1)  logical, difference basis for the input block
%       opts.pastSource       (1x1)  "measured" or "simulated"
%
%   THE HORIZON FLOOR: `minHorizon` is `nk + 1`. The applied sample
%   u(N+1) first appears in the regressor row at k = N+1+nk, so a horizon of nk
%   samples or fewer scores nothing that depends on it. Run assembly refuses a
%   shorter one; the loop asks the representation rather than reading `nk`.

    properties (SetAccess = immutable)
        nDim                % (1x1) double  = nInputs*nLagU + nLagY
        dependsOnSurrogate  % (1x1) logical ALWAYS false: built from raw signals
        nLagU              (1,1) double
        nLagY              (1,1) double
        nk                 (1,1) double
        nInputs            (1,1) double
        inputDifference    (1,1) logical
        nWarmup            (1,1) double   % samples with no complete regressor vector
        inputName          (1,:) string   % (1 x nInputs) base names "u1" ... "up"
        diffName           (1,:) string   % (1 x nInputs) base names "du1" ... "dup"
        coordInfo          (1,:) struct   % (1 x nDim) ONE ENTRY PER DIMENSION,
                                          % in coordinate order: .name .kind
                                          % .channel .lag. The single source of
                                          % both the name and its meaning
    end

    methods

        function obj = Narx(nLagU, nLagY, opts)
            arguments
                nLagU (1,1) double {mustBeInteger, mustBePositive}
                nLagY (1,1) double {mustBeInteger, mustBeNonnegative}
                opts.nk              (1,1) double {mustBeInteger, mustBeNonnegative} ...
                    = rid.defaultOptions().rep.nk
                opts.nInputs         (1,1) double {mustBeInteger, mustBePositive} ...
                    = rid.defaultOptions().rep.nInputs
                opts.inputDifference (1,1) logical ...
                    = rid.defaultOptions().rep.inputDifference
                opts.pastSource      (1,1) string ...
                    = rid.defaultOptions().rep.pastSource
            end

            obj@rid.rep.Representation(opts.pastSource);

            % The flag must not be a silent no-op: with a single input coordinate
            % there is no difference to take, and accepting it would let a
            % configuration mean something other than what it says.
            if opts.inputDifference && nLagU < 2
                error("rid:rep:invalidRepresentation", ...
                    "inputDifference = true needs nLagU >= 2; with nLagU = 1 " + ...
                    "the input block is a single level and there is no " + ...
                    "difference to form.");
            end

            obj.nLagU           = nLagU;
            obj.nLagY           = nLagY;
            obj.nk              = opts.nk;
            obj.nInputs         = opts.nInputs;
            obj.inputDifference = opts.inputDifference;

            obj.dependsOnSurrogate = false;
            obj.nDim               = opts.nInputs * nLagU + nLagY;
            obj.nWarmup            = max(opts.nk + nLagU - 1, nLagY);

            obj.inputName = "u" + string(1:opts.nInputs);
            if opts.inputDifference
                obj.diffName = "du" + string(1:opts.nInputs);
            else
                obj.diffName = strings(1, 0);
            end

            % The coordinates are built here, once, in exactly the order `build`
            % fills the columns of Zphys, so that coordinateInfo, bounds and build
            % cannot drift apart into three orderings. Each entry carries its own
            % KIND, CHANNEL and LAG, so no consumer has to recover them from the
            % name -- which is what made the lag `nk + iLag` unrecoverable
            % downstream, since nothing outside this constructor bounds it.
            info = repmat(rid.rep.Narx.coordinateEntry("x", "state", NaN, NaN), ...
                          1, obj.nDim);
            iDim = 0;
            for cIn = 1:opts.nInputs
                for iLag = 0:(nLagU - 1)
                    iDim = iDim + 1;
                    if opts.inputDifference && iLag > 0
                        thisLag    = opts.nk + iLag - 1;
                        info(iDim) = rid.rep.Narx.coordinateEntry( ...
                            rid.rep.Narx.coordinateName(obj.diffName(cIn), thisLag), ...
                            "inputDifference", cIn, thisLag);
                    else
                        thisLag    = opts.nk + iLag;
                        info(iDim) = rid.rep.Narx.coordinateEntry( ...
                            rid.rep.Narx.coordinateName(obj.inputName(cIn), thisLag), ...
                            "inputLevel", cIn, thisLag);
                    end
                end
            end
            for iLag = 1:nLagY
                iDim = iDim + 1;
                info(iDim) = rid.rep.Narx.coordinateEntry( ...
                    rid.rep.Narx.coordinateName("y", iLag), "output", NaN, iLag);
            end
            obj.coordInfo = info;
        end

        function info = coordinateInfo(obj)
            %COORDINATEINFO  What each coordinate IS, as data rather than as a name.
            %
            %   Exactly nDim entries in coordinate order, so `[info.name]` is the
            %   `dimName` that `bounds` returns and `requiredNames` (inherited,
            %   never overridden) is derived from it: one vocabulary, and one source
            %   for its meaning.
            %
            %   A NARX representation has no "state" coordinates; that kind exists
            %   for the StateSpace and is deliberately part of the shared closed
            %   set rather than of one subclass.
            %
            %   Outputs:
            %       info  (1 x nDim) struct  fields .name .kind .channel .lag
            info = obj.coordInfo;
        end

        function n = minHorizon(obj)
            %MINHORIZON  nk + 1: the first horizon row containing u(N+1).
            %
            %   The input block starts at u(k-nk), so the sample about to be applied
            %   first appears in the regressor row at k = N+1+nk -- the (nk+1)-th
            %   row of the horizon. With a shorter horizon `J` does not depend on
            %   the applied sample at all, and the loop would commit whatever the
            %   winning start happened to hold.
            %
            %   The LAG ORDERS do not enter. nLagU and nLagY deepen the warm-up,
            %   which is where the window starts, not where u(N+1) first appears.
            %
            %   Outputs:
            %       n  (1x1) double  = nk + 1
            n = obj.nk + 1;
        end

        function [box, dimName] = bounds(obj, roi)
            %BOUNDS  Ordered lookup of the per-dimension declared bounds.
            %
            %   Every dimension carries its own declared bounds and nothing is
            %   derived from anything else, so this is a lookup, not a rule. Two
            %   lags of the same channel may carry different ranges; a missing lag
            %   is an error, never a copy of its neighbour's.
            %
            %   Inputs:
            %       roi      (1x1) struct   one (1 x 2) range per requiredNames entry
            %   Outputs:
            %       box      (nDim x 2)     physical bounds, [lower upper]
            %       dimName  (1 x nDim)     the coordinate names, == requiredNames()
            arguments
                obj
                roi (1,1) struct
            end
            dimName = obj.requiredNames();
            box     = obj.lookupRanges(roi, dimName);
        end

        function [Zphys, kIndex] = build(obj, hist, srg)
            %BUILD  Physical regressor points from the signal history.
            %
            %   Row i is the regressor vector at time kIndex(i). The first nWarmup
            %   samples have no complete vector and are DROPPED, never padded
            %   (the warm-up row rule, see the class header).
            %
            %   Inputs:
            %       hist    (1x1) struct  .u .y .tS [.x0]; see the base class header
            %       srg     (1x1) object  surrogate; used only to obtain yHat when
            %                             pastSource == "simulated"
            %   Outputs:
            %       Zphys       (nRows x nDim)  PHYSICAL coordinates, samples in rows
            %       kIndex  (nRows x 1)     time index k of each row
            arguments
                obj
                hist (1,1) struct
                srg  = []
            end

            [u, y] = obj.signalsFromHistory(hist, srg, obj.nInputs);
            obj.reportFeedthrough(srg, obj.nk);

            nSteps = size(u, 1);
            nRows  = max(0, nSteps - obj.nWarmup);
            kIndex = (obj.nWarmup + 1 : nSteps).';
            Zphys      = zeros(nRows, obj.nDim);
            if nRows == 0
                return
            end

            iDim = 0;
            for cIn = 1:obj.nInputs
                for iLag = 0:(obj.nLagU - 1)
                    iDim = iDim + 1;
                    if obj.inputDifference && iLag > 0
                        % Successive difference at lag (nk + iLag - 1):
                        % u(k - nk - iLag + 1) - u(k - nk - iLag).
                        lagNew = obj.nk + iLag - 1;
                        Zphys(:, iDim) = u(kIndex - lagNew, cIn) ...
                                   - u(kIndex - lagNew - 1, cIn);
                    else
                        Zphys(:, iDim) = u(kIndex - obj.nk - iLag, cIn);
                    end
                end
            end
            for iLag = 1:obj.nLagY
                iDim = iDim + 1;
                Zphys(:, iDim) = y(kIndex - iLag);
            end
        end

        function s = describe(obj)
            %DESCRIBE  Human-readable summary, for logs and run records.
            %
            %   Outputs:
            %       s  (1x1) string
            basis = "lags";
            if obj.inputDifference
                basis = "level + differences";
            end
            s = sprintf(...
                "rid.rep.Narx  nDim = %d  (nInputs = %d, nLagU = %d, nLagY = %d, " + ...
                "nk = %d)\n  input basis : %s\n  pastSource  : %s\n" + ...
                "  warm-up     : %d sample(s) dropped\n  coordinates : %s\n" + ...
                "  dependsOnSurrogate = %d, requiresRebuildOnSurrogateUpdate = %d", ...
                obj.nDim, obj.nInputs, obj.nLagU, obj.nLagY, obj.nk, basis, ...
                obj.pastSource, obj.nWarmup, ...
                strjoin(obj.requiredNames(), ", "), ...
                obj.dependsOnSurrogate, obj.requiresRebuildOnSurrogateUpdate());
        end

    end

    methods (Static)
        function name = coordinateName(base, lag)
            %COORDINATENAME  The project's coordinate-naming convention, in one place.
            %
            %   `<base>_k` at lag 0, `<base>_km<n>` at lag n >= 1. Public and static
            %   because the convention is needed outside this class: +gen matches a
            %   `du<X>_k` admissible bound to the rate limit of channel X,
            %   and +viz renders these names for display. A second spelling of
            %   the same coordinate anywhere would violate the method description.
            %
            %   Inputs:
            %       base  (1x1) string  "u1", "du2", "y", ...
            %       lag   (1x1)         non-negative integer
            %   Outputs:
            %       name  (1x1) string
            arguments
                base (1,1) string
                lag  (1,1) double {mustBeInteger, mustBeNonnegative}
            end
            if lag == 0
                name = base + "_k";
            else
                name = base + "_km" + string(lag);
            end
        end
    end

end
