classdef StateSpace < rid.rep.Representation
%STATESPACE  Joint input--state-space representation: the input and the states at k.
%
%   The second of the two representation families the design commits to,
%   and neither is the reference case for the
%   other. The coordinates are
%
%       z(k) = [ u_1(k) ... u_p(k), x_1(k) ... x_n(k) ]
%
%   -- the input levels of every channel, then the states, all at lag 0. The row
%   with kIndex = k is [u(k), x(k)], where x(k) is the state BEFORE u(k) acts.
%
%   THERE IS NO `nk`, and an `nk` argument is refused rather than ignored. u(k) is
%   always present: it carries local persistent excitation, and the model
%   this space serves, x(k+1) = f(x(k), u(k)), needs it. `nk` belongs to the
%   external-dynamics input block -- the state-space coordinates are [u(k), x(k)]
%   whatever nk a NARX representation in the same study uses, so an nk = 1 NARX is
%   compared as [u(k-1), ..., y(k-1), ...] against [u(k), x(k)]. [u(k-1), x(k)]
%   would pair a state with an input it already contains. This is the reason
%   rid.rep.Narx refuses a difference flag it cannot honour: a setting that
%   silently does nothing says something false.
%
%   THE LOOK-BACK DEPTH IS ZERO. A point needs only u(k) and x(k), and x(1) is
%   known, so nWarmup = 0 and `build` returns ONE ROW PER SAMPLE with
%   kIndex = (1:nSteps)'. Two things follow:
%
%     - the one-row initial history of a run ALREADY yields a point, [u(1), x0] --
%       the rest state under the pre-run input, realized data that the history-seeding rule counts
%       before step 1. The xStart requirement's "a single row cannot form a regressor vector" is a
%       statement about NARX shapes with a lag and is false here;
%     - the candidate window is EMPTY, so the loop simulates a zero-length
%       input once per step, which the surrogate contract explicitly allows.
%
%   WHERE THE PAST STATES COME FROM is `pastSource`, exactly as for the past
%   outputs of rid.rep.Narx:
%
%       "simulated"   xTraj from srg.simulateStates(hist.u, hist.x0)
%       "measured"    hist.x, the states the process reported
%
%   `dependsOnSurrogate` is FALSE in both settings, as for Narx. That simulated
%   states move with the surrogate is already what pastSource says; the flag is
%   reserved for coordinates computed WITH the surrogate even though measurements
%   exist -- states estimated by an observer from measured y, which is not built.
%   `requiresRebuildOnSurrogateUpdate` therefore reduces to
%   pastSource == "simulated".
%
%   WHAT THIS REPRESENTATION ASSUMES AND CANNOT CHECK. The declared state
%   coordinates are the PROCESS's, and the surrogate's states are in the same
%   coordinates. Without that, a measured state is not an initial condition for
%   the surrogate and a simulated state is not comparable with a measured one.
%   `assertCompatible` checks the dimension, which is all that is checkable;
%   nothing can check the basis. This is why a state-space representation needs an
%   INTERPRETABLE model, and why no surrogate offers this class as its
%   `defaultRepresentation` -- a class cannot know whether its own states are
%   interpretable, so the representation is always constructed explicitly.
%
%   COORDINATE NAMES CARRY THE TIME INDEX: `x1_k`, not `x1` (the joint-space representation contract,
%   amending the coordinate-meaning rule). A state coordinate does have a time index, and x(k) versus
%   x(k+1) is exactly the silent off-by-one that writing the lag into the name
%   prevents for inputs. The names come from the single producer,
%   rid.rep.Narx.coordinateName, so there is no second naming rule to drift.
%
%   EVERY BOUND IS DECLARED, PER DIMENSION, in physical units:
%
%       roi = struct("u1_k", [0 5], "x1_k", [-1 1], "x2_k", [-3 3]);
%
%   Nothing is derived. Each state has its own physical range and its own units,
%   which is precisely the case for which no expansion rule exists.
%
%   Construction:
%       rep = rid.rep.StateSpace(nStates)
%       rep = rid.rep.StateSpace(nStates, "nInputs", 2, "pastSource", "measured")
%
%   Inputs:
%       nStates          (1x1)  state coordinates n, >= 1
%       opts.nInputs     (1x1)  input channels p, >= 1
%       opts.pastSource  (1x1)  "measured" or "simulated"
%
%   THE HORIZON FLOOR is the inherited 1: u(k) is in every row, so the
%   sample about to be committed is scored in the first horizon row.
%
%   THIS CLASS HAS NO DIFFERENCE COORDINATES AND NO STATE SUBSET. Either can be added
%   later without changing these names.

    properties (SetAccess = immutable)
        nDim                % (1x1) double  = nInputs + nStates
        dependsOnSurrogate  % (1x1) logical ALWAYS false; see the class header
        nStates            (1,1) double
        nInputs            (1,1) double
        nWarmup            (1,1) double   % 0; kept for `describe`, since the loop
                                          % reads the depth from kIndex
        coordInfo          (1,:) struct   % (1 x nDim) ONE ENTRY PER DIMENSION,
                                          % in coordinate order: .name .kind
                                          % .channel .lag
    end

    methods

        function obj = StateSpace(nStates, opts)
            arguments
                nStates (1,1) double {mustBeInteger, mustBePositive}
                opts.nInputs    (1,1) double {mustBeInteger, mustBePositive} ...
                    = rid.defaultOptions().rep.nInputs
                opts.pastSource (1,1) string ...
                    = rid.defaultOptions().rep.pastSource
            end

            obj@rid.rep.Representation(opts.pastSource);

            obj.nStates = nStates;
            obj.nInputs = opts.nInputs;

            obj.dependsOnSurrogate = false;
            obj.nDim               = opts.nInputs + nStates;
            obj.nWarmup            = 0;

            % Built here, once, in exactly the order `build` fills the columns,
            % so that coordinateInfo, bounds and build cannot drift into three
            % orderings. Each entry carries its own KIND, CHANNEL and LAG, so no
            % consumer recovers them from the name.
            info = repmat(rid.rep.StateSpace.coordinateEntry( ...
                       "x", "state", NaN, 0), 1, obj.nDim);
            for cIn = 1:opts.nInputs
                info(cIn) = rid.rep.StateSpace.coordinateEntry( ...
                    rid.rep.Narx.coordinateName("u" + cIn, 0), ...
                    "inputLevel", cIn, 0);
            end
            for iState = 1:nStates
                % Lag 0, and the SAME producer the inputs use: a state at k is
                % x1_k, named alike and declared alike.
                info(opts.nInputs + iState) = rid.rep.StateSpace.coordinateEntry( ...
                    rid.rep.Narx.coordinateName("x" + iState, 0), ...
                    "state", NaN, 0);
            end
            obj.coordInfo = info;
        end

        function info = coordinateInfo(obj)
            %COORDINATEINFO  What each coordinate IS, as data rather than as a name.
            %
            %   Exactly nDim entries in coordinate order, so `[info.name]` is the
            %   `dimName` that `bounds` returns and `requiredNames` (inherited,
            %   never overridden) is derived from it.
            %
            %   A state coordinate has kind "state", no channel, and LAG 0 -- not
            %   NaN. It is the state at k, and the instant it belongs to is part
            %   of what the coordinate is.
            %
            %   Outputs:
            %       info  (1 x nDim) struct  fields .name .kind .channel .lag
            info = obj.coordInfo;
        end

        function tf = needsMeasuredStates(obj)
            %NEEDSMEASUREDSTATES  True exactly when the past states are measured.
            %
            %   The loop reads this to decide whether to ask the process for a
            %   second output and whether to anchor each horizon at the last
            %   measured state. It is `pastSource` and nothing
            %   else: a second switch for the source of past STATES could
            %   contradict the one that already says where past predicted
            %   coordinates come from.
            %
            %   Outputs:
            %       tf  (1x1) logical
            tf = obj.pastSource == "measured";
        end

        function assertCompatible(obj, srg)
            %ASSERTCOMPATIBLE  The surrogate must have the states this space declares.
            %
            %   THE DIMENSION IS WHAT IS CHECKABLE. The real precondition is that
            %   the surrogate's states are the process's, in the declared
            %   coordinates -- without it a measured state is not an initial
            %   condition for the surrogate and a simulated state is not
            %   comparable with a measured one. No code can check a
            %   basis. Refusing the dimension mismatch at least makes the
            %   representable failures loud.
            %
            %   A MODEL THAT DECLARES NO nStates IS REFUSED TOO, and for the same
            %   reason it is refused by `simulateStates`: NaN is not a state count
            %   that happens to disagree, it is the absence of one, and NaN ~= n
            %   would otherwise decide the question by accident.
            %
            %   Called by rid.rhc.Setup at assembly and by `build` at every step,
            %   because a surrogate update may change the state dimension.
            %
            %   Inputs:
            %       srg  (1x1) object  the surrogate; duck-typed, as everywhere in
            %                          +rep (the package never names +surrogate)
            %   Outputs:
            %       none; throws rid:rep:stateDimensionMismatch
            if ~isobject(srg) || ~isprop(srg, "nStates")
                error("rid:rep:stateDimensionMismatch", ...
                    "This representation declares %d state coordinate(s) and " + ...
                    "needs a surrogate that states how many states it has. The " + ...
                    "object supplied (%s) has no nStates.", ...
                    obj.nStates, class(srg));
            end
            if ~isscalar(srg.nStates) || ~isfinite(srg.nStates) || ...
                    srg.nStates ~= obj.nStates
                error("rid:rep:stateDimensionMismatch", ...
                    "The surrogate has nStates = %s and this representation " + ...
                    "declares %d state coordinate(s). They must agree: the " + ...
                    "coordinates ARE the surrogate's states, and nothing can " + ...
                    "check that they are also the process's.", ...
                    mat2str(srg.nStates), obj.nStates);
            end
        end

        function [box, dimName] = bounds(obj, roi)
            %BOUNDS  Ordered lookup of the per-dimension declared bounds.
            %
            %   Every dimension carries its own declared bounds and nothing is
            %   derived. Here that is not merely a convention: each state has its
            %   own physical units, so there is no expansion pattern that could
            %   even be written down.
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
            %   Zphys = [hist.u, X] with X from `pastSource`, and
            %   kIndex = (1:nSteps)': the look-back depth is zero, so no row is
            %   dropped and every sample of the history is a point.
            %
            %   RE-SIMULATING FROM hist.x0 WITH THE CURRENT SURROGATE is the xStart re-derivation:
            %   the state trajectory is derived at every call and never carried
            %   forward, because a carried state is not merely stale after a
            %   surrogate update, it is type-incorrect. Carrying xTraj forward is
            %   possible and is adopted only on profiling evidence.
            %
            %   IN MEASURED MODE `srg` IS NEVER TOUCHED. Passing [] is legal and
            %   is how a recorded history is evaluated offline: measured
            %   state-space fill distance is rep.build over a history carrying
            %   hist.x, with no loop and no model.
            %
            %   Inputs:
            %       hist  (1x1) struct  .u .x0 .tS [.x]; see the base class header
            %       srg   (1x1) object  surrogate; used only when the states are
            %                           simulated
            %   Outputs:
            %       Zphys   (nSteps x nDim)  PHYSICAL coordinates, samples in rows
            %       kIndex  (nSteps x 1)     time index k of each row
            arguments
                obj
                hist (1,1) struct
                srg  = []
            end

            % STRUCTURAL FIRST: a malformed history must not reach the
            % pastSource guards, which would then say something confusing and
            % specific about a struct that is simply not a history.
            rid.assertHistory(hist);
            u = double(hist.u);
            if size(u, 2) ~= obj.nInputs
                error("rid:rep:dimensionMismatch", ...
                    "hist.u must have nInputs = %d columns (got %d). " + ...
                    "Samples are in rows.", obj.nInputs, size(u, 2));
            end
            nSteps = size(u, 1);

            switch obj.pastSource
                case "simulated"
                    % assertHistory already requires x0 of every history, so what
                    % remains is the case only pastSource can judge: a history
                    % that HAS an x0 with no states in it.
                    if isempty(hist.x0)
                        error("rid:rep:missingInitialState", ...
                            "pastSource = ""simulated"" needs hist.x0, the " + ...
                            "run's initial surrogate state. Set it explicitly, " + ...
                            "e.g. hist.x0 = srg.steadyState(hist.u(1,:)).");
                    end
                    % BEFORE simulating: a surrogate of the wrong order would
                    % otherwise fail inside simulateStates, or -- worse -- return
                    % a trajectory of the wrong width that only the shape check
                    % below would catch, with nothing to say about why.
                    obj.assertCompatible(srg);
                    [~, ~, X] = srg.simulateStates(u, hist.x0);
                    X = double(X);
                    if ~isequal(size(X), [nSteps obj.nStates])
                        error("rid:rep:dimensionMismatch", ...
                            "srg.simulateStates returned a %s trajectory; the " + ...
                            "contract is (nSteps x nStates) = (%d x %d).", ...
                            mat2str(size(X)), nSteps, obj.nStates);
                    end

                case "measured"
                    % SEMANTIC, and it stays here: whether a well-formed history
                    % must carry states is the representation's question, and the
                    % condition is pastSource.
                    if ~isfield(hist, "x")
                        error("rid:rep:missingStateMeasurement", ...
                            "pastSource = ""measured"" needs hist.x, the states " + ...
                            "the process reported. An absent field means no " + ...
                            "states were measured, and is never an empty " + ...
                            "placeholder -- a run whose states are not measured " + ...
                            "must set pastSource = ""simulated"" explicitly, " + ...
                            "since the mode is never inferred.");
                    end
                    X = double(hist.x);
                    % assertHistory has already matched x's row count to u's. How
                    % many states there are to measure is this class's question,
                    % so the WIDTH is checked here.
                    if size(X, 2) ~= obj.nStates
                        error("rid:rep:dimensionMismatch", ...
                            "hist.x must be (nSteps x nStates) = (%d x %d) " + ...
                            "(got %s).", nSteps, obj.nStates, mat2str(size(X)));
                    end
            end

            Zphys  = [u, X];
            kIndex = (1:nSteps).';
        end

        function s = describe(obj)
            %DESCRIBE  Human-readable summary, for logs and run records.
            %
            %   Outputs:
            %       s  (1x1) string
            if obj.needsMeasuredStates()
                source = "hist.x, as reported by the process";
            else
                source = "srg.simulateStates from hist.x0";
            end
            s = sprintf(...
                "rid.rep.StateSpace  nDim = %d  (nInputs = %d, nStates = %d)\n" + ...
                "  pastSource  : %s -- %s\n" + ...
                "  warm-up     : %d sample(s) dropped\n  coordinates : %s\n" + ...
                "  dependsOnSurrogate = %d, requiresRebuildOnSurrogateUpdate = %d", ...
                obj.nDim, obj.nInputs, obj.nStates, obj.pastSource, source, ...
                obj.nWarmup, strjoin(obj.requiredNames(), ", "), ...
                obj.dependsOnSurrogate, obj.requiresRebuildOnSurrogateUpdate());
        end

    end

end
