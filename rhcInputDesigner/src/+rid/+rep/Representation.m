classdef Representation
%REPRESENTATION  Abstract contract for a dynamics representation.
%
%   A REPRESENTATION IS A COORDINATE SYSTEM, NOT A MODEL.
%   The surrogate answers "where would the
%   process go if I applied this input?"; the representation answers "in which
%   coordinates is coverage measured?". They are related -- `defaultRepresentation`
%   offers the representation matching a model's structure -- but they are NOT
%   required to match, and any representation may be paired with any surrogate. An
%   industrial model with fifty states is perfectly usable behind a
%   three-dimensional regressor space, and forcing the two to agree would demand a
%   fifty-dimensional criterion.
%
%   "USES THE SURROGATE" IS NOT "DEPENDS ON THE SURROGATE". Every representation
%   uses the surrogate to predict the future output of a candidate. Only some
%   depend on it for the COORDINATES THEMSELVES -- an estimated state, an OBF pole,
%   a transformed axis. `dependsOnSurrogate` marks the second, and it is that flag
%   which drives the rebuild rule, never the first.
%
%   The contract every representation supplies:
%
%       names          = rep.requiredNames()      % (1 x nNames) string
%       [box, dimName] = rep.bounds(roi)          % (nDim x 2), (1 x nDim) string
%       [Zphys, kIndex]    = rep.build(hist, srg)     % (nRows x nDim) PHYSICAL, (nRows x 1)
%       nDim           = rep.nDim
%       flag           = rep.dependsOnSurrogate
%       n              = rep.minHorizon()         % (1x1) shortest SCORING horizon
%       tf             = rep.needsMeasuredStates()  % (1x1) logical
%       rep.assertCompatible(srg)                 % refuse a surrogate this space
%                                                 % cannot be built from
%       s              = rep.describe()
%
%   THE DECLARATION is a struct with one field per DIMENSION -- `requiredNames`
%   returns exactly nDim names, in coordinate order, and they are the same strings
%   `bounds` returns as `dimName`. One string serves as declaration key, coordinate
%   name and label:
%
%       roi = struct("u1_k", [0 5], "du1_k", [-0.5 0.5], "y_km1", [-1 2]);
%
%   EVERY BOUND IS DECLARED, PER DIMENSION; NOTHING IS DERIVED. There is no
%   inheritance across lags, no derivation from a rate limit, and no expansion
%   pattern -- which is the only rule that generalizes, since inheriting a range
%   across lags means something for NARX and nothing for a joint input--state-space
%   representation where each state has its own units. `bounds` is therefore an
%   ordered lookup: a missing declaration is an error naming what is missing
%   (rid:rep:missingRange), never a silently defaulted range.
%
%   The same rule governs the admissible region and the rate constraints,
%   and an inconsistent set of
%   declarations is caught by report, not by construction (the validate-never-correct rule).
%
%   THE HISTORY is the single source of truth. `build` receives it as a
%   struct -- the History object of +rhc will satisfy the same shape:
%
%       hist.u   (nSteps x nInputs)  applied input, PHYSICAL
%       hist.y   (nSteps x 1) or []  measured output, PHYSICAL; [] when offline
%       hist.tS  (1x1)               sample time, physical
%       hist.x0  (1 x nStates)       run's initial surrogate state; required only
%                                    when pastSource == "simulated"
%       hist.x   (nSteps x nStates)  measured states, present only when a process
%                                    reported them. Row k is the
%                                    state at the instant u(k) was applied
%
%   VALUE SEMANTICS. Representations are immutable value objects. Mutable state
%   exists in exactly two places in this design: rid.crit.CoverageState and the
%   update schedulers.
%
%   PACKAGE DEPENDENCY. `+rep` never names `rid.surrogate.*`: `srg` is used
%   duck-typed, so the documented dependency direction (rep -> roi only) holds and
%   a surrogate that is not a rid.surrogate.Surrogate subclass is usable.
%
%   COORDINATE NAMES ARE NOT A TYPE SYSTEM. `coordinateInfo` states what
%   each coordinate IS -- its kind, its input channel, its lag -- as data.
%   Consumers switch on `.kind`; they never infer meaning by regenerating a name
%   and matching the string. A name is for humans and for keying a declaration;
%   control flow reads the info struct. `coordinateName` remains the single
%   PRODUCER of the name; `coordinateInfo` is the single source of its
%   MEANING.

    properties (Abstract, SetAccess = immutable)
        % Declared without size/class validators: MATLAB does not permit them on a
        % property that implements an inherited abstract one, and repeating them in
        % every subclass is not an option either. Shapes:
        nDim                % (1x1) double   dimension of the regressor space
        dependsOnSurrogate  % (1x1) logical  do the COORDINATES depend on the surrogate?
    end

    properties (SetAccess = immutable)
        % Where the already-acquired output coordinates come from. It is
        % NEVER inferred from online/offline: the mode is not a property of the
        % coordinates.
        pastSource (1,1) string
    end

    methods

        function obj = Representation(pastSource)
            %REPRESENTATION  Base constructor; every subclass forwards pastSource here.
            arguments
                pastSource (1,1) string
            end
            if ~ismember(pastSource, ["measured", "simulated"])
                error("rid:rep:invalidPastSource", ...
                    "pastSource must be ""measured"" or ""simulated"" (got ""%s"").", ...
                    pastSource);
            end
            obj.pastSource = pastSource;
        end

        function tf = requiresRebuildOnSurrogateUpdate(obj)
            %REQUIRESREBUILDONSURROGATEUPDATE  The four-row rebuild table, in code.
            %
            %   A rebuild is required when the surrogate changed AND the historical
            %   coordinates depend on it. Two
            %   independent things make them depend on it, and the code reads BOTH
            %   -- never the online/offline mode:
            %
            %       pastSource == "simulated"    the past outputs are simulated
            %       dependsOnSurrogate           the coordinates themselves are
            %                                    reconstructed (estimated state, OBF
            %                                    pole, transformed axis)
            %
            %   These are deliberately NOT collapsed into one flag. Collapsing them
            %   would make rows 3 and 4 of that table indistinguishable and would
            %   force every future representation to re-OR them itself.
            %
            %   Outputs:
            %       tf  (1x1) logical
            tf = obj.dependsOnSurrogate || obj.pastSource == "simulated";
        end

        function n = minHorizon(~)
            %MINHORIZON  The shortest nHorizon that SCORES the applied input.
            %
            %   CONCRETE HERE, returning 1, and overridden only by a representation
            %   whose input block does not start at u(k). A representation that
            %   contains the current input scores the sample about to be committed
            %   in its very first horizon row, so 1 is the honest floor.
            %
            %   WHY THE REPRESENTATION STATES IT. With a NARX input delay the
            %   applied sample u(N+1) first appears in the regressor row at
            %   k = N+1+nk, so a horizon of nk samples or fewer scores nothing that
            %   depends on it: `J` is constant in the applied sample, the optimizer
            %   returns whatever its start held, and an arbitrary input is committed
            %   without an error. `nHorizon` has no default, so nothing else
            %   stopped a run from choosing 1 at nk = 1.
            %
            %   The loop must not branch on a representation class, and only the
            %   representation knows where its input block starts -- so `+rhc` reads
            %   this rather than `nk`, which is a property of one subclass.
            %
            %   WHAT IT DOES NOT GUARANTEE. For the INPUT coordinates the bound is
            %   exact. For the PREDICTED ones it is necessary but not sufficient: a
            %   surrogate with an internal dead time delays the applied input's
            %   effect on them further, and the loop cannot know that.
            %
            %   Outputs:
            %       n  (1x1) double  shortest admissible nHorizon, >= 1
            n = 1;
        end

        function tf = needsMeasuredStates(~)
            %NEEDSMEASUREDSTATES  Must the loop ask the process for states?
            %
            %   CONCRETE HERE, returning false, and overridden only by a
            %   representation whose coordinates are measured states --
            %   rid.rep.StateSpace under pastSource = "measured".
            %
            %   IT IS WHAT THE LOOP BRANCHES ON, and it exists so that the loop
            %   does not branch on a representation class. Two things follow from
            %   it: the process is called
            %   for a SECOND output, `[y, x] = process.apply(u)`, and each horizon
            %   is anchored one surrogate step from the last measured state rather
            %   than from a re-derived one. A process that measures no states
            %   keeps the one-output form it has always had.
            %
            %   IT IS NOT "hist.x EXISTS". Anchoring for every representation
            %   that happens to have measured states available was rejected: a
            %   NARX surrogate's state basis is arbitrary, so a measured state is
            %   not a valid initial condition for it.
            %
            %   Outputs:
            %       tf  (1x1) logical
            tf = false;
        end

        function assertCompatible(~, ~)
            %ASSERTCOMPATIBLE  Refuse a surrogate this space cannot be built from.
            %
            %   CONCRETE HERE AND A NO-OP, because for most representations there
            %   is nothing to refuse: the coordinates are arithmetic on recorded
            %   signals and any surrogate can supply the predicted ones. It is
            %   overridden where the pairing carries a real precondition --
            %   rid.rep.StateSpace needs a surrogate whose states ARE the states
            %   it declares.
            %
            %   WHY IT IS A SEPARATE MEMBER AND NOT A CHECK INSIDE `build`.
            %   rid.rhc.Setup calls it, so an incompatible pairing is refused at
            %   assembly rather than at step 1 of a run that has already started a
            %   process. `build` checks the same condition again, because a
            %   surrogate update may change the state dimension after assembly
            %   passed -- the two are not redundant, they answer the
            %   question at two different times.
            %
            %   WHAT IT CANNOT CHECK is the state BASIS. A dimension is
            %   comparable; a coordinate system is not. That assumption is stated
            %   in the joint-space representation contract and carried by the user.
            %
            %   Inputs:
            %       srg  (1x1) object  the surrogate the run will use
            %   Outputs:
            %       none; throws when the pairing is impossible
        end

        function names = requiredNames(obj)
            %REQUIREDNAMES  One name per dimension, in coordinate order.
            %
            %   CONCRETE, and derived from `coordinateInfo`, so that a name and its
            %   meaning cannot come from two places. A subclass implements
            %   `coordinateInfo` and gets this for free; it must not override it.
            %
            %   Outputs:
            %       names  (1 x nDim) string
            names = [obj.coordinateInfo().name];
        end

    end

    methods (Static, Access = protected)

        function e = coordinateEntry(name, kind, channel, lag)
            %COORDINATEENTRY  One well-formed `coordinateInfo` entry.
            %
            %   The field order, the field types and the closed set of `kind`
            %   values live HERE and nowhere else, so that the second
            %   representation (the StateSpace) cannot invent a fifth kind or a
            %   different field order and quietly break every consumer that
            %   switches on `.kind`.
            %
            %   `channel` and `lag` are NaN where the notion does not apply -- an
            %   output has no input channel, a state coordinate has no lag. NaN and
            %   not 0: zero is a valid channel-free lag, and a zero channel index
            %   would silently alias.
            %
            %   Inputs:
            %       name     (1x1) string   the coordinate name, from the naming
            %                               convention's single producer
            %       kind     (1x1) string   "inputLevel" | "inputDifference"
            %                               | "output" | "state"
            %       channel  (1x1) double   input channel index, NaN if not an input
            %       lag      (1x1) double   lag in samples, NaN if not applicable
            %   Outputs:
            %       e        (1x1) struct   fields .name .kind .channel .lag
            arguments
                name    (1,1) string
                kind    (1,1) string {mustBeMember(kind, ["inputLevel", ...
                            "inputDifference", "output", "state"])}
                channel (1,1) double
                lag     (1,1) double
            end
            e = struct("name", name, "kind", kind, ...
                       "channel", channel, "lag", lag);
        end

    end

    methods (Abstract)
        info           = coordinateInfo(obj)
        [box, dimName] = bounds(obj, roi)
        [Zphys, kIndex]    = build(obj, hist, srg)
        s              = describe(obj)
    end

    methods (Access = protected)

        function box = lookupRanges(~, roi, names)
            %LOOKUPRANGES  Ordered lookup of declared ranges. No expansion.
            %
            %   Every dimension carries its OWN declared bounds and nothing is
            %   derived from anything else, so this is a lookup and not a rule:
            %   there is no inheritance across lags, no derivation from a rate
            %   limit, and no expansion pattern to get wrong. `names` is in
            %   coordinate order and has exactly nDim entries, so the returned box
            %   is in coordinate order too.
            %
            %   THE OWNERSHIP SPLIT. This function decides which names must be
            %   present; it does NOT decide whether a range is a valid box. Ordering
            %   and finiteness belong to rid.roi.Normalization, which owns the box
            %   invariant and is already tested for it -- one owner, no drift. What
            %   is checked here is only the shape, so that a mis-shaped declaration
            %   fails legibly instead of as a cryptic assignment error.
            %
            %   Inputs:
            %       roi    (1x1) struct     one (1 x 2) [lower upper] field per name
            %       names  (1 x nNames)     string, in coordinate order
            %   Outputs:
            %       box    (nNames x 2)     physical bounds, in the order of `names`
            nNames = numel(names);
            box    = zeros(nNames, 2);

            declared = string(fieldnames(roi)).';
            missing  = setdiff(names, declared, "stable");
            if ~isempty(missing)
                error("rid:rep:missingRange", ...
                    "The declaration has no range for %s. Every dimension is " + ...
                    "declared separately and nothing is inherited from another " + ...
                    "dimension, so all %d names must be present: %s.", ...
                    strjoin("""" + missing + """", ", "), nNames, ...
                    strjoin("""" + names + """", ", "));
            end

            % A per-dimension declaration is longer than a per-signal one and so
            % more prone to typos, and an ignored field is exactly the silent
            % acceptance the project forbids. The typo shows up twice -- as a
            % missing name and as an unknown one -- and both are reported.
            unknown = setdiff(declared, names, "stable");
            if ~isempty(unknown)
                error("rid:rep:unknownRange", ...
                    "The declaration contains %s, which is not a dimension of " + ...
                    "this representation. Expected exactly: %s.", ...
                    strjoin("""" + unknown + """", ", "), ...
                    strjoin("""" + names + """", ", "));
            end

            for iName = 1:nNames
                r = roi.(names(iName));
                if ~isnumeric(r) || ~isreal(r) || ~isequal(size(r), [1 2])
                    error("rid:rep:invalidRange", ...
                        "The range declared for ""%s"" must be a real (1 x 2) " + ...
                        "row [lower upper] (got %s of class %s).", ...
                        names(iName), mat2str(size(r)), class(r));
                end
                box(iName, :) = double(r);
            end
        end

        function [u, y] = signalsFromHistory(obj, hist, srg, nInputs)
            %SIGNALSFROMHISTORY  The (u, y) pair the coordinates are built from.
            %
            %   Applies pastSource: "measured" takes the output from the
            %   process, "simulated" from the surrogate. The choice is explicit and
            %   is never inferred from whether a measurement happens to exist.
            %
            %   The surrogate is used ONLY to obtain the output trajectory;
            %   nothing about its order or
            %   internal structure reaches the coordinates, which is what lets any
            %   representation be paired with any surrogate.
            %
            %   Inputs:
            %       hist     (1x1) struct  see the class header
            %       srg      (1x1) object  surrogate; used only when simulating
            %       nInputs  (1x1)         input channels this representation expects
            %   Outputs:
            %       u  (nSteps x nInputs)  physical input
            %       y  (nSteps x 1)        physical output, measured or simulated

            % STRUCTURAL FIRST. A malformed history must not reach the
            % pastSource guards below, which would then say something confusing
            % and specific about a struct that is simply not a history.
            rid.assertHistory(hist);
            u = double(hist.u);
            if size(u, 2) ~= nInputs
                error("rid:rep:dimensionMismatch", ...
                    "hist.u must have nInputs = %d columns (got %d). " + ...
                    "Samples are in rows.", nInputs, size(u, 2));
            end
            nSteps = size(u, 1);

            switch obj.pastSource
                case "measured"
                    % SEMANTIC, and it stays here: the condition is pastSource,
                    % which only the representation knows.
                    if ~isfield(hist, "y") || isempty(hist.y)
                        error("rid:rep:missingMeasurement", ...
                            "pastSource = ""measured"" needs hist.y, which is " + ...
                            "empty. Offline runs must set pastSource = " + ...
                            """simulated"" explicitly -- the mode is never inferred.");
                    end
                    y = double(hist.y);
                    % assertHistory has already matched y's row count to u's; the
                    % SINGLE-OUTPUT convention is a property of the current model
                    % class, not of a history, so it is checked here.
                    if ~isequal(size(y), [nSteps 1])
                        error("rid:rep:dimensionMismatch", ...
                            "hist.y must be (nSteps x 1) = (%d x 1) (got %s).", ...
                            nSteps, mat2str(size(y)));
                    end

                case "simulated"
                    % `isfield` is NOT tested: assertHistory above requires x0 of
                    % every history, so an absent one cannot reach here and the
                    % pointer to steadyState for that case moved there with it
                    %. What remains is the case only pastSource can judge:
                    % a history that HAS an x0 with no states in it.
                    if isempty(hist.x0)
                        error("rid:rep:missingInitialState", ...
                            "pastSource = ""simulated"" needs hist.x0, the run's " + ...
                            "initial surrogate state. Set it explicitly, e.g. " + ...
                            "hist.x0 = srg.steadyState(hist.u(1,:)).");
                    end
                    y = srg.simulate(u, hist.x0);
                    y = double(y);
                    if ~isequal(size(y), [nSteps 1])
                        error("rid:rep:dimensionMismatch", ...
                            "srg.simulate returned %s; the contract is " + ...
                            "(nSteps x 1) = (%d x 1).", mat2str(size(y)), nSteps);
                    end
            end
        end

        function reportFeedthrough(~, srg, nk)
            %REPORTFEEDTHROUGH  Report, never refuse, the asymmetric mismatch.
            %
            %   A surrogate WITH a direct term against a representation WITHOUT u(k)
            %   leaves y(k) depending on a coordinate the space omits: the covered
            %   space under-determines the process, and a downstream model built on
            %   it is short of an input it needs. Nothing breaks computationally --
            %   the surrogate simulates from the signal, not from the regressor -- so
            %   this is a modelling inconsistency to be REPORTED, not refused.
            %
            %   The reverse combination is legitimate and silent: a strictly proper
            %   surrogate with u(k) in the space is the DEFAULT, because the input
            %   coordinates are decision variables, not predictions.
            %
            %   This is a correctness diagnostic, not a progress message, so it does
            %   not fall under the ban on console output in +rep (same reasoning as
            %   the floor warning). It cannot fire in the default configuration, where
            %   directTerm = 0.
            if nk == 0, return, end
            if ~isobject(srg) || ~isprop(srg, "directTerm"), return, end
            if any(srg.directTerm ~= 0)
                warning("rid:rep:feedthroughNotRepresented", ...
                    "The surrogate has a non-zero direct term, but this " + ...
                    "representation starts at u(k-%d) and so omits u(k). " + ...
                    "y(k) then depends on a coordinate the regressor space does " + ...
                    "not contain; the covered space under-determines the process " + ...
                    ". Set nk = 0 to represent it.", nk);
            end
        end

    end

end
