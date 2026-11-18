classdef Constraints
%CONSTRAINTS  The declared admissible region, classified into hard and soft limits.
%
%   THE ADMISSIBLE DECLARATION IS THE CONSTRAINT DECLARATION. Constraints are
%   enforced against the ADMISSIBLE region and never against the ROI:
%   the ROI is a ruler, not a fence. The
%   admissible region is declared per dimension, in exactly the ROI's format and
%   over exactly the same coordinate names, so this class consumes one
%   declaration struct and derives nothing from anything else. Every number it
%   holds is a number the user wrote down.
%
%   WHAT IS HARD AND WHAT IS SOFT. The split follows what the optimizer can
%   manipulate, and is a switch on the coordinate's declared KIND -- the table
%   verbatim, with no inference in between:
%
%       info(i).kind        class              enforcement
%       ------------------------------------------------------------------------
%       "inputLevel"        decision variable  HARD: lb / ub
%       "inputDifference"   decision variable  HARD: A*theta <= b
%       "output"            simulated          SOFT: penalty in the objective,
%       "state"             simulated                never a term in J
%
%   NAMES ARE NOT A TYPE SYSTEM. This class asks the representation what
%   each coordinate IS; it never regenerates a name and matches the string. The
%   earlier version did, which needed a lag bound that no consumer can know -- and
%   silently classified an input level as soft whenever `nk` pushed a lag past
%   `nDim`. `channel` comes from the same place, so `p` is never held here and a
%   representation and its constraints cannot disagree about how many inputs exist.
%
%   As a consequence `+gen` no longer names `rid.rep.*` at all: it consumes the
%   Representation interface, so the StateSpace works here unchanged.
%
%   A soft limit is a limit on a quantity the optimizer only PREDICTS -- an output,
%   a state, an output difference -- which is a nonlinear function of the decision
%   variables through the surrogate and therefore not expressible as a bound. It
%   enters the OBJECTIVE as a smooth penalty (rid.gen.Direct/penalty), never as a
%   term in J: J is a measurement of distribution quality and must stay comparable
%   across runs.
%
%   THE DEFAULT RATE LIMIT IS THE WIDEST RATE THE DECLARATION PERMITS: the
%   admissible bound on `du<c>_k` when the representation has that coordinate, and
%   otherwise +-(uMax_c - uMin_c). This invents nothing -- both numbers are already
%   declared -- and it is finite, which [-inf, inf] was not: an infinity propagates
%   into the random-walk construction and makes the rate-consistency check
%   vacuous. It is a DEFAULT, not a derivation: `opts.rateLimit` overrides it, and
%   rid.gen.reachabilityReport then compares the two. Under the default the rate
%   limit is non-binding, so a default run does not get the axis-aligned-box benefit
%   of the feedthrough choice until a real rate limit is declared.
%
%   CONTRADICTORY PER-LAG DECLARATIONS ARE AN ERROR, not an intersection.
%   `u1_k` and `u1_km1` are the same actuator at two instants; declaring them over
%   different boxes is a statement that contradicts itself, and silently taking one
%   of them -- or their intersection -- would be the silent acceptance the project's coding standard
%   forbids.
%
%   Construction:
%       cons = rid.gen.Constraints(rep, admissible)
%       cons = rid.gen.Constraints(rep, admissible, ...
%                                  "rateLimit", struct("u1", [-0.4 0.4]))
%
%   Inputs:
%       rep              (1x1) object  representation; supplies requiredNames and
%                                      validates the declaration through bounds
%       admissible       (1x1) struct  one (1 x 2) physical [lower upper] per
%                                      coordinate name, exactly as the ROI is
%                                      declared
%       opts.rateLimit   (1x1) struct  one (1 x 2) [duMin duMax] per CHANNEL base
%                                      name ("u1" ... "up"); overrides the direct-parametrization default
%                                      default for the channels it names
%
%   There is deliberately NO `nInputs` option: the channel count comes from the
%   representation's coordinates. A second declaration of `p` could disagree
%   with the first, and the disagreement was silent -- the surplus channels became
%   soft constraints, so an actuator lost its rate limit without any error.
%
%   PHYSICAL UNITS throughout: +gen works in physical units and normalization
%   happens at the representation boundary.

    properties (SetAccess = immutable)
        nDim          (1,1) double   % dimensions of the representation
        dimName       (1,:) string   % (1 x nDim) coordinate names, coordinate order
        admissibleBox (:,2) double   % (nDim x 2) the declaration, verbatim
        channelName   (1,:) string   % (1 x nChannels) "u1" ... "up", derived from
                                     % the representation's coordinates
        amplitude     (:,2) double   % (nChannels x 2) hard box on u, physical
        rate          (:,2) double   % (nChannels x 2) hard box on du, physical
        rateIsDeclared (1,:) logical % (1 x nChannels) did the user declare it?
        rateSource    (1,:) string   % (1 x nChannels) "declared" | "duCoordinate"
                                     %                | "amplitudeWidth"
        levelDim      (1,:) double   % (1 x nChannels) column of the regressor
                                     % array holding the SMALLEST-lag level
                                     % of each channel
        soft          (1,:) struct   % (1 x nSoft) .name .iDim .gMid .gHalf .s
    end

    methods

        function obj = Constraints(rep, admissible, opts)
            arguments
                rep
                admissible (1,1) struct
                opts.rateLimit (1,1) struct = struct()
            end

            % The declaration is validated by the representation that owns the
            % names: missing, unknown and mis-shaped entries fail with the declaration-validation rule
            % identifiers, so there is one validator for the ROI and the admissible
            % region alike.
            [box, names] = rep.bounds(admissible);

            % What each coordinate IS, as data. Everything below switches on
            % this; nothing below rebuilds a name to find out.
            info = rep.coordinateInfo();

            obj.nDim          = numel(names);
            obj.dimName       = names;
            obj.admissibleBox = box;

            % The admissible region is a declared specification, so an unbounded or
            % inverted interval is not declarable. This class owns the invariant for
            % the ADMISSIBLE box; rid.roi.Normalization owns it for the ROI, and the
            % two boxes are never the same object.
            for iDim = 1:obj.nDim
                if ~all(isfinite(box(iDim, :))) || box(iDim, 2) <= box(iDim, 1)
                    error("rid:gen:invalidAdmissibleBox", ...
                        "The admissible bounds declared for ""%s"" are [%g, %g]; " + ...
                        "they must be finite with the upper bound strictly " + ...
                        "greater than the lower.", names(iDim), ...
                        box(iDim, 1), box(iDim, 2));
                end
            end

            kind    = [info.kind];
            channel = [info.channel];
            lag     = [info.lag];

            % The channel count is READ, never declared here: `p` has one source.
            isInputCoord = kind == "inputLevel" | kind == "inputDifference";
            channels     = unique(channel(isInputCoord));
            nChannels    = numel(channels);
            if nChannels == 0
                error("rid:gen:missingAmplitude", ...
                    "The representation has no input coordinate, so there is " + ...
                    "nothing for the generator to constrain. Its coordinates " + ...
                    "are %s.", strjoin("""" + names + """", ", "));
            end

            obj.channelName    = "u" + string(channels);
            obj.amplitude      = zeros(nChannels, 2);
            obj.rate           = zeros(nChannels, 2);
            obj.rateIsDeclared = false(1, nChannels);
            obj.rateSource     = strings(1, nChannels);
            obj.levelDim       = nan(1, nChannels);

            for cIn = 1:nChannels
                thisChannel = channels(cIn);

                % Level coordinates of this channel, smallest lag first. Sorted
                % EXPLICITLY by the declared lag rather than trusting coordinate
                % order, so "the amplitude of the channel" always means the same
                % instant whatever order a representation emits.
                levelIdx = find(kind == "inputLevel" & channel == thisChannel);
                [~, ord] = sort(lag(levelIdx));
                levelIdx = levelIdx(ord);
                if isempty(levelIdx)
                    error("rid:gen:missingAmplitude", ...
                        "Channel ""%s"" has difference coordinates but no level " + ...
                        "coordinate, so its amplitude limit cannot be declared. " + ...
                        "The representation's coordinates are %s.", ...
                        obj.channelName(cIn), strjoin("""" + names + """", ", "));
                end
                obj.assertOneBox(box, names, levelIdx, ...
                    "rid:gen:inconsistentAmplitude", ...
                    "amplitude of channel """ + obj.channelName(cIn) + """");
                obj.levelDim(cIn)    = levelIdx(1);
                obj.amplitude(cIn,:) = box(levelIdx(1), :);

                % Difference coordinates of this channel. All lags of du are hard:
                % the band matrix constrains every consecutive pair, so one rate box
                % covers them all -- which is also why a disagreement between lags
                % is a contradiction rather than a refinement.
                diffIdx = find(kind == "inputDifference" & channel == thisChannel);
                if ~isempty(diffIdx)
                    obj.assertOneBox(box, names, diffIdx, ...
                        "rid:gen:inconsistentRate", ...
                        "rate of channel """ + obj.channelName(cIn) + """");
                end

                % declared limit, else the du_k admissible bound, else the
                % full amplitude width. Never [-inf, inf]. The lag-0 difference is
                % found by its DECLARED lag, not by rebuilding "du<c>_k".
                duZero = diffIdx(lag(diffIdx) == 0);
                if isfield(opts.rateLimit, obj.channelName(cIn))
                    r = opts.rateLimit.(obj.channelName(cIn));
                    obj.rate(cIn, :)        = obj.checkRateRow(r, obj.channelName(cIn));
                    obj.rateIsDeclared(cIn) = true;
                    obj.rateSource(cIn)     = "declared";
                elseif ~isempty(duZero)
                    obj.rate(cIn, :)    = box(duZero(1), :);
                    obj.rateSource(cIn) = "duCoordinate";
                else
                    width = obj.amplitude(cIn, 2) - obj.amplitude(cIn, 1);
                    obj.rate(cIn, :)    = [-width, width];
                    obj.rateSource(cIn) = "amplitudeWidth";
                end
                if obj.rate(cIn, 1) > 0 || obj.rate(cIn, 2) < 0
                    error("rid:gen:invalidRateLimit", ...
                        "The rate limit of channel ""%s"" is [%g, %g], which " + ...
                        "excludes a zero increment. Holding the input constant " + ...
                        "would then be infeasible, and the warm start of the warm-start initialization " + ...
                        "would need repair.", obj.channelName(cIn), ...
                        obj.rate(cIn, 1), obj.rate(cIn, 2));
                end
            end

            % Everything the optimizer only PREDICTS becomes a soft constraint, so
            % that every declared admissible bound is enforced somewhere: hard on
            % the decision variables, soft on the simulated quantities.
            % `kind` is a closed set, so this partition is total by construction --
            % a coordinate can no longer fall through into the wrong class.
            softIdx  = find(kind == "output" | kind == "state");
            obj.soft = repmat(struct("name", "", "iDim", 0, ...
                "gMid", 0, "gHalf", 0, "s", 0), 1, numel(softIdx));
            for iSoft = 1:numel(softIdx)
                iDim  = softIdx(iSoft);
                gHalf = (box(iDim, 2) - box(iDim, 1)) / 2;
                obj.soft(iSoft) = struct( ...
                    "name",  names(iDim), ...
                    "iDim",  iDim, ...
                    "gMid",  (box(iDim, 1) + box(iDim, 2)) / 2, ...
                    "gHalf", gHalf, ...
                    "s",     2 * gHalf);
            end

            % Unknown channels in the rate declaration are a typo, not a no-op.
            declaredChannels = string(fieldnames(opts.rateLimit)).';
            unknown = setdiff(declaredChannels, obj.channelName, "stable");
            if ~isempty(unknown)
                error("rid:gen:unknownRateLimit", ...
                    "rateLimit names %s, which is not an input channel. " + ...
                    "Expected one of %s.", ...
                    strjoin("""" + unknown + """", ", "), ...
                    strjoin("""" + obj.channelName + """", ", "));
            end
        end

        function s = describe(obj)
            %DESCRIBE  Human-readable summary, for logs and run records.
            %
            %   Outputs:
            %       s  (1x1) string
            lines = strings(numel(obj.channelName), 1);
            for cIn = 1:numel(obj.channelName)
                lines(cIn) = sprintf( ...
                    "  %-4s amplitude [%g, %g] from ""%s"" | rate [%g, %g] (%s)", ...
                    obj.channelName(cIn), obj.amplitude(cIn, 1), ...
                    obj.amplitude(cIn, 2), obj.dimName(obj.levelDim(cIn)), ...
                    obj.rate(cIn, 1), obj.rate(cIn, 2), obj.rateSource(cIn));
            end
            softLines = strings(numel(obj.soft), 1);
            for iSoft = 1:numel(obj.soft)
                softLines(iSoft) = sprintf("  %-16s [%g, %g]", obj.soft(iSoft).name, ...
                    obj.soft(iSoft).gMid - obj.soft(iSoft).gHalf, ...
                    obj.soft(iSoft).gMid + obj.soft(iSoft).gHalf);
            end
            if isempty(softLines), softLines = "  (none)"; end
            s = sprintf("rid.gen.Constraints  nDim = %d, channels = %d\n" + ...
                "hard:\n%s\nsoft (simulated quantities):\n%s", ...
                obj.nDim, numel(obj.channelName), strjoin(lines, newline), ...
                strjoin(softLines, newline));
        end

    end

    methods (Static, Access = private)

        function assertOneBox(box, names, idx, identifier, what)
            %ASSERTONEBOX  All lags of one signal must carry the same declared box.
            if numel(idx) < 2, return, end
            same = all(box(idx, :) == box(idx(1), :), 2);
            if ~all(same)
                bad = idx(find(~same, 1));
                error(identifier, ...
                    "The %s is declared twice and inconsistently: ""%s"" has " + ...
                    "[%g, %g] but ""%s"" has [%g, %g]. They are the same " + ...
                    "physical signal at two instants, so the declaration " + ...
                    "contradicts itself.", what, names(idx(1)), ...
                    box(idx(1), 1), box(idx(1), 2), names(bad), ...
                    box(bad, 1), box(bad, 2));
            end
        end

        function r = checkRateRow(r, channel)
            %CHECKRATEROW  Shape and validity of one declared rate limit.
            if ~isnumeric(r) || ~isreal(r) || ~isequal(size(r), [1 2]) ...
                    || ~all(isfinite(r)) || r(2) <= r(1)
                error("rid:gen:invalidRateLimit", ...
                    "The rate limit declared for channel ""%s"" must be a real, " + ...
                    "finite (1 x 2) row [duMin duMax] with duMax > duMin.", channel);
            end
            r = double(r);
        end

    end

end
