classdef CoordinateInfoContractTest < matlab.unittest.TestCase
%COORDINATEINFOCONTRACTTEST  The contract, over every representation shape.
%
%   ONE parameterized suite rather than a pile of per-class tests, so that adding a
%   representation means adding a PARAMETER and nothing else. Every invariant below
%   runs against every shape in `allRepresentations`, which is the point: `.kind`
%   has to be right across families and orders, not once for the default NARX.
%
%   The parameter set spans the two families the design commits to,
%   and BOTH are now real: external dynamics via
%   rid.rep.Narx, and the joint input--state-space family via rid.rep.StateSpace.
%   The shape double that stood in for the second is gone -- a double can only
%   ever agree with the contract it was written against, and this suite exists to
%   check the contract against what the code actually produces.
%
%   Invariants 1-8 are the representation's contract; invariant 9 is the consumer's
%   half of it -- that the hard/soft partition rid.gen.Constraints builds is TOTAL.
%   The two belong in one place because the coordinate-meaning rule exists precisely to keep them from
%   drifting apart.

    properties (Constant)
        Kinds = ["inputLevel" "inputDifference" "output" "state"]
    end

    properties (TestParameter)
        shape = CoordinateInfoContractTest.allRepresentations()
    end

    methods (Test)

        function nOneEntryPerDimension(tc, shape)
            % Invariant 1.
            rep = shape.make();
            tc.verifyEqual(numel(rep.coordinateInfo()), rep.nDim);
        end

        function requiredNamesIsExactlyTheInfoNamesInOrder(tc, shape)
            % Invariant 2 -- the no-drift invariant. requiredNames is CONCRETE on
            % the base class; a subclass that overrode it could disagree with its
            % own coordinateInfo, and every consumer would then pick a side.
            rep = shape.make();
            tc.verifyEqual(rep.requiredNames(), [rep.coordinateInfo().name]);
        end

        function everyKindIsInTheClosedSet(tc, shape)
            % Invariant 3, first half.
            rep = shape.make();
            tc.verifyTrue(all(ismember([rep.coordinateInfo().kind], tc.Kinds)));
        end

        function channelAndLagApplyExactlyWhereTheyShould(tc, shape)
            % Invariant 4. NaN and not 0 where the notion does not apply: a zero
            % channel index would silently alias, and lag 0 is a real lag.
            %
            % THE STATE CLAUSE CHANGED WITH THE STATE-SPACE REPRESENTATION. The coordinate-meaning rule expected
            % state names `x1`, `x2` with lag NaN; a state coordinate does carry
            % a time index, and x(k) versus x(k+1) is exactly the silent
            % off-by-one that writing the lag into the name prevents for inputs
            %. So a state's lag is finite and non-negative, and
            % invariant 5 below then covers `x1_k` as well.
            rep  = shape.make();
            info = rep.coordinateInfo();
            for iDim = 1:numel(info)
                e = info(iDim);
                switch e.kind
                    case {"inputLevel", "inputDifference"}
                        tc.verifyTrue(isfinite(e.channel), e.name);
                        tc.verifyTrue(ismember(e.channel, 1:shape.nInputs), e.name);
                        tc.verifyTrue(isfinite(e.lag) && e.lag >= 0, e.name);
                    case "output"
                        tc.verifyTrue(isnan(e.channel), e.name);
                    case "state"
                        tc.verifyTrue(isnan(e.channel), e.name);
                        tc.verifyTrue(isfinite(e.lag) && e.lag >= 0, e.name);
                end
            end
        end

        function theLagIsTheTrueLagAndIsNotDerivedFromPosition(tc, shape)
            % Invariant 5, generalized: wherever the name carries a lag, the
            % DECLARED lag agrees with it. The test may read the naming convention;
            % production code may not. This is what catches a lag recovered
            % from a coordinate's index, or bounded by nDim.
            rep  = shape.make();
            info = rep.coordinateInfo();
            for iDim = 1:numel(info)
                fromName = CoordinateInfoContractTest.lagInName(info(iDim).name);
                if ~isnan(fromName)
                    tc.verifyEqual(info(iDim).lag, fromName, ...
                        "Declared lag disagrees with the name " + info(iDim).name);
                end
            end
        end

        function namesAreUnique(tc, shape)
            % Invariant 6. One name per quantity; a duplicate
            % would make the declaration ambiguous rather than merely odd.
            rep   = shape.make();
            names = rep.requiredNames();
            tc.verifyEqual(numel(unique(names)), numel(names));
        end

        function boundsReturnsTheSameNamesInTheSameOrder(tc, shape)
            % Invariant 7: the declaration key and the coordinate name are one
            % vocabulary, so bounds cannot introduce a second ordering.
            rep = shape.make();
            [box, dimName] = rep.bounds( ...
                CoordinateInfoContractTest.declare(rep.requiredNames()));
            tc.verifyEqual(dimName, [rep.coordinateInfo().name]);
            tc.verifyEqual(size(box), [rep.nDim, 2]);
        end

        function everyChannelHasAtLeastOneInputLevel(tc, shape)
            % Invariant 8: without a level coordinate a channel has no declarable
            % amplitude, and the generator cannot bound the variable it commits.
            rep  = shape.make();
            info = rep.coordinateInfo();
            isLevel = [info.kind] == "inputLevel";
            for cIn = 1:shape.nInputs
                tc.verifyTrue(any(isLevel & [info.channel] == cIn), ...
                    sprintf("channel %d has no inputLevel coordinate", cIn));
            end
        end

        function theHardSoftPartitionIsTotal(tc, shape)
            % Invariant 9, the consumer's half. Every coordinate lands in exactly
            % one bucket. This is the assertion both the coordinate-meaning rule defects broke: a missed
            % input level showed up as an extra soft constraint and nothing else.
            rep  = shape.make();
            info = rep.coordinateInfo();
            cons = rid.gen.Constraints(rep, ...
                CoordinateInfoContractTest.declare(rep.requiredNames()));

            kind     = [info.kind];
            nHardExp = sum(kind == "inputLevel" | kind == "inputDifference");
            nSoftExp = sum(kind == "output" | kind == "state");

            tc.verifyEqual(nHardExp + nSoftExp, rep.nDim);       % closed set, no gap
            tc.verifyEqual(numel(cons.soft), nSoftExp);
            tc.verifyEqual(rep.nDim - numel(cons.soft), nHardExp);
            % ... and by name, not only by count.
            tc.verifyEqual(sort([cons.soft.name]), ...
                sort([info(kind == "output" | kind == "state").name]));
            tc.verifyEqual(numel(cons.channelName), shape.nInputs);
        end

    end

    methods (Test)   % suite-level: about the parameter set itself

        function noKindOutsideTheClosedSetCanBeConstructed(tc)
            % Invariant 3, second half: the set is closed where entries are MADE,
            % not merely where they are read, so a fifth kind cannot arrive by
            % accident. "outputDifference" is the one a reader is most
            % likely to reach for, the table having six rows to four kinds.
            tc.verifyError( ...
                @() SurrogateDependentRep.makeEntry("q", "outputDifference", NaN, 1), ...
                "MATLAB:validators:mustBeMember");
        end

        function theParameterSetContainsALagPastNDim(tc)
            % Guards invariant 5 against going vacuous again. Without a shape whose
            % largest lag exceeds nDim - 1, "the lag is not bounded by nDim" is
            % satisfied everywhere by accident and the invariant is decoration.
            % A future edit that drops narxLagPastNDim fails HERE, with the reason.
            shapes = CoordinateInfoContractTest.allRepresentations();
            found  = false;
            for name = string(fieldnames(shapes)).'
                rep = shapes.(name).make();
                lag = [rep.coordinateInfo().lag];
                found = found || any(lag(~isnan(lag)) > rep.nDim - 1);
            end
            tc.verifyTrue(found, ...
                "No shape has a lag past nDim - 1, so invariant 5 cannot fail. " + ...
                "Only nk raises a lag without raising nDim -- keep a shape with one.");
        end

        function everyKindOccursSomewhereInTheParameterSet(tc)
            % A kind that no parameter produces is a kind nobody tests. This fails
            % the day a kind is added without a representation that emits it.
            shapes = CoordinateInfoContractTest.allRepresentations();
            seen   = strings(1, 0);
            for name = string(fieldnames(shapes)).'
                rep  = shapes.(name).make();
                seen = [seen, [rep.coordinateInfo().kind]];      %#ok<AGROW>
            end
            for thisKind = tc.Kinds
                tc.verifyTrue(ismember(thisKind, seen), ...
                    "No parameter produces kind """ + thisKind + """.");
            end
        end

        function requiredNamesIsNeverOverriddenByASubclass(tc)
            % The structural half of invariant 2: it must stay INHERITED, or the
            % single source is only a convention.
            for className = ["rid.rep.Narx", "rid.rep.StateSpace"]
                mc      = meta.class.fromName(className);
                methods = string({mc.MethodList.Name});
                owner   = arrayfun(@(m) string(m.DefiningClass.Name), mc.MethodList);
                tc.verifyEqual(owner(methods == "requiredNames"), ...
                    "rid.rep.Representation", ...
                    className + " must not redefine requiredNames.");
            end
        end

    end

    methods (Static, Access = private)

        function s = allRepresentations()
            %ALLREPRESENTATIONS  The parameter set. ADD A REPRESENTATION HERE.
            %
            %   `.make` constructs it; `.nInputs` is the channel count the shape is
            %   expected to have, stated independently of the object so that a
            %   representation reporting the wrong one cannot agree with itself.
            mk = @(fh, p) struct("make", fh, "nInputs", p);

            % --- external dynamics (rid.rep.Narx), real ----------------------
            s.narxSisoDefault     = mk(@() rid.rep.Narx(2, 1, "nk", 0, ...
                                        "inputDifference", true), 1);
            s.narxSisoNoDiff      = mk(@() rid.rep.Narx(2, 1, "nk", 0, ...
                                        "inputDifference", false), 1);
            s.narxSisoDelayed     = mk(@() rid.rep.Narx(1, 1, "nk", 1, ...
                                        "inputDifference", false), 1);
            s.narxHigherOrder     = mk(@() rid.rep.Narx(3, 2, "nk", 0, ...
                                        "inputDifference", true), 1);
            s.narxMiso            = mk(@() rid.rep.Narx(2, 1, "nInputs", 2, ...
                                        "inputDifference", true), 2);
            s.narxMisoHigherOrder = mk(@() rid.rep.Narx(3, 2, "nInputs", 3, ...
                                        "inputDifference", false), 3);
            s.narxNfir            = mk(@() rid.rep.Narx(2, 0), 1);  % no output coord

            % The shape that makes invariant 5 BITE. Every other knob -- channels,
            % output lags, input lags -- grows nDim at least as fast as the largest
            % lag, so "the lag is not bounded by nDim" is vacuously true across all
            % of them. Only `nk` raises the lag without touching nDim. Here
            % nDim = 2 while the lags are 2 and 3.
            s.narxLagPastNDim     = mk(@() rid.rep.Narx(2, 0, "nk", 2, ...
                                        "inputDifference", false), 1);

            % --- joint input--state space (rid.rep.StateSpace), real ---------
            % BOTH pastSource settings are parameters, because the coordinates
            % are the same in either and a representation that reported
            % different ones would be describing where its values came from
            % instead of what they are.
            s.stateSpaceSimulated = mk(@() rid.rep.StateSpace(2, ...
                                        "pastSource", "simulated"), 1);
            s.stateSpaceMeasured  = mk(@() rid.rep.StateSpace(2, ...
                                        "pastSource", "measured"), 1);
            s.stateSpaceMiso      = mk(@() rid.rep.StateSpace(3, "nInputs", 2, ...
                                        "pastSource", "simulated"), 2);
            s.stateSpaceFirstOrder = mk(@() rid.rep.StateSpace(1), 1);
        end

        function lag = lagInName(name)
            %LAGINNAME  The lag a name encodes, or NaN if it encodes none.
            %
            %   `<base>_k` is lag 0 and `<base>_km<n>` is lag n. A TEST may read the
            %   convention this way; production code may not, which is the
            %   whole point of the invariant this serves.
            lag = NaN;
            tok = regexp(name, '_km(\d+)$', 'tokens', 'once');
            if ~isempty(tok)
                lag = str2double(tok{1});
            elseif endsWith(name, "_k")
                lag = 0;
            end
        end

        function roi = declare(names)
            %DECLARE  A valid declaration for an arbitrary name list.
            %
            %   The box includes zero so that a difference coordinate never trips
            %   the zero-increment check for reasons unrelated to the test.
            roi = struct();
            for iName = 1:numel(names)
                roi.(names(iName)) = [-2 3];
            end
        end

    end

end
