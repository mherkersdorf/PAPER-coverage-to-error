classdef StateSpaceTest < matlab.unittest.TestCase
%STATESPACETEST  rid.rep.StateSpace: the joint input--state-space coordinates.
%
%   The second representation family. Its coordinates are [u1_k ... up_k, x1_k ... xn_k]: the input
%   levels and the states at k, all at lag 0, with the row at kIndex = k reading
%   [u(k), x(k)] where x(k) is the state BEFORE u(k) acts.
%
%   THREE THINGS THIS SUITE IS FOR, in the order they would go wrong:
%
%   1. THE ROW PAIRING. x(k) against u(k), never x(k+1) against u(k). That
%      off-by-one produces plausible numbers and no error, so it is asserted
%      against hand-computed values rather than against the surrogate's own
%      trajectory, and `kIndex` is asserted alongside it: the first sample
%      already carries a complete point, and a warm-up of one would silently
%      drop it.
%   2. WHERE THE PAST STATES COME FROM. `pastSource` and nothing else (the joint-space representation contract) -- "simulated" re-simulates from hist.x0, "measured" reads hist.x
%      and never touches the surrogate. The empty surrogate in the measured
%      tests is what proves the second half.
%   3. THE ASSUMPTION IT CANNOT CHECK. The surrogate's states must be the
%      process's, in the same coordinates. `assertCompatible` checks the
%      dimension, which is all that is checkable, and refuses the rest loudly
%      rather than building a point cloud in two different bases.

    methods (Test)

        % --- the simulated build -----------------------

        function simulatedCoordinatesArePairedWithTheInputThatActsOnThem(tc)
            % HAND VALUES, not a second call to the surrogate: the whole risk
            % here is a trajectory shifted by one sample, and comparing against
            % srg.simulateStates would shift the reference with it.
            %
            % MinimalSurrogate is x(k+1) = x(k).*[0.5 0.8] + u(k), so from
            % x(1) = [1 2] and u = [0.5; -1; 2; 0; 3]:
            rep  = rid.rep.StateSpace(2, "pastSource", "simulated");
            u    = [0.5; -1; 2; 0; 3];
            hist = wellFormedHistory(u, "x0", [1 2]);
            srg  = MinimalSurrogate();

            expected = [ 0.5,  1.0,  2.0
                        -1.0,  1.0,  2.1
                         2.0, -0.5,  0.68
                         0.0,  1.75, 2.544
                         3.0,  0.875, 2.0352 ];

            [Zphys, kIndex] = rep.build(hist, srg);
            tc.verifyEqual(Zphys, expected, "AbsTol", 1e-12);
            tc.verifyEqual(kIndex, (1:5).');
        end

        function theFirstStateRowIsTheRunsInitialStateExactly(tc)
            % Row 1 is x0 BIT-EXACTLY, because x(1) is known and nothing has
            % acted on it yet. An implementation that recorded the state after
            % each step would put x(2) here and still look plausible.
            rep  = rid.rep.StateSpace(2, "pastSource", "simulated");
            x0   = [0.3 -0.7];
            Zphys = rep.build(wellFormedHistory([1; 2; 3], "x0", x0), ...
                              MinimalSurrogate());
            tc.verifyEqual(Zphys(1, 2:3), x0);
        end

        function aOneRowHistoryAlreadyYieldsAPoint(tc)
            % the look-back depth is zero, so the initial history
            % of every run contributes [u(1), x0] -- the rest state under the
            % pre-run input, realized data that the history-seeding rule counts before step 1.
            % The xStart requirement's "a single row cannot form a regressor vector" is a NARX
            % statement and is false here.
            rep = rid.rep.StateSpace(2, "pastSource", "simulated");
            tc.verifyEqual(rep.nWarmup, 0);

            [Zphys, kIndex] = rep.build( ...
                wellFormedHistory(1.5, "x0", [0.25 -0.5]), MinimalSurrogate());
            tc.verifyEqual(Zphys, [1.5, 0.25, -0.5]);
            tc.verifyEqual(kIndex, 1);
        end

        function aSimulatedBuildNeedsANonEmptyInitialState(tc)
            rep  = rid.rep.StateSpace(2, "pastSource", "simulated");
            hist = struct("u", [1; 2], "x0", zeros(1, 0), "tS", 0.1);
            tc.verifyError(@() rep.build(hist, MinimalSurrogate()), ...
                "rid:rep:missingInitialState");
        end

        % --- the measured build -------------------------------

        function measuredStatesAreReadFromTheHistoryAndTheSurrogateIsNeverCalled(tc)
            % The EMPTY surrogate is the assertion: in measured mode there is
            % nothing to simulate, and a build that reached for `srg` anyway
            % would fail here rather than quietly re-deriving the states it was
            % handed.
            rep  = rid.rep.StateSpace(2, "pastSource", "measured");
            u    = [0.5; -1; 2];
            x    = [1 2; -3 4; 0.5 0.25];
            hist = StateSpaceTest.measuredHistory(u, x);

            [Zphys, kIndex] = rep.build(hist, []);
            tc.verifyEqual(Zphys, [u, x]);
            tc.verifyEqual(kIndex, (1:3).');
        end

        function aMeasuredBuildWithoutMeasuredStatesIsRefused(tc)
            % Whether a well-formed history MUST carry states is the
            % representation's question, not the validator's: absent
            % `.x` is legal for rid.assertHistory and an error here.
            rep = rid.rep.StateSpace(2, "pastSource", "measured");
            tc.verifyError(@() rep.build(wellFormedHistory([1; 2]), []), ...
                "rid:rep:missingStateMeasurement");
        end

        function aMeasuredStateOfTheWrongWidthIsRefused(tc)
            % The width is checked HERE because only the representation knows
            % how many states it declares (rid.assertHistory deliberately does
            % not check it).
            rep  = rid.rep.StateSpace(3, "pastSource", "measured");
            hist = StateSpaceTest.measuredHistory([1; 2], [1 2; 3 4]);
            tc.verifyError(@() rep.build(hist, []), "rid:rep:dimensionMismatch");
        end

        function anInputOfTheWrongWidthIsRefused(tc)
            rep  = rid.rep.StateSpace(2, "nInputs", 2, "pastSource", "measured");
            hist = StateSpaceTest.measuredHistory([1; 2], [1 2; 3 4]);
            tc.verifyError(@() rep.build(hist, []), "rid:rep:dimensionMismatch");
        end

        % --- names and declarations -------------

        function theNamesAreTheInputLevelsThenTheStatesAtK(tc)
            % The joint-space representation contract amends the coordinate-meaning rule: a state coordinate DOES carry a time
            % index, and it comes from the naming convention's single producer.
            rep = rid.rep.StateSpace(2);
            tc.verifyEqual(rep.requiredNames(), ["u1_k", "x1_k", "x2_k"]);
            tc.verifyEqual(rep.nDim, 3);
        end

        function aMisoRepresentationListsEveryChannelThenTheStates(tc)
            rep = rid.rep.StateSpace(3, "nInputs", 2);
            tc.verifyEqual(rep.requiredNames(), ...
                ["u1_k", "u2_k", "x1_k", "x2_k", "x3_k"]);
            tc.verifyEqual(rep.nDim, 5);
        end

        function aStateDeclaredWithoutItsTimeIndexIsRefusedTwice(tc)
            % The typo shows up as a MISSING name and, once the rest is present,
            % as an UNKNOWN one. Both are reported, because an ignored field is
            % exactly the silent acceptance the declared-bounds rule forbids.
            rep = rid.rep.StateSpace(2);
            tc.verifyError(@() rep.bounds( ...
                struct("u1_k", [0 1], "x1", [0 1], "x2_k", [0 1])), ...
                "rid:rep:missingRange");
            tc.verifyError(@() rep.bounds( ...
                struct("u1_k", [0 1], "x1_k", [0 1], "x2_k", [0 1], ...
                       "x1", [0 1])), ...
                "rid:rep:unknownRange");
        end

        function boundsIsAnOrderedLookupOverTheDeclaration(tc)
            rep = rid.rep.StateSpace(2);
            [box, dimName] = rep.bounds( ...
                struct("u1_k", [0 5], "x1_k", [-1 1], "x2_k", [-3 3]));
            tc.verifyEqual(dimName, ["u1_k", "x1_k", "x2_k"]);
            tc.verifyEqual(box, [0 5; -1 1; -3 3]);
        end

        % --- flags ------------------------------

        function anInputDelayArgumentIsRefused(tc)
            % `nk` belongs to the external-dynamics input block: [u(k), x(k)] is
            % the space whatever nk a NARX representation in the same study
            % uses. A setting that silently did nothing would say something
            % false. MATLAB's own unknown-argument error is the
            % refusal; its identifier is not part of this contract.
            tc.verifyError(@() rid.rep.StateSpace(2, "nk", 1), ?MException);
        end

        function theFlagsTheRebuildRuleReadsAreTheOnesD107States(tc)
            for src = ["measured", "simulated"]
                rep = rid.rep.StateSpace(2, "pastSource", src);
                % FALSE in BOTH settings, as for Narx: that simulated states
                % move with the surrogate is what pastSource already says. The
                % flag is reserved for states ESTIMATED by an observer.
                tc.verifyFalse(rep.dependsOnSurrogate, src);
                tc.verifyEqual(rep.requiresRebuildOnSurrogateUpdate(), ...
                    src == "simulated", src);
                tc.verifyEqual(rep.needsMeasuredStates(), ...
                    src == "measured", src);
            end
        end

        function theHorizonFloorIsOneBecauseEveryRowHoldsTheAppliedInput(tc)
            % u(k) is in every row, so the sample about to be committed
            % is scored in the first horizon row and the inherited floor stands.
            tc.verifyEqual(rid.rep.StateSpace(2).minHorizon(), 1);
        end

        % --- compatibility ------------------------------------

        function aSurrogateWithADifferentStateCountIsRefused(tc)
            % MinimalSurrogate has two states.
            rep = rid.rep.StateSpace(3, "pastSource", "simulated");
            tc.verifyError(@() rep.assertCompatible(MinimalSurrogate()), ...
                "rid:rep:stateDimensionMismatch");
            % ... and `build` checks it again, because a surrogate update may
            % change the state dimension after assembly passed.
            tc.verifyError(@() rep.build( ...
                wellFormedHistory([1; 2], "x0", [1 2 3]), MinimalSurrogate()), ...
                "rid:rep:stateDimensionMismatch");
        end

        function aSurrogateThatDeclaresNoStateCountIsRefused(tc)
            % A black-box model that cannot say how many states it has is not
            % one this space can be built from -- and NaN must not compare
            % false-ish into acceptance.
            rep = rid.rep.StateSpace(2, "pastSource", "simulated");
            tc.verifyError(@() rep.assertCompatible(MinimalSurrogate(false)), ...
                "rid:rep:stateDimensionMismatch");
        end

        function aCompatibleSurrogateIsAccepted(tc)
            % The complement, so the refusals above are not vacuously satisfied
            % by a check that rejects everything.
            rep = rid.rep.StateSpace(2, "pastSource", "simulated");
            rep.assertCompatible(MinimalSurrogate());
            tc.verifyTrue(true);
        end

        % --- the rest --------------------------------------------------------

        function describeNamesTheCoordinatesAndTheSource(tc)
            s = rid.rep.StateSpace(2, "pastSource", "measured").describe();
            tc.verifyClass(s, "string");
            tc.verifySubstring(s, "rid.rep.StateSpace");
            tc.verifySubstring(s, "x2_k");
            tc.verifySubstring(s, "measured");
        end

        function anInvalidPastSourceIsRejectedAtConstruction(tc)
            tc.verifyError(@() rid.rep.StateSpace(2, "pastSource", "online"), ...
                "rid:rep:invalidPastSource");
        end

    end

    methods (Static, Access = private)

        function hist = measuredHistory(u, x)
            %MEASUREDHISTORY  A well-formed history carrying measured states.
            %
            %   `wellFormedHistory` has no `.x`, and adding one there would put
            %   the field in eighty call sites that do not want it. Row k of x
            %   is the state measured at the instant u(k) was applied.
            hist   = wellFormedHistory(u, "x0", x(1, :));
            hist.x = x;
        end

    end

end
