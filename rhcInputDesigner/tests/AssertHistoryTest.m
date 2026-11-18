classdef AssertHistoryTest < matlab.unittest.TestCase
%ASSERTHISTORYTEST  The one structural validator of the signal history.
%
%   The history is a plain struct and gains no methods, so the single thing a
%   class would have bought -- one place to validate -- has to be earned by this
%   function being the one everyone calls, and being called FIRST.
%
%   Two properties are asserted here that no individual check states:
%
%     - the ORDER. A malformed history must never reach a semantic guard, or the
%       guard says something confidently wrong about a struct that is not a
%       history at all.
%     - the SPLIT. What this function refuses and what it lets through are
%       complementary: emptiness is not malformation, and each package keeps the
%       question only it can answer.

    methods (Test)

        % --- what a well-formed history is --------------------------------------

        function theMinimalHistoryOfAnOfflineRunPasses(tc)
            % Offline is the ABSENCE of y, not an empty y: an empty one is a
            % sentinel, and sentinels invite isempty checks that conflate "no
            % measurement" with "not yet".
            tc.verifyWarningFree(@() rid.assertHistory( ...
                struct("u", [0; 0.5; 1], "x0", 0, "tS", 0.1)));
        end

        function anOnlineHistoryWithMeasurementsPasses(tc)
            tc.verifyWarningFree(@() rid.assertHistory( ...
                struct("u", [0; 0.5], "y", [1; 2], "x0", [0 1], "tS", 0.01)));
        end

        function measuredStatesArePresentOnlyWhenTheyWereMeasured(tc)
            % `.x` is the optional measured-state field, and it
            % is optional in the same sense `.y` is: a representation that does
            % not read states never creates it, and an empty placeholder would be
            % the sentinel that `.y` is deliberately not.
            %
            % Fails if: a structural validator starts REQUIRING `.x`. Whether a
            % history must carry it is the representation's question, and
            % rid.rep.StateSpace answers it with rid:rep:missingStateMeasurement.
            tc.verifyWarningFree(@() rid.assertHistory( ...
                struct("u", [0; 0.5], "y", [1; 2], "x", [0 1; 0.3 0.7], ...
                       "x0", [0 1], "tS", 0.01)));
            tc.verifyWarningFree(@() rid.assertHistory( ...
                struct("u", [0; 0.5], "y", [1; 2], "x0", [0 1], "tS", 0.01)));
        end

        function aHistoryWithNoStepsIsWellFormed(tc)
            % EMPTINESS IS NOT MALFORMATION. Whether a caller can work with zero
            % rows is that caller's question -- rid.gen.Direct answers it with a
            % message about the actuator's starting position that a generic
            % validator could not write.
            %
            % Fails if: assertHistory rejects an empty u, which would make
            % rid:gen:missingInitialInput unreachable and delete its explanation.
            tc.verifyWarningFree(@() rid.assertHistory( ...
                struct("u", zeros(0, 1), "x0", 0, "tS", 0.1)));
        end

        % --- what it refuses -----------------------------------------------------

        function everyRequiredFieldIsNamedWhenItIsMissing(tc)
            good = struct("u", [0; 1], "x0", 0, "tS", 0.1);
            for name = ["u", "x0", "tS"]
                tc.verifyError(@() rid.assertHistory(rmfield(good, name)), ...
                    "rid:missingHistoryField", "field " + name);
            end
        end

        function theMissingFieldMessageSaysWhatToDoAboutIt(tc)
            % The guidance that rid:rep:missingInitialState gave for an absent x0
            % must not vanish just because the check moved: the xStart requirement requires the user
            % to be told about the canonical producer, not merely that a field is
            % absent.
            %
            % Fails if: absorbing a semantic guard into a structural one keeps the
            % rejection and drops the advice.
            good = struct("u", [0; 1], "x0", 0, "tS", 0.1);
            tc.verifySubstring(tc.messageOf(@() rid.assertHistory( ...
                rmfield(good, "x0"))), "steadyState");
            tc.verifySubstring(tc.messageOf(@() rid.assertHistory( ...
                rmfield(good, "u"))), "unrecorded, not undefined");
        end

        function aNonStructOrAStructArrayIsRefused(tc)
            tc.verifyError(@() rid.assertHistory(42), "rid:invalidHistory");
            tc.verifyError(@() rid.assertHistory("hist"), "rid:invalidHistory");
            arr = struct("u", {1, 2}, "x0", {0, 0}, "tS", {0.1, 0.1});
            tc.verifyError(@() rid.assertHistory(arr), "rid:invalidHistory");
        end

        function aNonFiniteSignalIsMalformationAndNotMissingData(tc)
            % An absent measurement is an absent FIELD, so a NaN in a signal is a
            % corrupted entry. Letting it through would put NaN into dRaw, where it
            % propagates into J and every distance silently.
            %
            % Fails if: only shape is checked. NaN survives every size test.
            good = struct("u", [0; 1], "y", [1; 2], "x", [0 0; 1 1], "x0", 0, "tS", 0.1);
            for name = ["u", "y", "x", "x0"]
                bad = good;
                bad.(name)(1) = NaN;
                tc.verifyError(@() rid.assertHistory(bad), "rid:invalidHistory", name);
                bad.(name)(1) = Inf;
                tc.verifyError(@() rid.assertHistory(bad), "rid:invalidHistory", name);
            end
        end

        function theSampleTimeMustBeAPositiveFiniteScalar(tc)
            % tS was REQUIRED PRESENT by rid.rep.Representation and never checked
            % for value, so a tS of 0 or -1 passed everywhere until now.
            %
            % Fails if: presence is mistaken for validity.
            good = struct("u", [0; 1], "x0", 0, "tS", 0.1);
            for bad = {0, -1, Inf, NaN, [0.1 0.2], "0.1"}
                broken    = good;
                broken.tS = bad{1};
                tc.verifyError(@() rid.assertHistory(broken), "rid:invalidHistory", ...
                    mat2str(string(bad{1})));
            end
        end

        function theInitialStateIsARowBecauseThatIsHowSimulateTakesIt(tc)
            % A column x0 transposes silently into something that reads as several
            % single-state samples, which srg.simulate would accept and misuse.
            %
            % Fails if: only the element count is checked.
            good = struct("u", [0; 1], "x0", [0 1], "tS", 0.1);
            tc.verifyWarningFree(@() rid.assertHistory(good));
            broken    = good;
            broken.x0 = [0; 1];
            tc.verifyError(@() rid.assertHistory(broken), "rid:invalidHistory");
        end

        function theTwoSignalsMustAgreeAboutHowManyStepsHappened(tc)
            % Structural, not representational: u and y are sampled together, and
            % no representation has a better answer than any other. This is the
            % check that moved out of rid.rep and rid.surrogate.
            %
            % Fails if: the row counts are compared only where a representation
            % happens to look at them.
            bad = struct("u", [0; 1; 2], "y", [1; 2], "x0", 0, "tS", 0.1);
            tc.verifyError(@() rid.assertHistory(bad), "rid:invalidHistory");
        end

        function theMeasuredStatesCarryOneRowPerAppliedInput(tc)
            % Row k of `hist.x` is the state measured at the instant u(k) was
            % applied, so u and x count the same steps, for the
            % structural reason y does. The width is NOT checked here: how many
            % states there should be is the representation's question, and
            % rid.rep.StateSpace answers it with rid:rep:dimensionMismatch.
            %
            % Fails if: `.x` is let through unvalidated. rid.rep.StateSpace would
            % then index a short x against u and pair u(k) with the state of a
            % different instant -- or, with a long one, place points that were
            % never realized.
            good = struct("u", [0; 1; 2], "x", [0 0; 1 1; 2 2], "x0", [0 0], "tS", 0.1);
            tc.verifyWarningFree(@() rid.assertHistory(good));
            for nRows = [2, 4]
                broken   = good;
                broken.x = zeros(nRows, 2);
                tc.verifyError(@() rid.assertHistory(broken), "rid:invalidHistory", ...
                    "x with " + nRows + " rows");
            end
        end

        function aMeasuredStateThatIsNotARealMatrixIsRefused(tc)
            % The same judgement assertRealMatrix already makes about u and y: a
            % string, a complex value or a 3-D array is malformation, whatever
            % the row count says.
            %
            % Fails if: `.x` is checked only for its row count.
            good = struct("u", [0; 1], "x", [0 0; 1 1], "x0", [0 0], "tS", 0.1);
            for bad = {"0", {0}, [1i 0; 0 0], zeros(2, 2, 2)}
                broken   = good;
                broken.x = bad{1};
                tc.verifyError(@() rid.assertHistory(broken), "rid:invalidHistory", ...
                    class(bad{1}));
            end
        end

        % --- the order, and the split --------------------------------------------

        function everyConsumerRejectsAMalformedHistoryBeforeItsOwnGuard(tc)
            % THE "CALLED FIRST" RULE, asserted on all four consumers at once. A
            % struct with no u reaches each of them; each must answer
            % "that is not a history" rather than its own, more specific and here
            % entirely misleading, complaint.
            %
            % Fails if: any consumer runs its semantic check first -- +gen would
            % then talk about rate limits, +rep about pastSource and +surrogate
            % about there being nothing to refit on, none of which is the problem.
            notAHistory = struct("y", 1);
            gen = rid.gen.Direct(rid.gen.Constraints(rid.rep.Narx(2, 1), ...
                      struct("u1_k", [-2 3], "du1_k", [-0.4 0.4], "y_km1", [-2 3])), 4);
            rep = rid.rep.Narx(2, 1, "inputDifference", false);
            srg = rid.surrogate.LinearFirstOrder(0.5, 0.1, 1);
            sch = rid.surrogate.ImmediateScheduler("updateInterval", 1);

            tc.verifyError(@() gen.feasibleSet(notAHistory), "rid:missingHistoryField");
            tc.verifyError(@() rep.build(notAHistory, srg), "rid:missingHistoryField");
            tc.verifyError(@() srg.update(notAHistory), "rid:missingHistoryField");
            tc.verifyError(@() sch.poll(notAHistory, srg), "rid:missingHistoryField");
        end

        function aWellFormedHistoryStillReachesEachPackagesOwnGuard(tc)
            % The other half of the split, and the one that would silently rot: if
            % assertHistory grew to answer the semantic questions too, these
            % identifiers would disappear and their messages with them.
            %
            % Fails if: a structural validator starts refusing things that are
            % merely not what THIS caller needs.
            gen = rid.gen.Direct(rid.gen.Constraints(rid.rep.Narx(2, 1), ...
                      struct("u1_k", [-2 3], "du1_k", [-0.4 0.4], "y_km1", [-2 3])), 4);
            srg = rid.surrogate.LinearFirstOrder(0.5, 0.1, 1);
            repMeasured = rid.rep.Narx(2, 1, "inputDifference", false);

            % +gen: well formed, but no last row to take u(k-1) from.
            tc.verifyError(@() gen.feasibleSet(wellFormedHistory(zeros(0, 1))), ...
                "rid:gen:missingInitialInput");

            % +rep: well formed and offline, but pastSource says "measured".
            tc.verifyError(@() repMeasured.build(wellFormedHistory((1:8).'), srg), ...
                "rid:rep:missingMeasurement");

            % +surrogate: well formed and offline, so there is nothing to refit on.
            tc.verifyError(@() srg.update(wellFormedHistory((1:8).')), ...
                "rid:surrogate:insufficientData");
        end

    end

    methods (Static, Access = private)

        function msg = messageOf(fcn)
            try
                fcn();
                msg = "";
            catch err
                msg = string(err.message);
            end
        end

    end

end
