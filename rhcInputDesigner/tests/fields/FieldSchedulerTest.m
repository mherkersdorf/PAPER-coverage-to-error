classdef FieldSchedulerTest < matlab.unittest.TestCase
%FIELDSCHEDULERTEST  The refinement-field update seam.
%
%   Each test names the input that would make it fail. Three of them exist because
%   the interesting failures here are all SILENT:
%
%     - a scheduler that fires by default would change v1 behaviour away from pure
%       space-filling, and nothing would raise;
%     - a scheduler that estimated `relevance`, or that dropped it, would replace a
%       declared field with ones and yield a worse experiment forever;
%     - an estimator that raises must not end a run, but the rejection must be
%       visible in the event log rather than absorbed.

    methods (Test)

        function theDefaultNeverFires(tc)
            % the MECHANISM exists from the start, the POLICY does
            % not. With updateInterval = Inf nothing is re-estimated and the run is
            % the pure space-filling design.
            %
            % Fails if: the default interval becomes finite, or the first poll is
            % treated as "no update yet, so update now" -- an off-by-one that would
            % fire once at step 1 on every run.
            sch = rid.fields.ImmediateScheduler();
            tc.verifyEqual(sch.updateInterval, Inf);

            [fieldsNew, info] = FieldSchedulerTest.pollWith(sch, 50);
            tc.verifyEmpty(fieldsNew);
            tc.verifyEqual(info.status, "none");
            tc.verifyFalse(info.updated);
            tc.verifyEqual(sch.nUpdates, 0);
        end

        function itFiresOnTheIntervalAndOverridesOnlyTheEstimatedMembers(tc)
            % The load-bearing test of the amended the refinement-field scheduler. `relevance` is declared a
            % priori and passed through untouched; `complexity` and `noiseLevel` are
            % replaced. The declared relevance here is NOT nominal, so a scheduler
            % that dropped it, held its own copy, or estimated it would return ones
            % and fail -- with nominal relevance every wrong implementation passes.
            %
            % Fails if: the scheduler holds a relevance of its own, or rebuilds the
            % struct from scratch instead of overriding two members of the one it
            % was given.
            sch = rid.fields.ImmediateScheduler( ...
                "updateInterval", 4, ...
                "complexity", ConstantFieldEstimator(3), ...
                "noiseLevel", ConstantFieldEstimator(7));

            declared            = rid.fields.nominalFields(5);
            declared.relevance  = (1:5).';       % deliberately NOT ones

            [none, info] = FieldSchedulerTest.pollWith(sch, 3, declared);
            tc.verifyEmpty(none, "3 steps have not reached an interval of 4");
            tc.verifyEqual(info.status, "none");

            [fieldsNew, info] = FieldSchedulerTest.pollWith(sch, 4, declared);
            tc.verifyEqual(info.status, "applied");
            tc.verifyTrue(info.updated);
            tc.verifyEqual(info.estimated, ["complexity", "noiseLevel"]);
            tc.verifyEqual(info.k0, 4);
            tc.verifyEqual(info.k1, info.k0, ...
                "estimation is blocking, so the snapshot and the swap coincide");

            tc.verifyEqual(fieldsNew.relevance,  declared.relevance);
            tc.verifyEqual(fieldsNew.complexity, 3 * ones(5, 1));
            tc.verifyEqual(fieldsNew.noiseLevel, 7 * ones(5, 1));
            tc.verifyEqual(sch.nUpdates, 1);
            tc.verifyEqual(sch.lastUpdateStep, 4);

            % The interval runs from the last ACCEPTED update, not from step zero.
            tc.verifyEmpty(FieldSchedulerTest.pollWith(sch, 7, declared));
            tc.verifyNotEmpty(FieldSchedulerTest.pollWith(sch, 8, declared));
        end

        function anIncompleteFieldSetIsRefused(tc)
            % An absent member defaults to NOMINAL in rid.crit.demandedResolution
            %, which is right there and wrong here: a scheduler handed a set
            % without `relevance` would return one without it, and a declared
            % relevance would silently become ones. It is the member the scheduler
            % must not estimate, so it is the member whose loss nothing downstream
            % can see.
            %
            % Fails if: the completeness check is dropped, or is applied AFTER the
            % interval test -- an incomplete set would then be accepted silently on
            % every step the scheduler does not fire, which is most of them.
            sch     = rid.fields.ImmediateScheduler("updateInterval", 1);
            partial = rmfield(rid.fields.nominalFields(5), "relevance");
            tc.verifyError(@() FieldSchedulerTest.pollWith(sch, 4, partial), ...
                "rid:fields:incompleteFields");

            % Also refused on a step that would NOT have fired.
            quiet = rid.fields.ImmediateScheduler();
            tc.verifyError(@() FieldSchedulerTest.pollWith(quiet, 4, partial), ...
                "rid:fields:incompleteFields");
        end

        function aRejectedEstimateKeepsTheFieldsAndIsLogged(tc)
            % The rejected-update handling's precedent: one bad estimate must not end a run, and the
            % refusal must be visible in the event log -- +fields may not write to
            % the console.
            %
            % Fails if: the estimator's error propagates and ends the run, or is
            % swallowed with status "none", which would be indistinguishable from
            % the policy simply not firing.
            sch = rid.fields.ImmediateScheduler("updateInterval", 1, ...
                "complexity", FailingFieldEstimator());

            [fieldsNew, info] = FieldSchedulerTest.pollWith(sch, 4);
            tc.verifyEmpty(fieldsNew);
            tc.verifyEqual(info.status, "rejected");
            tc.verifyFalse(info.updated);
            tc.verifySubstring(info.message, "rid:fields:estimatorRefused");
            tc.verifyEqual(sch.nRejected, 1);
            tc.verifyEqual(sch.nUpdates, 0);
            tc.verifyEqual(sch.lastUpdateStep, 0, ...
                "a refused attempt must not reset the interval clock");
        end

        function aMalformedHistoryIsRefusedBeforeAnyEstimatorSeesIt(tc)
            % structural validation first, so a malformed history never
            % reaches a semantic check that would report a confusing symptom.
            %
            % Fails if: assertHistory is called after the interval test, or not at
            % all -- the estimator would then receive a struct that is not a history.
            sch = rid.fields.ImmediateScheduler("updateInterval", 1);
            tc.verifyError(@() sch.poll(struct("u", "not a signal"), 1, 1, 1, ...
                zeros(5, 1), rid.fields.nominalFields(5)), ...
                "rid:missingHistoryField");
        end

        function theSchedulerIsAHandleAndDescribesItself(tc)
            % Schedulers are one of the two intended homes of mutable state
            % so the update counters must survive being passed around
            % by value; and the run record stores a human-readable summary.
            sch = rid.fields.ImmediateScheduler("updateInterval", 1);
            tc.verifyTrue(isa(sch, "rid.fields.Scheduler"));
            tc.verifyTrue(isa(sch, "handle"));

            FieldSchedulerTest.pollWith(sch, 4);
            copy = sch;
            FieldSchedulerTest.pollWith(copy, 5);
            tc.verifyEqual(sch.nUpdates, 2, "a value copy would report 1");

            tc.verifySubstring(sch.describe(), "ImmediateScheduler");
            tc.verifySubstring(sch.describe(), "Nominal");
        end

    end

    methods (Static, Access = private)

        function [fieldsNew, info] = pollWith(sch, nSteps, fields)
            %POLLWITH  Poll at a given history length; the surrogate is unused here.
            %
            %   rid.fields.Nominal ignores every argument by design, so the
            %   placeholders below are the point rather than a shortcut: what this
            %   test class exercises is the SEAM, and the estimators behind it are
            %   tested in NominalTest.
            if nargin < 3
                fields = rid.fields.nominalFields(5);
            end
            hist = wellFormedHistory(zeros(nSteps, 1));
            [fieldsNew, info] = sch.poll(hist, 1, 1, 1, zeros(5, 1), fields);
        end

    end

end
