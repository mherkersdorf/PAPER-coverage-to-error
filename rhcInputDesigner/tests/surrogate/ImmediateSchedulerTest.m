classdef ImmediateSchedulerTest < matlab.unittest.TestCase
%IMMEDIATESCHEDULERTEST  The seam: blocking updates behind a scheduler.
%
%   The scheduler is the one place besides rid.crit.CoverageState where mutable
%   state is intended, so it is a handle and the tests check that its counters
%   advance. What it must NOT do is decide anything the model owns: the fit is a
%   method of the model, and the scheduler only decides WHEN to ask and what to
%   record.
%
%   `info` is asserted field by field because AsyncScheduler fills the same
%   fields with something else, and the event log must not change shape when it
%   arrives.

    properties (Constant)
        TS = 0.1
    end

    methods (Access = private)
        function [hist, truth] = fittableHistory(tc, nSteps, seed)
            truth = rid.surrogate.LinearFirstOrder(0.7, tc.TS, 2.5, "offset", -0.3);
            rs    = RandStream("twister", "Seed", seed);
            u     = 2 * rand(rs, nSteps, 1) - 1;
            hist  = wellFormedHistory(u, ...
                "y", truth.simulate(u, truth.steadyState(u(1, :))), "tS", tc.TS);
        end
    end

    methods (Test)

        function schedulersAreHandleObjects(tc)
            % the project's coding standard, Architecture: mutable state lives in exactly two places.
            tc.verifyTrue(isa(rid.surrogate.ImmediateScheduler(), "handle"));
        end

        function theDefaultPolicyNeverUpdates(tc)
            % the update MECHANISM exists from the start, the POLICY does
            % not. A default that silently refits would presuppose one.
            tc.verifyEqual(rid.defaultOptions().surrogate.updateInterval, Inf);

            sch  = rid.surrogate.ImmediateScheduler();
            srg  = rid.surrogate.LinearFirstOrder(2, tc.TS, 3);
            hist = tc.fittableHistory(500, 31);

            [srgNew, info] = sch.poll(hist, srg);
            tc.verifyEmpty(srgNew);
            tc.verifyEqual(info.status, "none");
            tc.verifyFalse(info.updated);
            tc.verifyEqual(sch.nUpdates, 0);
        end

        function aFiniteIntervalUpdatesAndReportsAZeroModelAge(tc)
            % Blocking: the snapshot and the swap are the same step, so k0 == k1
            % and the model age is 0. The fields are reported anyway because
            % AsyncScheduler fills them differently.
            sch  = rid.surrogate.ImmediateScheduler("updateInterval", 100);
            srg  = rid.surrogate.LinearFirstOrder(2, tc.TS, 3);
            [hist, truth] = tc.fittableHistory(300, 32);

            [srgNew, info] = sch.poll(hist, srg);

            tc.verifyClass(srgNew, "rid.surrogate.LinearFirstOrder");
            tc.verifyEqual(srgNew.timeConstant, truth.timeConstant, "RelTol", 1e-9);
            tc.verifyEqual(info.status, "applied");
            tc.verifyTrue(info.updated);
            tc.verifyEqual(info.k0, 300);
            tc.verifyEqual(info.k1, 300);
            tc.verifyEqual(info.modelAge, 0);
            tc.verifyEqual(sch.nUpdates, 1);
            tc.verifyEqual(sch.lastUpdateStep, 300);
        end

        function theIntervalIsMeasuredFromTheLastAcceptedSwap(tc)
            sch = rid.surrogate.ImmediateScheduler("updateInterval", 100);
            srg = rid.surrogate.LinearFirstOrder(2, tc.TS, 3);

            [srgNew, ~] = sch.poll(tc.fittableHistory(120, 33), srg);
            tc.verifyNotEmpty(srgNew);                    % 120 - 0 >= 100

            [srgNew, info] = sch.poll(tc.fittableHistory(180, 33), srgNew);
            tc.verifyEmpty(srgNew);                       % 180 - 120 < 100
            tc.verifyEqual(info.status, "none");

            [srgNew, ~] = sch.poll(tc.fittableHistory(240, 33), srg);
            tc.verifyNotEmpty(srgNew);                    % 240 - 120 >= 100
            tc.verifyEqual(sch.nUpdates, 2);
        end

        function aModelThatCannotRefitIsNotPolledAndSaysWhy(tc)
            % "unsupported", not "none": the event log must record WHY the
            % surrogate never changes.
            sch  = rid.surrogate.ImmediateScheduler("updateInterval", 1);
            srg  = MinimalSurrogate();
            hist = wellFormedHistory(ones(50, 1), "y", ones(50, 1), "tS", tc.TS);

            [srgNew, info] = sch.poll(hist, srg);
            tc.verifyEmpty(srgNew);
            tc.verifyEqual(info.status, "unsupported");
            tc.verifySubstring(char(info.message), "canUpdate");
            tc.verifyEqual(sch.nUpdates, 0);
        end

        function aRejectedFitIsReportedAndTheRunSurvives(tc)
            % one stretch of history that happens not to excite the
            % parameters must not end a run. The previous surrogate is kept and the
            % rejection is recorded for the event log rather than printed.
            sch  = rid.surrogate.ImmediateScheduler("updateInterval", 1);
            srg  = rid.surrogate.LinearFirstOrder(2, tc.TS, 3);
            hist = wellFormedHistory(ones(200, 1), "tS", tc.TS, ...
                "y", srg.simulate(ones(200, 1), 0));      % constant input

            [srgNew, info] = sch.poll(hist, srg);

            tc.verifyEmpty(srgNew);
            tc.verifyEqual(info.status, "rejected");
            tc.verifySubstring(char(info.message), "rid:surrogate:insufficientData");
            tc.verifyEqual(sch.nRejected, 1);
            tc.verifyEqual(sch.nUpdates, 0);
            tc.verifyEqual(sch.lastUpdateStep, 0);        % a rejection is not a swap
        end

        function anUnstableFitIsRejectedThroughTheSameePath(tc)
            sch = rid.surrogate.ImmediateScheduler("updateInterval", 1);
            srg = rid.surrogate.LinearFirstOrder(2, tc.TS, 3);

            rs     = RandStream("twister", "Seed", 34);
            nSteps = 150;
            u = rand(rs, nSteps, 1);
            y = zeros(nSteps, 1);
            for kStep = 2:nSteps
                y(kStep) = 1.05 * y(kStep - 1) + 0.5 * u(kStep - 1);
            end
            hist = wellFormedHistory(u, "y", y, "tS", tc.TS);

            [srgNew, info] = sch.poll(hist, srg);
            tc.verifyEmpty(srgNew);
            tc.verifyEqual(info.status, "rejected");
            tc.verifySubstring(char(info.message), "rid:surrogate:unstableUpdate");
        end

        function everyInfoCarriesTheEventLogFields(tc)
            expected = ["status" "updated" "k0" "k1" "modelAge" "message"];
            sch  = rid.surrogate.ImmediateScheduler("updateInterval", 100);
            srg  = rid.surrogate.LinearFirstOrder(2, tc.TS, 3);
            [~, info] = sch.poll(tc.fittableHistory(300, 35), srg);
            tc.verifyEqual(sort(string(fieldnames(info)).'), sort(expected));
        end

    end

end
