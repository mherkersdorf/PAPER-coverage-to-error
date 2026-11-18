classdef ReplayTest < matlab.unittest.TestCase
%REPLAYTEST  Invariant 5: a run is replayable from (seed, opts, hist, eventLog).
%
%   The central case is deliberately a run CONTAINING A SURROGATE UPDATE. A replay
%   of a run without one exercises nothing: with the model fixed, every
%   step is an ordinary deterministic optimization, and the parts that make replay
%   hard -- the swap instant, the cache rebuild that follows it, the affine map
%   rebuilt across it -- never happen. A test suite that only ever replayed static
%   runs would report green while the interesting half was untested.
%
%   Bit-for-bit, never to a tolerance: a replay that agreed only approximately
%   would hide exactly the reorderings a canary exists to catch (the lesson,
%   one layer up).
%
%   A MEASURED-STATE RUN IS THE OTHER CASE THAT CANNOT BE FAKED.
%   Its horizons are anchored on `hist.x(end,:)`, so the record must carry the
%   states and the playback must hand them back: a replay that lost them would
%   design a different signal from step one, and there would be nothing left to
%   compare against.

    methods (Test)

        function anOfflineRunReplaysBitForBit(tc)
            % The baseline. Fails if anything in the loop reads state the record
            % does not carry -- the global RNG, the wall clock, a persistent
            % variable.
            rec        = ReplayTest.recordOf("nSteps", 5);
            [ok, diff] = rid.io.replay(rec);

            tc.verifyTrue(diff.replayable);
            tc.verifyTrue(ok, ReplayTest.explain(diff));
            tc.verifyEqual(diff.result.hist.u, rec.outputs.u);
            tc.verifyEqual(diff.result.status, rec.outputs.status);
        end

        function aRunContainingASurrogateUpdateReplaysBitForBit(tc)
            % THE CASE THAT MATTERS. The surrogate is refit part-way
            % through, which moves the model the candidates are simulated with, and
            % -- because pastSource is "measured" here and the representation is
            % built from raw signals -- is followed by the cache and affine-map
            % rebuilds the update triggers.
            %
            % Fails if: the update instant is re-derived from wall-clock rather
            % than reproduced by the policy, if any nondeterminism inside the fit
            % is not surfaced, or if the schedulers are reused instead of reset --
            % a carried interval clock would make the replay update at a different
            % step, which is precisely what the event-log comparison catches.
            rec = ReplayTest.recordOf("nSteps", 10, "online", true, ...
                "pastSource", "measured", "updateInterval", 4);

            kinds = string({rec.outputs.eventLog.kind});
            tc.assertTrue(any(kinds == "surrogateUpdate"), ...
                "fixture: the run must contain a surrogate update, or this " + ...
                "test is the static case again");
            % Filtered by kind: the log may also carry coverage resets, whose
            % detail is a different struct -- which is the point of a CLOSED kind
            % set that code branches on rather than a free-form message.
            swaps   = rec.outputs.eventLog(kinds == "surrogateUpdate");
            applied = arrayfun(@(e) e.detail.updated, swaps);
            tc.assertTrue(any(applied), ...
                "fixture: an update must have been ACCEPTED, not merely attempted");

            [ok, diff] = rid.io.replay(rec);
            tc.verifyTrue(ok, ReplayTest.explain(diff));
            tc.verifyEqual(diff.result.hist.u, rec.outputs.u);
            tc.verifyEqual(diff.result.hist.y, rec.outputs.y);

            % The swap happened at the same step, which is what invariant 5 is
            % about -- not merely that the histories agree.
            replayKinds = string({diff.result.eventLog.kind});
            tc.verifyEqual(replayKinds, kinds);
            tc.verifyEqual([diff.result.eventLog.step], ...
                           [rec.outputs.eventLog.step]);
        end

        function aMeasuredStateRunReplaysBitForBitIncludingItsStates(tc)
            % without the states in the record, a measured-state
            % run is not replayable AT ALL. The anchor each horizon starts from
            % is one surrogate step from `hist.x(end,:)`, so a replay that could
            % not see the measured states would design a different signal from
            % the first step -- and there would be nothing to compare it with.
            %
            % Fails if: `rec.outputs.x` is not stored; `hist0.x` is not seeded
            % from its first row; the playback is built without the states.
            rec = ReplayTest.measuredStateRecord();

            tc.assertTrue(isfield(rec.outputs, "x"), ...
                "the record must carry the measured states");
            tc.assertEqual(size(rec.outputs.x, 1), size(rec.outputs.u, 1), ...
                "fixture: one measured state per applied sample");

            [ok, diff] = rid.io.replay(rec);
            tc.verifyTrue(diff.replayable);
            tc.verifyTrue(ok, ReplayTest.explain(diff));
            tc.verifyEqual(diff.result.hist.u, rec.outputs.u);
            tc.verifyEqual(diff.result.hist.y, rec.outputs.y);
            tc.verifyEqual(diff.result.hist.x, rec.outputs.x);
        end

        function aPrimedOfflineRunReplaysBitForBit(tc)
            % A run may start from a history longer than the one pre-run row
            %: the NARX warm-up, or a priming signal. Those rows are the
            % initial condition, so the replay must hand ALL of them to the loop
            % -- the first size(u,1) - nSteps rows of the stored u -- and not
            % only row 1.
            %
            % Fails if: hist0 is rebuilt from u(1,:) alone; the history then
            % differs from the second row on.
            rec = ReplayTest.recordOf("nSteps", 4, "nPrime", 5);
            tc.assertEqual(size(rec.outputs.u, 1), 5 + 4, ...
                "fixture: the record must hold the priming rows and the designed ones");

            [ok, diff] = rid.io.replay(rec);
            tc.verifyTrue(ok, ReplayTest.explain(diff));
            tc.verifyEqual(diff.result.hist.u, rec.outputs.u);
        end

        function aPrimedOnlineRunPlaysBackFromItsFirstDesignedStep(tc)
            % The online case of the same: the priming rows carry measurements,
            % and playback must start AFTER them, at the first designed step.
            % Starting it at row 2 would feed the loop a priming measurement as
            % the response to its first designed input.
            %
            % Fails if: hist0.y is the first row only, or the playback is offset
            % by one row instead of by the priming length.
            rec = ReplayTest.recordOf("nSteps", 6, "nPrime", 5, "online", true, ...
                "pastSource", "measured", "updateInterval", 4);
            tc.assertEqual(size(rec.outputs.y, 1), 5 + 6, ...
                "fixture: one measurement per applied sample, priming included");

            [ok, diff] = rid.io.replay(rec);
            tc.verifyTrue(ok, ReplayTest.explain(diff));
            tc.verifyEqual(diff.result.hist.u, rec.outputs.u);
            tc.verifyEqual(diff.result.hist.y, rec.outputs.y);
        end

        function anOfflineRecordSaysItMeasuredNoStates(tc)
            % The counterpart, and the reason the field is an EMPTY matrix rather
            % than absent: a record is a serialized snapshot, so `isempty` there
            % means "no states were measured" exactly as it does for `y`. An
            % offline replay must then build a process-free run as before.
            %
            % Fails if: the record invents a state block for a run that had none.
            rec = ReplayTest.recordOf("nSteps", 3);
            tc.verifyEqual(rec.outputs.x, zeros(0, 0));

            [ok, diff] = rid.io.replay(rec);
            tc.verifyTrue(ok, ReplayTest.explain(diff));
            tc.verifyFalse(isfield(diff.result.hist, "x"));
        end

        function aPlaybackWithNoStatesRefusesToInventThem(tc)
            % A playback built from an offline record has no states to hand back,
            % and asking it for them is a contract error rather than a zero: an
            % invented state would anchor every horizon at the origin and the
            % divergence would surface many steps later as a puzzling `u`.
            %
            % Fails if: the second output is silently zeros(nApplied, 0) or [].
            playback = rid.io.PlaybackProcess([1; 2; 3]);
            tc.verifyEqual(playback.apply(0), 1, ...
                "the one-output form of the process contract is untouched");
            tc.verifyError(@() ReplayTest.bothOutputs(playback, 0), ...
                "rid:io:noStatePlayback");
        end

        function aPlaybackHandsBackItsStoredStatesInStep(tc)
            % The states run alongside the outputs, one row per applied sample,
            % and they run out together.
            playback = rid.io.PlaybackProcess([1; 2; 3], [10 11; 20 21; 30 31]);

            [y, x] = playback.apply([0; 0]);
            tc.verifyEqual(y, [1; 2]);
            tc.verifyEqual(x, [10 11; 20 21]);

            [y, x] = playback.apply(0);
            tc.verifyEqual(y, 3);
            tc.verifyEqual(x, [30 31]);
            tc.verifyError(@() playback.apply(0), "rid:io:playbackExhausted");
        end

        function aPlaybackWhoseStatesDoNotMatchItsOutputsIsRefused(tc)
            % The two are measurements of the SAME applied samples, so a state
            % block of a different length cannot be a record of the same run.
            % Refused at construction, where the mismatch is still legible.
            tc.verifyError(@() rid.io.PlaybackProcess([1; 2; 3], [10; 20]), ...
                "rid:io:dimensionMismatch");
        end

        function aRecordWithoutQIsReportedNonReplayableRatherThanRegenerated(tc)
            % Regeneration from the
            % seed is NOT a reproducibility mechanism: the generator belongs to a
            % toolbox and may change between releases, so a regenerated Q would
            % produce a mismatch that looks exactly like a code regression.
            %
            % Fails if: replay regenerates Q, or raises instead of reporting.
            rec = ReplayTest.recordOf("nSteps", 3, "storeQ", false);
            tc.assertFalse(rec.problem.storedQ);

            [ok, diff] = rid.io.replay(rec);
            tc.verifyFalse(ok);
            tc.verifyFalse(diff.replayable);
            tc.verifySubstring(diff.reason, "storeQ");
        end

        function aDivergenceIsNamedByStepAndQuantity(tc)
            % "The histories differ" is useless to someone deciding whether a
            % change was intended; "u diverged at step 4" says where to look.
            %
            % The divergence is manufactured by editing the STORED history, which
            % is the one thing a test may do to a record it owns -- it stands in
            % for the code having changed underneath an old record, which is the
            % situation replay exists to detect.
            rec = ReplayTest.recordOf("nSteps", 5);
            rec.outputs.u(4, :) = rec.outputs.u(4, :) + 0.05;

            [ok, diff] = rid.io.replay(rec);
            tc.verifyFalse(ok);
            tc.verifyTrue(diff.replayable);
            tc.verifyEqual(diff.quantity, "u");
            tc.verifyEqual(diff.step, 4);
            tc.verifyNotEmpty(diff.expected);
            tc.verifyNotEmpty(diff.actual);
        end

        function replayNeverTouchesAProcess(tc)
            % a replay must never have a process. An online run is
            % reproduced from the record's own measurements through a playback
            % source that runs out rather than inventing anything.
            %
            % Fails if: replay constructs a plant, or silently pads when the stored
            % measurements are exhausted -- a replay that had already diverged
            % would then continue and report a different, later divergence.
            rec = ReplayTest.recordOf("nSteps", 4, "online", true, ...
                                      "pastSource", "measured");
            [ok, diff] = rid.io.replay(rec);
            tc.verifyTrue(ok, ReplayTest.explain(diff));

            % Truncating the stored measurements must surface as an exhausted
            % playback, not as a quietly shorter run.
            short = rec;
            short.outputs.y = short.outputs.y(1:3, :);
            tc.verifyError(@() rid.io.replay(short), "rid:io:playbackExhausted");
        end

        function aRunStoppedByItsWallClockGuardReplaysToItsRecordedStep(tc)
            % The step a wall-clock guard stops at depends on the machine,
            % so a replay that kept `maxWallTime` stopped wherever its own clock
            % said: a real guard-stopped record (arm i, 30 steps, 60 s) diverged
            % at its last step. The fixture stops after 3 steps and then stores a
            % clock guard of 1e-9 s in place of the step guard -- a run stopped by
            % its clock, replayed on a machine infinitely slower. That makes the
            % defect deterministic: the old replay stopped before step 1.
            %
            % Fails if: the replay keeps `maxWallTime`, or does not stop at the
            % recorded step with the recorded status.
            rec = ReplayTest.recordOf("nSteps", 6, "maxSteps", 3);
            tc.assertEqual(rec.outputs.status, "guardTripped", ...
                "fixture: the run must be stopped by its guard");
            tc.assertEqual(rec.outputs.nSteps, 3);
            rec.configuration.effective.rhc.maxSteps    = Inf;
            rec.configuration.effective.rhc.maxWallTime = 1e-9;

            [ok, diff] = rid.io.replay(rec);
            tc.verifyTrue(ok, ReplayTest.explain(diff));
            tc.verifyEqual(diff.result.nSteps, 3);
            tc.verifyEqual(diff.result.status, "guardTripped");
        end

        function aBudgetRunIsReplayedWithoutItsWallClockGuard(tc)
            % The replay step-count rule, for a run that was NOT stopped by its guard: a finite
            % `maxWallTime` that the run never reached must not stop the replay
            % either -- a slower machine is not a different design.
            %
            % Fails if: the replay keeps `maxWallTime`.
            rec = ReplayTest.recordOf("nSteps", 4);
            tc.assertEqual(rec.outputs.status, "budgetReached");
            rec.configuration.effective.rhc.maxWallTime = 1e-9;

            [ok, diff] = rid.io.replay(rec);
            tc.verifyTrue(ok, ReplayTest.explain(diff));
            tc.verifyEqual(diff.result.status, "budgetReached");
        end

        function aRecordThatWentThroughDiskReplaysToo(tc)
            % The record replay actually gets is usually one that was loaded, and
            % load reconstructs every handle in it. If replay depended on sharing a
            % handle with the original run it would pass in-memory and fail here.
            root = string(tempname());
            mkdir(root);
            tc.addTeardown(@() rmdir(root, "s"));

            rec  = ReplayTest.recordOf("nSteps", 4);
            path = rid.io.saveRun(rec, "", "root", root);
            back = rid.io.loadRun(path);

            [ok, diff] = rid.io.replay(back);
            tc.verifyTrue(ok, ReplayTest.explain(diff));
        end

        function aMeasuredStateRecordSurvivesDiskAndReplays(tc)
            % THE ROUND TRIP AN EXPERIMENT ACTUALLY MAKES. The
            % in-memory test above replays the record the run returned; a record
            % a figure cites months later is one that went out through `saveRun`
            % and came back through `loadRun`. `run.mat` is where `rec.outputs.x`
            % has to survive, and the measured states are the one field a
            % measured-state record has that no other run kind does -- so a writer
            % or a reader that quietly dropped them would stay green on every
            % record in this suite, and the loss would surface as a run that
            % simply could not be replayed.
            %
            % Fails if: the states do not survive the write or the read, or the
            % replay depends on sharing a handle with the original run.
            root = string(tempname());
            mkdir(root);
            tc.addTeardown(@() rmdir(root, "s"));

            rec  = ReplayTest.measuredStateRecord();
            tc.assertNotEmpty(rec.outputs.x, ...
                "fixture: this record must carry measured states");

            path = rid.io.saveRun(rec, "", "root", root);
            back = rid.io.loadRun(path);

            tc.assertTrue(isfield(back.outputs, "x"), ...
                "the loaded record lost its measured states");
            tc.verifyEqual(back.outputs.x, rec.outputs.x, ...
                "the states must come back bit-identical");

            [ok, diff] = rid.io.replay(back);
            tc.verifyTrue(diff.replayable);
            tc.verifyTrue(ok, ReplayTest.explain(diff));
            tc.verifyEqual(diff.result.hist.x, rec.outputs.x);
        end

    end

    methods (Static)

        function rec = recordOf(varargin)
            %RECORDOF  A record from a short run, built by running the loop.
            args = struct(varargin{:});
            opt  = @(name, fallback) ReplayTest.pick(args, name, fallback);

            pastSource = opt("pastSource", "simulated");
            f = SetupTest.fixture("nSup", 10, "epsReq0", 0.6, ...
                                  "pastSource", pastSource);
            if opt("updateInterval", Inf) < Inf
                % A surrogate that CAN refit, on a schedule that fires inside the
                % run: rid.surrogate.LinearFirstOrder.update needs measurements, so
                % this pairs with an online run.
                f.sch = rid.surrogate.ImmediateScheduler( ...
                    "updateInterval", opt("updateInterval", Inf), ...
                    "minSamples", 3);
            end

            ov = RunTest.overrides("nSteps", opt("nSteps", 5));
            if isfield(args, "maxSteps")
                ov.rhc.maxSteps = args.maxSteps;
            end
            ov.crit.epsReq0 = 0.6;
            ov.io.storeQ    = opt("storeQ", true);

            process = [];
            online  = opt("online", false);
            if online
                process = FirstOrderProcess(0.8, 0.1, 1.4);
            end

            hist0 = RunTest.history("srg", f.srg, "online", online);
            if opt("nPrime", 1) > 1
                % A priming history: the pre-run row, then data the run starts
                % with; online, with a measurement for every row.
                hist0 = RunTest.primingHistory(f.srg);
                hist0.u = hist0.u(1:opt("nPrime", 1));
                if online
                    hist0.y = 0.5 * hist0.u + 0.05;
                end
            end
            r     = rid.rhc.run(process, SetupTest.setupFrom(f), hist0, ov);
            rec = r.record;
        end

        function rec = measuredStateRecord(opts)
            %MEASUREDSTATERECORD  A record of an ONLINE run whose states are measured.
            %
            %   The joint input--state space, `pastSource = "measured"`,
            %   against a plant that shares the surrogate's state basis but not
            %   its dynamics -- so the anchor of part 6 actually moves the design
            %   and a replay that lost the states could not reproduce it.
            arguments
                opts.nSteps  (1,1) double = 4
                opts.nStates (1,1) double = 2
            end
            f  = SetupTest.stateSpaceFixture(opts.nStates, ...
                     "pastSource", "measured", "epsReq0", RunTest.StateEpsReq0);
            ov = RunTest.stateOverrides(opts.nSteps);
            ov.io.storeQ = true;

            r = rid.rhc.run(RunTest.statePlant(opts.nStates), ...
                    SetupTest.setupFrom(f), ...
                    RunTest.measuredHistory(f.srg, zeros(1, opts.nStates)), ov);
            rec = r.record;
        end

        function [y, x] = bothOutputs(process, uApplied)
            %BOTHOUTPUTS  Ask for two outputs, so that a test can verify the refusal.
            [y, x] = process.apply(uApplied);
        end

        function value = pick(args, name, fallback)
            value = fallback;
            if isfield(args, name)
                value = args.(name);
            end
        end

        function text = explain(diff)
            %EXPLAIN  The divergence, for a failure message worth reading.
            text = "replay diverged: " + diff.reason;
        end

    end

end
