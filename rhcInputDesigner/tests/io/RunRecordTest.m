classdef RunRecordTest < matlab.unittest.TestCase
%RUNRECORDTEST  Run records: what they contain, and that they survive the disk.
%
%   A plot made six months ago is reproducible if and only if its record is
%   complete. The tests below check the two properties that
%   make that true and that would otherwise fail silently:
%
%     - the record stores the EFFECTIVE options, so a later change of default
%       cannot alter what an old run means. Asserted by loading a record whose
%       stored value DIFFERS from the current default, which a loader that merged
%       defaults would quietly overwrite;
%     - records are IMMUTABLE, so a second run never overwrites a result a figure
%       already cites.
%
%   Records are written into a temporary folder, never into the project's `runs/`:
%   a test suite that leaves experiment directories behind is indistinguishable
%   from an experiment.

    properties
        Root (1,1) string
    end

    methods (TestMethodSetup)
        function makeRoot(tc)
            tc.Root = string(tempname());
            mkdir(tc.Root);
            tc.addTeardown(@() rmdir(tc.Root, "s"));
        end
    end

    methods (Test)

        function saveThenLoadReturnsAnIdenticalRecord(tc)
            % Round trip of the record.
            %
            % Fails if: a field is dropped on the way to disk, or the loader
            % reconstructs anything rather than reading what was stored.
            rec  = RunRecordTest.record();
            path = rid.io.saveRun(rec, "a note", "root", tc.Root);
            tc.verifyTrue(isfolder(path));
            tc.verifyTrue(isfile(fullfile(path, "run.mat")));
            tc.verifyTrue(isfile(fullfile(path, "run.json")));

            back = rid.io.loadRun(path);
            expected      = rec;
            expected.note = "a note";
            tc.verifyEqual(back.identity.runId, expected.identity.runId);
            tc.verifyEqual(back.note, expected.note);
            tc.verifyEqual(back.outputs.u, expected.outputs.u);
            tc.verifyEqual(back.outputs.status, expected.outputs.status);
            tc.verifyEqual(back.configuration.effective, expected.configuration.effective);
            tc.verifyEqual(back.setup.epsReqInitial, expected.setup.epsReqInitial);
        end

        function theEffectiveOptionsSurviveAChangeOfDefault(tc)
            % THE LOAD-BEARING ONE (section 1). A record storing only "defaults
            % plus epsReq0 = 0.05" changes meaning the day a default changes, and
            % defaults will change while the psiMax = 1 default and the betaNoise default are open.
            %
            % The stored psiMax below differs from the current default, so a loader
            % that merged `rid.defaultOptions()` under the overrides would return
            % the default and this would fail. Simulating a changed default by
            % storing a non-default value is the only way to test this without
            % editing rid.defaultOptions from a test.
            rec = RunRecordTest.record();
            tc.assertNotEqual(rec.configuration.effective.crit.psiMax, ...
                rid.defaultOptions().crit.psiMax, ...
                "fixture: the stored value must differ from today's default");

            path = rid.io.saveRun(rec, "", "root", tc.Root);
            back = rid.io.loadRun(path);
            tc.verifyEqual(back.configuration.effective.crit.psiMax, ...
                rec.configuration.effective.crit.psiMax);
            % The overrides are kept ALONGSIDE, because they preserve the intent.
            tc.verifyEqual(back.configuration.overrides, ...
                rec.configuration.overrides);
        end

        function epsReqInitialIsTheStartingResolutionAndSetupNoLongerHoldsIt(tc)
            % THE REASON THE SNAPSHOT EXISTS (the single-setup rule).
            % `setup.state` is a handle the loop ADVANCED, so at the end of a run
            % `setup` no longer holds the configuration the run started from.
            %
            % A coverage reset alone does NOT demonstrate this: a reset clears the
            % cache and changes no `epsReq` at all -- which is exactly why it is
            % its own event kind and not a mislabelled "resolutionChange". So the
            % run below carries BOTH, and each half is asserted against the thing
            % it actually moved: a field update moves `epsReq`, and the reset moves
            % the cache.
            %
            % Fails if: the record stores only `setup`, or stores `epsReqInitial`
            % by reading `setup.state.epsReq` at the end rather than at the start.
            [rec, epsReqAtStart] = RunRecordTest.resetAndRefineRecord();

            kinds = string({rec.outputs.eventLog.kind});
            tc.assertTrue(any(kinds == "coverageReset"), ...
                "fixture: a reset must have fired");
            tc.assertTrue(any(kinds == "resolutionChange"), ...
                "fixture: the resolution must have changed, or epsReqInitial " + ...
                "and the final epsReq cannot differ at all");

            finalState = rec.setup.setup.state;
            tc.verifyEqual(rec.setup.epsReqInitial, epsReqAtStart);
            tc.verifyNotEqual(rec.setup.epsReqInitial, finalState.epsReq, ...
                "the stored snapshot must not be the resolution the run ended on");

            % And the snapshot is what a replay rebuilds from, so it must still
                % describe a usable state.
            tc.verifyEqual(size(rec.setup.epsReqInitial), size(finalState.epsReq));
            tc.verifyTrue(all(rec.setup.epsReqInitial > 0));
        end

        function consecutiveResetsAreReportedAsAWarningWithoutANewStatus(tc)
            % The coverage reset's degenerate case: if a single horizon covers the whole ROI the
            % cache is empty again at every step, the problem posed is identical
            % each time, and the input freezes. That is a statement about the
            % CONFIGURATION and the user needs it.
            %
            % It travels in the event detail and in `result.warnings`, never
            % through `warning` -- +rhc may not write to the console,
            % and diagnostics travel in returned structs. And NO fifth termination
            % status: the four stay, and this run ends on its budget like any other.
            %
            % Fails if: the first reset is already called a warning (one reset is
            % ordinary), or a status is added, or the freeze goes unreported.
            r = RunRecordTest.frozenRun();

            tc.verifyEqual(r.status, "budgetReached", ...
                "the four statuses stay; a degenerate configuration is not a fifth");
            resets = r.eventLog(string({r.eventLog.kind}) == "coverageReset");
            tc.assertGreaterThan(numel(resets), 2, ...
                "fixture: the run must reset on consecutive steps");

            tc.verifyEqual(resets(1).detail.consecutive, 1);
            tc.verifyEqual(resets(1).detail.severity, "info", ...
                "one reset is the ordinary the coverage reset case, not a warning");
            tc.verifyEqual(resets(2).detail.severity, "warning");
            tc.verifyNotEmpty(r.warnings);
            tc.verifySubstring(r.warnings(1), "consecutive coverage");
        end

        function anIsolatedResetIsNotAWarningAndAConsecutiveOneIs(tc)
            % THE OTHER HALF, and it lives inside ONE run rather than in a second
            % fixture. An isolated reset is the intended coverage-reset behaviour -- the
            % cache empties, the design re-covers the ROI over several steps, a
            % second sweep -- and must stay silent; only a reset immediately
            % following another says the configuration is degenerate.
            %
            % Testing both branches in one run also removes a fragility: a run
            % configured to reset in isolation can drift into the consecutive
            % regime as coverage gets easier, which is exactly what this fixture
            % does after step 10. What is asserted is the DISCRIMINATION, not a
            % property of the whole run.
            %
            % Fails if: every reset is a warning, or none is.
            r      = RunTest.offlineRun("nSup", 10, "epsReq0", 0.75, "nSteps", 12);
            resets = r.eventLog(string({r.eventLog.kind}) == "coverageReset");
            steps  = [resets.step];
            tc.assertGreaterThan(numel(resets), 1, "fixture: several resets");

            isFirstOfARun = [true, diff(steps) > 1];
            tc.assertTrue(any(isFirstOfARun), "fixture: an isolated reset");
            tc.assertTrue(any(~isFirstOfARun), "fixture: a consecutive reset");

            for iReset = 1:numel(resets)
                expected = "warning";
                if isFirstOfARun(iReset)
                    expected = "info";
                end
                tc.verifyEqual(resets(iReset).detail.severity, expected, ...
                    sprintf("reset at step %d", steps(iReset)));
            end
            tc.verifyEqual(numel(r.warnings), sum(~isFirstOfARun), ...
                "one warning per consecutive reset, and none for an isolated one");
        end

        function aRecordIsNeverOverwritten(tc)
            % Section 2 rule 2. Re-running produces a new record; a collision is an
            % error rather than a silent replacement of a result already cited.
            rec = RunRecordTest.record();
            rid.io.saveRun(rec, "", "root", tc.Root);
            tc.verifyError(@() rid.io.saveRun(rec, "", "root", tc.Root), ...
                "rid:io:recordExists");
        end

        function aDirtyTreeEmbedsItsDiffAndRequireCleanTreeRefusesIt(tc)
            % Section 2 rule 4. The diff is EMBEDDED so the record stays
            % self-contained, and written out as dirty.patch so it can be applied
            % without MATLAB.
            rec              = RunRecordTest.record();
            rec.code.isDirty = true;
            rec.code.diff    = "diff --git a/x b/x" + newline + "+changed";
            path             = rid.io.saveRun(rec, "", "root", tc.Root);
            tc.verifyTrue(isfile(fullfile(path, "dirty.patch")));

            % And the refusal happens BEFORE a run, not after it: the
            % pre-flight capture is what makes that possible.
            f  = SetupTest.fixture("epsReq0", RunTest.EpsReq0);
            ov = RunTest.overrides("nSteps", 2);
            ov.io.requireCleanTree = true;
            if rid.io.gitState(RunRecordTest.projectRoot()).isDirty
                tc.verifyError(@() rid.rhc.run([], SetupTest.setupFrom(f), ...
                    RunTest.history("srg", f.srg), ov), ...
                    "rid:io:dirtyWorkingTree");
            else
                % A clean tree must NOT be refused -- the other half.
                r = rid.rhc.run([], SetupTest.setupFrom(f), ...
                    RunTest.history("srg", f.srg), ov);
                tc.verifyEqual(r.nSteps, 2);
            end
        end

        function anAbsoluteRunsDirRoundTripsThroughWriteAndRead(tc)
            % THE BUG THIS COMMIT CLOSES. The resolution of opts.io.runsDir lived
            % in three places: rid.rhc.run handled an absolute path, saveRun and
            % loadRun each joined the project root on unconditionally. With an
            % absolute runsDir the run wrote to the right directory while loadRun
            % looked under <projectRoot>/<abs/path> -- WRITE AND READ DISAGREED,
            % and the symptom was a record reported missing moments after it was
            % written, which reads as a missing record rather than as a path bug.
            %
            % Fails if: any of the three recomputes the root for itself again, or
            % if fullfile is used where an absolute-path test is needed -- joining
            % a root onto an absolute path does not error, it silently produces a
            % third directory that is neither.
            root = string(tempname());
            tc.addTeardown(@() RunTest.removeIfPresent(root));

            f  = SetupTest.fixture("nSup", 8, "epsReq0", RunTest.EpsReq0);
            ov = RunTest.overrides("nSteps", 2);
            ov.io.saveRun = true;
            ov.io.runsDir = root;      % ABSOLUTE

            r = rid.rhc.run([], SetupTest.setupFrom(f), ...
                            RunTest.history("srg", f.srg), ov);

            % 1. the run wrote where the option said,
            tc.verifyTrue(startsWith(r.recordPath, root), ...
                sprintf("wrote to %s, not under the overridden runsDir", ...
                        r.recordPath));
            tc.verifyTrue(isfile(fullfile(r.recordPath, "run.mat")));

            % 2. the record NAMES the directory it actually wrote to,
            tc.verifyEqual(r.record.configuration.effective.io.runsDir, root);

            % 3. and the reader finds it -- by runId alone, through the same
            %    resolution, which is the half that used to fail.
            back = rid.io.loadRun(r.record.identity.runId, "root", root);
            tc.verifyEqual(back.identity.runId, r.record.identity.runId);
            tc.verifyEqual(rid.io.runsRoot(root), root, ...
                "an absolute runsDir resolves to itself");
        end

        function oneResolutionServesTheWriterTheReaderAndTheLoop(tc)
            % Not "three functions agree today" but "there is one function": a
            % relative runsDir resolves against the project root, an absolute one
            % to itself, and rid.io.projectRoot is the single definition of where
            % the project is -- it was three copies, one per file.
            tc.verifyEqual(rid.io.runsRoot("runs"), ...
                string(fullfile(rid.io.projectRoot(), "runs")));
            tc.verifyEqual(rid.io.runsRoot(), ...
                rid.io.runsRoot(rid.defaultOptions().io.runsDir), ...
                "the omitted argument must mean the documented default");
            tc.verifyTrue(isfolder(fullfile(rid.io.projectRoot(), "src")), ...
                "the project root must be the folder holding src/");
        end

        function theProjectRootChecksItsOwnInvariantRatherThanAssumingIt(tc)
            % The derivation is a fixed number of `fileparts` calls, so it
            % is silently wrong the moment the file moves -- one package deeper and
            % it returns `src/`, which EXISTS, is readable, and is not the project
            % root. Nothing downstream would notice: runsRoot would resolve
            % `src/runs`, saveRun would create it, and the failure would surface
            % much later as a record reported missing right after it was written.
            %
            % The relocation is genuine rather than simulated: the same file, byte
            % for byte, is placed one package deeper in a differently-named package
            % (so it shadows nothing), over a tree where the CORRECT root would
            % pass. Its derivation then lands on `<temp>/src`, which has neither
            % marker, and the error must name what was computed and what was
            % missing.
            %
            % Fails if: the check is dropped, or it tests only that the path is
            % non-empty, or it uses one marker -- a single common name can be
            % matched by an unrelated ancestor.
            temp = string(tempname());
            tc.addTeardown(@() RunTest.removeIfPresent(temp));

            % A tree where the correct root is `temp`: it has both markers.
            mkdir(fullfile(temp, "src"));
            mkdir(fullfile(temp, "tests"));

            % The same file, one package too deep, in a package named differently
            % so that nothing on the real path is shadowed.
            deep = fullfile(temp, "src", "+ridProbe", "+io", "+extra");
            mkdir(deep);
            copyfile(which("rid.io.projectRoot"), fullfile(deep, "projectRoot.m"));
            addpath(fullfile(temp, "src"));
            tc.addTeardown(@() rmpath(fullfile(temp, "src")));

            % Its derivation reaches <temp>/src -- a real folder, and the wrong one.
            % Caught rather than only asserted, because the MESSAGE is half the
            % point: an error saying merely "invalid project root" would leave the
            % reader no better off than the missing record would have.
            raised = false;
            try
                ridProbe.io.extra.projectRoot();
            catch err
                raised = true;
                tc.verifyEqual(string(err.identifier), "rid:io:invalidProjectRoot");
                tc.verifySubstring(err.message, fullfile(temp, "src"), ...
                    "the error must name the path it computed");
                tc.verifySubstring(err.message, "src/", ...
                    "and what was missing from it");
            end
            tc.verifyTrue(raised, ...
                "the relocated copy derived a rootless path and did not object");

            % THE OTHER HALF: the real one passes, and its answer is the folder
            % that holds src/ and tests/. A check that always raised would pass the
            % assertion above on its own.
            root = rid.io.projectRoot();
            tc.verifyTrue(isfolder(fullfile(root, "src")));
            tc.verifyTrue(isfolder(fullfile(root, "tests")));
        end

        function theAffineProbeSettingsReachTheRecordAndTheProbe(tc)
            % Both were literals in rid.gen.affineMap's arguments block, so a
            % record's effective options could not say what the probe did.
            %
            % They are not the same kind of number, and the two halves below are
            % the difference. `affineProbeStep` is CONDITIONING: the probe
            % evaluates the linear part rather than approximating it, so h cannot
            % change the map or the verdict. `affineCheckTol` decides the VERDICT:
            % tightening it past what any real composition can meet must drop the
            % run off the analytic path.
            %
            % THE CONDITIONING HALF IS ASSERTED TO A ROUNDING TOLERANCE, NOT
            % BIT-FOR-BIT, and the distinction is real rather than a hedge. "Exact
            % for any h" is exact in exact arithmetic; in floating point the
            % recovered A differs in its last ulps, and a multi-start gradient
            % search amplifies that into a designed signal differing at ~1e-13
            % relative. What must not change is the MAP and the PATH, and 1e-9 is
            % far below any difference a real change of behaviour would produce.
            %
            % Fails if: run keeps calling affineMap without passing them, in which
            % case the override is recorded and ignored -- and the second half
            % would still report "analytic".
            baseline = RunTest.offlineRun("nSteps", 3);
            tc.verifyEqual(baseline.record.configuration.effective.gen.affineProbeStep, ...
                1e-3);
            tc.verifyEqual(baseline.record.configuration.effective.gen.affineCheckTol, ...
                1e-9);
            tc.assertEqual(baseline.steps(1).gradientPath, "analytic");

            f  = SetupTest.fixture("nSup", 12, "epsReq0", RunTest.EpsReq0);
            ov = RunTest.overrides("nSteps", 3);
            ov.gen = struct("affineProbeStep", 1e-1);
            conditioned = rid.rhc.run([], SetupTest.setupFrom(f), ...
                                      RunTest.history("srg", f.srg), ov);
            tc.verifyEqual(conditioned.hist.u, baseline.hist.u, "RelTol", 1e-9, ...
                "h is conditioning only: a hundredfold change of probe step must " + ...
                "not move the design beyond rounding");
            tc.verifyEqual(conditioned.steps(1).gradientPath, "analytic", ...
                "and must not change which gradient path the run takes");
            tc.verifySubstring(conditioned.steps(1).affineTier, "affine");

            f2  = SetupTest.fixture("nSup", 12, "epsReq0", RunTest.EpsReq0);
            ov2 = RunTest.overrides("nSteps", 2);
            ov2.gen = struct("affineCheckTol", 1e-30);
            strict = rid.rhc.run([], SetupTest.setupFrom(f2), ...
                                 RunTest.history("srg", f2.srg), ov2);
            tc.verifyEqual(strict.steps(1).gradientPath, "finiteDifference", ...
                "checkTol decides the verdict, so it must reach the probe");
            tc.verifySubstring(strict.steps(1).affineTier, "notAffine");
        end

        % --- the record reports what was USED ------------------------

        function theOptimizerSeedIsTheOptimizersOwnAndNotTheGenerators(tc)
            % The used-values record rule. `rec.seeds.optimizer` was assigned `opts.gen.seed` -- the
            % generator's seed under a second name -- so the record asserted an
            % equality that is not a fact and stayed plausible for exactly as long
            % as the two happened to coincide. Every fixture in this suite makes
            % them coincide, which is why nothing caught it.
            %
            % Fails if: the field is copied from an option rather than read from
            % the object that acted.
            f   = SetupTest.fixture("nSup", 8, "epsReq0", RunTest.EpsReq0);
            f.sol = rid.opt.MultiStart("seed", 7);
            tc.assertNotEqual(f.sol.seed, rid.defaultOptions().gen.seed, ...
                "fixture: the two seeds must DIFFER, or the copy passes");

            r = rid.rhc.run([], SetupTest.setupFrom(f), ...
                            RunTest.history("srg", f.srg), ...
                            RunTest.overrides("nSteps", 2));

            tc.verifyEqual(r.record.seeds.optimizer, f.sol.seed, ...
                "the record must name the seed the optimizer actually used");
            tc.verifyNotEqual(r.record.seeds.optimizer, r.record.seeds.generator, ...
                "two fields under two names must not be one value");
        end

        function anOptimizerWithNoSeedRecordsAnExplicitAbsence(tc)
            % rid.opt.GridSearch is exhaustive and has no starts at all, so it has
            % no seed (the affineness-check sampling records that it reports nStarts = 0). NaN, never 0:
            % zero would be a seed the optimizer never used, which is the kind of
            % plausible-but-false entry the used-values record rule exists to remove.
            f     = SetupTest.fixture("nSup", 6, "epsReq0", 1.0, "nHorizon", 2);
            f.sol = rid.opt.GridSearch("nLevels", 3);
            r     = rid.rhc.run([], SetupTest.setupFrom(f), ...
                                RunTest.history("srg", f.srg), ...
                                RunTest.overrides("nSteps", 2, "epsReq0", 1.0));

            tc.verifyTrue(isnan(r.record.seeds.optimizer), ...
                "an optimizer with no seed must record an absence, not a number");
        end

        function aSuppliedQRecordsItsOriginAndNoInventedSeed(tc)
            % Until now `rid.roi.supportingPoints` was never called from
            % `src/`, so a Q built by hand -- as this one is -- had its provenance
            % discarded by every caller: the record asserted `opts.roi.seed` as
            % the seed of a generator that never ran. qInfo is the fix.
            %
            % Fails if: qOrigin is inferred from opts.roi.* rather than read from
            % qInfo, or rec.seeds.supportingPoints falls back to opts.roi.seed.
            f = SetupTest.fixture("nSup", 8, "epsReq0", RunTest.EpsReq0);
            qInfo = struct("origin", "supplied", "method", "", "seed", NaN, ...
                           "hash", "");
            setup = rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, f.fields, ...
                                  f.gen, f.sol, f.state, f.sch, f.fieldSch, qInfo);

            ov = RunTest.overrides("nSteps", 2);
            % A DIFFERENT seed than any generator used for this Q -- the fixture's
            % own f.Q was seeded with opts.seed = 0 inside SetupTest.fixture, and
            % no generator ran here at all. If the record ever fell back to this
            % option it would be caught here rather than passing by coincidence.
            ov.roi = struct("seed", 42);

            r = rid.rhc.run([], setup, RunTest.history("srg", f.srg), ov);

            tc.verifyEqual(r.record.problem.qOrigin, "supplied");
            tc.verifyTrue(isnan(r.record.seeds.supportingPoints), ...
                "no generator produced this Q, so no seed may be recorded");
            tc.verifyNotEqual(r.record.seeds.supportingPoints, 42, ...
                "opts.roi.seed must not leak into the record as if it had " + ...
                "generated Q");
        end

        function theRecordedUpdateIntervalIsTheSchedulersNotTheOptions(tc)
            % The schedulers read `updateInterval` from rid.defaultOptions at
            % CONSTRUCTION, so an already-built scheduler keeps its own value while
            % a run's override sits in the effective options unread: the record
            % could name an interval the run did not use.
            %
            % Fails if: the interval is copied from opts rather than read from the
            % scheduler objects that polled.
            f     = SetupTest.fixture("nSup", 8, "epsReq0", RunTest.EpsReq0);
            sch   = rid.surrogate.ImmediateScheduler("updateInterval", 3);
            fldSch = rid.fields.ImmediateScheduler("updateInterval", 5);
            setup = rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, f.fields, ...
                                  f.gen, f.sol, f.state, sch, fldSch, f.qInfo);

            ov = RunTest.overrides("nSteps", 2);
            % The option says something ELSE, which is the whole point: before
            % the used-values record rule this is what the record would have carried.
            ov.surrogate = struct("updateInterval", 99);
            ov.fields    = struct("updateInterval", 77);

            r = rid.rhc.run([], setup, RunTest.history("srg", f.srg), ov);

            tc.verifyEqual(r.record.problem.updateInterval.surrogate, 3);
            tc.verifyEqual(r.record.problem.updateInterval.fields, 5);
        end

        function theMeasuredStatesAreInTheOutputsAndSurviveTheDisk(tc)
            % Section 1, the `outputs` row: `x` when the run measured states
            %. The history is the source of truth WITHIN a run; across
            % runs the record is, so a measured-state run that stored no states
            % could never be replayed or re-plotted -- and the omission would be
            % invisible, because every other quantity would still be there.
            %
            % Fails if: `x` is not written, or is written as a copy of `y`.
            rec = ReplayTest.measuredStateRecord("nSteps", 3);

            tc.verifyEqual(size(rec.outputs.x), ...
                           [size(rec.outputs.u, 1), 2]);
            tc.verifyTrue(all(isfinite(rec.outputs.x), "all"));

            path = rid.io.saveRun(rec, "", "root", tc.Root);
            back = rid.io.loadRun(path);
            tc.verifyEqual(back.outputs.x, rec.outputs.x);
        end

        function aMissingRecordIsReportedRatherThanReturnedEmpty(tc)
            tc.verifyError(@() rid.io.loadRun("no-such-run", "root", tc.Root), ...
                "rid:io:recordNotFound");
        end

        function theGitCaptureNeverBlocksAndNeverThrows(tc)
            % A run must not hang or die because of git. This is not hypothetical:
            % MATLAB's `system` hands the child a terminal, so `git diff` starts
            % its pager and waits for a keypress that never comes -- the symptom
            % being a run that produces no output at all, before the first step.
            %
            % Fails if: --no-pager is dropped, or a non-repository path raises
            % instead of degrading to "unknown".
            code = rid.io.gitState(RunRecordTest.projectRoot());
            tc.verifyClass(code.commit, "string");
            tc.verifyClass(code.isDirty, "logical");
            tc.verifyGreaterThan(strlength(code.srcHash), 0);

            outside = rid.io.gitState(string(tempdir()));
            tc.verifyEqual(outside.commit, "unknown");
            tc.verifyFalse(outside.gitAvailable);
        end

    end

    methods (Static)

        function rec = record()
            %RECORD  A real record from a short run, with a non-default psiMax.
            %
            %   Built by running the loop rather than assembled by hand: a record
            %   whose shape the tests invented would drift from the one rid.rhc.run
            %   actually writes, and the drift would be invisible.
            f     = SetupTest.fixture("nSup", 8, "epsReq0", RunTest.EpsReq0);
            state = rid.crit.CoverageState(f.Q, f.epsReq, "psiMax", 2);
            setup = rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, f.fields, ...
                                  f.gen, f.sol, state, f.sch, f.fieldSch, f.qInfo);
            ov = RunTest.overrides("nSteps", 3);
            ov.crit.psiMax = 2;          % deliberately not today's default
            ov.io.storeQ   = true;
            r   = rid.rhc.run([], setup, RunTest.history("srg", f.srg), ov);
            rec = r.record;
        end

        function [rec, epsReqAtStart] = resetAndRefineRecord()
            %RESETANDREFINERECORD  A run carrying BOTH a reset and a field update.
            %
            %   The reset needs a configuration the design can cover; the field
            %   update needs a scheduler that fires with values differing from the
            %   ones the state was built with, so that epsReq actually moves.
            f = SetupTest.fixture("nSup", 10, "epsReq0", 0.75);
            f.fieldSch = rid.fields.ImmediateScheduler("updateInterval", 4, ...
                "complexity", ConstantFieldEstimator(4));
            ov = RunTest.overrides("nSteps", 10);
            ov.crit.epsReq0 = 0.75;
            ov.io.storeQ    = true;

            epsReqAtStart = f.state.epsReq;
            r   = rid.rhc.run([], SetupTest.setupFrom(f), ...
                              RunTest.history("srg", f.srg), ov);
            rec = r.record;
        end

        function r = frozenRun()
            %FROZENRUN  THE COVERAGE RESET's degenerate configuration: one horizon covers the ROI.
            %
            %   Four supporting points and a very loose epsReq0, which is the case
            %   the coverage reset records as observed. The run certifies at every step, so the
            %   resets are consecutive and the input freezes.
            r = RunTest.offlineRun("nSup", 4, "epsReq0", 1.4, "nSteps", 8);
        end

        function root = projectRoot()
            %PROJECTROOT  The folder holding src/, tests/ and runs/.
            root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
        end

    end

end
