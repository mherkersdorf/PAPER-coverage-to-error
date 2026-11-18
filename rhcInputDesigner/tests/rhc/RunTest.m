classdef RunTest < matlab.unittest.TestCase
%RUNTEST  The receding-horizon loop of the method description.
%
%   Every test names the input that would make it fail, and several exist in TWO
%   halves because the one-sided version is vacuous: a check that always refuses
%   passes a refusal test on its own, and a loop that resets its coverage and then
%   does nothing passes a reset test on its own.
%
%   The fixture is `SetupTest.fixture` -- one coherent SISO problem, first-order
%   surrogate, NARX space -- so a change to the problem shows up in both test
%   classes rather than in one of two copies.

    properties (Constant)
        EpsReq0 = 0.30
        % The joint input--state space is wider and its fixture carries fewer
        % supporting points, so it has its own baseline. It must agree with
        % `SetupTest.stateSpaceFixture`'s default: the loop refuses a state whose
        % epsReq is not what these fields and this epsReq0 produce.
        StateEpsReq0 = 0.50
    end

    methods (Test)

        % --- the loop runs, and what it produces is what it committed ---------

        function anOfflineRunDesignsASignalAndReportsItsBudget(tc)
            % The basic case: an offline run produces a signal, J falls, and a
            % record is built as part of the result.
            %
            % Fails if: the loop commits nothing, or stops before its budget.
            r = RunTest.offlineRun("nSteps", 6);

            tc.verifyEqual(r.status, "budgetReached");
            tc.verifyEqual(r.nSteps, 6);
            % One pre-run actuator position plus six designed samples.
            tc.verifyEqual(size(r.hist.u), [7 1]);
            tc.verifyEqual(numel(r.steps), 6);
            tc.verifyLessThan(r.J, r.Jtrace(1));
            tc.verifyNotEmpty(fieldnames(r.record));
            tc.verifyEqual(r.recordPath, "", ...
                "saveRun was off, so nothing may have been written");
        end

        function onlyTheAppliedUnitIsCommittedAndTheHorizonIsDiscarded(tc)
            % Invariant 1 and the applied-unit declaration. The horizon is FIVE samples and the applied
            % unit is ONE, so four samples per step are look-ahead applied to
            % nothing. A test at nHorizon = 1 would pass on a loop that commits the
            % whole trajectory, which is why the horizon is deliberately long.
            %
            % Fails if: the loop commits best.uSeq, or assumes "the first sample"
            % rather than asking the parametrization what its applied unit is.
            r = RunTest.offlineRun("nSteps", 4, "nHorizon", 5);

            tc.assertEqual(size(r.hist.u), [5 1], ...
                "four steps must commit four rows, not four horizons");

            anyLookAheadDiffered = false;
            for kStep = 1:r.nSteps
                uSeq      = reshape(r.steps(kStep).theta, [], 1);
                committed = r.hist.u(kStep + 1, :);
                tc.assertEqual(size(uSeq), [5 1]);
                tc.verifyEqual(size(r.steps(kStep).uApplied), [1 1]);
                tc.verifyEqual(committed, r.steps(kStep).uApplied, ...
                    "the committed row must be exactly what `applied` returned");
                tc.verifyEqual(committed, uSeq(1), ...
                    "for rid.gen.Direct the applied unit is the first sample");
                anyLookAheadDiffered = anyLookAheadDiffered || ...
                    any(abs(uSeq(2:end) - uSeq(1)) > 1e-9);
            end

            % THE DISCARD, MADE VISIBLE. If the look-ahead were merely a repeat of
            % the applied sample, committing the horizon would be indistinguishable
            % from committing one sample, and this test would prove nothing.
            tc.verifyTrue(anyLookAheadDiffered, ...
                "fixture is degenerate: every winning trajectory was constant");
        end

        % --- the four termination statuses (section 5) ------------------------

        function theBudgetStatusIsReachedWhenTheLengthIs(tc)
            r = RunTest.offlineRun("nSteps", 3);
            tc.verifyEqual(r.status, "budgetReached");
            tc.verifyEqual(r.nSteps, 3);
        end

        function theCertifiedStatusIsReachedInCertificationMode(tc)
            % A coarse epsReq0 over few supporting points, so certification is
            % attainable inside the guard.
            %
            % Fails if: the loop tests J == 0 itself instead of asking the
            % criterion, or never checks at all.
            r = RunTest.offlineRun("mode", "certification", "nSup", 6, ...
                "epsReq0", 1.0, "maxSteps", 60);
            tc.verifyEqual(r.status, "certified");
            tc.verifyTrue(r.certification.isSatisfied);
            tc.verifyEqual(r.J, 0);
            % Section 5 requires the achieved bound to be reported with it.
            tc.verifyEqual(r.certification.bound, ...
                r.certification.epsReqMax + r.certification.fillDistanceQ);
        end

        function theGuardStatusIsReachedWhenMaxStepsIsBelowTheBudget(tc)
            % maxSteps is a GUARD and nSteps is the BUDGET: two different rules
            % with two different statuses, and the tighter one must win.
            %
            % Fails if: the two are conflated into one counter, which would report
            % "budgetReached" here.
            r = RunTest.offlineRun("nSteps", 20, "maxSteps", 3);
            tc.verifyEqual(r.status, "guardTripped");
            tc.verifyEqual(r.nSteps, 3);
        end

        function theStalledStatusIsReachedWhenJStopsImproving(tc)
            % A stall is a legitimate outcome to report, not an error to hide: it
            % usually means the remaining supporting points are dynamically
            % unreachable. Here they are made so -- epsReq0 is far finer than the
            % design can achieve, so J flattens well inside the budget.
            %
            % Fails if: the stall counter is never incremented, or is reset by
            % anything other than an improvement or a logged discontinuity.
            r = RunTest.offlineRun("nSteps", 80, "epsReq0", 0.02, ...
                "stallSteps", 3, "nSup", 40);
            tc.verifyEqual(r.status, "stalled");
            tc.verifyLessThan(r.nSteps, 80, ...
                "a stalled run must stop before its budget, or nothing stalled");
            tail = r.Jtrace(end - 2:end);
            tc.verifyEqual(tail, repmat(tail(1), 3, 1), "AbsTol", 1e-12, ...
                "the last steps must genuinely have improved nothing");
        end

        % --- the monotonicity assert (invariant 3) ---------------------------

        function jIsMonotoneNonIncreasingWithinASegment(tc)
            % Adding samples can only shrink nearest-neighbour distances, so J is
            % non-increasing between updates. The assert is on by default, so an
            % ordinary run passing it is half the statement.
            r = RunTest.offlineRun("nSteps", 8);
            tc.verifyTrue(all(diff(r.Jtrace) <= 8 * eps(1)), ...
                sprintf("Jtrace = %s", mat2str(r.Jtrace.', 6)));
        end

        function theMonotonicityAssertFiresOnAStateWhoseJRises(tc)
            % THE HALF THAT MATTERS: an assert no input can trip is decoration.
            %
            % The violation is produced by a coverage state that reports a RISING
            % J as samples are committed -- the symptom of a stale cache, or of
            % rows committed that were never applied. It drives the loop's own
            % assert through the loop's own entry point; reproducing the comparison
            % inside the test would assert against the test's copy of it.
            %
            % Fails if: the assert is removed, compares in the wrong direction, or
            % carries a tolerance wide enough to swallow a rise of 0.1.
            tc.verifyError(@() RunTest.risingRun(true, false), ...
                "rid:rhc:monotonicityViolated");
        end

        function monotonicityAssertIsSilentWhenOffOrAcrossADiscontinuity(tc)
            % Two more halves, and each pins a different way of getting it wrong.
            %
            % Switched off, the same run must complete: an assert that cannot be
            % disabled is not a debug assert.
            r = RunTest.risingRun(false, false);
            tc.verifyEqual(r.nSteps, 4);

            % And with a resolution change logged at EVERY step, the assert must
            % not fire even though J rises, because J is not comparable across one.
            %   A loop that compared against a
            % baseline spanning the discontinuity would fail here; one that never
            % compared at all would fail the test above.
            r = RunTest.risingRun(true, true);
            tc.verifyEqual(r.nSteps, 4);
            kinds = string({r.eventLog.kind});
            tc.verifyTrue(all(kinds == "resolutionChange"));
            tc.verifyEqual(numel(kinds), 4);
        end

        % --- the budget-mode coverage reset --------------------------

        function budgetModeResetsAtCertificationAndKeepsDesigning(tc)
            % The coverage reset in both halves. A run that certifies early must NOT stop -- the
            % length is the commitment -- and it must keep DESIGNING afterwards.
            %
            % Without the second half this passes on a loop that resets and then
            % produces arbitrary input, which is the exact failure the reset exists
            % to prevent: with every point resolved the active set is empty, the
            % gradient is zero everywhere, and every start is dead on arrival.
            %
            % Fails if: the reset is not performed, or is performed without being
            % logged (the assert then fires across it), or the loop stops at J == 0
            % in budget mode.
            [r, f] = RunTest.offlineRun("nSup", 10, "epsReq0", 0.75, "nSteps", 12);

            tc.assertEqual(r.status, "budgetReached");
            tc.assertEqual(r.nSteps, 12);

            resetAt = find([r.steps.coverageReset]);
            tc.assertNotEmpty(resetAt, ...
                "fixture: the run must certify inside its budget, or the coverage reset is " + ...
                "never exercised");
            kReset = resetAt(1);
            tc.assertLessThan(kReset, 12, ...
                "fixture: the reset must leave budget to spend afterwards");

            kinds = string({r.eventLog.kind});
            tc.verifyTrue(any(kinds == "coverageReset"));
            tc.verifyEqual(r.eventLog(find(kinds == "coverageReset", 1)).step, kReset);

            % IT KEEPS DESIGNING, and this is the half that a loop which resets and
            % then produces arbitrary input would fail.
            after = r.steps(kReset + 1:end);
            tc.verifyGreaterThan(after(1).nActive, 0, ...
                "the reset must return supporting points to the active set");
            tc.verifyTrue(any([after.JBefore] > [after.JAfter]), ...
                "at least one step after the reset must IMPROVE coverage");
            uAfter = r.hist.u(kReset + 2:end, :);
            tc.verifyGreaterThan(numel(uniquetol(uAfter, 1e-9)), 1, ...
                "the input must keep moving after the reset");

            % AND THE REPORTED COVERAGE IS HONEST. Not the reset's fiction: the
            % final J must equal what a state that never reset would report over
            % the same realized history.
            fresh = rid.crit.CoverageState(f.Q, f.epsReq);
            fresh.rebuild(f.nrm.toNormalized(f.rep.build(r.hist, f.srg)));
            tc.verifyEqual(r.J, fresh.value(), "AbsTol", 1e-12, ...
                "the final rebuild must recover the coverage the run achieved");
        end

        function theResetCanBeSwitchedOff(tc)
            % The coverage reset is on by default in budget mode and switchable off. With it
            % off, the same configuration reaches J == 0 and simply stays there --
            % which is the behaviour the reset exists to replace, and seeing it
            % confirms the earlier test measured the reset rather than the fixture.
            r = RunTest.offlineRun("nSup", 10, "epsReq0", 0.75, "nSteps", 12, ...
                "resetOnCertification", false);
            tc.verifyEqual(r.status, "budgetReached");
            tc.verifyEmpty(find([r.steps.coverageReset], 1));
            tc.verifyEqual(r.Jtrace(end), 0);
        end

        % --- a reset survives a rebuild in the same step -------------

        function aResetIsNotUndoneByASurrogateRebuildInTheSameStep(tc)
            % The loop reset first and rebuilt second, so a step that did
            % both re-added the whole history and restored full coverage: the
            % reset was annulled inside the step that performed it. From then on
            % the run sat permanently in the state the coverage reset was written to prevent
            % -- empty active set, zero gradient everywhere, every start dead on
            % arrival -- while resetting again at every step.
            %
            % `requiresRebuildOnSurrogateUpdate` is true for ANY simulated-past
            % run, so this is ordinary offline design with updates enabled, not a
            % corner.
            r = RunTest.updatingRun("nSup", 10, "epsReq0", 0.55, "nSteps", 12, ...
                                    "updateInterval", 1);

            tc.assertEqual(r.status, "budgetReached");
            resetAt = find([r.steps.coverageReset]);
            tc.assertNotEmpty(resetAt, ...
                "fixture: the run must reset at least once, or the persistent reset boundary is never " + ...
                "exercised");
            kinds = string({r.eventLog.kind});
            tc.assertTrue(any(kinds == "surrogateUpdate"), ...
                "fixture: the surrogate must actually update, or the reset " + ...
                "has no rebuild to survive");

            % THE OUTCOME THE PERSISTENT RESET BOUNDARY REQUIRES. After a step that reset, the supporting
            % points are back in play. Zero here is the rebuild having undone it.
            for kStep = resetAt
                tc.verifyGreaterThan(r.steps(kStep).nActive, 0, sprintf( ...
                    "step %d reset and then rebuilt to full coverage: the " + ...
                    "reset was annulled inside its own step", kStep));
            end
        end

        function theConsecutiveResetWarningDoesNotFireForARebuild(tc)
            % THE SECOND HARM, AND THE WORSE ONE. With the reset annulled every
            % step, the run emitted the coverage reset's consecutive-reset warning naming
            % nSup, epsReq0 and nHorizon -- a correct diagnostic pointing at an
            % innocent configuration. A warning that misattributes is worse than
            % no warning: it sends the reader to change settings that were never
            % the cause.
            %
            % Fails if: the reset is undone, which manufactures the consecutive
            % run of resets that trips the warning.
            % THE FIXTURE IS NOT DEGENERATE, and that is what makes the
            % assertion meaningful. epsReq0 = 0.55 over ten supporting points
            % certifies several times in twelve steps but needs several steps to
            % re-cover after each reset, so the resets are genuinely isolated.
            % At 0.75 a single sample re-covers everything and the warning is
            % CORRECT -- the input really does freeze into a repeat -- which is
            % why that setting cannot be used to test misattribution.
            r = RunTest.updatingRun("nSup", 10, "epsReq0", 0.55, "nSteps", 12, ...
                                    "updateInterval", 1);

            kinds  = string({r.eventLog.kind});
            resets = r.eventLog(kinds == "coverageReset");
            tc.assertGreaterThan(numel(resets), 1, ...
                "fixture: several resets are needed, or ''not consecutive'' " + ...
                "is satisfied by there being nothing to be consecutive");
            u = r.hist.u(2:end, :);
            tc.assertGreaterThan(numel(uniquetol(u(end-2:end), 1e-9)), 1, ...
                "fixture: the input froze, so this configuration IS degenerate " + ...
                "and the warning would be right");

            severities = arrayfun(@(e) string(e.detail.severity), resets);
            tc.verifyEqual(severities, repmat("info", size(severities)), ...
                "consecutive resets here are the rebuild annulling each one, " + ...
                "not a degenerate configuration");
            tc.verifyEmpty(r.warnings, sprintf( ...
                "the run blamed an innocent configuration: %s", ...
                strjoin(r.warnings, " ")));
        end

        % --- the initial history is counted --------------------------

        function theCoverageStateCountsAPrimingHistory(tc)
            % The loop admitted only rows created by the step just
            % committed (`kIndex > kNow - nApplied`), so the rows `rep.build`
            % already yields from `hist0` were added by NOTHING: the run
            % overstated the work remaining and spent budget re-covering ground
            % the supplied data already covered.
            %
            % THE FIXTURE IS THE TEST. A one-row hist0 yields ZERO regressor
            % rows, so the incremental path and a rebuild agree by accident --
            % which is why this went unseen. The priming history below is long
            % enough to yield several rows, and the assert on that count is what
            % stops this test decaying into the vacuous one.
            f     = SetupTest.fixture("nSup", 12, "epsReq0", RunTest.EpsReq0);
            hist0 = RunTest.primingHistory(f.srg);

            nPrimed = size(f.rep.build(hist0, f.srg), 1);
            tc.assertGreaterThan(nPrimed, 1, ...
                "fixture: the priming history must yield several regressor " + ...
                "rows, or it cannot distinguish the two counting paths");

            r = rid.rhc.run([], SetupTest.setupFrom(f), hist0, ...
                            RunTest.overrides("nSteps", 6));

            % BEFORE STEP 1 IS WHERE THIS IS SHARP. The history-seeding rule says the priming rows
            % enter the state BEFORE the first step, and `steps(1).JBefore` is
            % that state's value. Comparing only the FINAL J is too weak: the
            % designed signal may later cover the priming rows' neighbourhood
            % anyway, and the two paths then agree for a reason that has nothing
            % to do with seeding.
            primed = rid.crit.CoverageState(f.Q, f.epsReq);
            primed.rebuild(f.nrm.toNormalized(f.rep.build(hist0, f.srg)));
            tc.assertLessThan(primed.value(), 1, ...
                "fixture: the priming history must actually resolve something, " + ...
                "or an unseeded state would score the same");
            tc.verifyEqual(r.steps(1).JBefore, primed.value(), "AbsTol", 1e-12, ...
                "the run began as though the priming history did not exist");

            % And the invariant as the history-seeding rule states it, over the whole run: the
            % incremental path and a full rebuild agree for ANY hist0.
            fresh = rid.crit.CoverageState(f.Q, f.epsReq);
            fresh.rebuild(f.nrm.toNormalized(f.rep.build(r.hist, f.srg)));
            tc.verifyEqual(r.J, fresh.value(), "AbsTol", 1e-12, ...
                "the loop's incremental count must equal a rebuild over the " + ...
                "same realized history, priming rows included");
        end

        function theTwoCountingPathsAgreeForAOneRowHistoryToo(tc)
            % THE CONTROL. This case passed before the history-seeding rule was fixed and must keep
            % passing: seeding the state must not double-count the trivial
            % history, whose regressor row count is zero.
            f = SetupTest.fixture("nSup", 12, "epsReq0", RunTest.EpsReq0);
            r = RunTest.offlineRun("nSteps", 6);

            fresh = rid.crit.CoverageState(f.Q, f.epsReq);
            fresh.rebuild(f.nrm.toNormalized(f.rep.build(r.hist, f.srg)));
            tc.verifyEqual(r.J, fresh.value(), "AbsTol", 1e-12);
        end

        % --- the guided target levels -------------------------

        function theGuidedTargetLevelsAreBuiltForADefaultRun(tc)
            % The setup-time target levels. `targetLevelsFor` was called only in the surrogate-update
            % branch, and `opts.surrogate.updateInterval` is Inf by default
            % -- so a DEFAULT run performed no update and the vector was
            % never built. Every run initialized unguided, which is the fallback
            % the surrogate-guided initialization exists to replace.
            %
            % The run below is the ordinary one: nothing here switches updates
            % on. That is the whole point -- a test that enabled updates would
            % pass on the broken loop.
            %
            % Fails if: the levels are built only after an update.
            r = RunTest.offlineRun("nSteps", 4);

            tc.assertEmpty(r.eventLog, ...
                "fixture: this run must perform NO surrogate update, or it " + ...
                "cannot distinguish the setup-time target levels from the update branch");

            levels = r.record.problem.targetLevels;
            tc.verifyTrue(levels.available, ...
                "the levels must be built from the surrogate Setup already holds");
            tc.verifyNotEqual(levels.reason, "not built yet", ...
                "a finished run must not report the levels as pending: " + ...
                """not built yet"" reads as LATER and meant NEVER");
            tc.verifyEqual(numel(levels.yAchieved), ...
                rid.defaultOptions().rhc.nTargetLevels, ...
                "opts.rhc.nTargetLevels must actually govern something");

            % A NARX space HAS an output coordinate, so it keeps tiling the
            % output range and the state-locus branch of the joint-space representation contract must not
            % capture it. This is the guard on the branch, and the canaries are
            % the guard on the values.
            tc.verifyEqual(levels.basis, "output");
            tc.verifyEmpty(levels.xAchieved, ...
                "an output-range tiling achieves no STATE values");
        end

        function theGuidedLevelsReachInitialGuesses(tc)
            % THE HALF THAT MATTERS. Building the vector and not passing it on
            % would satisfy the test above while every start stayed unguided.
            % `startInfo.guided` is the generator's own report of which rows it
            % drew from targetLevels, so this reads the loop's diagnostics rather
            % than re-deriving them.
            %
            % The surrogate-guided initialization keeps the LAST row unguided as a deliberate hedge, so the
            % assertion is "some but not all", never "all".
            r = RunTest.offlineRun("nSteps", 4);

            for kStep = 1:r.nSteps
                guided = r.steps(kStep).optDiag.startInfo.guided;
                tc.assertNotEmpty(guided, ...
                    sprintf("step %d reported no starts at all", kStep));
                tc.verifyTrue(any(guided), sprintf( ...
                    "step %d initialized entirely unguided: targetLevels did " + ...
                    "not reach initialGuesses", kStep));
                tc.verifyFalse(guided(end), ...
                    "the last row stays blind, by the surrogate-guided initialization");
            end
        end

        % --- the guided levels on the steady-state locus -------

        function guidedLevelsTileTheSteadyStateLocusWithoutAnOutputCoordinate(tc)
            % rid.rep.StateSpace has NO output coordinate, so the output-range
            % tiling of the surrogate-guided initialization finds nothing to tile and every state-space run
            % would start unguided -- which would confound the comparison of the two spaces, since the
            % comparison would then differ in its initialization as well as in
            % its space. The levels are spread along the steady-state LOCUS
            % instead.
            %
            % THE FIXTURE SEPARATES THE TWO TILINGS. Its rest state is x1 = u^3,
            % so levels equally spaced along the locus are clustered in `u` and
            % the other way round; for a linear model the two sets coincide and
            % either implementation would pass. The assertion is on `xAchieved`,
            % what the levels ACHIEVED, rather than on what was aimed at.
            %
            % Fails if: the levels are tiled in `u`, or in an arc length measured
            % along `u` rather than along the locus.
            r = RunTest.stateSpaceRun("nSteps", 2, "nStates", 1, ...
                         "srg", RunTest.cubicLocusSurrogate());
            info    = r.record.problem.targetLevels;
            nLevels = rid.defaultOptions().rhc.nTargetLevels;

            tc.assertTrue(info.available, ...
                "the locus levels were not built: " + info.reason);
            tc.verifyEqual(info.basis, "stateLocus");
            tc.assertEqual(size(info.xAchieved), [nLevels 1]);

            % The probe grid is the loop's own, so the rounding a nearest-probe
            % choice can introduce is KNOWN rather than guessed: at most half the
            % widest gap between neighbouring rest states.
            uProbe  = linspace(info.uRange(1), info.uRange(2), ...
                               rid.defaultOptions().rhc.nStaticProbes).';
            x1Probe = uProbe .^ 3;
            tolX    = max(diff(x1Probe));
            tc.verifyEqual(info.xAchieved, ...
                linspace(x1Probe(1), x1Probe(end), nLevels).', "AbsTol", tolX);

            % AND THE OTHER HALF, without which a tiling in `u` could still pass
            % on a nearly straight locus: the inputs that achieve those states are
            % strongly NON-uniform, which is the whole point of aiming at the
            % locus instead of at the input range.
            dU = diff(nthroot(info.xAchieved, 3));
            tc.verifyGreaterThan(max(dU) / min(dU), 3, ...
                "equally spaced along this locus means unequally spaced in u");
        end

        function theLocusIsTiledOnlyWhereItLiesInsideTheDeclaredStateBox(tc)
            % Step 3 of the construction: the probes are kept
            % inside the DECLARED state box. Here the locus leaves it -- the model
            % rests at [u, 0] with `u` declared over [-1 1] and both states over
            % [-0.5 0.5] -- so the levels may span only the intersection. Tiling
            % the whole locus would aim half the starts at rest states the ROI
            % excludes, which is the inconsistency the validate-never-correct rule reports rather than
            % absorbs.
            %
            % The rest locus of a position/velocity model is x2 = 0, so tiling it
            % in arc length reduces to tiling x1 and no rule is needed for which
            % state to tile.
            %
            % Fails if: the inside-the-box filter is dropped.
            xRange = [-0.5 0.5];
            r = RunTest.stateSpaceRun("nSteps", 2, "nStates", 2, ...
                         "srg", RunTest.unitGainSurrogate(), "stateRange", xRange);
            info    = r.record.problem.targetLevels;
            nLevels = rid.defaultOptions().rhc.nTargetLevels;

            tc.assertTrue(info.available, ...
                "the locus levels were not built: " + info.reason);
            tc.verifyEqual(info.basis, "stateLocus");

            % x1 = u along this locus, so the probe spacing in u IS the achievable
            % resolution in x1. Two of them: one for the nearest-probe choice, one
            % because the probe that sits exactly ON the box boundary may fall
            % either side of it by a rounding error in the rest-state solve.
            uProbe = linspace(info.uRange(1), info.uRange(2), ...
                              rid.defaultOptions().rhc.nStaticProbes).';
            tolX   = 2 * max(diff(uProbe));
            tc.verifyEqual(info.xAchieved(:, 1), ...
                linspace(xRange(1), xRange(2), nLevels).', "AbsTol", tolX, ...
                "the levels must tile the part of the locus inside the box");
            tc.verifyEqual(info.xAchieved(:, 2), zeros(nLevels, 1), ...
                "AbsTol", 1e-12, "the rest velocity is zero all along this locus");
        end

        function aLocusEntirelyOutsideTheDeclaredBoxIsReportedUnavailable(tc)
            % The other end of the same clip: with fewer than two
            % probes inside the box there is no arc length to measure, and the
            % levels are REPORTED unavailable rather than invented from the probes
            % outside it. The run then initializes unguided -- the fallback,
            % the same code path with a different vector -- and must still finish,
            % because an ROI that excludes the model's whole rest locus is a
            % modelling inconsistency for `reachabilityReport` to point at,
            % not a reason for the loop to fail.
            %
            % Fails if: the filter is dropped, or an empty locus is tiled anyway.
            r = RunTest.stateSpaceRun("nSteps", 2, "nStates", 2, ...
                         "srg", RunTest.unitGainSurrogate(), "stateRange", [-3 -2]);
            info = r.record.problem.targetLevels;

            tc.verifyEqual(r.status, "budgetReached");
            tc.assertEqual(r.nSteps, 2, "the run must have stepped at all");
            tc.verifyFalse(info.available);
            tc.verifyEqual(info.basis, "stateLocus", ...
                "the branch that could not deliver must still name itself");
            tc.verifyNotEmpty(char(info.reason), ...
                "an unavailable vector is reported WITH a reason");
            tc.verifyEmpty(info.xAchieved);

            for kStep = 1:r.nSteps
                guided = r.steps(kStep).optDiag.startInfo.guided;
                tc.assertNotEmpty(guided, ...
                    sprintf("step %d reported no starts at all", kStep));
                tc.verifyFalse(any(guided), ...
                    sprintf("step %d drew from levels that do not exist", kStep));
            end
        end

        function theChordAcrossAGapInTheLocusIsNotArcLength(tc)
            % THE CLAUSE THAT HAD NO TEST: arc length accumulates
            % over CONSECUTIVE probes only. Where the locus leaves the declared
            % state box and comes back, the chord between the last probe before
            % and the first probe after spans a piece of locus that was never
            % probed inside the box -- so it is not arc length, and bridging it
            % would stretch the measure and move every level after the gap.
            %
            % THE FIXTURE IS A TWO-STATE ONE FOR A REASON. A one-state locus
            % cannot leave an interval and re-enter it anywhere but where it
            % left, so its bridging chord is zero and both rules agree: the
            % clause is invisible on a one-state fixture, which is how it came to
            % be implemented untested. Here `x1` goes on moving while `x2` is
            % outside the box, so the exit and the re-entry are different points.
            %
            % The expectation is re-derived from the fixture's own locus under
            % BOTH rules, and the two are asserted to disagree before either is
            % compared against: a mirror that could not tell them apart would
            % prove nothing.
            %
            % Fails if: the gap is bridged.
            xRange  = [-0.5 0.5];
            srg     = RunTest.gapLocusSurrogate();
            r       = RunTest.stateSpaceRun("nSteps", 2, "nStates", 2, ...
                          "srg", srg, "stateRange", xRange);
            info    = r.record.problem.targetLevels;
            nLevels = rid.defaultOptions().rhc.nTargetLevels;

            tc.assertTrue(info.available, ...
                "the locus levels were not built: " + info.reason);
            tc.verifyEqual(info.basis, "stateLocus");

            % The loop's own probe grid, and the locus through it.
            uProbe = linspace(info.uRange(1), info.uRange(2), ...
                              rid.defaultOptions().rhc.nStaticProbes).';
            p      = srg.params;
            xProbe = [p.g1 * uProbe, p.c - p.b * uProbe.^2];
            iIn    = find(all(xProbe >= xRange(1) & xProbe <= xRange(2), 2));

            % FIXTURE, PART ONE: the locus must actually leave the box and return.
            tc.assertGreaterThan(numel(iIn), 3, ...
                "fixture: too few probes inside the box to tile anything");
            tc.assertTrue(any(diff(iIn) > 1), ...
                "fixture: the locus must LEAVE the declared box and re-enter " + ...
                "it, or there is no gap for the clause to refuse");

            xNorm    = (xProbe - xRange(1)) / (xRange(2) - xRange(1));
            unbridged = RunTest.locusLevels(xProbe, xNorm, iIn, nLevels, false);
            bridged   = RunTest.locusLevels(xProbe, xNorm, iIn, nLevels, true);

            % FIXTURE, PART TWO: the two rules must DISAGREE on this locus, or
            % the comparison below is satisfied by an implementation that bridges.
            tc.assertNotEqual(bridged, unbridged, ...
                "fixture: bridging the gap must change the chosen levels");

            tc.verifyEqual(info.xAchieved, unbridged, "AbsTol", 1e-12);
        end

        function aLocusOfZeroArcLengthIsReportedUnavailable(tc)
            % THE THIRD WAY THE LOCUS CANNOT BE TILED, and the only one that was
            % not a refusal but a silent success (the joint-space representation contract). Every probed rest state is the SAME state, so there is
            % no arc to spread levels along and the nearest probe to every target
            % is the first one: the construction returned nine copies of one
            % input and reported them available. That is a guided initialization
            % in name only -- nine starts where the surrogate-guided initialization asked for a spread -- and it
            % is worse than the unguided fallback, which at least spreads its
            % starts blindly.
            %
            % `Bd = 0` is the honest form of the case: an input that never reaches
            % the state. Nothing DISAGREES here, as it does when the locus leaves
            % the declared box; the locus is simply a point, and a point has no
            % arc length.
            %
            % Fails if: a zero arc length is tiled anyway.
            srg  = RunTest.flatLocusSurrogate();
            r    = RunTest.stateSpaceRun("nSteps", 2, "nStates", 1, "srg", srg);
            info = r.record.problem.targetLevels;

            % THE FIXTURE SHOWS THE DEFECT. If the rest states at
            % the two ends of the u range differed at all, the locus would have an
            % arc length and this would be the ordinary available case.
            tc.assertEqual(srg.steadyState(info.uRange(1)), ...
                           srg.steadyState(info.uRange(2)), ...
                "fixture: every input must rest at the same state");

            tc.verifyEqual(r.status, "budgetReached");
            tc.assertEqual(r.nSteps, 2, "the run must have stepped at all");
            tc.verifyFalse(info.available);
            tc.verifyEqual(info.basis, "stateLocus", ...
                "the branch that could not deliver must still name itself");
            tc.verifyNotEmpty(char(info.reason), ...
                "an unavailable vector is reported WITH a reason");
            tc.verifyEmpty(info.xAchieved);

            for kStep = 1:r.nSteps
                guided = r.steps(kStep).optDiag.startInfo.guided;
                tc.assertNotEmpty(guided, ...
                    sprintf("step %d reported no starts at all", kStep));
                tc.verifyFalse(any(guided), sprintf( ...
                    "step %d drew from levels that carry no spread", kStep));
            end
        end

        function bothOfflineStateSpaceRunsBuildLocusLevelsThatReachTheStarts(tc)
            % THE PAIR OF OFFLINE STATE-SPACE RUNS: the linear model and the smoothly nonlinear one
            % must BOTH initialize guided, on the ordinary fixture rather than on
            % one built to make a locus visible.
            %
            % The tests above measure what `targetLevelsFor` COMPUTES, on fixtures
            % chosen to separate a tiling of the locus from a tiling of the input
            % range. What none of them reads is whether the vector ARRIVES. That is
            % the defect exactly -- a vector built and never passed on, with
            % every run initializing unguided and nothing saying so -- and the
            % state-locus basis is a second place it can be dropped. It would
            % confound the comparison, whose two arms must differ in their SPACE and not in
            % their initialization.
            %
            % Fails if: the locus branch reports levels it did not produce, or the
            % vector does not reach `initialGuesses`.
            fixtures = struct("name", {"LinearStateSpace", "StateSpaceFunction"}, ...
                              "srg",  {[], RunTest.duffingSurrogate()});
            nLevels  = rid.defaultOptions().rhc.nTargetLevels;

            for iCase = 1:numel(fixtures)
                what = fixtures(iCase).name;
                r    = RunTest.stateSpaceRun("nSteps", 2, "srg", fixtures(iCase).srg);
                info = r.record.problem.targetLevels;

                tc.assertTrue(info.available, what + ": " + info.reason);
                tc.verifyEqual(info.basis, "stateLocus", what + ...
                    ": a space with no output coordinate tiles the rest locus");
                tc.verifyEqual(size(info.xAchieved), [nLevels 2], ...
                    what + ": one achieved rest state per level");

                for kStep = 1:r.nSteps
                    guided = r.steps(kStep).optDiag.startInfo.guided;
                    tc.assertNotEmpty(guided, sprintf( ...
                        "%s step %d reported no starts at all", what, kStep));
                    tc.verifyTrue(any(guided), sprintf( ...
                        "%s step %d initialized entirely unguided: the locus " + ...
                        "levels did not reach initialGuesses", what, kStep));
                    tc.verifyFalse(guided(end), ...
                        what + ": the last row stays blind, by the surrogate-guided initialization");
                end
            end
        end

        % --- online, offline, and the incoherent pair ----------------

        function offlineWithMeasuredPastSourceIsRefused(tc)
            % An offline run has nothing to measure. Two declarations the user
            % made, compared and neither corrected.
            f = SetupTest.fixture("pastSource", "measured", ...
                                  "epsReq0", RunTest.EpsReq0);
            tc.verifyError(@() rid.rhc.run([], SetupTest.setupFrom(f), ...
                RunTest.history("srg", f.srg), RunTest.overrides("nSteps", 2)), ...
                "rid:rhc:incoherentDeclaration");
        end

        function aHorizonTooShortToScoreTheAppliedInputIsRefused(tc)
            % The shortest-horizon declaration, at assembly. With nk = 1 the applied sample u(N+1) first
            % appears in the regressor row at k = N+2, so a one-sample horizon
            % scores nothing that depends on it: J is constant in the sample about
            % to be committed, every start evaluates the same, and the loop commits
            % whatever the winning start happened to hold -- silently.
            %
            % Fails if: the refusal is dropped, or reads nk instead of asking the
            % representation for its floor.
            tc.verifyError(@() RunTest.delayedRun(1, 1), ...
                "rid:rhc:horizonTooShort");

            % The other half: the floor itself must run. A check that always
            % refused would pass the line above on its own.
            r = RunTest.delayedRun(1, 2);
            tc.verifyEqual(r.status, "budgetReached");
            tc.verifyEqual(r.nSteps, 2);
        end

        function theHorizonFloorDoesNotBiteAtTheDefaultInputDelay(tc)
            % nk = 0 puts u(k) in every row, so the floor is 1 and no ordinary run
            % is affected. The fixture uses nHorizon = 4; the sharp case is
            % nHorizon = 1, which must still run.
            r = RunTest.delayedRun(0, 1);
            tc.verifyEqual(r.status, "budgetReached");
        end

        function bothCoherentDeclarationsRun(tc)
            % THE OTHER HALF. A check that always refuses would pass the test above
            % on its own, so both coherent pairings must actually run: offline with
            % "simulated", and online with "measured".
            offline = RunTest.offlineRun("nSteps", 2);
            tc.verifyEqual(offline.nSteps, 2);

            online = RunTest.onlineRun("nSteps", 2, "pastSource", "measured");
            tc.verifyEqual(online.nSteps, 2);
            tc.verifyEqual(size(online.hist.y), [3 1]);
        end

        function anOnlineRunRecordsWhatTheProcessMeasured(tc)
            % The process contract's shape: the applied unit in, its measurement out. The plant
            % differs from the surrogate on purpose, so a loop that stored its own
            % prediction instead would agree with the model rather than the plant.
            [r, plant] = RunTest.onlineRun("nSteps", 5, "pastSource", "simulated");
            tc.verifyEqual(plant.nApplications, 5, ...
                "the process must be applied exactly once per step");
            tc.verifyEqual(size(r.hist.y), [6 1]);
            tc.verifyTrue(all(isfinite(r.hist.y)));

            % The measurement is the PLANT's, not the surrogate's prediction.
            predicted = r.record.setup.setup.srg.simulate(r.hist.u, r.hist.x0);
            tc.verifyGreaterThan(max(abs(predicted - r.hist.y)), 1e-6, ...
                "fixture: the plant and the surrogate must differ, or this " + ...
                "test cannot tell a measurement from a prediction");
        end

        function aMisshapedMeasurementRaisesAndIsNotCaught(tc)
            % Recovery is policy and stays out, so the error must reach the
            % caller rather than being absorbed into a status.
            for mode = ["tooManyRows", "notFinite"]
                f = SetupTest.fixture("pastSource", "measured", ...
                                      "epsReq0", RunTest.EpsReq0);
                ov = RunTest.overrides("nSteps", 2);
                tc.verifyError(@() rid.rhc.run(MisshapedProcess(mode), ...
                    SetupTest.setupFrom(f), ...
                    RunTest.history("srg", f.srg, "online", true), ov), ...
                    "rid:rhc:invalidMeasurement", "mode = " + mode);
            end
        end

        % --- the joint input--state space in the loop ----------------

        function anOfflineStateSpaceRunDesignsASignalAndCountsItsRestPoint(tc)
            % The loop over rid.rep.StateSpace, offline, with the simulated past
            %. Three things at once, because they are the three
            % that a representation of look-back depth ZERO changes:
            %
            %   - it runs at all: the window is EMPTY here, so every step
            %     asks the surrogate for a zero-length simulation, which the
            %     contract of the joint-space representation contract had to be extended to allow;
            %   - the analytic path is still taken: the state's
            %     Markov parameters replace the output's, and `ed` and the
            %     starting state enter only the constant term, so the map from
            %     theta to the normalized points stays affine;
            %   - the REST POINT is already a point.
            %
            % Fails if: the zero-length window raises; the affine check is
            % skipped; or the one-row history is treated as yielding nothing.
            [r, f] = RunTest.stateSpaceRun("nSteps", 6);

            tc.verifyEqual(r.status, "budgetReached");
            tc.verifyEqual(r.nSteps, 6);
            tc.verifyEqual(size(r.hist.u), [7 1]);
            tc.verifyFalse(isfield(r.hist, "x"), ...
                "an offline run measures nothing, so no hist.x may appear");
            tc.verifyEqual(r.steps(1).windowDepth, 0, ...
                "the look-back depth of this space is zero");

            tc.verifyEqual(r.steps(1).gradientPath, "analytic");
            tc.verifySubstring(r.steps(1).affineTier, "affine");

            % THE REST POINT IS ALREADY A POINT. The
            % one-row hist0 yields [u(1), x0] -- the rest state under the pre-run
            % input, realized data that must be in the coverage state BEFORE step
            % 1. For every NARX shape with a lag that row does not exist, which is
            % what the xStart requirement's "a single row cannot form a regressor vector" was about
            % and what makes this case new rather than a repeat.
            hist0 = RunTest.history("srg", f.srg);
            tc.assertEqual(size(f.rep.build(hist0, f.srg), 1), 1, ...
                "fixture: the one-row history must yield exactly the rest point");
            primed = rid.crit.CoverageState(f.Q, f.epsReq);
            primed.rebuild(f.nrm.toNormalized(f.rep.build(hist0, f.srg)));
            tc.assertLessThan(primed.value(), 1, ...
                "fixture: the rest point must resolve something, or a state " + ...
                "that never saw it would score the same");
            tc.verifyEqual(r.steps(1).JBefore, primed.value(), "AbsTol", 1e-12, ...
                "the run began as though the rest point did not exist");

            % And the history-seeding rule's invariant over the whole run: the incremental count and
            % a full rebuild agree.
            fresh = rid.crit.CoverageState(f.Q, f.epsReq);
            fresh.rebuild(f.nrm.toNormalized(f.rep.build(r.hist, f.srg)));
            tc.verifyEqual(r.J, fresh.value(), "AbsTol", 1e-12);
        end

        function aNonlinearStateSpaceSurrogateTakesTheValueOnlyPath(tc)
            % THE OTHER HALF of the tier detection, and the reason the test
            % process is SMOOTHLY nonlinear. The affine check must
            % genuinely fail here and the loop must drop to the value-only
            % objective -- not raise, and not quietly keep using a Jacobian
            % measured once at the first step.
            %
            % A piecewise-linear process would PASS the check whenever its probe
            % and check points share one linear region, and this test would then
            % report the analytic tier for a nonlinear model. The Duffing cubic
            % has no linear region to hide in.
            %
            % Fails if: the affine verdict is assumed rather than measured, or
            % `useAnalytic` is not lowered when the check fails.
            r = RunTest.stateSpaceRun("nSteps", 4, "srg", RunTest.duffingSurrogate());

            tc.verifyEqual(r.status, "budgetReached");
            tc.verifyEqual(r.steps(1).gradientPath, "finiteDifference");
            tc.verifySubstring(r.steps(1).affineTier, "notAffine");
        end

        function anOnlineRunCollectsMeasuredStatesAndAnchorsEachHorizonOnThem(tc)
            % The joint-space representation contract, both halves of it.
            %
            % The loop asks a state-measuring process for its second output and
            % appends it to `hist.x`; and each horizon starts ONE SURROGATE STEP
            % FROM THE LAST MEASURED STATE rather than from a state re-simulated
            % through the whole history. Without the anchor a measured state would
            % change the realized points and never the design -- the run would
            % collect measurements and ignore them, which is silent.
            %
            % Fails if: the second output is never requested; the states are not
            % appended; or `xNow` is re-derived from hist.x0 as it is for NARX.
            [r, f, plant] = RunTest.measuredStateRun("nSteps", 5);

            tc.verifyEqual(r.nSteps, 5);
            tc.verifyEqual(plant.nApplications, 5, ...
                "the process must be applied exactly once per step");
            tc.verifyEqual(plant.lastNargout, 2, ...
                "a representation that reads measured states must make the " + ...
                "loop ask for them");
            tc.verifyEqual(size(r.hist.x, 1), size(r.hist.u, 1), ...
                "hist.x carries one row per applied sample, as hist.y does");
            tc.verifyEqual(size(r.hist.x, 2), 2);
            tc.verifyTrue(all(isfinite(r.hist.x), "all"));

            % THE ANCHOR, COMPUTED INDEPENDENTLY. At the start of step k the
            % history ended at row k, so the horizon started one surrogate step
            % from the state measured there.
            for kStep = 1:r.nSteps
                [~, anchored] = f.srg.simulate(r.hist.u(kStep, :), ...
                                               r.hist.x(kStep, :));
                tc.verifyEqual(r.steps(kStep).xNow, anchored, "AbsTol", 1e-12, ...
                    sprintf("step %d did not start from the measured state", kStep));
            end

            % THE FIXTURE SHOWS THE DEFECT. The alternative --
            % `xNow` re-simulated from hist.x0 through the whole input history --
            % has to be a DIFFERENT number, or the assertion above is satisfied by
            % a loop that never anchored anything. It is different exactly because
            % the plant is not the surrogate.
            [~, reSimulated] = f.srg.simulate(r.hist.u(1:r.nSteps, :), r.hist.x0);
            tc.assertGreaterThan(max(abs(reSimulated - r.steps(end).xNow)), 1e-3, ...
                "fixture: the plant and the surrogate must disagree, or the " + ...
                "anchored and the re-simulated state cannot be told apart");
        end

        function aProcessThatCannotMeasureItsStatesIsRefused(tc)
            % Three ways to fail the two-output contract, all of them silent
            % otherwise: a process that returns one output only would stop the run
            % with MATLAB's own "Too many output arguments", which names the loop
            % rather than the declaration that asked; and a mis-shaped or
            % non-finite `x` would be appended to a history whose width then
            % disagrees with the representation's, surfacing one `build` later.
            %
            % Fails if: the translation or either shape check is removed.
            tc.verifyError(@() RunTest.measuredStateRun("nSteps", 2, ...
                "process", FirstOrderProcess(0.8, 0.1, 1.4)), ...
                "rid:rhc:invalidMeasurement");

            for defect = ["wrongWidth", "notFinite"]
                tc.verifyError(@() RunTest.measuredStateRun("nSteps", 2, ...
                    "defect", defect), "rid:rhc:invalidMeasurement", ...
                    "defect = " + defect);
            end

            % The other half: a check that always refused would pass the three
            % above on its own.
            r = RunTest.measuredStateRun("nSteps", 2);
            tc.verifyEqual(r.nSteps, 2);
        end

        function aMeasuredStateRunWithoutInitialStatesIsRefused(tc)
            % The fifth hist0 check of the method description.
            % The loop appends one state per step from the first, so an `x` with
            % fewer rows than `u` is a malformed history one append later -- and
            % an absent `x` is how a history says NO states were measured, which
            % must not be silently turned into a run of a different kind.
            %
            % Fails if: the check is dropped, or the mode is inferred from the
            % presence of the field instead of being asked of the representation.
            f     = SetupTest.stateSpaceFixture(2, "pastSource", "measured");
            hist0 = RunTest.history("srg", f.srg, "online", true);   % no `x`
            tc.verifyError(@() rid.rhc.run(RunTest.statePlant(2), ...
                SetupTest.setupFrom(f), hist0, RunTest.stateOverrides(2)), ...
                "rid:rhc:invalidInitialHistory");
        end

        function anOfflineRunThatReadsMeasuredStatesIsRefused(tc)
            % The process contract's incoherent pair, in the space where it now has a second
            % meaning: an offline run has nothing to measure, states included.
            % Two declarations the user made, compared and neither corrected.
            f = SetupTest.stateSpaceFixture(2, "pastSource", "measured");
            tc.verifyError(@() rid.rhc.run([], SetupTest.setupFrom(f), ...
                RunTest.measuredHistory(f.srg, zeros(1, 2)), ...
                RunTest.stateOverrides(2)), "rid:rhc:incoherentDeclaration");
        end

        function aNarxRunAsksAStateMeasuringProcessForItsOutputOnly(tc)
            % THE COMPATIBILITY HALF. A process that CAN report
            % its states must not be asked for them by a representation that does
            % not read them: the one-output form is unchanged, and every
            % existing adapter goes on working. The process reports how many
            % outputs it was asked for, because nothing else can show that the
            % loop did not ask.
            %
            % Fails if: the loop requests two outputs whenever the process can
            % supply them, which would also put an unused hist.x into the history
            % and into every record built from it.
            f     = SetupTest.fixture("epsReq0", RunTest.EpsReq0);
            plant = StateSpaceProcess(0.75, 0.2, 1, 0);
            r     = rid.rhc.run(plant, SetupTest.setupFrom(f), ...
                        RunTest.history("srg", f.srg, "online", true), ...
                        RunTest.overrides("nSteps", 3));

            tc.verifyEqual(plant.lastNargout, 1, ...
                "a NARX representation reads no measured states");
            tc.verifyFalse(isfield(r.hist, "x"), ...
                "no hist.x may appear merely because the process could supply one");
            tc.verifyEqual(size(r.hist.y), [4 1]);
        end

        % --- options, and the commitments that have no default ---------------

        function theCommitmentsWithoutDefaultsAreRequired(tc)
            f     = SetupTest.fixture("epsReq0", RunTest.EpsReq0);
            setup = SetupTest.setupFrom(f);

            hist0 = RunTest.history("srg", f.srg);

            noEps      = RunTest.overrides("nSteps", 2);
            noEps.crit = rmfield(noEps.crit, "epsReq0");
            tc.verifyError(@() rid.rhc.run([], setup, hist0, noEps), ...
                "rid:rhc:missingOption");

            tc.verifyError(@() rid.rhc.run([], setup, hist0, ...
                RunTest.overrides()), "rid:rhc:missingOption");
        end

        function theInitialConditionArrivesAsAHistoryAndIsValidated(tc)
            % an initial condition is DATA, and the object for it
            % already exists. `rid.assertHistory` is its one structural validator
            %, so a malformed hist0 is refused by the same function every
            % other consumer uses rather than by a check invented here.
            %
            % Fails if: run assembles a history from options again, or validates
            % hist0 with its own partial copy of assertHistory.
            f     = SetupTest.fixture("epsReq0", RunTest.EpsReq0);
            setup = SetupTest.setupFrom(f);
            ov    = RunTest.overrides("nSteps", 2);

            tc.verifyError(@() rid.rhc.run([], setup, ...
                struct("u", 0, "tS", 0.1), ov), "rid:missingHistoryField");
            tc.verifyError(@() rid.rhc.run([], setup, ...
                struct("u", 0, "x0", 0, "tS", -1), ov), "rid:invalidHistory");

            % And the run's first row IS the supplied one: the pre-run actuator
            % position is the caller's declaration, not something run invents.
            hist0 = RunTest.history("srg", f.srg, "u0", 0.4);
            r     = rid.rhc.run([], setup, hist0, ov);
            tc.verifyEqual(r.hist.u(1, :), 0.4);
            tc.verifyEqual(r.hist.x0, hist0.x0);
            tc.verifyEqual(r.hist.tS, hist0.tS);
        end

        function aHistoryThatDoesNotMatchTheProblemIsRefused(tc)
            % Three mismatches, each silent otherwise. The channel count would
            % transpose into a differently shaped design; an online run without
            % measurements would produce a y shorter than its u, which is a
            % malformed history one append later; and a surrogate discretized at a
            % different tS advances every simulated horizon by a different amount
            % of physical time than the data it is fitted to.
            f     = SetupTest.fixture("epsReq0", RunTest.EpsReq0);
            setup = SetupTest.setupFrom(f);
            ov    = RunTest.overrides("nSteps", 2);

            wide = RunTest.history("srg", f.srg);
            wide.u = [0 0];
            tc.verifyError(@() rid.rhc.run([], setup, wide, ov), ...
                "rid:rhc:invalidInitialHistory");

            offlineShaped = RunTest.history("srg", f.srg);   % no `y`
            tc.verifyError(@() rid.rhc.run(FirstOrderProcess(0.8, 0.1, 1.4), ...
                setup, offlineShaped, ov), "rid:rhc:invalidInitialHistory");

            wrongRate    = RunTest.history("srg", f.srg);
            wrongRate.tS = f.srg.tS * 2;
            tc.verifyError(@() rid.rhc.run([], setup, wrongRate, ov), ...
                "rid:rhc:inconsistentSampleTime");
        end

        function aRunWithNoFiniteStoppingRuleIsRefused(tc)
            % assembly refuses a configuration in which
            % no stopping rule is finite, rather than inventing a cap. Certification
            % may never be reached -- the remaining supporting points can be
            % dynamically unreachable, which is what "stalled" exists to report --
            % so every guard at Inf would loop forever.
            %
            % Fails if: a default cap is introduced, or the check reads only some
            % of the four rules.
            f  = SetupTest.fixture("epsReq0", RunTest.EpsReq0, ...
                                   "mode", "certification");
            ov = RunTest.overrides();
            ov.crit.mode = "certification";
            hist0 = RunTest.history("srg", f.srg);
            tc.verifyError(@() rid.rhc.run([], SetupTest.setupFrom(f), hist0, ov), ...
                "rid:rhc:noStoppingRule");

            % THE OTHER HALF, and one case per guard: ANY of the four being finite
            % is enough. A check that refused unless a particular one was set would
            % pass the test above and fail three of these.
            %
            % epsReq0 is coarser here and nSup smaller, because a CERTIFICATION run
            % must first be certifiable at all: resolutionReport refuses a baseline
            % finer than the discretization can resolve, and that refusal
            % would mask the one being tested.
            for guard = ["maxSteps", "maxWallTime", "stallSteps"]
                certifiable = SetupTest.fixture("nSup", 6, "epsReq0", 1.0, ...
                                                "mode", "certification");
                relaxed              = RunTest.overrides();
                relaxed.crit.mode    = "certification";
                relaxed.crit.epsReq0 = 1.0;
                relaxed.rhc.(guard)  = 3;
                r = rid.rhc.run([], SetupTest.setupFrom(certifiable), ...
                                RunTest.history("srg", certifiable.srg), relaxed);
                tc.verifyNotEmpty(r.status, "guard = " + guard);
            end

            % And in budget mode nSteps alone suffices, since it IS a stopping rule.
            budget = RunTest.offlineRun("nSteps", 2);
            tc.verifyEqual(budget.status, "budgetReached");
        end

        function everyStoppingRuleOptionComesFromDefaultOptions(tc)
            % the option tables live in
            % rid.defaultOptions and nowhere else, so the effective options struct
            % a record stores says what the run actually did.
            %
            % Fails if: a value is hard-coded at its use site, in which case it is
            % absent here and the record's effective struct is incomplete; or if a
            % documented option is missing, in which case an override of it would be
            % refused as an unknown name.
            opts = rid.defaultOptions();

            tc.verifyEqual(opts.rhc.maxSteps,             Inf);
            tc.verifyEqual(opts.rhc.maxWallTime,          Inf);
            tc.verifyEqual(opts.rhc.stallSteps,           Inf);
            tc.verifyEqual(opts.rhc.stallTolerance,       0);
            tc.verifyTrue(opts.rhc.assertMonotonicity);
            tc.verifyTrue(opts.rhc.resetOnCertification);
            tc.verifyEqual(opts.rhc.nTargetLevels,        9);
            tc.verifyEqual(opts.rhc.nStaticProbes,        201);

            tc.verifyTrue(opts.io.saveRun);
            tc.verifyFalse(opts.io.storeQ);
            tc.verifyFalse(opts.io.requireCleanTree);
            tc.verifyEqual(opts.io.runsDir, "runs");
            tc.verifyEqual(opts.gen.affineProbeStep, 1e-3);
            tc.verifyEqual(opts.gen.affineCheckTol,  1e-9);
            tc.verifyEqual(opts.io.slug,    "run");

            % nSteps has NO default, on the epsReq0 precedent: the budget
            % length is the commitment in budget mode. Its absence is what makes
            % rid.rhc.run able to refuse a run that omits it.
            tc.verifyFalse(isfield(opts.rhc, "nSteps"));
            % And the initial condition is not an option at all any more.
            tc.verifyFalse(isfield(opts.rhc, "u0"));
            tc.verifyFalse(isfield(opts.rhc, "y0"));
            tc.verifyFalse(isfield(opts.rhc, "fields"));
        end

        function theAnalyticPathIsTakenAndReported(tc)
            % The affine check must PASS for the default configuration, and the
            % path taken must be visible: a fall-back to finite differences cannot
            % masquerade as the fast path.
            r = RunTest.offlineRun("nSteps", 3);
            tc.verifyEqual(r.steps(1).gradientPath, "analytic");
            tc.verifySubstring(r.steps(1).affineTier, "affine");
        end

        function pNormInfFallsBackToTheGradientFreePathWithoutFailing(tc)
            % pNorm = Inf is a supported mode of the criterion whose GRADIENT is
            % not implemented. The loop must discover that by CALLING and
            % drop to the tier-3 path, not raise -- and must report which ran.
            %
            % Fails if: the loop reads opts.crit.pNorm and branches on it, or
            % catches every error rather than exactly the one identifier.
            f     = SetupTest.fixture("epsReq0", RunTest.EpsReq0);
            state = rid.crit.CoverageState(f.Q, f.epsReq, "pNorm", Inf);
            setup = rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, f.fields, ...
                                  f.gen, f.sol, state, f.sch, f.fieldSch, f.qInfo);
            ov = RunTest.overrides("nSteps", 2);
            ov.crit.pNorm = Inf;
            r = rid.rhc.run([], setup, RunTest.history("srg", f.srg), ov);
            tc.verifyEqual(r.steps(1).gradientPath, "finiteDifference");
            tc.verifyEqual(r.nSteps, 2);
        end

    end

    methods (Static)

        function hist0 = history(varargin)
            %HISTORY  The experiment's initial condition, as the history struct.
            %
            %   `u(1,:)` is the actuator's position before the run began -- one
            %   row, physical -- and `x0` comes from the canonical producer,
            %   `srg.steadyState`. `y` is present only for an ONLINE run,
            %   because an absent y is how an offline history says it has none.
            args  = struct(varargin{:});
            u0    = RunTest.pick(args, "u0", 0);
            srg   = RunTest.pick(args, "srg", rid.surrogate.LinearFirstOrder(0.5, 0.1, 1));
            hist0 = struct("u", u0, "x0", srg.steadyState(u0), "tS", srg.tS);
            if RunTest.pick(args, "online", false)
                hist0.y = RunTest.pick(args, "y0", 0);
            end
        end

        function hist0 = primingHistory(srg)
            %PRIMINGHISTORY  An initial condition that already yields regressors.
            %
            %   `u(1,:)` is still the pre-run actuator position; the rows
            %   after it are data the experiment starts with. Long enough that
            %   `rep.build` returns several rows, which is what the history-seeding rule's invariant
            %   needs in order to be testable at all.
            uPrime = [0; 0.4; -0.3; 0.6; -0.5];
            hist0  = struct("u", uPrime, "x0", srg.steadyState(uPrime(1)), ...
                            "tS", srg.tS);
        end

        function removeIfPresent(folder)
            %REMOVEIFPRESENT  Teardown that survives a test failing before it wrote.
            if isfolder(folder)
                rmdir(folder, "s");
            end
        end

        function value = pick(args, name, fallback)
            value = fallback;
            if isfield(args, name)
                value = args.(name);
            end
        end

        function ov = overrides(varargin)
            %OVERRIDES  The commitments every run needs, plus what a test varies.
            args    = struct(varargin{:});
            ov      = struct();
            ov.crit = struct("epsReq0", RunTest.EpsReq0);
            ov.rhc  = struct();
            ov.io   = struct("saveRun", false);
            for name = string(fieldnames(args)).'
                switch name
                    case {"epsReq0", "mode", "pNorm"}
                        ov.crit.(name) = args.(name);
                    otherwise
                        ov.rhc.(name) = args.(name);
                end
            end
        end

        function [r, f] = offlineRun(varargin)
            %OFFLINERUN  A run with no process; pastSource is "simulated".
            %
            %   The fixture is returned as well, for the tests that must rebuild a
            %   quantity independently rather than read the loop's own report.
            args    = struct(varargin{:});
            [f, ov] = RunTest.problem(args, "simulated");
            r       = rid.rhc.run([], SetupTest.setupFrom(f), ...
                                  RunTest.history("srg", f.srg), ov);
        end

        function [r, plant] = onlineRun(varargin)
            %ONLINERUN  A run against a plant that is NOT the surrogate.
            args   = struct(varargin{:});
            source = "simulated";
            if isfield(args, "pastSource")
                source = args.pastSource;
            end
            [f, ov] = RunTest.problem(args, source);
            plant   = FirstOrderProcess(0.8, 0.1, 1.4);
            r       = rid.rhc.run(plant, SetupTest.setupFrom(f), ...
                                  RunTest.history("srg", f.srg, "online", true), ov);
        end

        function r = delayedRun(nk, nHorizon)
            %DELAYEDRUN  An offline run over a NARX space with an input delay.
            %
            %   The priming history is TWO rows rather than the usual one, because
            %   at nk = 1 the space discards two warm-up samples and a
            %   one-row history cannot fill the first horizon at all. That raises
            %   `rid:rhc:historyTooShort` -- a different refusal from the one under
            %   test, and one that would mask it.
            f     = SetupTest.fixture("nk", nk, "nHorizon", nHorizon, ...
                        "pastSource", "simulated", "epsReq0", RunTest.EpsReq0);
            hist0 = struct("u", [0; 0.2], "x0", f.srg.steadyState(0), ...
                           "tS", f.srg.tS);
            r     = rid.rhc.run([], SetupTest.setupFrom(f), hist0, ...
                                RunTest.overrides("nSteps", 2));
        end

        function hist0 = measuredHistory(srg, xMeasured0)
            %MEASUREDHISTORY  hist0 for a run whose representation reads states.
            %
            %   `x(1,:)` is the state measured before the run began (the joint-space representation contract), and it is the PLANT's: `hist.x` is what the process reported,
            %   while `hist.x0` is the surrogate's own initial state. The two are
            %   different quantities and the loop keeps them apart.
            hist0   = RunTest.history("srg", srg, "online", true);
            hist0.x = xMeasured0;
        end

        function ov = stateOverrides(nSteps)
            %STATEOVERRIDES  The commitments of a state-space run, at its baseline.
            ov = RunTest.overrides("epsReq0", RunTest.StateEpsReq0, ...
                                   "nSteps", nSteps);
        end

        function plant = statePlant(nStates, opts)
            %STATEPLANT  A plant in the surrogate's state basis but not its model.
            %
            %   `SetupTest.stateSpaceSurrogate` is diag(0.5 + 0.1*(0:n-1)) with
            %   Bd = 0.1*(1:n)'. Both the poles and the gain move here, so a
            %   horizon anchored on a measurement and one re-simulated through the
            %   surrogate diverge -- which is what makes the joint-space representation contract testable.
            %   The BASIS is deliberately shared: that is the assumption a
            %   measured-state run rests on and no code can check (part 4).
            arguments
                nStates     (1,1) double {mustBeInteger, mustBePositive}
                opts.defect (1,1) string = "none"
            end
            Ad    = diag(0.62 + 0.08 * (0:nStates - 1));
            Bd    = 0.17 * (1:nStates).';
            Cd    = [1, zeros(1, nStates - 1)];
            plant = StateSpaceProcess(Ad, Bd, Cd, zeros(1, nStates), ...
                                      "defect", opts.defect);
        end

        function srg = duffingSurrogate()
            %DUFFINGSURROGATE  The smoothly nonlinear model.
            %
            %   A damped, driven Duffing oscillator, discretized by fixed-step RK4
            %   at `tS`. Damped and with a positive linear stiffness, so it has a
            %   single rest state at every input in the ROI and the base-class
            %   `steadyState` fallback settles.
            %
            %   THE HORIZON MUST REACH THE NONLINEAR REGION, or the fixture is a
            %   linear one in disguise. The affine check compares an absolute
            %   deviation in NORMALIZED coordinates against 1e-9, and at a fine
            %   sample time a four-sample horizon from rest leaves `x1` at a few
            %   thousandths, where `beta*x1^3` is of the order of that tolerance
            %   and the verdict is decided by rounding. `tS` and the gain are set
            %   so the states travel a visible fraction of their declared box
            %   within the horizon, where the cubic term is a few tens of percent
            %   of the linear one. Still well inside the box, so the soft limits
            %   of the hard/soft constraint split do not dominate the objective instead.
            params = struct("tS", 0.2, "delta", 0.3, "alpha", 1, ...
                            "beta", 2, "gain", 2);
            srg    = rid.surrogate.StateSpaceFunction(@duffingStepFcn, ...
                         @duffingOutputFcn, params, 2, 1, params.tS);
        end

        function srg = cubicLocusSurrogate()
            %CUBICLOCUSSURROGATE  A one-state model whose rest state is u^3.
            %
            %   The fixture of the locus tiling: a curved
            %   steady-state locus, so tiling it and tiling the input range are
            %   different sets of levels. `linearOutputFcn` supplies the output,
            %   which the locus construction never reads -- it is the STATES this
            %   space is tiled over.
            params = struct("a", 0.6, "Cd", 1, "directTerm", 0, "fd", 0);
            srg    = rid.surrogate.StateSpaceFunction(@cubicLocusStepFcn, ...
                         @linearOutputFcn, params, 1, 1, 0.1);
        end

        function srg = gapLocusSurrogate()
            %GAPLOCUSSURROGATE  A two-state model whose rest locus leaves the box and returns elsewhere.
            %
            %   The fixture of the gap clause; `gapLocusStepFcn`
            %   carries the geometry and the reason it needs two states. `Cd`
            %   reads `x1`, which the locus construction never looks at -- it is
            %   the STATES this space is tiled over.
            params = struct("a", 0.5, "g1", 0.45, "c", 0.75, "b", 3, ...
                            "Cd", [1 0], "directTerm", 0, "fd", 0);
            srg    = rid.surrogate.StateSpaceFunction(@gapLocusStepFcn, ...
                         @linearOutputFcn, params, 2, 1, 0.1);
        end

        function xAchieved = locusLevels(xProbe, xNorm, iIn, nLevels, bridge)
            %LOCUSLEVELS  The arc-length tiling, bridging the gaps or not.
            %
            %   Both rules from one place, so the test that compares them cannot
            %   differ from itself in some second way. `bridge` true adds the
            %   chord across a gap -- the defect -- and false adds it only
            %   between consecutive probes, which is the specification.
            %
            %   Inputs:
            %       xProbe   (nProbes x nState)  rest states, physical units
            %       xNorm    (nProbes x nState)  the same, scaled by the box
            %       iIn      (nInside x 1)       which probes lie inside the box
            %       nLevels  (1x1)
            %       bridge   (1x1) logical
            %   Outputs:
            %       xAchieved (nLevels x nState)
            s = zeros(numel(iIn), 1);
            for ii = 2:numel(iIn)
                chord = 0;
                if bridge || iIn(ii) == iIn(ii - 1) + 1
                    chord = norm(xNorm(iIn(ii), :) - xNorm(iIn(ii - 1), :));
                end
                s(ii) = s(ii - 1) + chord;
            end

            sTarget   = linspace(s(1), s(end), nLevels).';
            xAchieved = zeros(nLevels, size(xProbe, 2));
            for iLevel = 1:nLevels
                [~, iNear]           = min(abs(s - sTarget(iLevel)));
                xAchieved(iLevel, :) = xProbe(iIn(iNear), :);
            end
        end

        function srg = flatLocusSurrogate()
            %FLATLOCUSSURROGATE  A one-state model whose rest state does not move with u.
            %
            %   `Bd = 0`: the input never reaches the state. Every input therefore
            %   rests at the same state, so the steady-state locus is a single
            %   POINT and has no arc length -- the degenerate case
            %   part 9. Stable, so the constructor accepts it and `steadyState`
            %   is exact rather than a converged fallback.
            srg = rid.surrogate.LinearStateSpace(0.5, 0, 1, 0.1);
        end

        function srg = unitGainSurrogate()
            %UNITGAINSURROGATE  A second-order unit-gain model resting at [u, 0].
            %
            %   Position and velocity, so the rest locus is exactly the line
            %   x2 = 0 that the joint-space representation contract names: tiling it in arc length reduces
            %   to tiling x1, and no rule is needed for WHICH state to tile.
            w0   = 2;
            zeta = 0.3;
            Ac   = [0, 1; -w0^2, -2 * zeta * w0];
            Bc   = [0; w0^2];
            Cc   = [1, 0];
            srg  = rid.surrogate.LinearStateSpace.fromContinuous(Ac, Bc, Cc, 0.1);
        end

        function [r, f] = stateSpaceRun(opts)
            %STATESPACERUN  An OFFLINE run over rid.rep.StateSpace, states simulated.
            arguments
                opts.nSteps     (1,1) double = 6
                opts.nStates    (1,1) double = 2
                opts.srg                     = []
                opts.stateRange (1,:) double = []
            end
            f = SetupTest.stateSpaceFixture(opts.nStates, ...
                    "epsReq0", RunTest.StateEpsReq0, "srg", opts.srg, ...
                    "stateRange", opts.stateRange);
            r = rid.rhc.run([], SetupTest.setupFrom(f), ...
                    RunTest.history("srg", f.srg), ...
                    RunTest.stateOverrides(opts.nSteps));
        end

        function [r, f, plant] = measuredStateRun(opts)
            %MEASUREDSTATERUN  An ONLINE run whose past states are MEASURED.
            %
            %   `opts.process` substitutes a process that breaks the contract; the
            %   default plant keeps it.
            arguments
                opts.nSteps  (1,1) double = 5
                opts.nStates (1,1) double = 2
                opts.defect  (1,1) string = "none"
                opts.process              = []
            end
            f = SetupTest.stateSpaceFixture(opts.nStates, ...
                    "pastSource", "measured", "epsReq0", RunTest.StateEpsReq0);
            plant = opts.process;
            if isempty(plant)
                plant = RunTest.statePlant(opts.nStates, "defect", opts.defect);
            end
            r = rid.rhc.run(plant, SetupTest.setupFrom(f), ...
                    RunTest.measuredHistory(f.srg, zeros(1, opts.nStates)), ...
                    RunTest.stateOverrides(opts.nSteps));
        end

        function r = updatingRun(varargin)
            %UPDATINGRUN  An ONLINE run whose surrogate actually refits.
            %
            %   Online, because a refit needs measurements; `pastSource` stays
            %   "simulated", which is what makes
            %   `requiresRebuildOnSurrogateUpdate` true and puts a rebuild in the
            %   same step as a reset. The interval is a CONSTRUCTION argument of
            %   the scheduler, not a run override.
            args     = struct(varargin{:});
            interval = RunTest.pick(args, "updateInterval", 1);
            args     = rmfield(args, intersect(fieldnames(args), {'updateInterval'}));
            [f, ov]  = RunTest.problem(args, "simulated");
            sch      = rid.surrogate.ImmediateScheduler( ...
                           "updateInterval", interval, "minSamples", 3);
            setup    = rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, f.fields, ...
                                     f.gen, f.sol, f.state, sch, f.fieldSch, f.qInfo);
            plant    = FirstOrderProcess(0.8, 0.1, 1.4);
            r        = rid.rhc.run(plant, setup, ...
                           RunTest.history("srg", f.srg, "online", true), ov);
        end

        function r = risingRun(assertOn, logEveryStep)
            %RISINGRUN  A run whose J rises, for the monotonicity assert.
            %
            %   The rise comes from RisingCoverageState; `logEveryStep` swaps in a
            %   field scheduler that fires every step, so the same rise happens
            %   across a logged discontinuity instead of within a segment.
            f     = SetupTest.fixture("nSup", 8, "epsReq0", RunTest.EpsReq0);
            state = RisingCoverageState(f.Q, f.epsReq, 0.1);
            fieldSch = f.fieldSch;
            if logEveryStep
                fieldSch = EagerFieldScheduler();
            end
            setup = rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, f.fields, ...
                                  f.gen, f.sol, state, f.sch, fieldSch, f.qInfo);
            ov = RunTest.overrides("nSteps", 4);
            ov.rhc.assertMonotonicity = assertOn;
            r = rid.rhc.run([], setup, RunTest.history("srg", f.srg), ov);
        end

        function [f, ov] = problem(args, pastSource)
            %PROBLEM  One fixture and one override set, from a test's own keywords.
            mode = "budget";
            if isfield(args, "mode")
                mode = args.mode;
            end
            fixtureArgs = {"pastSource", pastSource, "mode", mode, ...
                           "epsReq0", RunTest.epsReq0Of(args)};
            for name = ["nHorizon", "nSup"]
                if isfield(args, name)
                    fixtureArgs = [fixtureArgs, {name, args.(name)}]; %#ok<AGROW>
                end
            end
            f = SetupTest.fixture(fixtureArgs{:});

            ov      = struct();
            ov.crit = struct("epsReq0", RunTest.epsReq0Of(args), "mode", mode);
            ov.rhc  = struct();
            ov.io   = struct("saveRun", false);
            for name = string(fieldnames(args)).'
                switch name
                    case {"nHorizon", "nSup", "epsReq0", "pastSource", "mode", ...
                          "online", "u0", "y0"}
                        % consumed by the fixture or by RunTest.history
                    otherwise
                        ov.rhc.(name) = args.(name);
                end
            end
        end

        function value = epsReq0Of(args)
            value = RunTest.pick(args, "epsReq0", RunTest.EpsReq0);
        end

    end

end
