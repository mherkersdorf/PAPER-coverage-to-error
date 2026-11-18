classdef ProgressTest < matlab.unittest.TestCase
%PROGRESSTEST  rid.viz.progress: one mark per event, and the segments J lives in.
%
%   THE TWO CLAIMS THIS PLOT MAKES. That every event in the log is marked --
%   all seven kinds, not the one kind that happens to be interesting -- and
%   that the trace is broken exactly where `rid.rhc.run` breaks its own
%   monotonicity baseline. Both are countable, which is what makes them
%   assertable when the picture is not.
%
%   THE FIXTURE CARRIES A REJECTED SURROGATE UPDATE. That is the whole reason
%   `busy` runs a scheduler that refits: a `"surrogateUpdate"` is logged
%   whenever the scheduler had something to say, and only `"applied"` replaced
%   the model. A segment rule keyed on the KIND rather than on `detail.status`
%   would invent a discontinuity the run does not have, and a fixture with no
%   rejected fit could not tell the two rules apart.

    properties (Access = private)
        Quiet   % a run with no events at all
        Busy    % the coverage reset's degenerate regime: resets, refits, one rejected fit
        Refit   % the surrogate moves; coverage never resets
    end

    methods (TestClassSetup)

        function buildRuns(tc)
            tc.Quiet = VizFixture.quiet(3);
            tc.Busy  = VizFixture.busy();
            tc.Refit = VizFixture.refitting();
        end

    end

    methods (Test)

        function aRunWithNoEventsIsOneUnbrokenSegment(tc)
            % Fails if: a boundary is invented where the loop never reset its
            % baseline, which would make a perfectly comparable trace look
            % discontinuous.
            tc.assertEmpty(tc.Quiet.eventLog, "fixture: no events");

            [fig, info] = rid.viz.progress(tc.Quiet, "visible", "off");
            tc.addTeardown(@close, fig);

            tc.verifyEqual(info.nSegments, 1);
            tc.verifyEqual(info.nMarks,    0);
            tc.verifyEqual(info.segments(1).first, 1);
            tc.verifyEqual(info.segments(1).last,  info.nSteps);
            tc.verifyEqual(numel(findobj(fig, "Tag", "segment")), 1, ...
                "one drawn segment line");
        end

        function everyEventIsMarkedExactlyOnce(tc)
            % Fails if: marks are drawn only for the kinds that break a segment,
            % or only for surrogateUpdate. A jump in J with no mark beside it is
            % the plot nobody can read six months later.
            [fig, info] = rid.viz.progress(tc.Busy, "visible", "off");
            tc.addTeardown(@close, fig);

            nEvents = numel(tc.Busy.eventLog);
            tc.assertGreaterThan(nEvents, 0, "fixture: the run must log events");
            tc.verifyEqual(info.nMarks, nEvents);

            drawn = 0;
            for kind = info.kinds
                drawn = drawn + numel(findobj(fig, "Tag", char(kind)));
            end
            tc.verifyEqual(drawn, nEvents, "one drawn mark per logged event");
        end

        function everyKindInTheClosedSetHasItsOwnMark(tc)
            % Fails if: a kind is dropped from the styling, which is how a
            % newly-emitted kind would become invisible.
            [fig, info] = rid.viz.progress(tc.Busy, "visible", "off");
            tc.addTeardown(@close, fig);
            tc.verifyEqual(info.kinds, ["surrogateUpdate", "resolutionChange", ...
                "coverageReset", "supportChange", "normalizationChange", ...
                "modeChange", "seed"], ...
                "the closed set of the method description, in its order");
        end

        function theSegmentsBreakWhereTheLoopBreaksItsBaseline(tc)
            % Fails if: the boundary set is keyed on the kind rather than on
            % detail.status, or is placed before the step that logged the event
            % rather than after it.
            % ON `Refit`, NOT ON `Busy`. In `busy` a coverage reset fires at
            % nearly every step, so every segment is one step long whatever the
            % surrogate rule says and the two candidate rules agree by
            % construction. Here the only boundaries are
            % applied surrogate updates.
            tc.assertEmpty(ProgressTest.kindsOf(tc.Refit, "coverageReset"), ...
                "fixture: no reset may fire, or the segment rule is untestable");

            [fig, info] = rid.viz.progress(tc.Refit, "visible", "off");
            tc.addTeardown(@close, fig);

            expected = ProgressTest.expectedSegments(tc.Refit);
            tc.verifyEqual([info.segments.first], [expected.first]);
            tc.verifyEqual([info.segments.last],  [expected.last]);
            tc.verifyEqual(numel(findobj(fig, "Tag", "segment")), info.nSegments);
        end

        function theSegmentsPartitionEveryStepExactlyOnce(tc)
            % Fails if: a step is dropped from the plot, or drawn twice.
            for r = {tc.Quiet, tc.Busy}
                [fig, info] = rid.viz.progress(r{1}, "visible", "off");
                tc.addTeardown(@close, fig);
                covered = [];
                for iSeg = 1:info.nSegments
                    covered = [covered, info.segments(iSeg).first:info.segments(iSeg).last]; %#ok<AGROW>
                end
                tc.verifyEqual(sort(covered), 1:info.nSteps);
            end
        end

        function aRejectedSurrogateUpdateIsMarkedButBreaksNothing(tc)
            % THE SHARP ONE. The model never moved, so the values on either side
            % are still comparable. Marking it and breaking on it are different
            % decisions, and only the first is right.
            marks = ProgressTest.marksOf(tc, tc.Refit);
            rejected = marks(arrayfun(@(m) m.kind == "surrogateUpdate" && ...
                                            m.status == "rejected", marks));
            tc.assertNotEmpty(rejected, ...
                "fixture: the run must carry a rejected surrogate update, or " + ...
                "a kind-keyed rule and a status-keyed rule agree here");
            tc.verifyFalse(any([rejected.breaks]), ...
                "a rejected fit left the model in place: nothing to break");

            applied = marks(arrayfun(@(m) m.kind == "surrogateUpdate" && ...
                                           m.status == "applied", marks));
            tc.assertNotEmpty(applied, "fixture: an applied update too");
            tc.verifyTrue(all([applied.breaks]));
        end

        function aConsecutiveResetIsDrawnAsTheWarningItIs(tc)
            % Fails if: severity is dropped. Consecutive resets are a statement
            % about the CONFIGURATION, not a routine event, and they are
            % the one thing on this axes a reader must not skim past.
            marks    = ProgressTest.marksOf(tc, tc.Busy);
            warnings = marks(arrayfun(@(m) m.severity == "warning", marks));
            tc.assertNotEmpty(warnings, ...
                "fixture: the run must enter consecutive resets");
            tc.verifyTrue(all(arrayfun(@(m) m.kind == "coverageReset", warnings)));
        end

        function aRecordDrawsTheSamePlotAsItsResult(tc)
            [figResult, fromResult] = rid.viz.progress(tc.Busy,        "visible", "off");
            [figRecord, fromRecord] = rid.viz.progress(tc.Busy.record, "visible", "off");
            tc.addTeardown(@close, figResult);
            tc.addTeardown(@close, figRecord);

            tc.verifyEqual(fromRecord.nSegments, fromResult.nSegments);
            tc.verifyEqual(fromRecord.nMarks,    fromResult.nMarks);
            tc.verifyEqual(fromRecord.runId,     fromResult.runId);
        end

        function aRecordFromDiskDrawsItToo(tc)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            root = tc.applyFixture(TemporaryFolderFixture()).Folder;

            path = rid.io.saveRun(tc.Busy.record, "", "root", string(root));
            rec  = rid.io.loadRun(path);

            [figDisk, onDisk] = rid.viz.progress(rec, "visible", "off");
            tc.addTeardown(@close, figDisk);
            [figLive, live]   = rid.viz.progress(tc.Busy, "visible", "off");
            tc.addTeardown(@close, figLive);

            tc.verifyEqual(onDisk.nSegments, live.nSegments);
            tc.verifyEqual(onDisk.nMarks,    live.nMarks);
        end

        function everyDimensionIsDrawn(tc)
            % The plot does not depend on nDim, but the seam it reads through
            % does, so it is exercised at each one all the same.
            for d = [2 3 5]
                [fig, info] = rid.viz.progress(VizFixture.quiet(d), "visible", "off");
                tc.addTeardown(@close, fig);
                tc.verifyGreaterThan(info.nSteps, 0, sprintf("nDim = %d", d));
            end
        end

        function aSourceThatIsNeitherIsRefused(tc)
            tc.verifyError(@() rid.viz.progress(struct("u", 1), "visible", "off"), ...
                "rid:viz:unknownSource");
        end

        function theFigureIsInvisibleOnlyWhenAsked(tc)
            hidden = rid.viz.progress(tc.Quiet, "visible", "off");
            tc.addTeardown(@close, hidden);
            tc.verifyEqual(hidden.Visible, matlab.lang.OnOffSwitchState("off"));

            shown = rid.viz.progress(tc.Quiet);
            tc.addTeardown(@close, shown);
            tc.verifyEqual(shown.Visible, matlab.lang.OnOffSwitchState("on"));
        end

    end

    methods (Static, Access = private)

        function seg = expectedSegments(r)
            %EXPECTEDSEGMENTS  The boundaries, worked out from the log by hand.
            %
            %   Written out here rather than read from the code under test: the
            %   loop stores Jtrace(k) BEFORE polling the schedulers, so an event
            %   at step k separates Jtrace(k) from Jtrace(k+1).
            nSteps  = numel(r.Jtrace);
            breaks  = false(nSteps, 1);
            for iEvent = 1:numel(r.eventLog)
                e = r.eventLog(iEvent);
                switch string(e.kind)
                    case "coverageReset"
                        doesBreak = true;
                    case {"surrogateUpdate", "resolutionChange"}
                        doesBreak = isfield(e.detail, "status") && ...
                                    string(e.detail.status) == "applied";
                    otherwise
                        doesBreak = false;
                end
                if doesBreak && e.step >= 1 && e.step <= nSteps
                    breaks(e.step) = true;
                end
            end

            seg   = struct("first", {}, "last", {});
            first = 1;
            for kStep = 1:nSteps
                if breaks(kStep) || kStep == nSteps
                    seg(end + 1) = struct("first", first, "last", kStep); %#ok<AGROW>
                    first = kStep + 1;
                end
            end
        end

        function events = kindsOf(r, kind)
            %KINDSOF  Every logged event of one kind.
            keep   = arrayfun(@(e) string(e.kind) == kind, r.eventLog);
            events = r.eventLog(keep);
        end

        function marks = marksOf(tc, r)
            [fig, info] = rid.viz.progress(r, "visible", "off");
            tc.addTeardown(@close, fig);
            marks = info.marks;
        end

    end

end
