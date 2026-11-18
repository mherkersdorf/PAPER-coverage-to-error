classdef CoverageTest < matlab.unittest.TestCase
%COVERAGETEST  rid.viz.coverage: the mapping, since the picture cannot be asserted.
%
%   WHAT IS ASSERTED. Not that the figure looks right -- the figure is accepted when the
%   picture is recognizable as the behaviour from the earlier papers, and that
%   is not a test's judgement. What is asserted is the MAPPING: that every
%   supporting point is drawn exactly once, that each of the three states holds
%   the points it should, that a higher-dimensional space is projected rather
%   than truncated, and that a record which cannot supply Q is refused instead
%   of reconstructed.
%
%   THE STATE COUNTS ARE CHECKED AGAINST THE ORACLE, not against the cache the
%   figure read. `rid.crit.evaluate` is the permanent independent recomputation,
%   so recomputing dTilde from the run's own realized points and
%   classifying it here shares no code with the path the figure takes --
%   `rid.crit.coverageReport` reads the cache. Reading the cache twice would
%   agree with itself whatever it held.

    properties (TestParameter)
        nDim = {2, 3, 5}
    end

    properties (Access = private)
        Quiet   % struct: nDim -> a run result with Q stored
    end

    methods (TestClassSetup)

        function buildRuns(tc)
            tc.Quiet = struct();
            for d = [2 3 5]
                tc.Quiet.("d" + d) = VizFixture.quiet(d);
            end
        end

    end

    methods (Test)

        function everySupportingPointIsDrawnExactlyOnceInEveryPanel(tc, nDim)
            % Fails if: a state series is dropped, double-drawn, or a panel is
            % filled from a different selection than the other panels.
            r   = tc.runOf(nDim);
            fig = tc.figureOf(r);

            drawn = tc.perPanelCounts(fig);
            tc.verifyEqual(sum(drawn, 2), repmat(r.record.problem.nSup, size(drawn, 1), 1), ...
                "every panel must draw nSup supporting points, in total");
        end

        function theThreeStateCountsMatchTheOracle(tc, nDim)
            % Fails if: the saturated state is folded into the active one, the
            % boundary is put on the wrong side, or the figure reads a cache
            % that is not the run's honest final one.
            r        = tc.runOf(nDim);
            expected = CoverageTest.oracleCounts(r);

            [fig, info] = rid.viz.coverage(r, "visible", "off");
            tc.addTeardown(@close, fig);

            tc.verifyEqual(info.label, ["resolved", "active", "saturated"]);
            tc.verifyEqual(info.count, expected, ...
                "resolved / active / saturated, against rid.crit.evaluate");
            tc.verifyEqual(tc.perPanelCounts(fig), ...
                repmat(expected, tc.nPanels(nDim), 1), ...
                "the drawn series must hold what info reports");
        end

        function aFixtureInWhichAllThreeStatesAreOccupied(tc, nDim)
            % NOT A TEST OF THE CODE -- a test of the fixture.
            % With no saturated point, a two-colour resolved/unresolved plot
            % agrees with a three-state one by construction, and the counts
            % above cannot fail. This is the assertion that keeps the fixture
            % one in which the defect is visible, at EVERY dimension.
            counts = CoverageTest.oracleCounts(tc.runOf(nDim));
            tc.verifyEqual(nnz(counts), 3, ...
                "fixture: resolved, active AND saturated must all be occupied, " + ...
                "or the three-way split is untestable at this dimension");
        end

        function theRealizedSamplesAreDrawnOverTheSupportingPoints(tc, nDim)
            % Fails if: the realized points are omitted, or drawn beneath the
            % supporting points, where they are invisible.
            r   = tc.runOf(nDim);
            fig = tc.figureOf(r);

            realized = findobj(fig, "Tag", "realized");
            tc.assertNotEmpty(realized, "the realized samples must be drawn");
            tc.verifyEqual(numel(realized), tc.nPanels(nDim));
            tc.verifyGreaterThan(numel(realized(1).XData), 0);

            ax       = ancestor(realized(1), "axes");
            children = ax.Children;                     % last drawn is first
            tc.verifyEqual(children(1).Tag, 'realized', ...
                "the realized samples must be the topmost series");
        end

        function aHigherDimensionalSpaceIsProjectedAndSaidToBe(tc)
            % Fails if: nDim > 3 silently plots the first two coordinates, or a
            % panel is drawn without saying it is a shadow of the full space.
            r           = tc.runOf(5);
            [fig, info] = rid.viz.coverage(r, "visible", "off");
            tc.addTeardown(@close, fig);

            tc.verifyTrue(info.isProjection);
            tc.verifyEqual(info.nAxes, 10, "all 10 pairs of a 5-D space");
            tc.verifyEqual(sort(info.pairs, 1), sort(nchoosek(1:5, 2), 1));

            titles = string(arrayfun(@(a) string(a.Title.String), ...
                findobj(fig, "Type", "axes"), "UniformOutput", false));
            tc.verifyTrue(all(startsWith(titles, "projection")), ...
                "every panel must be labelled as a projection");

            layout = findobj(fig, "Type", "tiledlayout");
            tc.verifyTrue(contains(string(layout.Subtitle.String), "PROJECTION"), ...
                "the reader must be told resolution is a property of the FULL space");
        end

        function twoAndThreeDimensionsAreDrawnWhole(tc)
            % Fails if: a whole space is projected anyway, which would hide the
            % one case where the picture IS the truth.
            for d = [2 3]
                [fig, info] = rid.viz.coverage(tc.runOf(d), "visible", "off");
                tc.addTeardown(@close, fig);
                tc.verifyFalse(info.isProjection, sprintf("nDim = %d", d));
                tc.verifyEqual(info.nAxes, 1, sprintf("nDim = %d", d));
                tc.verifyEqual(size(info.pairs, 2), d, sprintf("nDim = %d", d));
            end
        end

        function theAxesAreFramedOnTheRoiAndNotOnTheData(tc, nDim)
            % Fails if: the axes are cropped to the points drawn. Coverage is a
            % statement about a REGION, so a design that reached a third of the
            % space would otherwise look exactly like one that filled it -- the
            % plot would rescale away the very thing it is drawn to show.
            r   = tc.runOf(nDim);
            fig = tc.figureOf(r);
            box = r.record.setup.setup.nrm.normalizedBox();

            [fig2, info] = rid.viz.coverage(r, "visible", "off");
            tc.addTeardown(@close, fig2);

            axList = flipud(findobj(fig, "Type", "axes"));
            for iAx = 1:numel(axList)
                pair = info.pairs(iAx, :);
                tc.verifyLessThanOrEqual(axList(iAx).XLim(1), box(pair(1), 1));
                tc.verifyGreaterThanOrEqual(axList(iAx).XLim(2), box(pair(1), 2));
                tc.verifyLessThanOrEqual(axList(iAx).YLim(1), box(pair(2), 1));
                tc.verifyGreaterThanOrEqual(axList(iAx).YLim(2), box(pair(2), 2));
            end
        end

        function theFrameIsWiderThanTheDesignActuallyReached(tc)
            % The fixture half of the claim: a run that happened
            % to fill its ROI would frame identically either way, so the limits
            % must be checked against a design that did NOT.
            r   = tc.runOf(3);
            fig = tc.figureOf(r);
            box = r.record.setup.setup.nrm.normalizedBox();

            realized = findobj(fig, "Tag", "realized");
            tc.assertNotEmpty(realized);
            reach = max(realized(1).XData) - min(realized(1).XData);
            tc.assertLessThan(reach, 0.9 * (box(1, 2) - box(1, 1)), ...
                "fixture: the design must leave part of the ROI unreached");

            ax = findobj(fig, "Type", "axes");
            tc.verifyGreaterThan(ax(1).XLim(2) - ax(1).XLim(1), reach);
        end

        function aRefittedRunStillDrawsItsRealizedSamples(tc)
            % THE POINT OF THE INITIAL-VALUE RULE. Z is derived and never stored, so the
            % realized points are rebuilt through rep.build -- and after a refit
            % the loop remaps every past coordinate under the NEW model
            % (invariant 3). With only the initial surrogate in the record, the
            % positions that existed at the end corresponded to a model the
            % record did not contain, so this layer had to be withheld.
            %
            % Fails if: the final surrogate is not stored, or the rebuild
            % silently uses the initial one -- which would draw points the run
            % never visited, the worse failure because it is a plausible one.
            r = VizFixture.refitting();
            tc.assertTrue(CoverageTest.surrogateMoved(r), ...
                "fixture: the surrogate must actually be refitted, or the " + ...
                "initial and final models are the same object and this " + ...
                "cannot be tested at all");
            tc.assertTrue(r.record.setup.setup.rep.requiresRebuildOnSurrogateUpdate(), ...
                "fixture: and the coordinates must follow the surrogate, or " + ...
                "which model is used makes no difference");

            [fig, info] = rid.viz.coverage(r.record, "visible", "off");
            tc.addTeardown(@close, fig);

            tc.verifyGreaterThan(info.nRealized, 0, ...
                "the realized layer must now be drawn");
            tc.verifyNotEmpty(findobj(fig, "Tag", "realized"));
        end

        function theRebuildUsesTheFinalSurrogateAndNotTheInitialOne(tc)
            % THE SHARP ANCHOR. That points are drawn does not show WHICH model
            % produced them; drawing them under the initial surrogate would also
            % be non-empty. The two coordinate sets must actually differ, and the
            % figure must match the final one.
            r     = VizFixture.refitting();
            setup = r.record.setup.setup;

            [Zinitial, ~] = setup.rep.build(r.hist, setup.srg);
            [Zfinal,   ~] = setup.rep.build(r.hist, r.record.outputs.srg);
            tc.assertFalse(isequal(Zinitial, Zfinal), ...
                "fixture: the two models must give different coordinates, or " + ...
                "reading either would pass");

            [fig, info] = rid.viz.coverage(r.record, "visible", "off");
            tc.addTeardown(@close, fig);
            drawn = findobj(fig, "Tag", "realized");

            expected = setup.nrm.toNormalized(Zfinal);
            tc.verifyEqual(info.nRealized, size(expected, 1));
            tc.verifyEqual(drawn(1).XData(:), expected(:, 1), "RelTol", 1e-12, ...
                "the drawn layer must be the FINAL model's coordinates");
        end

        function aMeasuredStateRunsRealizedPointsAreRebuiltFromItsStates(tc)
            % THE VIEW MUST CARRY `hist.x`. `realizedPoints` calls
            % `rep.build(v.hist, v.srg)`, and a "measured" state-space
            % representation reads its past states from the history rather than
            % simulating them -- so a view that dropped the field would not draw
            % the wrong points, it would raise `rid:rep:missingStateMeasurement`
            % and the run could not be plotted at all.
            %
            % Fails if: viewOf builds its history from `u`, `y`, `x0` and `tS`
            % alone, as it did before the measured states existed.
            rec = ReplayTest.measuredStateRecord("nSteps", 4);
            rep = rec.setup.setup.rep;
            tc.assertTrue(rep.needsMeasuredStates(), ...
                "fixture: the representation must READ the measured states, " + ...
                "or simulating them would pass just as well");

            [fig, info] = rid.viz.coverage(rec, "visible", "off");
            tc.addTeardown(@close, fig);

            % Against an independent rebuild, not against what the figure read.
            [Zphys, ~] = rep.build(struct("u",  rec.outputs.u, ...
                                          "y",  rec.outputs.y, ...
                                          "x",  rec.outputs.x, ...
                                          "x0", rec.outputs.x0, ...
                                          "tS", rec.outputs.tS), rec.outputs.srg);
            expected = rec.setup.setup.nrm.toNormalized(Zphys);

            drawn = findobj(fig, "Tag", "realized");
            tc.assertNotEmpty(drawn);
            tc.verifyEqual(info.nRealized, size(expected, 1));
            tc.verifyEqual(drawn(1).XData(:), expected(:, 1), "RelTol", 1e-12);
        end

        function theStoredFinalSurrogateSurvivesTheDisk(tc)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            root = tc.applyFixture(TemporaryFolderFixture()).Folder;

            r    = VizFixture.refitting();
            path = rid.io.saveRun(r.record, "", "root", string(root));
            rec  = rid.io.loadRun(path);

            tc.verifyEqual(rec.outputs.srg.describe(), ...
                           r.record.outputs.srg.describe());
            [figDisk, onDisk] = rid.viz.coverage(rec, "visible", "off");
            tc.addTeardown(@close, figDisk);
            tc.verifyGreaterThan(onDisk.nRealized, 0);
        end

        function aRecordDrawsTheSameFigureAsItsResult(tc, nDim)
            % Fails if: the record path reads a different field layout, or the
            % result path reads something the record does not carry.
            r = tc.runOf(nDim);

            [figResult, fromResult] = rid.viz.coverage(r,        "visible", "off");
            [figRecord, fromRecord] = rid.viz.coverage(r.record, "visible", "off");
            tc.addTeardown(@close, figResult);
            tc.addTeardown(@close, figRecord);

            tc.verifyEqual(fromRecord.count,     fromResult.count);
            tc.verifyEqual(fromRecord.nRealized, fromResult.nRealized);
            tc.verifyEqual(fromRecord.runId,     fromResult.runId);
        end

        function aRecordFromDiskDrawsTheSameFigureToo(tc)
            % Fails if: anything the figure needs does not survive saveRun ->
            % loadRun. This is the path +io exists for; an in-memory record
            % would not exercise it.
            import matlab.unittest.fixtures.TemporaryFolderFixture
            root = tc.applyFixture(TemporaryFolderFixture()).Folder;

            r    = tc.runOf(3);
            path = rid.io.saveRun(r.record, "", "root", string(root));
            rec  = rid.io.loadRun(path);

            [figLive, live]   = rid.viz.coverage(r,   "visible", "off");
            [figDisk, onDisk] = rid.viz.coverage(rec, "visible", "off");
            tc.addTeardown(@close, figLive);
            tc.addTeardown(@close, figDisk);

            tc.verifyEqual(onDisk.count,     live.count);
            tc.verifyEqual(onDisk.nRealized, live.nRealized);
        end

        function aRecordWithoutQIsRefusedAndSaysWhat(tc)
            % Fails if: Q is reconstructed from its seed. That would draw a
            % picture of points the run never used, which is what section 7
            % refuses.
            %
            % THE BACK DOOR IS GONE. This test used to assert the
            % opposite -- that `setup.setup.Q` held the points whatever `storeQ`
            % said -- because "the setup is stored in full" kept two further
            % copies and the refusal here had to be shown to ignore them. The
            % record now contains Q nowhere, so what is asserted is that, and
            % the refusal on top of it.
            r = VizFixture.quiet(2, struct("storeQ", false));
            tc.assertFalse(logical(r.record.problem.storedQ), ...
                "fixture: the record must genuinely lack Q");
            tc.assertEmpty(r.record.setup.setup.Q, ...
                "the stored Setup carries no Q either");
            tc.assertEmpty(r.record.setup.setup.state.Q, ...
                "and neither does its coverage state");

            tc.verifyError(@() rid.viz.coverage(r,        "visible", "off"), ...
                "rid:viz:recordWithoutQ");
            tc.verifyError(@() rid.viz.coverage(r.record, "visible", "off"), ...
                "rid:viz:recordWithoutQ");
        end

        function aSourceThatIsNeitherIsRefused(tc)
            tc.verifyError(@() rid.viz.coverage(struct("u", 1), "visible", "off"), ...
                "rid:viz:unknownSource");
        end

        function theFigureIsInvisibleOnlyWhenAsked(tc)
            % Fails if: `visible` is ignored, which is what would leave a test
            % run littered with windows -- and, on a headless batch, failing.
            hidden = rid.viz.coverage(tc.runOf(2), "visible", "off");
            tc.addTeardown(@close, hidden);
            tc.verifyEqual(hidden.Visible, matlab.lang.OnOffSwitchState("off"));

            shown = rid.viz.coverage(tc.runOf(2));
            tc.addTeardown(@close, shown);
            tc.verifyEqual(shown.Visible, matlab.lang.OnOffSwitchState("on"), ...
                "a plotting function shows its figure unless told not to");
        end

    end

    methods (Access = private)

        function r = runOf(tc, nDim)
            r = tc.Quiet.("d" + nDim);
        end

        function fig = figureOf(tc, r)
            fig = rid.viz.coverage(r, "visible", "off");
            tc.addTeardown(@close, fig);
        end

        function n = nPanels(~, nDim)
            n = 1;
            if nDim > 3
                n = nchoosek(nDim, 2);
            end
        end

        function counts = perPanelCounts(~, fig)
            %PERPANELCOUNTS  (nPanels x 3) points drawn per state, per panel.
            axList = flipud(findobj(fig, "Type", "axes"));
            labels = ["resolved", "active", "saturated"];
            counts = zeros(numel(axList), 3);
            for iAx = 1:numel(axList)
                for iState = 1:3
                    h = findobj(axList(iAx), "Tag", char(labels(iState)));
                    counts(iAx, iState) = numel(h.XData);
                end
            end
        end

    end

    methods (Static, Access = private)

        function tf = surrogateMoved(r)
            %SURROGATEMOVED  Did a logged surrogate update actually replace the model?
            tf = false;
            for iEvent = 1:numel(r.eventLog)
                e = r.eventLog(iEvent);
                if e.kind == "surrogateUpdate" && isfield(e.detail, "status") ...
                        && string(e.detail.status) == "applied"
                    tf = true;
                    return
                end
            end
        end

        function counts = oracleCounts(r)
            %ORACLECOUNTS  resolved / active / saturated, from rid.crit.evaluate.
            %
            %   Independent of the cache the figure reads: `evaluate` shares no
            %   code with `CoverageState` on purpose, so this is
            %   a genuine second opinion rather than the same number twice.
            setup      = r.record.setup.setup;
            [Zphys, ~] = setup.rep.build(r.hist, setup.srg);
            Znorm      = setup.nrm.toNormalized(Zphys);

            [~, info] = rid.crit.evaluate(Znorm, setup.Q, setup.state.epsReq, ...
                "psiMax", setup.state.psiMax, "pNorm", setup.state.pNorm);

            resolved  = info.dTilde <= 1;
            saturated = info.dTilde >= 1 + setup.state.psiMax;
            counts    = [sum(resolved), sum(~resolved & ~saturated), sum(saturated)];
        end

    end

end
