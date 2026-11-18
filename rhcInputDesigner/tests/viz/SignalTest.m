classdef SignalTest < matlab.unittest.TestCase
%SIGNALTEST  rid.viz.signal: the trajectory, its output, and its bounds.
%
%   WHAT MAKES THIS PLOT WORTH DRAWING is that a saturated run and a stalled
%   run look identical without the bounds on the axes, and they call for
%   opposite responses. So the assertions are about the bounds and about the
%   saturation count, not about the line: that the box comes from the
%   generator's own constraints, that the marked samples are exactly the ones
%   sitting on it, and that samples belonging to the initial condition are not
%   counted as the design's doing.

    properties (Access = private)
        Offline   % no process: no measured output
        Online    % a plant: y is present
    end

    methods (TestClassSetup)

        function buildRuns(tc)
            tc.Offline = VizFixture.quiet(3);
            tc.Online  = VizFixture.refitting();
        end

    end

    methods (Test)

        function theInputIsDrawnAsHeldSamplesWithItsBox(tc)
            % Fails if: the bounds are omitted, so a run pinned at its limit is
            % indistinguishable from one that stopped moving.
            [fig, info] = rid.viz.signal(tc.Offline, "visible", "off");
            tc.addTeardown(@close, fig);

            tc.verifyEqual(numel(findobj(fig, "Tag", "u")), info.nInputs, ...
                "one held-sample trace per input channel");
            tc.verifyEqual(numel(findobj(fig, "Tag", "uBound")), 2 * info.nInputs, ...
                "an upper and a lower amplitude bound per channel");
            tc.verifyEqual(numel(findobj(fig, "Tag", "duBound")), 2 * info.nInputs, ...
                "and the rate box too: a rate-limited signal is flat for a " + ...
                "reason the amplitude axes does not show");
        end

        function theBoxIsTheGeneratorsOwnAndNotAnInvention(tc)
            % Fails if: the drawn bounds come from the ROI, from the
            % normalization, or from anywhere but the constraints the run
            % actually optimized under.
            fig = tc.figureOf(tc.Offline);
            con = tc.Offline.record.setup.setup.gen.constraints;

            drawn = sort(arrayfun(@(h) h.Value, findobj(fig, "Tag", "uBound")));
            tc.verifyEqual(drawn, sort(con.amplitude(:)), "AbsTol", 0, ...
                "the amplitude bounds, verbatim");

            drawnRate = sort(arrayfun(@(h) h.Value, findobj(fig, "Tag", "duBound")));
            tc.verifyEqual(drawnRate, sort(con.rate(:)), "AbsTol", 0);
        end

        function aSaturatedSampleIsMarkedAndCounted(tc)
            % THE SHARP ONE, and it needs a run that actually saturates. The
            % marks and the count must agree, and both must be the samples on
            % the bound -- not all of them, and not none.
            r = VizFixture.saturating();
            [fig, info] = rid.viz.signal(r, "visible", "off");
            tc.addTeardown(@close, fig);

            expected = SignalTest.countOnAmplitudeBound(r, info.nPriming);
            tc.assertGreaterThan(sum(expected), 0, ...
                "fixture: the run must reach its amplitude bound, or " + ...
                "marking and not marking agree here");
            tc.assertLessThan(sum(expected), info.nDesigned * info.nInputs, ...
                "fixture: and must not sit on it at every sample, or " + ...
                "marking everything would pass too");

            tc.verifyEqual(info.nAtAmplitude, expected);
            tc.verifyEqual(numel(findobj(fig, "Tag", "uSaturated")), ...
                nnz(expected), "one marked series per channel that saturates");
        end

        function theInitialConditionIsNotCountedAsTheDesignsDoing(tc)
            % Fails if: the priming samples are folded into the designed ones.
            % `u(1,:)` is the pre-run actuator position and the rows
            % after it are data the experiment started with; neither is a
            % decision the loop made.
            r = tc.Offline;
            [fig, info] = rid.viz.signal(r, "visible", "off");
            tc.addTeardown(@close, fig);

            committed = 0;
            for kStep = 1:numel(r.steps)
                committed = committed + size(r.steps(kStep).uApplied, 1);
            end
            tc.verifyEqual(info.nDesigned, committed, ...
                "only the APPLIED UNIT of each step is committed");
            tc.verifyEqual(info.nPriming, size(r.hist.u, 1) - committed);
            tc.verifyGreaterThan(info.nPriming, 0, ...
                "fixture: there must be an initial condition to exclude");
            tc.verifyNotEmpty(findobj(fig, "Tag", "designStart"), ...
                "the boundary must be visible, not merely accounted for");
        end

        function anInitialSampleOnTheBoundIsNotCountedAsSaturation(tc)
            % THE SHARP ONE FOR THE SPLIT. A sample the experiment STARTED at
            % is not the design saturating; only what the loop committed is
            %. `saturating` cannot show this to be wrong -- its
            % initial condition never touches a bound, so counting those rows
            % and not counting them give the same number.
            r = VizFixture.saturatingWithPrimedBound();
            [fig, info] = rid.viz.signal(r, "visible", "off");
            tc.addTeardown(@close, fig);

            con   = r.record.setup.setup.gen.constraints;
            u     = r.hist.u;
            width = con.amplitude(1, 2) - con.amplitude(1, 1);
            onBound = abs(u - con.amplitude(1, 1)) <= 1e-6 * width | ...
                      abs(u - con.amplitude(1, 2)) <= 1e-6 * width;

            tc.assertGreaterThan(sum(onBound(1:info.nPriming)), 0, ...
                "fixture: the initial condition must reach a bound, or " + ...
                "excluding it and not excluding it agree here");
            tc.verifyEqual(info.nAtAmplitude, sum(onBound(info.nPriming + 1:end)), ...
                "only the samples the loop committed");
            tc.verifyLessThan(info.nAtAmplitude, sum(onBound), ...
                "and strictly fewer than every sample on the bound");
        end

        function aMeasuredOutputIsDrawnWhenThereIsOne(tc)
            % Fails if: y is invented for an offline run, or dropped for an
            % online one. An absent y is how a history says it has none.
            [figOn, online] = rid.viz.signal(tc.Online, "visible", "off");
            tc.addTeardown(@close, figOn);
            tc.assertTrue(isfield(tc.Online.hist, "y"), "fixture: online");
            tc.verifyTrue(online.hasOutput);
            tc.verifyEqual(online.nAxes, 3);
            tc.verifyNotEmpty(findobj(figOn, "Tag", "y"));

            [figOff, offline] = rid.viz.signal(tc.Offline, "visible", "off");
            tc.addTeardown(@close, figOff);
            tc.verifyFalse(offline.hasOutput);
            tc.verifyEqual(offline.nAxes, 2);
            tc.verifyEmpty(findobj(figOff, "Tag", "y"));
            tc.verifyTrue(any(contains(offline.missing, "no measured output")), ...
                "the absence must be stated, not silent");
        end

        function theBoundsComeFromTheContractAndAreNeverOptional(tc)
            % The constraints contract moved `constraints` onto the rid.gen.Generator contract, so
            % this plot no longer asks whether a box exists. A generator that
            % cannot state its actuator limits fails at CONSTRUCTION -- see
            % GeneratorContractTest -- and never reaches a figure, which is why
            % the duck-typed stand-in this test used to build is gone.
            %
            % Fails if: the presence check comes back, in which case a run pinned
            % against its limits could again be drawn as one that never
            % approached them.
            for r = {tc.Offline, tc.Online, VizFixture.saturating()}
                [fig, info] = rid.viz.signal(r{1}, "visible", "off");
                tc.addTeardown(@close, fig);
                tc.verifyEqual(numel(findobj(fig, "Tag", "uBound")), ...
                    2 * info.nInputs, "bounds are unconditional");
                tc.verifyEmpty(info.missing(contains(info.missing, "bounds")), ...
                    "there is no longer any reason to report bounds missing");
            end
        end

        function everyDimensionIsDrawn(tc)
            for d = [2 3 5]
                [fig, info] = rid.viz.signal(VizFixture.quiet(d), "visible", "off");
                tc.addTeardown(@close, fig);
                tc.verifyGreaterThan(info.nDesigned, 0, sprintf("nDim = %d", d));
            end
        end

        function everyPanelSharesOneTimeAxis(tc)
            % Fails if: the panels are left to scale themselves. u starts at the
            % experiment's first sample and du only at its second, so without
            % linking, the same horizontal position means two different instants
            % and an input excursion appears to line up with the wrong
            % measurement.
            [fig, info] = rid.viz.signal(tc.Online, "visible", "off");
            tc.addTeardown(@close, fig);
            tc.assertEqual(info.nAxes, 3, "fixture: three panels to align");

            axList = findobj(fig, "Type", "axes");
            limits = vertcat(axList.XLim);
            tc.verifyEqual(limits, repmat(limits(1, :), size(limits, 1), 1), ...
                "every panel must share one time axis");
        end

        function aRecordDrawsTheSameFigureAsItsResult(tc)
            [figResult, fromResult] = rid.viz.signal(tc.Online,        "visible", "off");
            [figRecord, fromRecord] = rid.viz.signal(tc.Online.record, "visible", "off");
            tc.addTeardown(@close, figResult);
            tc.addTeardown(@close, figRecord);

            tc.verifyEqual(fromRecord.nDesigned,    fromResult.nDesigned);
            tc.verifyEqual(fromRecord.nAtAmplitude, fromResult.nAtAmplitude);
            tc.verifyEqual(fromRecord.hasOutput,    fromResult.hasOutput);
        end

        function aRecordFromDiskDrawsItToo(tc)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            root = tc.applyFixture(TemporaryFolderFixture()).Folder;

            path = rid.io.saveRun(tc.Online.record, "", "root", string(root));
            rec  = rid.io.loadRun(path);

            [figDisk, onDisk] = rid.viz.signal(rec, "visible", "off");
            tc.addTeardown(@close, figDisk);
            [figLive, live]   = rid.viz.signal(tc.Online, "visible", "off");
            tc.addTeardown(@close, figLive);

            tc.verifyEqual(onDisk.nDesigned,    live.nDesigned);
            tc.verifyEqual(onDisk.nAtAmplitude, live.nAtAmplitude);
            tc.verifyEqual(onDisk.hasOutput,    live.hasOutput);
        end

        function aRecordWithoutQStillDrawsTheSignal(tc)
            % The signal is in `outputs`, not in `problem`, so storeQ has
            % nothing to do with it. Only the scatter needs Q.
            r = VizFixture.quiet(2, struct("storeQ", false));
            [fig, info] = rid.viz.signal(r, "visible", "off");
            tc.addTeardown(@close, fig);
            tc.verifyGreaterThan(info.nDesigned, 0);
        end

        function aSourceThatIsNeitherIsRefused(tc)
            tc.verifyError(@() rid.viz.signal(struct("u", 1), "visible", "off"), ...
                "rid:viz:unknownSource");
        end

        function theFigureIsInvisibleOnlyWhenAsked(tc)
            hidden = rid.viz.signal(tc.Offline, "visible", "off");
            tc.addTeardown(@close, hidden);
            tc.verifyEqual(hidden.Visible, matlab.lang.OnOffSwitchState("off"));

            shown = rid.viz.signal(tc.Offline);
            tc.addTeardown(@close, shown);
            tc.verifyEqual(shown.Visible, matlab.lang.OnOffSwitchState("on"));
        end

    end

    methods (Access = private)

        function fig = figureOf(tc, r)
            fig = rid.viz.signal(r, "visible", "off");
            tc.addTeardown(@close, fig);
        end

    end

    methods (Static, Access = private)

        function n = countOnAmplitudeBound(r, nPriming)
            %COUNTONAMPLITUDEBOUND  The designed samples on a bound, worked out here.
            %
            %   Independent of the figure: it reads the constraints and the
            %   history, not the counts the plot reported.
            con = r.record.setup.setup.gen.constraints;
            u   = r.hist.u;
            n   = zeros(1, size(u, 2));
            for iCh = 1:size(u, 2)
                width = con.amplitude(iCh, 2) - con.amplitude(iCh, 1);
                at    = abs(u(:, iCh) - con.amplitude(iCh, 1)) <= 1e-6 * width | ...
                        abs(u(:, iCh) - con.amplitude(iCh, 2)) <= 1e-6 * width;
                at(1:nPriming) = false;
                n(iCh) = sum(at);
            end
        end

    end

end
