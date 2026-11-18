classdef EventKindsTest < matlab.unittest.TestCase
%EVENTKINDSTEST  rid.eventKinds: the closed set has one owner.
%
%   THE SET WAS WRITTEN OUT TWICE -- once in `rid.rhc.run`'s local `event()` as
%   the guard that closes it, once in `+viz/private/eventKinds` because `+rhc` is
%   not reachable from `+viz` (`+viz` depends on results only). Two copies of one
%   normative list are two things that can drift, and the coordinate-meaning rule is that a
%   closed set has a single source its consumers read rather than restate.
%
%   WHY IT SITS AT `rid.eventKinds` AND NOT IN A PACKAGE. It is read by `+rhc`,
%   which emits the events, and by `+viz`, which draws them; putting it in either
%   would make one package depend on the other, and `+viz` depending on `+rhc`
%   is the dependency the architecture forbids outright. `rid.defaultOptions` and
%   `rid.assertHistory` already establish the top level as where something both
%   sides need lives.
%
%   WHAT THIS DOES NOT CHANGE. An unrecognized kind still draws, in grey. The
%   list decides colours, not membership, so a kind emitted before the list knows
%   about it is still visible -- which is the property that kept the second copy
%   harmless while it existed, and is worth keeping now that it does not.

    methods (Test)

        function theSetIsTheOneDocsRunRecordsSectionOneDeclares(tc)
            % Fails if: a kind is added, removed or reordered here without the
            % specification moving first. The order is section 1's order, and
            % consumers that colour by index depend on it.
            tc.verifyEqual(rid.eventKinds(), ...
                ["surrogateUpdate", "resolutionChange", "coverageReset", ...
                 "supportChange", "normalizationChange", "modeChange", "seed"]);
        end

        function itIsAStringRowOfSeven(tc)
            kinds = rid.eventKinds();
            tc.verifyClass(kinds, "string");
            tc.verifySize(kinds, [1 7]);
            tc.verifyEqual(numel(unique(kinds)), 7, "no duplicates");
        end

        function everyKindInTheSetGetsItsOwnColour(tc)
            % `+viz` colours by INDEX into this set, so a reordering or a
            % removal silently repaints every kind after it. Drawing one event of
            % each kind and counting distinct colours is what catches that.
            rec = EventKindsTest.recordWithKinds(rid.eventKinds());
            [fig, info] = rid.viz.progress(rec, "visible", "off");
            tc.addTeardown(@close, fig);

            tc.verifyEqual(info.nMarks, 7);
            colours = unique(vertcat(findobj(fig, "Type", "line").Color), "rows");
            tc.verifyEqual(size(colours, 1), 8, ...
                "seven kinds plus the J trace, all distinct");
        end

        function AKindOutsideTheSetIsStillDrawn(tc)
            % THE PROPERTY THAT KEPT THE SECOND COPY HARMLESS, and the reason it
            % is worth keeping now there is only one. Nothing branches on
            % MEMBERSHIP: the set decides colours, so a kind emitted before this
            % list knows about it appears in grey rather than vanishing.
            %
            % Fails if: `+viz` starts treating the set as a whitelist, in which
            % case a newly-emitted kind would be invisible at exactly the moment
            % somebody needs to see it.
            rec = EventKindsTest.recordWithKinds("somethingNobodyHasAddedYet");
            [fig, info] = rid.viz.progress(rec, "visible", "off");
            tc.addTeardown(@close, fig);

            tc.verifyEqual(info.nMarks, 1);
            drawn = findobj(fig, "Tag", "somethingNobodyHasAddedYet");
            tc.verifyNotEmpty(drawn, "an unknown kind must still be drawn");
            tc.verifyEqual(drawn.Color, [0.40 0.40 0.40], "in grey");
        end

        function thereIsExactlyOneCopyOfTheSetInTheSource(tc)
            % The coordinate-meaning rule made checkable. Fails if a second copy comes back --
            % which is how the first one arrived: not by anyone deciding to
            % duplicate the list, but by a consumer being unable to reach it.
            root  = fullfile(rid.io.projectRoot(), "src");
            files = dir(fullfile(root, "**", "*.m"));

            holders = strings(1, 0);
            for iFile = 1:numel(files)
                path = fullfile(files(iFile).folder, files(iFile).name);
                text = string(fileread(path));
                % The marker is the LAST two kinds adjacent: no consumer needs
                % to name those together except a copy of the whole set.
                if contains(text, """modeChange""") && contains(text, """seed""")
                    holders(end + 1) = string(files(iFile).name); %#ok<AGROW>
                end
            end
            tc.verifyEqual(holders, "eventKinds.m", ...
                "the closed set belongs to rid.eventKinds and nowhere else");
        end

    end

    methods (Static, Access = private)

        function rec = recordWithKinds(kinds)
            %RECORDWITHKINDS  A real record, with one event of each kind spliced in.
            %
            %   The record comes from an actual run so every field `+viz` reads
            %   is the one `rid.rhc.run` produces; only the event log is
            %   replaced. The loop emits three of the seven kinds today, so a
            %   run alone cannot exercise the other four -- and a fixture in
            %   which four of the cases are absent is not a fixture for them
            %  .
            rec = VizFixture.quiet(2).record;
            log = struct("step", {}, "kind", {}, "detail", {});
            for iKind = 1:numel(kinds)
                log(end + 1, 1) = struct( ...
                    "step",   min(iKind, numel(rec.outputs.Jtrace)), ...
                    "kind",   string(kinds(iKind)), ...
                    "detail", struct()); %#ok<AGROW>
            end
            rec.outputs.eventLog = log;
        end

    end

end
