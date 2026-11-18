classdef RecordFormatTest < matlab.unittest.TestCase
%RECORDFORMATTEST  `run.mat` is -v7, -v7.3 above 2 GB, and both round-trip.
%
%   Each format is checked end to end -- save, load, replay -- because a record
%   that loads but no longer replays is a record that has lost what makes it worth
%   keeping. The format is read from the file's own header, not inferred from its
%   size: "MATLAB 5.0 MAT-file" for -v7, "MATLAB 7.3 MAT-file" for -v7.3.
%
%   Records are written into a temporary folder, never into the project's `runs/`.

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

        function aRecordIsWrittenAsV7AndReplaysAfterTheRoundTrip(tc)
            % Fails if: saveRun writes -v7.3, or a -v7 record loses anything
            % replay needs (the handle members of the stored Setup above all).
            rec  = ReplayTest.recordOf("nSteps", 4);
            path = rid.io.saveRun(rec, "", "root", tc.Root);
            tc.verifyEqual(RecordFormatTest.formatOf(path), "-v7");

            back          = rid.io.loadRun(path);
            expected      = rec;
            expected.note = "";
            tc.verifyEqual(back, expected);

            [ok, diff] = rid.io.replay(back);
            tc.verifyTrue(ok, "replay diverged: " + diff.reason);
        end

        function anExistingV73RecordStaysReadableAndReplays(tc)
            % Every record written before the record file format is -v7.3, and loadRun is unchanged:
            % it must still read them, and they must still replay.
            rec  = ReplayTest.recordOf("nSteps", 4);
            path = fullfile(tc.Root, rec.identity.runId);
            mkdir(path);
            save(fullfile(path, "run.mat"), "rec", "-v7.3");
            tc.assertEqual(RecordFormatTest.formatOf(path), "-v7.3", ...
                "fixture: the record must be -v7.3");

            back = rid.io.loadRun(path);
            tc.verifyEqual(back, rec);

            [ok, diff] = rid.io.replay(back);
            tc.verifyTrue(ok, "replay diverged: " + diff.reason);
        end

    end

    methods (Test, TestTags = {'LargeMemory'})

        function aRecordOverTwoGigabytesFallsBackToV73AndIsNotLost(tc)
            % `save -v7` of a variable of 2 GB or more writes a file WITHOUT it
            % and only warns. A real record of 2 GB + 8 B is used rather than a
            % lowered threshold, so the test exercises MATLAB's own behaviour
            % (about 4 GB of memory at its peak, a few seconds).
            %
            % Fails if: saveRun keeps the empty -v7 file (loadRun then finds no
            % `rec`), or falls back by any route other than the file's contents.
            rec = ReplayTest.recordOf("nSteps", 2);
            rec.outputs.padding = zeros(2^28 + 1, 1);
            tc.assertGreaterThanOrEqual(8 * numel(rec.outputs.padding), 2^31, ...
                "fixture: the record must exceed what -v7 can hold");

            path = rid.io.saveRun(rec, "", "root", tc.Root);
            tc.verifyEqual(RecordFormatTest.formatOf(path), "-v7.3");

            back = rid.io.loadRun(path);
            rec.note = "";
            % isequaln, not verifyEqual: a diagnostic over a 2 GB field is useless.
            tc.verifyTrue(isequaln(back, rec), "the record did not survive the fallback");
        end

    end

    methods (Static)

        function format = formatOf(path)
            %FORMATOF  "-v7" or "-v7.3", from the text header of run.mat.
            fid = fopen(fullfile(path, "run.mat"), "r");
            header = string(fread(fid, [1 19], "*char"));
            fclose(fid);
            switch header
                case "MATLAB 5.0 MAT-file"
                    format = "-v7";
                case "MATLAB 7.3 MAT-file"
                    format = "-v7.3";
                otherwise
                    format = "unknown: " + header;
            end
        end

    end

end
