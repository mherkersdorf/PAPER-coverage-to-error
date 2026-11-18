classdef StoredQTest < matlab.unittest.TestCase
%STOREDQTEST  `storeQ = false` actually removes `Q`, and says so truthfully.
%
%   THE DEFECT THIS CLOSES. The flag turned off `rec.problem.Q` alone, while
%   the single-setup rule's "the setup is stored in full" kept `Q` in two further places --
%   `setup.Q` and `setup.state`. The record therefore carried the points AND
%   announced that it did not, and `rid.io.replay` refused on the flag rather
%   than on the contents: stricter than the record required, for a reason that
%   was not true. Not a record misreporting the run, but a record misreporting
%   itself.
%
%   THE FLAG EXISTS FOR SIZE, so "the file is
%   measurably smaller" is not a nice-to-have here -- it is the whole purpose,
%   and a definition under which the file never shrinks would leave the option
%   doing nothing at all. That is why the size assertion is in this class and
%   not treated as a performance detail.
%
%   REPLAY REFUSES ON CONTENTS, NOT ON THE FLAG, so the two can no longer
%   disagree: whichever way a record was produced, what decides is whether the
%   points are there.

    properties (Access = private)
        Root
    end

    methods (TestMethodSetup)

        function useATemporaryRunsRoot(tc)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            tc.Root = string(tc.applyFixture(TemporaryFolderFixture()).Folder);
        end

    end

    methods (Test)

        function withoutStoreQTheRecordContainsQNowhere(tc)
            % THE SHARP ONE. All three places, named individually, because the
            % defect was precisely that one of them was cleared and two were not.
            rec = ReplayTest.recordOf("nSteps", 3, "storeQ", false);

            tc.verifyFalse(logical(rec.problem.storedQ), "the flag says no");
            tc.verifyEmpty(rec.problem.Q,               "problem.Q");
            tc.verifyEmpty(rec.setup.setup.Q,           "setup.Q");
            tc.verifyEmpty(rec.setup.setup.state.Q,     "setup.state.Q");
        end

        function withStoreQTheRecordCarriesQInTheDeclaredPlace(tc)
            rec = ReplayTest.recordOf("nSteps", 3, "storeQ", true);
            tc.verifyTrue(logical(rec.problem.storedQ));
            tc.verifySize(rec.problem.Q, [10, 3]);
            tc.verifyEqual(rec.setup.setup.Q, rec.problem.Q, ...
                "and the stored Setup agrees with it");
        end

        function theFlagStaysTrueOfTheRecordAfterARoundTrip(tc)
            % Fails if: the strip happens on the way to disk but not in the
            % record the caller holds, or the reverse. `result.record` and the
            % record read back must make the same statement.
            path = rid.io.saveRun(ReplayTest.recordOf("nSteps", 3, "storeQ", false), ...
                       "", "root", tc.Root);
            rec  = rid.io.loadRun(path);

            tc.verifyEmpty(rec.problem.Q);
            tc.verifyEmpty(rec.setup.setup.Q);
            tc.verifyEmpty(rec.setup.setup.state.Q);
            tc.verifyFalse(logical(rec.problem.storedQ));
        end

        function replayDeclinesARecordWithoutQ(tc)
            [ok, diff] = rid.io.replay( ...
                ReplayTest.recordOf("nSteps", 3, "storeQ", false));

            tc.verifyFalse(ok);
            tc.verifyFalse(diff.replayable);
            tc.verifySubstring(diff.reason, "Q");
        end

        function replayReproducesARecordWithQ(tc)
            % The other direction, or the refusal above proves nothing: a replay
            % that never succeeds would pass that test too.
            [ok, diff] = rid.io.replay( ...
                ReplayTest.recordOf("nSteps", 3, "storeQ", true));

            tc.verifyTrue(diff.replayable);
            tc.verifyTrue(ok, "a record WITH Q must reproduce its own history");
        end

        function replayRefusesOnTheContentsAndNotOnTheFlag(tc)
            % THE ANCHOR FOR "the two cannot disagree". A record whose flag lies
            % -- Q genuinely absent, storedQ forced true -- must still be
            % refused, because what decides is whether the points are there.
            rec = ReplayTest.recordOf("nSteps", 3, "storeQ", false);
            rec.problem.storedQ = true;               % the flag now lies

            [ok, diff] = rid.io.replay(rec);
            tc.verifyFalse(ok);
            tc.verifyFalse(diff.replayable, ...
                "the contents decide, not the announcement");
        end

        function theSavedFileIsMeasurablySmallerWithoutQ(tc)
            % WHAT THE FLAG IS FOR (section 7). While Q sat inside the stored
            % Setup the file did not actually shrink, so the option did nothing
            % -- which is the argument the storeQ rule gives for fixing the flag rather
            % than redefining it.
            withQ    = StoredQTest.bytesOf(tc, StoredQTest.runWith(true),  "withQ");
            withoutQ = StoredQTest.bytesOf(tc, StoredQTest.runWith(false), "withoutQ");

            % THREE copies of an (nSup x 3) double were stored -- problem.Q,
            % setup.Q and setup.state.Q -- and the flag used to drop only the
            % first. So the threshold is TWO copies: at one copy the assertion
            % would already have passed before the fix and could not have failed
            % (the seen-to-fail rule, anchored where the statement is sharp).
            oneCopy = StoredQTest.NSup * 3 * 8;
            tc.verifyGreaterThan(withQ - withoutQ, 2 * oneCopy, ...
                sprintf("withQ = %d bytes, withoutQ = %d, one copy of Q = %d", ...
                        withQ, withoutQ, oneCopy));
        end

    end

    properties (Constant, Access = private)
        NSup = 2000   % large enough that Q dominates the difference in file size
    end

    methods (Static, Access = private)

        function r = runWith(storeQ)
            %RUNWITH  One short run with a LARGE Q, for the file-size test only.
            %
            %   The assertions about presence and about replay use
            %   `ReplayTest.recordOf`, which is the fixture already known to
            %   reproduce bit-for-bit. This one exists because the size claim
            %   needs a Q big enough to dominate the file, and nothing else.
            rep = rid.rep.Narx(2, 1, "pastSource", "simulated");
            roi = VizFixture.region(rep);
            [box, dimName] = rep.bounds(roi);
            nrm = rid.roi.Normalization(box, dimName);
            srg = rid.surrogate.LinearFirstOrder(0.5, 0.1, 1);

            nSup   = StoredQTest.NSup;
            [Q, ~, qInfo] = rid.roi.supportingPoints(nrm, nSup, "seed", 1, ...
                                                      "progress", "off");
            fields = rid.fields.nominalFields(nSup);
            fillQ  = rid.roi.estimateFillDistance(Q, nrm);
            epsReq = rid.crit.demandedResolution(fields, "epsReq0", 0.3, ...
                         "mode", "budget", "fillDistanceQ", fillQ);
            state  = rid.crit.CoverageState(Q, epsReq);

            setup = rid.rhc.Setup(srg, rep, nrm, roi, Q, fields, ...
                        rid.gen.Direct(rid.gen.Constraints(rep, roi), 5), ...
                        rid.opt.GridSearch(), state, ...
                        rid.surrogate.ImmediateScheduler(), ...
                        rid.fields.ImmediateScheduler(), ...
                        qInfo);

            ov = struct("crit", struct("epsReq0", 0.3), ...
                        "rhc",  struct("nSteps", 2), ...
                        "io",   struct("saveRun", false, "storeQ", storeQ));
            uPrime = [0; 0.4; -0.3; 0.6; -0.5];
            hist0  = struct("u", uPrime, "x0", srg.steadyState(uPrime(1)), ...
                            "tS", srg.tS);
            r = rid.rhc.run([], setup, hist0, ov);
        end

        function bytes = bytesOf(tc, r, slug)
            %BYTESOF  Size of the written run.mat.
            r.record.identity.runId = r.record.identity.runId + "_" + slug;
            path = rid.io.saveRun(r.record, "", "root", tc.Root);
            info = dir(fullfile(path, "run.mat"));
            bytes = info.bytes;
        end

    end

end
