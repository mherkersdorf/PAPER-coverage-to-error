classdef SupportingSetTest < matlab.unittest.TestCase
%SUPPORTINGSETTEST  Stored supporting-point sets: round trip and refusals.
%
%   rid.roi.saveSupportingSet writes what opts.roi.method = "stored" loads back
%   through rid.roi.supportingPoints. Three things are checked on load and
%   corrected nowhere (the pattern): nSup, nDim and the normalized box
%   extent. The stored fillDistanceQ is a CHECKSUM, re-measured and compared by
%   the probe brackets overlap, never trusted on its own; identity is the CONTENT
%   HASH, not the path.
%
%   Every file this class writes is under the real, git-ignored
%   supportingSets/ -- the stored-Q reuse deliberately gives rid.roi.saveSupportingSet no
%   root override, unlike rid.io.saveRun, so there is nowhere else to redirect
%   it to. Each test uses a unique name and deletes its own file in teardown.
%
%   THE PROBE COUNT IS HELD FIXED AT 400, ON BOTH SIDES OF EVERY ROUND TRIP.
%   With matching nProbes and the default probeSeed the Monte-Carlo estimate
%   is reproducible bit-for-bit (same Q, same settings), so the bracket-overlap
%   check passes by construction whenever Q is unchanged and only fails when it
%   is meant to -- the damaged-points test collapses half the design instead of
%   relying on probe-count noise to manufacture a disagreement.

    properties (Constant)
        Seed    = 20260914
        NProbes = 400
    end

    methods (Access = private)

        function nrm = square(~)
            nrm = rid.roi.Normalization([0 1; 0 1], ["u" "y"]);
        end

        function nrm = weightedSquare(~)
            nrm = rid.roi.Normalization([0 1; 0 1], ["u" "y"], ...
                "dimWeights", [2 1]);
        end

        function [Q, fillDistanceQ, info] = generatedOn(tc, nrm, nSup)
            [Q, fillDistanceQ, info] = rid.roi.supportingPoints( ...
                nrm, nSup, "seed", tc.Seed, "nProbes", tc.NProbes, ...
                "progress", "off");
        end

        function [Q, fillDistanceQ, info] = generated(tc, nSup)
            [Q, fillDistanceQ, info] = tc.generatedOn(tc.square(), nSup);
        end

        function [Qback, fillDistanceQ, info] = loaded(tc, nrm, nSup, name)
            [Qback, fillDistanceQ, info] = rid.roi.supportingPoints(nrm, nSup, ...
                "method", "stored", "supportingSet", name, ...
                "nProbes", tc.NProbes, "progress", "off");
        end

        function [name, path] = save(tc, Q, nrm, fillDistanceQ, info)
            name = string(java.util.UUID.randomUUID()) + ".mat";
            path = rid.roi.saveSupportingSet(Q, nrm, fillDistanceQ, info, name);
            tc.addTeardown(@() SupportingSetTest.removeIfPresent(path));
        end

    end

    methods (Static, Access = private)

        function removeIfPresent(path)
            if isfile(path)
                delete(path);
            end
        end

    end

    methods (Test)

        % --- the round trip ---------------------------------------------------

        function roundTripPreservesQExactlyAndNamesTheHash(tc)
            % THE CENTRAL CASE. Fails if: loading rescales, regenerates, or
            % otherwise returns anything but the exact bytes that were saved --
            % or if the hash it reports is not reproducible across reloads of
            % the same file, which is what "identity is the content" means.
            nSup = 30;
            [Q, fillDistanceQ, info] = tc.generated(nSup);
            [name, path] = tc.save(Q, tc.square(), fillDistanceQ, info);
            tc.assertTrue(isfile(path), "fixture: the file must exist");

            [Qback, ~, info1] = tc.loaded(tc.square(), nSup, name);

            tc.verifyEqual(Qback, Q, "the loaded Q must be bit-identical");
            tc.verifyEqual(info1.origin, "stored");
            tc.verifyEqual(info1.method, "stored");
            tc.verifyTrue(isnan(info1.seed), "no generator ran THIS call");
            tc.verifyNotEmpty(info1.hash);
            tc.verifyNotEqual(info1.hash, "unknown");

            [~, ~, info2] = tc.loaded(tc.square(), nSup, name);
            tc.verifyEqual(info2.hash, info1.hash, ...
                "the hash is a function of the content: reloading must reproduce it");
        end

        function theRunRecordNamesTheOriginAndTheHashOfAStoredSet(tc)
            % The Q provenance and the stored-Q reuse together: a run built from a "stored" Q must
            % report so in the record, by hash, never by inferring it from
            % opts.roi.* -- and storeQ's promised default for a non-generated
            % Q must have actually fired.
            nSup = 10;
            f = SetupTest.fixture("nSup", nSup, "epsReq0", RunTest.EpsReq0, ...
                "seed", tc.Seed);

            [Q, fillDistanceQ, info] = tc.generatedOn(f.nrm, nSup);
            [name, ~] = tc.save(Q, f.nrm, fillDistanceQ, info);
            [Qback, ~, qInfo] = tc.loaded(f.nrm, nSup, name);

            fillQ  = rid.roi.estimateFillDistance(Qback, f.nrm);
            epsReq = rid.crit.demandedResolution(f.fields, ...
                "epsReq0", RunTest.EpsReq0, "mode", "budget", ...
                "fillDistanceQ", fillQ);
            state  = rid.crit.CoverageState(Qback, epsReq);
            setup  = rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, Qback, f.fields, ...
                f.gen, f.sol, state, f.sch, f.fieldSch, qInfo);

            ov = RunTest.overrides("nSteps", 2);
            r  = rid.rhc.run([], setup, RunTest.history("srg", f.srg), ov);

            tc.verifyEqual(r.record.problem.qOrigin, "stored");
            tc.verifyEqual(r.record.problem.qHash, qInfo.hash);
            tc.verifyTrue(logical(r.record.problem.storedQ), ...
                "storeQ must default to true for a non-generated Q");
        end

        % --- save-time refusals -------------------------------------------------

        function saveRefusesToOverwriteAnExistingFile(tc)
            [Q, fillDistanceQ, info] = tc.generated(15);
            [name, ~] = tc.save(Q, tc.square(), fillDistanceQ, info);
            tc.verifyError(@() rid.roi.saveSupportingSet(Q, tc.square(), ...
                fillDistanceQ, info, name), "rid:roi:supportingSetExists");
        end

        function saveRefusesAPointOutsideTheDeclaredBox(tc)
            [Q, fillDistanceQ, info] = tc.generated(15);
            bad = Q;
            bad(1, 1) = 5;   % outside [0, 1]
            name = string(java.util.UUID.randomUUID()) + ".mat";
            tc.verifyError(@() rid.roi.saveSupportingSet(bad, tc.square(), ...
                fillDistanceQ, info, name), "rid:roi:supportingSetOutOfBounds");
        end

        function saveRefusesAnAnalyticFillDistanceProvenance(tc)
            % Only a Monte-Carlo fillDistanceQ has the bracket a later
            % load re-checks against; "grid" is cheap to regenerate anyway.
            [Q, fillDistanceQ, info] = rid.roi.supportingPoints(tc.square(), ...
                16, "method", "grid", "progress", "off");
            name = string(java.util.UUID.randomUUID()) + ".mat";
            tc.verifyError(@() rid.roi.saveSupportingSet(Q, tc.square(), ...
                fillDistanceQ, info, name), ...
                "rid:roi:unsupportedFillDistanceProvenance");
        end

        % --- load-time refusals: nothing is corrected ------------------

        function loadRefusesAMissingFile(tc)
            tc.verifyError(@() tc.loaded(tc.square(), 10, ...
                string(java.util.UUID.randomUUID()) + ".mat"), ...
                "rid:roi:supportingSetNotFound");
        end

        function loadRefusesANSupMismatch(tc)
            [Q, fillDistanceQ, info] = tc.generated(20);
            [name, ~] = tc.save(Q, tc.square(), fillDistanceQ, info);
            tc.verifyError(@() tc.loaded(tc.square(), 21, name), ...
                "rid:roi:supportingSetSizeMismatch");
        end

        function loadRefusesADimensionMismatch(tc)
            [Q, fillDistanceQ, info] = tc.generated(20);
            [name, ~] = tc.save(Q, tc.square(), fillDistanceQ, info);
            nrm3 = rid.roi.Normalization([0 1; 0 1; 0 1], ["u" "y" "z"]);
            tc.verifyError(@() tc.loaded(nrm3, 20, name), ...
                "rid:roi:supportingSetDimensionMismatch");
        end

        function loadRefusesABoxExtentMismatch(tc)
            % A stretched box is the anisotropy failure exactly: a
            % maximin design scaled per dimension is no longer maximin.
            [Q, fillDistanceQ, info] = tc.generatedOn(tc.square(), 25);
            [name, ~] = tc.save(Q, tc.square(), fillDistanceQ, info);
            tc.verifyError(@() tc.loaded(tc.weightedSquare(), 25, name), ...
                "rid:roi:supportingSetBoxMismatch");
        end

        function loadRefusesDamagedPoints(tc)
            % The stored fillDistanceQ is a CHECKSUM: collapsing half the
            % design onto its own centroid stays inside the box (so the
            % bounds check alone would not catch it) but drastically raises
            % the true fill distance, which re-measurement must catch.
            [Q, fillDistanceQ, info] = tc.generated(60);
            [name, path] = tc.save(Q, tc.square(), fillDistanceQ, info);

            raw     = load(path, "stored");
            damaged = raw.stored;
            half    = 1:floor(size(damaged.Q, 1) / 2);
            damaged.Q(half, :) = repmat(mean(damaged.Q, 1), numel(half), 1);
            stored  = damaged; %#ok<NASGU>
            save(path, "stored");

            tc.verifyError(@() tc.loaded(tc.square(), 60, name), ...
                "rid:roi:supportingSetFillDistanceMismatch");
        end

        % --- the option coupling ---------------------------------------

        function supportingSetIsRequiredByStoredAndRefusedOtherwise(tc)
            tc.verifyError(@() rid.roi.supportingPoints(tc.square(), 10, ...
                "method", "stored", "progress", "off"), ...
                "rid:roi:missingSupportingSet");
            tc.verifyError(@() rid.roi.supportingPoints(tc.square(), 10, ...
                "method", "lhs", "supportingSet", "x.mat", "progress", "off"), ...
                "rid:roi:supportingSetNotApplicable");
        end

    end

end
