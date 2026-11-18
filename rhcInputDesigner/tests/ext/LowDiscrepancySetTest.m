classdef LowDiscrepancySetTest < matlab.unittest.TestCase
%LOWDISCREPANCYSETTEST  Analytic and property tests for rid.ext.lowDiscrepancySet.
%
%   Layers 1 and 2 of the project's test plan. The point of the function is COVERAGE, not
%   equidistribution in the abstract: rid.roi.estimateFillDistance takes a maximum,
%   and a maximum is decided by where the probe set leaves gaps. The
%   coverage test below is therefore the contract, not decoration.

    properties (Constant)
        Seed      = 20260910
        OtherSeed = 20260911
    end

    properties (Access = private)
        SuiteStream
    end

    methods (TestMethodSetup)
        function captureGlobalStream(tc)
            tc.SuiteStream = RandStream.getGlobalStream();
        end
    end

    methods (TestMethodTeardown)
        function restoreGlobalStream(tc)
            RandStream.setGlobalStream(tc.SuiteStream);
        end
    end

    methods (Access = private)
        function d = seedDiagnostic(tc)
            d = sprintf("seed = %d", tc.Seed);
        end

        function h = coveringRadius(~, X, probes)
            % Worst distance from a probe to its nearest point of X.
            h = max(min(rid.ext.pairwiseDistance(probes, X), [], 2));
        end
    end

    methods (Test)

        function hasRequestedShapeAndLiesInUnitCube(tc)
            X = rid.ext.lowDiscrepancySet(500, 4, tc.Seed);
            tc.verifySize(X, [500 4], tc.seedDiagnostic());
            tc.verifyGreaterThanOrEqual(X, 0, tc.seedDiagnostic());
            tc.verifyLessThan(X, 1, tc.seedDiagnostic());
        end

        function coversTheCubeBetterThanAnIidSample(tc)
            % The property the estimator buys it for. Measured against a common,
            % independent reference set so the comparison is of the two point
            % sets and not of two measurements.
            nDim   = 4;
            nPoints = 2000;
            reference = rand(RandStream("twister", "Seed", tc.OtherSeed), 20000, nDim);

            lowDisc = rid.ext.lowDiscrepancySet(nPoints, nDim, tc.Seed);
            iid     = rand(RandStream("twister", "Seed", tc.Seed), nPoints, nDim);

            tc.verifyLessThan(tc.coveringRadius(lowDisc, reference), ...
                tc.coveringRadius(iid, reference), tc.seedDiagnostic());
        end

        function aShorterSetIsAPrefixOfALongerOne(tc)
            % What makes the convergence trace of rid.roi.estimateFillDistance
            % comparable across probe counts: a shorter run is a prefix, never a
            % resample.
            long  = rid.ext.lowDiscrepancySet(800, 3, tc.Seed);
            short = rid.ext.lowDiscrepancySet(200, 3, tc.Seed);
            tc.verifyEqual(long(1:200, :), short, tc.seedDiagnostic());
        end

        function deterministicUnderAFixedSeedAndDiffersAcrossSeeds(tc)
            a = rid.ext.lowDiscrepancySet(300, 3, tc.Seed);
            b = rid.ext.lowDiscrepancySet(300, 3, tc.Seed);
            c = rid.ext.lowDiscrepancySet(300, 3, tc.OtherSeed);
            tc.verifyEqual(a, b, tc.seedDiagnostic());
            tc.verifyNotEqual(a, c, ...
                sprintf("seeds = %d, %d", tc.Seed, tc.OtherSeed));
        end

        function consecutiveDrawsFromOneStreamAreIndependentSets(tc)
            % rid.roi.estimateFillDistance takes two sets from one stream: the
            % probes, and the independent set that measures the probes' own fill
            % distance. A second draw that repeated the first would make that
            % measurement zero.
            stream = RandStream("twister", "Seed", tc.Seed);
            first  = rid.ext.lowDiscrepancySet(300, 3, stream);
            second = rid.ext.lowDiscrepancySet(300, 3, stream);
            tc.verifyNotEqual(first, second, tc.seedDiagnostic());
        end

        function theDigitalShiftKeepsTheSetInsideTheCube(tc)
            % An XOR on the binary expansion is a permutation of the dyadic boxes,
            % so it maps the cube onto itself. A Cranley-Patterson shift would
            % instead wrap around the boundary -- which is the half of the ROI
            % that matters most here.
            for seed = [0 1 7 12345]
                X = rid.ext.lowDiscrepancySet(1000, 5, seed);
                tc.verifyGreaterThanOrEqual(X, 0, sprintf("seed = %d", seed));
                tc.verifyLessThan(X, 1, sprintf("seed = %d", seed));
            end
        end

        function doesNotTouchTheGlobalStream(tc)
            RandStream.setGlobalStream(RandStream("philox", "Seed", 909));
            before      = RandStream.getGlobalStream();
            stateBefore = before.State;

            rid.ext.lowDiscrepancySet(1000, 3, tc.Seed);

            after = RandStream.getGlobalStream();
            tc.verifySameHandle(after, before, tc.seedDiagnostic());
            tc.verifyEqual(after.State, stateBefore, tc.seedDiagnostic());
        end

        function rejectsAnInvalidSeed(tc)
            tc.verifyError(@() rid.ext.lowDiscrepancySet(10, 2, -1), ...
                "rid:ext:invalidSeed");
            tc.verifyError(@() rid.ext.lowDiscrepancySet(10, 2, 1.5), ...
                "rid:ext:invalidSeed");
        end

    end

end
