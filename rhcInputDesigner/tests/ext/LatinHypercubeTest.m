classdef LatinHypercubeTest < matlab.unittest.TestCase
%LATINHYPERCUBETEST  Analytic and property tests for rid.ext.latinHypercube.
%
%   Layers 1 (analytic cases) and 2 (property tests) of the project's test plan. The seeds are
%   literals, reported in the failure diagnostics. Several tests deliberately set
%   the global stream in order to observe what the function does to it;
%   the fixture below restores whatever stream the suite was run with.

    properties (Constant)
        Seed      = 20260828
        OtherSeed = 20260829
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
    end

    methods (Test)

        function hasRequestedShapeAndLiesInUnitCube(tc)
            X = rid.ext.latinHypercube(12, 3, tc.Seed);
            tc.verifySize(X, [12 3], tc.seedDiagnostic());
            tc.verifyGreaterThan(X, 0, tc.seedDiagnostic());
            tc.verifyLessThan(X, 1, tc.seedDiagnostic());
        end

        function everyStratumHoldsExactlyOnePoint(tc)
            % The defining Latin-hypercube property, per dimension.
            nPoints = 12;
            X = rid.ext.latinHypercube(nPoints, 3, tc.Seed);
            for kDim = 1:size(X, 2)
                stratum = floor(X(:, kDim) * nPoints);
                tc.verifyEqual(sort(stratum), (0:nPoints - 1).', ...
                    tc.seedDiagnostic());
            end
        end

        function unsmoothedSampleSitsOnStratumMidpoints(tc)
            % Analytic: with smoothing off and no criterion iteration, each column
            % is a permutation of (0.5 : nPoints-0.5) / nPoints, exactly.
            nPoints = 4;
            X = rid.ext.latinHypercube(nPoints, 2, tc.Seed, ...
                "smooth", false, "criterion", "none");
            expected = ((1:nPoints).' - 0.5) / nPoints;
            for kDim = 1:size(X, 2)
                tc.verifyEqual(sort(X(:, kDim)), expected, tc.seedDiagnostic());
            end
        end

        function isDeterministicUnderAFixedSeed(tc)
            % Within a MATLAB release -- the Latin-hypercube seed handling; run records store Q rather than
            % relying on regeneration.
            X1 = rid.ext.latinHypercube(9, 2, tc.Seed);
            X2 = rid.ext.latinHypercube(9, 2, tc.Seed);
            tc.verifyEqual(X1, X2, tc.seedDiagnostic());
        end

        function differentSeedsGiveDifferentDesigns(tc)
            X1 = rid.ext.latinHypercube(9, 2, tc.Seed);
            X2 = rid.ext.latinHypercube(9, 2, tc.OtherSeed);
            tc.verifyNotEqual(X1, X2, ...
                sprintf("seeds = %d, %d", tc.Seed, tc.OtherSeed));
        end

        function leavesTheGlobalStreamUntouched(tc)
            % The decisive test for the Latin-hypercube seed handling: the caller's next draw is unaffected.
            rng(42, "twister");
            expected = rand;

            rng(42, "twister");
            streamBefore = RandStream.getGlobalStream();
            stateBefore  = streamBefore.State;

            rid.ext.latinHypercube(9, 2, tc.Seed);

            streamAfter = RandStream.getGlobalStream();
            tc.verifySameHandle(streamAfter, streamBefore);
            tc.verifyEqual(streamAfter.State, stateBefore);
            tc.verifyEqual(rand, expected);
        end

        function restoresTheGlobalStreamWhenTheCallFails(tc)
            % onCleanup, not the happy path: lhsdesign rejects the criterion after
            % the borrowed stream has been installed.
            rng(42, "twister");
            streamBefore = RandStream.getGlobalStream();
            stateBefore  = streamBefore.State;

            tc.verifyError(@() rid.ext.latinHypercube(9, 2, tc.Seed, ...
                "criterion", "notACriterion"), ?MException);

            streamAfter = RandStream.getGlobalStream();
            tc.verifySameHandle(streamAfter, streamBefore);
            tc.verifyEqual(streamAfter.State, stateBefore);
        end

        function doesNotDependOnTheAmbientGlobalState(tc)
            rng(1, "twister");
            X1 = rid.ext.latinHypercube(9, 2, tc.Seed);

            rng(999, "twister");
            rand(17, 1);                      % advance the global stream arbitrarily
            X2 = rid.ext.latinHypercube(9, 2, tc.Seed);

            tc.verifyEqual(X1, X2, tc.seedDiagnostic());
        end

        function acceptsARandStreamInsteadOfASeed(tc)
            % The form +roi needs for substreams inside parfor.
            X1 = rid.ext.latinHypercube(9, 2, RandStream("twister", "Seed", tc.Seed));
            X2 = rid.ext.latinHypercube(9, 2, tc.Seed);
            tc.verifyEqual(X1, X2, tc.seedDiagnostic());
        end

        function handlesSinglePointAndSingleDimension(tc)
            X = rid.ext.latinHypercube(1, 3, tc.Seed);
            tc.verifySize(X, [1 3], tc.seedDiagnostic());
            tc.verifyGreaterThan(X, 0, tc.seedDiagnostic());
            tc.verifyLessThan(X, 1, tc.seedDiagnostic());

            nPoints = 6;
            X = rid.ext.latinHypercube(nPoints, 1, tc.Seed);
            tc.verifySize(X, [nPoints 1], tc.seedDiagnostic());
            tc.verifyEqual(sort(floor(X * nPoints)), (0:nPoints - 1).', ...
                tc.seedDiagnostic());
        end

        function rejectsInvalidSizes(tc)
            tc.verifyError(@() rid.ext.latinHypercube(2.5, 2, tc.Seed), ...
                "MATLAB:validators:mustBeInteger");
            tc.verifyError(@() rid.ext.latinHypercube(0, 2, tc.Seed), ...
                "MATLAB:validators:mustBePositive");
            tc.verifyError(@() rid.ext.latinHypercube(9, -1, tc.Seed), ...
                "MATLAB:validators:mustBePositive");
        end

        function rejectsAnInvalidSeed(tc)
            tc.verifyError(@() rid.ext.latinHypercube(9, 2, -1), ...
                "rid:ext:invalidSeed");
            tc.verifyError(@() rid.ext.latinHypercube(9, 2, 1.5), ...
                "rid:ext:invalidSeed");
            tc.verifyError(@() rid.ext.latinHypercube(9, 2, "x"), ...
                "rid:ext:invalidSeed");
        end

        function rejectsAnUnknownCriterion(tc)
            % The value list is lhsdesign's, so the identifier is not pinned here.
            tc.verifyError(@() rid.ext.latinHypercube(9, 2, tc.Seed, ...
                "criterion", "nope"), ?MException);
        end

    end

end
