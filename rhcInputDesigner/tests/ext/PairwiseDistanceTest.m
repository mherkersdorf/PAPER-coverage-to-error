classdef PairwiseDistanceTest < matlab.unittest.TestCase
%PAIRWISEDISTANCETEST  Analytic and property tests for rid.ext.pairwiseDistance.
%
%   Layers 1 (analytic cases) and 2 (property tests) of the project's test plan. Randomized
%   cases draw from a private RandStream built from the literal seed below, which
%   is reported in every failure diagnostic; the global stream is never touched.

    properties (Constant)
        Seed = 20260828
    end

    methods (Access = private)

        function s = newStream(tc)
            s = RandStream("twister", "Seed", tc.Seed);
        end

        function d = seedDiagnostic(tc)
            d = sprintf("seed = %d", tc.Seed);
        end

    end

    methods (Test)

        function analyticSmallCase(tc)
            % Hand-computable: two points against three, unit square.
            A = [0 0; 1 0];
            B = [0 0; 0 1; 1 1];
            expected = [0 1 sqrt(2); 1 sqrt(2) 1];
            D = rid.ext.pairwiseDistance(A, B);
            tc.verifyEqual(D, expected, "AbsTol", 1e-15);
        end

        function analyticOneDimensional(tc)
            % In 1-D the Euclidean distance is the absolute difference, exactly.
            a = [0; 1; 2.5];
            b = [-1; 3];
            D = rid.ext.pairwiseDistance(a, b);
            tc.verifyEqual(D, abs(a - b.'), "AbsTol", 1e-15);
        end

        function shapeIsNumPointsAByNumPointsB(tc)
            s = tc.newStream();
            A = rand(s, 7, 3);
            B = rand(s, 5, 3);
            D = rid.ext.pairwiseDistance(A, B);
            tc.verifySize(D, [7 5], tc.seedDiagnostic());
        end

        function selfDistanceIsZeroAndMatrixSymmetric(tc)
            s = tc.newStream();
            A = rand(s, 6, 4);
            D = rid.ext.pairwiseDistance(A, A);
            tc.verifyEqual(diag(D), zeros(6, 1), tc.seedDiagnostic());
            tc.verifyEqual(D, D.', "AbsTol", 1e-15, tc.seedDiagnostic());
        end

        function swappingArgumentsTransposesResult(tc)
            s = tc.newStream();
            A = rand(s, 7, 3);
            B = rand(s, 5, 3);
            tc.verifyEqual(rid.ext.pairwiseDistance(A, B), ...
                rid.ext.pairwiseDistance(B, A).', ...
                "AbsTol", 1e-15, tc.seedDiagnostic());
        end

        function satisfiesTriangleInequality(tc)
            s = tc.newStream();
            P = rand(s, 8, 3);
            D = rid.ext.pairwiseDistance(P, P);
            for jSup = 1:size(P, 1)
                slack = D(:, jSup) + D(jSup, :) - D;   % (8 x 8) by expansion
                tc.verifyGreaterThanOrEqual(slack, -1e-12, tc.seedDiagnostic());
            end
        end

        function invariantUnderTranslation(tc)
            s = tc.newStream();
            A = rand(s, 7, 3);
            B = rand(s, 5, 3);
            shift = 10 * rand(s, 1, 3);
            tc.verifyEqual(rid.ext.pairwiseDistance(A + shift, B + shift), ...
                rid.ext.pairwiseDistance(A, B), ...
                "AbsTol", 1e-12, tc.seedDiagnostic());
        end

        function scalesWithTheCoordinates(tc)
            % Scaling epsReq0 by alpha is equivalent to scaling
            % all distances by 1/alpha -- rests on this homogeneity.
            s = tc.newStream();
            A = rand(s, 7, 3);
            B = rand(s, 5, 3);
            alpha = 3.7;
            tc.verifyEqual(rid.ext.pairwiseDistance(alpha * A, alpha * B), ...
                alpha * rid.ext.pairwiseDistance(A, B), ...
                "RelTol", 1e-12, tc.seedDiagnostic());
        end

        function nearestNeighbourReductionMatchesNaiveLoop(tc)
            % The exact expression +crit uses, against an
            % obviously-correct double loop.
            s = tc.newStream();
            Q = rand(s, 10, 3);
            Znorm = rand(s, 13, 3);

            expected = zeros(10, 1);
            for jSup = 1:10
                best = Inf;
                for iSample = 1:13
                    best = min(best, norm(Q(jSup, :) - Znorm(iSample, :), 2));
                end
                expected(jSup) = best;
            end

            dRaw = min(rid.ext.pairwiseDistance(Q, Znorm), [], 2);
            tc.verifyEqual(dRaw, expected, "AbsTol", 1e-12, tc.seedDiagnostic());
        end

        function emptyPointSetsGiveEmptyResult(tc)
            s = tc.newStream();
            B = rand(s, 5, 3);
            A = rand(s, 4, 3);
            tc.verifySize(rid.ext.pairwiseDistance(zeros(0, 3), B), [0 5]);
            tc.verifySize(rid.ext.pairwiseDistance(A, zeros(0, 3)), [4 0]);
        end

        function rejectsMismatchedDimension(tc)
            tc.verifyError(@() rid.ext.pairwiseDistance([1 2], [1 2 3]), ...
                "rid:ext:dimensionMismatch");
        end

        function rejectsNonFiniteInput(tc)
            tc.verifyError(@() rid.ext.pairwiseDistance([NaN 2], [1 2]), ...
                "MATLAB:validators:mustBeFinite");
            tc.verifyError(@() rid.ext.pairwiseDistance([1 2], [Inf 2]), ...
                "MATLAB:validators:mustBeFinite");
        end

        function rejectsComplexInput(tc)
            tc.verifyError(@() rid.ext.pairwiseDistance([1i 2], [1 2]), ...
                "MATLAB:validators:mustBeReal");
        end

        function rejectsThreeDimensionalArray(tc)
            tc.verifyError(@() rid.ext.pairwiseDistance(cat(3, eye(2), eye(2)), [1 2]), ...
                "MATLAB:validation:IncompatibleSize");
        end

    end

end
