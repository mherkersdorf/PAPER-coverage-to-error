classdef NormalizationTest < matlab.unittest.TestCase
%NORMALIZATIONTEST  Analytic and property tests for rid.roi.Normalization.
%
%   Layers 1 (analytic cases) and 2 (property tests) of the project's test plan. The transform is
%   deterministic, so no seed is involved. The load-bearing tests are
%   roundTripIsExactOutsideTheUnitCube and transformNeverClamps: the ROI is a
%   ruler, not a fence, so normalized values
%   outside [0,1] are legitimate and must survive both directions.

    properties (Constant)
        % Deliberately asymmetric and of different magnitudes per dimension, so
        % that a swapped offset or scale cannot pass unnoticed.
        Box     = [-2 6; 100 300; -0.5 -0.1]
        DimName = ["u" "y" "dy"]
        Weights = [3 1 0.25]
        RelTol  = 1e-12
    end

    methods (Access = private)
        function nrm = unitWeighted(tc)
            nrm = rid.roi.Normalization(tc.Box, tc.DimName);
        end
        function nrm = weighted(tc)
            nrm = rid.roi.Normalization(tc.Box, tc.DimName, "dimWeights", tc.Weights);
        end
    end

    methods (Test)

        % --- construction ----------------------------------------------------

        function constructorStoresBoxNamesAndWeights(tc)
            nrm = tc.weighted();
            tc.verifyEqual(nrm.nDim, 3);
            tc.verifyEqual(nrm.box, tc.Box);
            tc.verifyEqual(nrm.dimName, tc.DimName);
            tc.verifyEqual(nrm.dimWeights, tc.Weights);
            tc.verifyEqual(nrm.offset, tc.Box(:, 1).');
            tc.verifyEqual(nrm.span, (tc.Box(:, 2) - tc.Box(:, 1)).');
        end

        function dimWeightsDefaultToOne(tc)
            nrm = tc.unitWeighted();
            tc.verifyEqual(nrm.dimWeights, ones(1, 3));
        end

        function isValueClassNotHandle(tc)
            % Value semantics: everything outside CoverageState and the update
            % schedulers is pure.
            tc.verifyFalse(isa(tc.unitWeighted(), "handle"));
        end

        % --- the map ---------------------------------------------------------

        function unitWeightsMapBoxOntoUnitCube(tc)
            % Analytic: lower -> 0, midpoint -> 0.5, upper -> 1, exactly.
            nrm = tc.unitWeighted();
            Zphys = [tc.Box(:, 1).'; mean(tc.Box, 2).'; tc.Box(:, 2).'];
            tc.verifyEqual(nrm.toNormalized(Zphys), ...
                [0 0 0; 0.5 0.5 0.5; 1 1 1], "RelTol", tc.RelTol, "AbsTol", tc.RelTol);
        end

        function nonUnitWeightsMapBoxOntoWeightedBox(tc)
            % The normalized ROI is [0,w_1] x ... x [0,w_nDim].
            nrm = tc.weighted();
            Zphys = [tc.Box(:, 1).'; tc.Box(:, 2).'];
            tc.verifyEqual(nrm.toNormalized(Zphys), [0 0 0; tc.Weights], ...
                "RelTol", tc.RelTol, "AbsTol", tc.RelTol);
        end

        function normalizedBoxIsZeroToWeights(tc)
            nrm = tc.weighted();
            tc.verifyEqual(nrm.normalizedBox(), [zeros(3, 1), tc.Weights.']);
        end

        function dimensionsAreIndependent(tc)
            % The transform is affine PER DIMENSION: perturbing one physical
            % coordinate may move only its own normalized column.
            nrm = tc.weighted();
            base = mean(tc.Box, 2).';
            for kDim = 1:3
                perturbed = base;
                perturbed(kDim) = perturbed(kDim) + 0.1 * nrm.span(kDim);
                delta = nrm.toNormalized(perturbed) - nrm.toNormalized(base);
                tc.verifyEqual(delta([1:kDim - 1, kDim + 1:end]), ...
                    zeros(1, 2), "AbsTol", tc.RelTol, ...
                    sprintf("perturbing dimension %d moved another column", kDim));
                tc.verifyGreaterThan(abs(delta(kDim)), 0);
            end
        end

        % --- the two rules that make the ROI a ruler, not a fence -------------

        function roundTripIsExactOutsideTheUnitCube(tc)
            % Points far below the lower bound and far above the upper bound must
            % round-trip: the admissible region may be wider than the ROI.
            nrm = tc.weighted();
            lower = tc.Box(:, 1).';
            span  = (tc.Box(:, 2) - tc.Box(:, 1)).';
            Zphys = [lower - 3 * span;          % normalized ~ -3 * w
                 lower - 0.25 * span;
                 lower + 0.5 * span;        % inside
                 lower + 1.75 * span;
                 lower + 10 * span];        % normalized ~ 10 * w
            tc.verifyEqual(nrm.toPhysical(nrm.toNormalized(Zphys)), Zphys, "RelTol", tc.RelTol);
        end

        function roundTripIsExactInBothDirections(tc)
            % ... and starting from normalized coordinates, again outside [0,1].
            nrm = tc.weighted();
            Znorm = [-2.5 -0.3 0.4; 0.1 0.9 0.2; 7.0 4.5 -1.25] .* tc.Weights;
            tc.verifyEqual(nrm.toNormalized(nrm.toPhysical(Znorm)), Znorm, "RelTol", tc.RelTol);
        end

        function transformNeverClamps(tc)
            % Values below 0 and above the weighted upper edge are returned as they
            % are, in both directions. A clamp here would silently move samples.
            nrm = tc.weighted();
            below = tc.Box(:, 1).' - (tc.Box(:, 2) - tc.Box(:, 1)).';
            above = tc.Box(:, 2).' + (tc.Box(:, 2) - tc.Box(:, 1)).';
            Znorm = nrm.toNormalized([below; above]);
            tc.verifyEqual(Znorm(1, :), -tc.Weights, "RelTol", tc.RelTol);
            tc.verifyEqual(Znorm(2, :), 2 * tc.Weights, "RelTol", tc.RelTol);

            Zphys = nrm.toPhysical([-5 -5 -5; 5 5 5]);
            tc.verifyLessThan(Zphys(1, :), tc.Box(:, 1).');
            tc.verifyGreaterThan(Zphys(2, :), tc.Box(:, 2).');
        end

        function emptyPointSetRoundTrips(tc)
            % The acquired set is empty before the first sample is committed.
            nrm = tc.weighted();
            Zphys = zeros(0, 3);
            tc.verifySize(nrm.toNormalized(Zphys), [0 3]);
            tc.verifySize(nrm.toPhysical(Zphys), [0 3]);
        end

        function describeNamesEveryDimension(tc)
            s = tc.weighted().describe();
            tc.verifyClass(s, "string");
            for kDim = 1:3
                tc.verifyTrue(contains(s, tc.DimName(kDim)), ...
                    sprintf("describe() omits dimension name %s", tc.DimName(kDim)));
            end
        end

        % --- rejections ------------------------------------------------------

        function rejectsUpperBoundNotAboveLower(tc)
            tc.verifyError(@() rid.roi.Normalization([0 1; 2 2], ["a" "b"]), ...
                "rid:roi:invalidBox");
            tc.verifyError(@() rid.roi.Normalization([0 1; 3 2], ["a" "b"]), ...
                "rid:roi:invalidBox");
        end

        function rejectsNonFiniteBox(tc)
            tc.verifyError(@() rid.roi.Normalization([0 Inf], "a"), ...
                "rid:roi:invalidBox");
            tc.verifyError(@() rid.roi.Normalization([NaN 1], "a"), ...
                "rid:roi:invalidBox");
        end

        function rejectsEmptyBox(tc)
            tc.verifyError(@() rid.roi.Normalization(zeros(0, 2), strings(1, 0)), ...
                "rid:roi:invalidBox");
        end

        function rejectsNameCountMismatch(tc)
            tc.verifyError(@() rid.roi.Normalization(tc.Box, ["u" "y"]), ...
                "rid:roi:dimensionMismatch");
        end

        function rejectsWeightCountMismatch(tc)
            tc.verifyError(@() rid.roi.Normalization(tc.Box, tc.DimName, ...
                "dimWeights", [1 2]), "rid:roi:dimensionMismatch");
        end

        function rejectsNonPositiveOrNonFiniteWeights(tc)
            tc.verifyError(@() rid.roi.Normalization(tc.Box, tc.DimName, ...
                "dimWeights", [1 0 1]), "rid:roi:invalidWeights");
            tc.verifyError(@() rid.roi.Normalization(tc.Box, tc.DimName, ...
                "dimWeights", [1 -1 1]), "rid:roi:invalidWeights");
            tc.verifyError(@() rid.roi.Normalization(tc.Box, tc.DimName, ...
                "dimWeights", [1 Inf 1]), "rid:roi:invalidWeights");
        end

        function rejectsPointsWithWrongColumnCount(tc)
            nrm = tc.weighted();
            tc.verifyError(@() nrm.toNormalized([1 2]), "rid:roi:dimensionMismatch");
            tc.verifyError(@() nrm.toPhysical([1 2 3 4]), "rid:roi:dimensionMismatch");
        end

    end

end
