classdef SupportingPointsTest < matlab.unittest.TestCase
%SUPPORTINGPOINTSTEST  Analytic and property tests for rid.roi.supportingPoints.
%
%   Layers 1 (analytic cases) and 2 (property tests) of the project's test plan. Seeds are
%   literals and are reported in the failure diagnostics.
%
%   The grid tests pin the CELL-CENTRED convention of the cell-centred grid with exact values: a
%   grid including the endpoints would cover the ROI corners exactly and give a
%   different fill distance, so the convention has to be asserted, not assumed.

    properties (Constant)
        Seed      = 20260829
        OtherSeed = 20260830
        RelTol    = 1e-12
    end

    methods (Access = private)

        function d = seedDiagnostic(tc)
            d = sprintf("seed = %d", tc.Seed);
        end

        function nrm = unitInterval(~)
            nrm = rid.roi.Normalization([-4 12], "u");
        end

        function nrm = unitSquare(~)
            nrm = rid.roi.Normalization([0 1; -100 300], ["u" "y"]);
        end

        function nrm = weightedSquare(~)
            % Normalized ROI [0,2] x [0,1].
            nrm = rid.roi.Normalization([0 1; -100 300], ["u" "y"], ...
                "dimWeights", [2 1]);
        end

    end

    methods (Test)

        % --- grid: the convention, asserted exactly --------------------------

        function gridIsCellCentredIn1D(tc)
            % Analytic: nSup = 4 on [0,1] gives the cell centres 1/8, 3/8, 5/8, 7/8.
            % An endpoint-inclusive grid would give 0, 1/3, 2/3, 1.
            [Q, ~, info] = rid.roi.supportingPoints(tc.unitInterval(), 4, ...
                "method", "grid", ...
                "progress", "off");
            tc.verifyEqual(Q, [0.125; 0.375; 0.625; 0.875], "RelTol", tc.RelTol);
            tc.verifyEqual(info.nPerDim, 4);
        end

        function gridIsCellCentredIn2D(tc)
            % nSup = 9 in 2-D: the tensor product of {1/6, 1/2, 5/6} with itself.
            Q = rid.roi.supportingPoints(tc.unitSquare(), 9, "method", "grid", ...
                "progress", "off");
            tc.verifySize(Q, [9 2]);
            expectedAxis = [1/6; 1/2; 5/6];
            for kDim = 1:2
                tc.verifyEqual(unique(Q(:, kDim)), expectedAxis, "RelTol", tc.RelTol);
            end
            tc.verifyEqual(sortrows(Q), sortrows(combvec2(expectedAxis)), ...
                "RelTol", tc.RelTol);
        end

        function gridFillDistanceMatchesTheConvention(tc)
            % Analytic: sqrt(nDim) * h / 2 with h = 1/3, the half-diagonal of a cell.
            [~, fillDistanceQ] = rid.roi.supportingPoints(tc.unitSquare(), 9, ...
                "method", "grid", ...
                "progress", "off");
            tc.verifyEqual(fillDistanceQ, sqrt(2) / 6, "RelTol", tc.RelTol);
        end

        function gridFillDistanceWithDimWeightsIsTheAnisotropicHalfDiagonal(tc)
            % ROI [0,2] x [0,1], nSup = 4 -> 2 per dimension -> spacing [1, 0.5].
            [~, fillDistanceQ, info] = rid.roi.supportingPoints(tc.weightedSquare(), ...
                4, "method", "grid", ...
                "progress", "off");
            tc.verifyEqual(info.fillDistance.spacing, [1 0.5], "RelTol", tc.RelTol);
            tc.verifyEqual(fillDistanceQ, 0.5 * norm([1 0.5]), "RelTol", tc.RelTol);
        end

        function gridLiesInsideTheNormalizedRoi(tc)
            nrm = tc.weightedSquare();
            Q = rid.roi.supportingPoints(nrm, 25, "method", "grid", "progress", "off");
            tc.verifyGreaterThan(Q, 0);
            tc.verifyLessThan(Q, nrm.dimWeights);
        end

        function gridReportsAchievedCountWhenNSupIsNotAPerfectPower(tc)
            % nSup = 2000 in 3-D gives 12 per dimension, hence 1728 points. Both
            % counts are reported; the count is never silently adjusted, because a
            % smaller Q changes the measured fillDistanceQ.
            nrm = rid.roi.Normalization([0 1; 0 1; 0 1], ["u" "y" "dy"]);
            [Q, fillDistanceQ, info] = rid.roi.supportingPoints(nrm, 2000, ...
                "method", "grid", ...
                "progress", "off");
            tc.verifyEqual(info.nSupRequested, 2000);
            tc.verifyEqual(info.nPerDim, 12);
            tc.verifyEqual(info.nSup, 1728);
            tc.verifySize(Q, [1728 3]);
            tc.verifyEqual(fillDistanceQ, 0.5 * sqrt(3) / 12, "RelTol", tc.RelTol);
        end

        function gridNeverOverspendsTheRequestedBudget(tc)
            % Property: nSup is a budget. The achieved count is the largest
            % nPerDim^nDim that fits inside it, never more.
            nrm3 = rid.roi.Normalization([0 1; 0 1; 0 1], ["u" "y" "dy"]);
            nrm2 = tc.unitSquare();
            for nSup = [1 2 5 7 8 26 27 63 64 999 1000 1001 2000]
                [Q3, ~, info3] = rid.roi.supportingPoints(nrm3, nSup, "method", "grid", ...
                    "progress", "off");
                tc.verifyLessThanOrEqual(info3.nSup, nSup, ...
                    sprintf("3-D grid overspent the budget at nSup = %d", nSup));
                tc.verifyEqual(info3.nSup, size(Q3, 1));
                tc.verifyEqual(info3.nSup, info3.nPerDim^3);
                tc.verifyGreaterThan((info3.nPerDim + 1)^3, nSup, ...
                    sprintf("3-D grid left a whole layer unspent at nSup = %d", nSup));

                [~, ~, info2] = rid.roi.supportingPoints(nrm2, nSup, "method", "grid", ...
                    "progress", "off");
                tc.verifyLessThanOrEqual(info2.nSup, nSup, ...
                    sprintf("2-D grid overspent the budget at nSup = %d", nSup));
                tc.verifyGreaterThan((info2.nPerDim + 1)^2, nSup, ...
                    sprintf("2-D grid left a whole layer unspent at nSup = %d", nSup));
            end
        end

        function gridUsesTheWholeBudgetOnAPerfectPower(tc)
            % nPerDim comes from an exact integer test, not floor(nSup^(1/nDim)):
            % 1000^(1/3) evaluates to 9.999999999999998 and would lose a layer.
            nrm = rid.roi.Normalization([0 1; 0 1; 0 1], ["u" "y" "dy"]);
            [Q, ~, info] = rid.roi.supportingPoints(nrm, 1000, "method", "grid", ...
                "progress", "off");
            tc.verifyEqual(info.nPerDim, 10);
            tc.verifyEqual(info.nSup, 1000);
            tc.verifySize(Q, [1000 3]);
        end

        function gridFillDistanceHalvesWhenTheResolutionDoubles(tc)
            % Property: the fill distance of a cell-centred grid scales with 1/n.
            nrm = tc.unitSquare();
            [~, coarse] = rid.roi.supportingPoints(nrm, 4, "method", "grid", ...
                "progress", "off");
            [~, fine]   = rid.roi.supportingPoints(nrm, 16, "method", "grid", ...
                "progress", "off");
            tc.verifyEqual(fine, coarse / 2, "RelTol", tc.RelTol);
        end

        function gridIsSeedIndependent(tc)
            nrm = tc.unitSquare();
            a = rid.roi.supportingPoints(nrm, 16, "method", "grid", "seed", tc.Seed, ...
                "progress", "off");
            b = rid.roi.supportingPoints(nrm, 16, "method", "grid", "seed", tc.OtherSeed, ...
                "progress", "off");
            tc.verifyEqual(a, b);
        end

        function gridFillDistanceIsReportedAsAnalytic(tc)
            [~, ~, info] = rid.roi.supportingPoints(tc.unitSquare(), 9, ...
                "method", "grid", ...
                "progress", "off");
            tc.verifyEqual(info.fillDistance.method, "gridAnalytic");
        end

        % --- lhs -------------------------------------------------------------

        function lhsIsTheDefaultMethod(tc)
            [~, ~, info] = rid.roi.supportingPoints(tc.unitSquare(), 20, ...
                "nProbes", 500, ...
                "progress", "off");
            tc.verifyEqual(info.method, "lhs");
        end

        function lhsHasRequestedShapeAndLiesInsideTheNormalizedRoi(tc)
            nrm = tc.weightedSquare();
            [Q, ~, info] = rid.roi.supportingPoints(nrm, 30, ...
                "seed", tc.Seed, "nProbes", 500, ...
                "progress", "off");
            tc.verifySize(Q, [30 2], tc.seedDiagnostic());
            tc.verifyEqual(info.nSup, 30);
            tc.verifyEqual(info.nSupRequested, 30);
            tc.verifyGreaterThan(Q, 0, tc.seedDiagnostic());
            tc.verifyLessThan(Q, nrm.dimWeights, tc.seedDiagnostic());
        end

        function lhsPreservesTheLatinPropertyAfterWeighting(tc)
            % Scaling by dimWeights is per dimension, so each column still holds
            % exactly one point per stratum.
            nSup = 12;
            nrm  = tc.weightedSquare();
            Q = rid.roi.supportingPoints(nrm, nSup, "seed", tc.Seed, "nProbes", 500, ...
                "progress", "off");
            for kDim = 1:2
                stratum = floor(Q(:, kDim) / nrm.dimWeights(kDim) * nSup);
                tc.verifyEqual(sort(stratum), (0:nSup - 1).', tc.seedDiagnostic());
            end
        end

        function lhsIsDeterministicUnderAFixedSeed(tc)
            nrm = tc.unitSquare();
            a = rid.roi.supportingPoints(nrm, 25, "seed", tc.Seed, "nProbes", 500, ...
                "progress", "off");
            b = rid.roi.supportingPoints(nrm, 25, "seed", tc.Seed, "nProbes", 500, ...
                "progress", "off");
            tc.verifyEqual(a, b, tc.seedDiagnostic());
        end

        function lhsDiffersAcrossSeeds(tc)
            nrm = tc.unitSquare();
            a = rid.roi.supportingPoints(nrm, 25, "seed", tc.Seed, "nProbes", 500, ...
                "progress", "off");
            b = rid.roi.supportingPoints(nrm, 25, "seed", tc.OtherSeed, "nProbes", 500, ...
                "progress", "off");
            tc.verifyNotEqual(a, b, ...
                sprintf("seeds = %d, %d", tc.Seed, tc.OtherSeed));
        end

        function lhsFillDistanceEqualsADirectEstimatorCall(tc)
            % The measurement is delegated, not reimplemented: same probe seed and
            % probe count must give the identical number.
            nrm = tc.weightedSquare();
            [Q, fillDistanceQ, info] = rid.roi.supportingPoints(nrm, 40, ...
                "seed", tc.Seed, "probeSeed", 11, "nProbes", 5000, ...
                "progress", "off");
            direct = rid.roi.estimateFillDistance(Q, nrm, ...
                "seed", 11, "nProbes", 5000);
            tc.verifyEqual(fillDistanceQ, direct, tc.seedDiagnostic());
            tc.verifyEqual(info.fillDistance.method, "monteCarlo");
            tc.verifyEqual(info.fillDistance.seed, 11);
            tc.verifyEqual(info.fillDistance.nProbes, 5000);
        end

        function lhsFillDistanceCarriesTheConvergenceTrace(tc)
            % The trace counts STOCHASTIC probes: under the
            % seeded default the 2^nDim + 2*nDim deterministic probes are excluded,
            % because they sit at the deepest holes and would flatten the trace by
            % construction. In 2-D that is 8 of the 4000.
            [~, ~, info] = rid.roi.supportingPoints(tc.unitSquare(), 40, ...
                "seed", tc.Seed, "nProbes", 4000, ...
                "progress", "off");
            nStochastic = 4000 - info.fillDistance.nSeededProbes;
            tc.verifyEqual(info.fillDistance.nSeededProbes, 8);
            tc.verifyEqual(info.fillDistance.trace.nProbes, ...
                [ceil(nStochastic / 4) ceil(nStochastic / 2) nStochastic]);
            tc.verifySize(info.fillDistance.trace.value, [1 3]);
        end

        function moreSupportingPointsGiveASmallerFillDistance(tc)
            % Property, with the probe set held fixed so the comparison is of the
            % designs and not of two different probe draws.
            nrm = tc.unitSquare();
            [~, coarse] = rid.roi.supportingPoints(nrm, 20, ...
                "seed", tc.Seed, "probeSeed", 7, "nProbes", 20000, ...
                "progress", "off");
            [~, fine]   = rid.roi.supportingPoints(nrm, 200, ...
                "seed", tc.Seed, "probeSeed", 7, "nProbes", 20000, ...
                "progress", "off");
            tc.verifyLessThan(fine, coarse, tc.seedDiagnostic());
        end

        % --- the improving methods -----------------------------------

        function improvedMethodsPermuteTheLhsDesignWithinEachDimension(tc)
            % What makes it legitimate to improve the ALREADY-SCALED design: a
            % swap exchanges two points' coordinates within ONE dimension, so
            % every column keeps its multiset of values. The Latin property and
            % the box therefore survive untouched.
            nrm  = tc.weightedSquare();
            nSup = 60;
            base = rid.roi.supportingPoints(nrm, nSup, "seed", tc.Seed, ...
                "nProbes", 500, "progress", "off");
            for method = ["lhsImproved" "lhsImprovedFull"]
                Q = rid.roi.supportingPoints(nrm, nSup, "method", method, ...
                    "seed", tc.Seed, "nProbes", 500, "progress", "off");
                tc.verifySize(Q, [nSup 2], tc.seedDiagnostic());
                tc.verifyEqual(sort(Q), sort(base), ...
                    sprintf("%s: columns are not a permutation of the lhs design", method));
                tc.verifyGreaterThan(Q, 0, tc.seedDiagnostic());
                tc.verifyLessThan(Q, nrm.dimWeights, tc.seedDiagnostic());
            end
        end

        function improvedMethodsRaiseTheMinimumInterPointDistance(tc)
            % The criterion the external function actually maximizes (maximin /
            % packing), measured in the dimWeights-scaled coordinates it was
            % given. This is the contract, and deliberately NOT an ordering on
            % fillDistanceQ: that is a covering radius, a different functional,
            % and the chain through "lhsImproved" was measured to be non-monotone
            % in it (the Latin-hypercube method option). Do not "restore" the
            % fill-distance ordering here -- it has counterexamples at nSup 50
            % and 100 in this very fixture family.
            nrm  = tc.weightedSquare();
            nSup = 100;
            minDist = @(method) min(pdist(rid.roi.supportingPoints(nrm, nSup, ...
                "method", method, "seed", tc.Seed, "nProbes", 500, ...
                "progress", "off")));
            plain    = minDist("lhs");
            improved = minDist("lhsImproved");
            full     = minDist("lhsImprovedFull");
            tc.verifyGreaterThan(improved, plain, tc.seedDiagnostic());
            tc.verifyGreaterThanOrEqual(full, improved, tc.seedDiagnostic());
        end

        function improvementMaximizesTheDimWeightedMetricNotTheUnitCube(tc)
            % The Latin-hypercube method option's anisotropy question, settled empirically and then pinned:
            % improving the scaled design beats improving the unit cube and
            % scaling afterwards, in the metric +crit measures in.
            nrm  = tc.weightedSquare();
            nSup = 60;
            Q = rid.roi.supportingPoints(nrm, nSup, "method", "lhsImproved", ...
                "seed", tc.Seed, "nProbes", 500, "progress", "off");

            unitCube = rid.ext.latinHypercube(nSup, 2, tc.Seed);
            isotropic = rid.ext.improveLatinHypercube(unitCube, ...
                'plot', 0, 'display', false, 'normalize', false, ...
                'maxCalculationTime', 60) .* nrm.dimWeights;

            tc.verifyGreaterThan(min(pdist(Q)), min(pdist(isotropic)), ...
                tc.seedDiagnostic());
        end

        function improvedMethodsAreDeterministic(tc)
            % The improvement contains no randomness, so the Latin-hypercube seed handling's seed discipline
            % is untouched: same seed, same X1, same Q.
            nrm = tc.unitSquare();
            args = {"method", "lhsImprovedFull", "seed", tc.Seed, ...
                    "nProbes", 500, "progress", "off"};
            a = rid.roi.supportingPoints(nrm, 40, args{:}, "progress", "off");
            b = rid.roi.supportingPoints(nrm, 40, args{:}, "progress", "off");
            tc.verifyEqual(a, b, tc.seedDiagnostic());
        end

        function infoReportsTheImprovementAndItsCap(tc)
            [~, ~, info] = rid.roi.supportingPoints(tc.unitSquare(), 60, ...
                "method", "lhsImproved", "seed", tc.Seed, "nProbes", 500, ...
                "progress", "off");
            tc.verifyEqual(info.method, "lhsImproved");
            tc.verifyGreaterThan(info.improveElapsed, 0);
            tc.verifyClass(info.timeCapHit, "logical");
            tc.verifyFalse(info.timeCapHit, ...
                "the default cap must not bite on a 60-point 2-D design");
            tc.verifyEqual(info.improveMetric, "dimWeighted");
        end

        function infoReportsNoImprovementForTheNonImprovingMethods(tc)
            % The fields are always present, so a run record never has to branch
            % on which method built Q.
            for method = ["lhs" "grid"]
                [~, ~, info] = rid.roi.supportingPoints(tc.unitSquare(), 25, ...
                    "method", method, "seed", tc.Seed, "nProbes", 500, ...
                    "progress", "off");
                tc.verifyEqual(info.improveElapsed, 0, method);
                tc.verifyFalse(info.timeCapHit, method);
                tc.verifyEqual(info.improveMetric, "none", method);
            end
        end

        function aTimeCapThatBitesIsReportedNeverSilent(tc)
            % A cap that bites turns a deterministic design into a
            % machine-dependent one, so it must show up in info.
            [Q, ~, info] = rid.roi.supportingPoints(tc.weightedSquare(), 200, ...
                "method", "lhsImprovedFull", "seed", tc.Seed, "nProbes", 500, ...
                "improveMaxTime", 1e-6, "progress", "off");
            tc.verifyTrue(info.timeCapHit, tc.seedDiagnostic());
            tc.verifySize(Q, [200 2]);
            tc.verifyGreaterThan(Q, 0, tc.seedDiagnostic());
            tc.verifyLessThan(Q, tc.weightedSquare().dimWeights, tc.seedDiagnostic());
        end

        function theDefaultMethodStillBuildsTheUnimprovedLhsDesign(tc)
            % the addition is purely additive. Nothing frozen moves only
            % if the default and an explicit "lhs" are bit-identical.
            nrm = tc.weightedSquare();
            tc.verifyEqual(rid.defaultOptions().roi.method, "lhs");
            [Qdefault, fDefault] = rid.roi.supportingPoints(nrm, 40, ...
                "seed", tc.Seed, "nProbes", 500, "progress", "off");
            [Qlhs, fLhs] = rid.roi.supportingPoints(nrm, 40, "method", "lhs", ...
                "seed", tc.Seed, "nProbes", 500, "progress", "off");
            tc.verifyEqual(Qdefault, Qlhs, tc.seedDiagnostic());
            tc.verifyEqual(fDefault, fLhs, tc.seedDiagnostic());
        end

        % --- the liveness signal --------------------------------------

        function progressAutoIsSilentOnTheFastMethods(tc)
            % "lhs" and "grid" are fast enough that a message there is noise.
            % Captured, not eyeballed.
            printed = evalc(['rid.roi.supportingPoints(tc.unitSquare(), 40, ' ...
                '"method", "lhs", "seed", tc.Seed, "nProbes", 500, ' ...
                '"progress", "auto");']);
            tc.verifyEmpty(printed, "method ""lhs"" must print nothing");
            printed = evalc(['rid.roi.supportingPoints(tc.unitSquare(), 25, ' ...
                '"method", "grid", "progress", "auto");']);
            tc.verifyEmpty(printed, "method ""grid"" must print nothing");
        end

        function progressAutoPrintsForTheImprovingMethods(tc)
            printed = evalc(['rid.roi.supportingPoints(tc.unitSquare(), 60, ' ...
                '"method", "lhsImprovedFull", "seed", tc.Seed, ' ...
                '"nProbes", 500, "improveMaxTime", 60, "progress", "auto");']);
            tc.verifyNotEmpty(printed);
            % What the line has to say, because its purpose fixes its content.
            tc.verifySubstring(printed, "supporting points");
            tc.verifySubstring(printed, "not a stall");
            tc.verifySubstring(printed, "BEFORE the first input step");
            tc.verifySubstring(printed, "nSup = 60");
            tc.verifySubstring(printed, "nDim = 2");
            tc.verifySubstring(printed, "lhsImprovedFull");
            tc.verifySubstring(printed, "time cap 60 s");
            tc.verifySubstring(printed, "lhsImproved"" or ""lhs");
            tc.verifySubstring(printed, "fillDistanceQ =");
            tc.verifySubstring(printed, "the time cap was not hit");
        end

        function progressAutoReportsACapThatBit(tc)
            printed = evalc(['rid.roi.supportingPoints(tc.weightedSquare(), 200, ' ...
                '"method", "lhsImprovedFull", "seed", tc.Seed, ' ...
                '"nProbes", 500, "improveMaxTime", 1e-6, "progress", "auto");']);
            tc.verifySubstring(printed, "TIME CAP WAS HIT");
            tc.verifySubstring(printed, "machine-dependent");
        end

        function progressOffSilencesTheImprovingMethods(tc)
            printed = evalc(['rid.roi.supportingPoints(tc.unitSquare(), 40, ' ...
                '"method", "lhsImproved", "seed", tc.Seed, "nProbes", 500, ' ...
                '"progress", "off");']);
            tc.verifyEmpty(printed);
        end

        function rejectsUnknownProgressSetting(tc)
            tc.verifyError(@() rid.roi.supportingPoints(tc.unitSquare(), 10, ...
                "progress", "loud"), "rid:roi:invalidProgress");
        end

        % --- provenance and boundaries ---------------------------------------

        function infoRecordsGeneratorSeedAndDimension(tc)
            [~, ~, info] = rid.roi.supportingPoints(tc.weightedSquare(), 20, ...
                "seed", tc.Seed, "nProbes", 500, ...
                "progress", "off");
            tc.verifyEqual(info.seed, tc.Seed);
            tc.verifyEqual(info.nDim, 2);
        end

        function supportingPointsMapBackInsideThePhysicalBox(tc)
            % Q discretizes the ROI, so it must land inside the declared box.
            nrm = tc.weightedSquare();
            Q = rid.roi.supportingPoints(nrm, 50, "seed", tc.Seed, "nProbes", 500, ...
                "progress", "off");
            Zphys = nrm.toPhysical(Q);
            tc.verifyGreaterThan(Zphys, nrm.box(:, 1).', tc.seedDiagnostic());
            tc.verifyLessThan(Zphys, nrm.box(:, 2).', tc.seedDiagnostic());
        end

        function singleSupportingPointIsLegal(tc)
            [Q, fillDistanceQ] = rid.roi.supportingPoints(tc.unitInterval(), 1, ...
                "method", "grid", ...
                "progress", "off");
            tc.verifyEqual(Q, 0.5, "RelTol", tc.RelTol);
            tc.verifyEqual(fillDistanceQ, 0.5, "RelTol", tc.RelTol);
        end

        function rejectsNonPositiveOrFractionalNSup(tc)
            nrm = tc.unitSquare();
            tc.verifyError(@() rid.roi.supportingPoints(nrm, 0, "progress", "off"), ...
                "rid:roi:invalidSupportingPointCount");
            tc.verifyError(@() rid.roi.supportingPoints(nrm, -5, "progress", "off"), ...
                "rid:roi:invalidSupportingPointCount");
            tc.verifyError(@() rid.roi.supportingPoints(nrm, 12.5, "progress", "off"), ...
                "rid:roi:invalidSupportingPointCount");
        end

        function rejectsUnknownMethod(tc)
            tc.verifyError(@() rid.roi.supportingPoints(tc.unitSquare(), 10, ...
                "method", "sobol", ...
                "progress", "off"), "rid:roi:unknownMethod");
        end

    end

end

% -----------------------------------------------------------------------------

function P = combvec2(axisCoords)
%COMBVEC2  Tensor product of one axis with itself, as (n^2 x 2) points.
[A, B] = ndgrid(axisCoords, axisCoords);
P = [A(:), B(:)];
end
