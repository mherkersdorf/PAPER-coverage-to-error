classdef EstimateFillDistanceTest < matlab.unittest.TestCase
%ESTIMATEFILLDISTANCETEST  Analytic and property tests for rid.roi.estimateFillDistance.
%
%   Layers 1 (analytic cases) and 2 (property tests) of the project's test plan. Seeds are
%   literals and are reported in the failure diagnostics.
%
%   The estimator's defining property is that it is BIASED LOW: a maximum
%   over finitely many probes of a supremum can only under-report. Two tests assert
%   that direction on cases whose true value is known analytically -- an
%   over-reporting estimator would be a bug, an under-reporting one is the
%   documented behaviour.

    properties (Constant)
        Seed      = 20260829
        OtherSeed = 20260830
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

        function nrm = unitInterval(~)
            % One dimension, normalized ROI [0,1].
            nrm = rid.roi.Normalization([0 1], "u");
        end

        function nrm = unitSquare(~)
            nrm = rid.roi.Normalization([0 1; 0 1], ["u" "y"]);
        end

        function nrm = weightedSquare(~)
            % Normalized ROI [0,2] x [0,1].
            nrm = rid.roi.Normalization([0 1; 0 1], ["u" "y"], "dimWeights", [2 1]);
        end

        function [Q, h] = cellCentredGrid(~, nPerDim, extent)
            % Cell-centred grid on [0,extent(1)] x ... ; the convention.
            nDim = numel(extent);
            axisCoords = cell(1, nDim);
            for kDim = 1:nDim
                axisCoords{kDim} = extent(kDim) * ((0.5:nPerDim - 0.5).' / nPerDim);
            end
            cols = cell(1, nDim);
            [cols{:}] = ndgrid(axisCoords{:});
            Q = zeros(nPerDim^nDim, nDim);
            for kDim = 1:nDim
                Q(:, kDim) = cols{kDim}(:);
            end
            h = extent ./ nPerDim;
        end

    end

    methods (Test)

        % --- analytic cases --------------------------------------------------

        function singlePointInUnitIntervalApproachesOneHalf(tc)
            % Analytic: with P = 0.5 on [0,1] the supremum is 0.5, attained at both
            % ends. The estimate must approach it from below and never exceed it.
            nrm = tc.unitInterval();
            est = rid.roi.estimateFillDistance(0.5, nrm, ...
                "nProbes", 1e5, "seed", tc.Seed);
            tc.verifyLessThanOrEqual(est, 0.5, tc.seedDiagnostic());
            tc.verifyGreaterThan(est, 0.4995, tc.seedDiagnostic());
        end

        function cornerPointRecoversTheNormalizedExtent(tc)
            % P = {0} on [0,2]: the supremum is 2, at the far end. This is what
            % fails if the probes are drawn on the wrong extent -- the unit cube
            % instead of the weighted box, or the physical bounds.
            nrm = rid.roi.Normalization([-7 13], "u", "dimWeights", 2);
            est = rid.roi.estimateFillDistance(0, nrm, ...
                "nProbes", 1e5, "seed", tc.Seed);
            tc.verifyLessThanOrEqual(est, 2, tc.seedDiagnostic());
            tc.verifyGreaterThan(est, 2 - 1e-3, tc.seedDiagnostic());
        end

        function neverExceedsTheAnalyticGridValue(tc)
            % The bias is one-directional, in 1-D and in 2-D. sqrt(nDim)*h/2 is the
            % half-diagonal of a cell, the exact supremum for a cell-centred grid.
            %
            % A FEW ULPS OF SLACK, and why it is not a weakened contract. Under the
            % seeded default a probe now sits EXACTLY at the maximizer, so
            % the estimate attains the supremum instead of approaching it, and the
            % comparison becomes one between two floating-point spellings of the
            % same real number: pdist2's sqrt(x^2 + y^2) against 0.5*sqrt(2)/6,
            % which differ by 3 ulps here. The one-directional claim is unchanged;
            % only exact-equality-in-binary was never part of it.
            nrm1  = tc.unitInterval();
            truth1 = 0.5 * (1 / 8);
            est1  = rid.roi.estimateFillDistance(tc.cellCentredGrid(8, 1), nrm1, ...
                "nProbes", 5e4, "seed", tc.Seed);
            tc.verifyLessThanOrEqual(est1, truth1 + 8 * eps(truth1), ...
                tc.seedDiagnostic());

            nrm2   = tc.unitSquare();
            truth2 = 0.5 * sqrt(2) / 6;
            est2   = rid.roi.estimateFillDistance(tc.cellCentredGrid(6, [1 1]), ...
                nrm2, "nProbes", 5e4, "seed", tc.Seed);
            tc.verifyLessThanOrEqual(est2, truth2 + 8 * eps(truth2), ...
                tc.seedDiagnostic());
        end

        function twoDimensionalGridEstimateIsCloseToTheAnalyticValue(tc)
            % Biased low, but not by much once the probe set is large.
            nrm   = tc.unitSquare();
            Q     = tc.cellCentredGrid(5, [1 1]);
            truth = 0.5 * sqrt(2) / 5;
            est   = rid.roi.estimateFillDistance(Q, nrm, ...
                "nProbes", 2e5, "seed", tc.Seed);
            tc.verifyGreaterThan(est, 0.97 * truth, tc.seedDiagnostic());
        end

        function respectsDimWeights(tc)
            % Anisotropic ROI [0,2] x [0,1] with anisotropic spacing: the truth is
            % the half-diagonal of a cell, 0.5 * norm([2 1]./4).
            nrm     = tc.weightedSquare();
            [Q, h]  = tc.cellCentredGrid(4, [2 1]);
            truth   = 0.5 * norm(h);
            est     = rid.roi.estimateFillDistance(Q, nrm, ...
                "nProbes", 2e5, "seed", tc.Seed);
            tc.verifyLessThanOrEqual(est, truth, tc.seedDiagnostic());
            tc.verifyGreaterThan(est, 0.97 * truth, tc.seedDiagnostic());
        end

        % --- the convergence trace -------------------------------------------

        function traceIsNonDecreasing(tc)
            % A running maximum over a growing probe set cannot fall.
            nrm = tc.unitSquare();
            [~, info] = rid.roi.estimateFillDistance(tc.cellCentredGrid(4, [1 1]), ...
                nrm, "nProbes", 4000, "seed", tc.Seed, "probeMethod", "sobol");
            tc.verifyEqual(sort(info.trace.value), info.trace.value, ...
                tc.seedDiagnostic());
        end

        function traceLastEntryEqualsTheReturnedValue(tc)
            nrm = tc.unitSquare();
            [est, info] = rid.roi.estimateFillDistance(tc.cellCentredGrid(4, [1 1]), ...
                nrm, "nProbes", 4000, "seed", tc.Seed, "probeMethod", "sobol");
            tc.verifyEqual(info.trace.value(end), est, tc.seedDiagnostic());
        end

        function traceProbeCountsAreQuarterHalfAndAll(tc)
            nrm = tc.unitSquare();
            [~, info] = rid.roi.estimateFillDistance([0.5 0.5], nrm, ...
                "nProbes", 4000, "seed", tc.Seed, "probeMethod", "sobol");
            tc.verifyEqual(info.trace.nProbes, [1000 2000 4000]);
            tc.verifySize(info.trace.value, [1 3]);
        end

        function traceEntriesMatchIndependentRunsOfTheSameLength(tc)
            % The quarter- and half-set entries are the values a shorter run of the
            % same seed would report: the probe sequence is a prefix, not a resample.
            nrm = tc.unitSquare();
            P   = tc.cellCentredGrid(4, [1 1]);
            [~, info] = rid.roi.estimateFillDistance(P, nrm, ...
                "nProbes", 4000, "seed", tc.Seed, "probeMethod", "sobol");
            quarter = rid.roi.estimateFillDistance(P, nrm, ...
                "nProbes", 1000, "seed", tc.Seed, "probeMethod", "sobol");
            half    = rid.roi.estimateFillDistance(P, nrm, ...
                "nProbes", 2000, "seed", tc.Seed, "probeMethod", "sobol");
            tc.verifyEqual(info.trace.value(1), quarter, tc.seedDiagnostic());
            tc.verifyEqual(info.trace.value(2), half, tc.seedDiagnostic());
        end

        % --- the probe set, and the bracket ---------------------------

        function theDefaultProbeSetCoversBetterThanIidUniform(tc)
            % What the probe bracket actually buys, measured on the quantity it is bought for:
            % the probe set's OWN fill distance, which is the width of the bracket.
            % A gap in the probe set is where the deepest hole of P hides.
            nrm = rid.roi.Normalization([0 1; 0 1; 0 1; 0 1], ["a" "b" "c" "d"]);
            Q   = tc.cellCentredGrid(5, [1 1 1 1]);
            [~, sobol]   = rid.roi.estimateFillDistance(Q, nrm, ...
                "nProbes", 4000, "seed", tc.Seed, "probeMethod", "sobol");
            [~, uniform] = rid.roi.estimateFillDistance(Q, nrm, ...
                "nProbes", 4000, "seed", tc.Seed, "probeMethod", "uniform");
            tc.verifyEqual(sobol.probeMethod, "sobol");
            tc.verifyLessThan(sobol.probeFillDistance, ...
                uniform.probeFillDistance, tc.seedDiagnostic());
        end

        function theBracketIsTheEstimateAndTheProbeSetsOwnFillDistance(tc)
            nrm = tc.unitSquare();
            [est, info] = rid.roi.estimateFillDistance(tc.cellCentredGrid(4, [1 1]), ...
                nrm, "nProbes", 4000, "seed", tc.Seed);
            tc.verifyGreaterThan(info.probeFillDistance, 0, tc.seedDiagnostic());
            tc.verifyEqual(info.bracket, ...
                [est, est + info.probeFillDistance], tc.seedDiagnostic());
        end

        function theBracketContainsTheAnalyticTruth(tc)
            % The 1-Lipschitz argument: a probe set of fill distance
            % hProbe cannot miss the supremum by more than hProbe. Asserted on the
            % cases whose supremum is known exactly.
            nrm3 = rid.roi.Normalization([0 1; 0 1; 0 1], ["a" "b" "c"]);
            [~, info3] = rid.roi.estimateFillDistance(tc.cellCentredGrid(6, [1 1 1]), ...
                nrm3, "nProbes", 1e4, "seed", tc.Seed);
            tc.verifyGreaterThanOrEqual(info3.bracket(2), 0.5 * sqrt(3) / 6, ...
                tc.seedDiagnostic());

            nrm2 = tc.unitSquare();
            [~, info2] = rid.roi.estimateFillDistance(tc.cellCentredGrid(5, [1 1]), ...
                nrm2, "nProbes", 1e4, "seed", tc.Seed);
            tc.verifyGreaterThanOrEqual(info2.bracket(2), 0.5 * sqrt(2) / 5, ...
                tc.seedDiagnostic());
        end

        function theUniformPathIsTheProbeSequenceItClaimsToBe(tc)
            % The i.i.d. path is still selectable and still draws the probes in
            % probe order, one probe's coordinates at a time. The oracle is built
            % here, independently of the implementation's blocking.
            nrm     = tc.weightedSquare();
            Q       = tc.cellCentredGrid(4, [2 1]);
            nProbes = 3000;

            stream = RandStream("twister", "Seed", tc.Seed);
            probes = rand(stream, 2, nProbes).' .* [2 1];
            expected = max(min(rid.ext.pairwiseDistance(probes, Q), [], 2));

            est = rid.roi.estimateFillDistance(Q, nrm, "nProbes", nProbes, ...
                "seed", tc.Seed, "probeMethod", "uniform");
            tc.verifyEqual(est, expected, tc.seedDiagnostic());
        end

        % --- corner seeding, as apparatus for corner seeding's measurement --------------

        function cornerSeedingIsExactOnACellCentredGridAndThatIsTautological(tc)
            % The supremum of a cell-centred grid is attained AT A CELL CORNER, and
            % the ROI's own vertex (0,...,0) is one: it sits at the cell
            % half-diagonal from the first cell's centre. A probe set containing
            % the vertices therefore returns the analytic value to machine
            % precision -- by construction, not by covering better. This test
            % exists to pin that fact where a reader of the measurement will
            % meet it, because it is what stops "100 % of the truth at every
            % dimension" from being read as evidence of coverage.
            nDim = 5;
            nrm  = rid.roi.Normalization([zeros(nDim, 1) ones(nDim, 1)], ...
                "d" + string(1:nDim));
            Q    = tc.cellCentredGrid(3, ones(1, nDim));
            est  = rid.roi.estimateFillDistance(Q, nrm, ...
                "nProbes", 4000, "seed", tc.Seed, "probeMethod", "sobolCorners");
            tc.verifyEqual(est, sqrt(nDim) / 6, "RelTol", 1e-12, ...
                tc.seedDiagnostic());
        end

        function cornerSeedingFindsWhatTheOtherTwoLeaveOnAnLhsDesign(tc)
            % The case that is NOT tautological: an lhs design places no point at
            % a vertex, so the corners are genuinely deep holes and no probe set
            % that misses them can report them. The estimator is biased low in one
            % direction, so a larger estimate is a better one.
            nDim = 4;
            nrm  = rid.roi.Normalization([zeros(nDim, 1) ones(nDim, 1)], ...
                "d" + string(1:nDim));
            Q    = rid.roi.supportingPoints(nrm, 500, "method", "lhs", "seed", 0);

            plain   = rid.roi.estimateFillDistance(Q, nrm, "nProbes", 4000, ...
                "seed", tc.Seed, "probeMethod", "sobol");
            seeded  = rid.roi.estimateFillDistance(Q, nrm, "nProbes", 4000, ...
                "seed", tc.Seed, "probeMethod", "sobolCorners");
            tc.verifyGreaterThan(seeded, plain, tc.seedDiagnostic());
        end

        function theCeilingIsReportedAndFallsBackToSobol(tc)
            % 2^nDim corners are free to about nDim = 12 and unaffordable beyond.
            % Above nProbes the method falls back and SAYS SO -- the cell-centred grid
            % precedent, where a count is adjusted, both values are returned and
            % nothing is silent. Refusing would make the default unusable in high
            % dimension for a reason the caller cannot act on.
            nDim = 8;                                  % 2^8 + 16 = 272 seeds
            nrm  = rid.roi.Normalization([zeros(nDim, 1) ones(nDim, 1)], ...
                "d" + string(1:nDim));
            Q    = rid.roi.supportingPoints(nrm, 200, "method", "lhs", "seed", 0);

            [est, info] = rid.roi.estimateFillDistance(Q, nrm, ...
                "nProbes", 100, "seed", tc.Seed, "probeMethod", "sobolCorners");

            tc.verifyTrue(info.probeSeedingFellBack, tc.seedDiagnostic());
            tc.verifyEqual(info.probeMethodRequested, "sobolCorners");
            tc.verifyEqual(info.probeMethod, "sobol");
            tc.verifyEqual(info.nSeededProbes, 0);
            tc.verifyEqual(info.nStochasticProbes, 100);

            % NOT TRUNCATED: the fallback is exactly the "sobol" probe set, not
            % the first 100 of a 272-corner list. An arbitrary subset of the
            % vertices is a different probe set, not a smaller one.
            plainSobol = rid.roi.estimateFillDistance(Q, nrm, ...
                "nProbes", 100, "seed", tc.Seed, "probeMethod", "sobol");
            tc.verifyEqual(est, plainSobol, tc.seedDiagnostic());
        end

        function belowTheCeilingNothingFallsBack(tc)
            % The other half of the report: when the seeds fit, both names agree
            % and the count of them is stated rather than implied.
            nDim = 3;                                   % 2^3 + 6 = 14 seeds
            nrm  = rid.roi.Normalization([zeros(nDim, 1) ones(nDim, 1)], ...
                "d" + string(1:nDim));
            [~, info] = rid.roi.estimateFillDistance(tc.cellCentredGrid(4, ones(1, nDim)), ...
                nrm, "nProbes", 1000, "seed", tc.Seed, "probeMethod", "sobolCorners");

            tc.verifyFalse(info.probeSeedingFellBack, tc.seedDiagnostic());
            tc.verifyEqual(info.probeMethod, "sobolCorners");
            tc.verifyEqual(info.probeMethodRequested, "sobolCorners");
            tc.verifyEqual(info.nSeededProbes, 14);
            tc.verifyEqual(info.nStochasticProbes, 986);
        end

        function cornerSeedingIsTheDefault(tc)
            % Corner seeding: adopted at every dimension.
            tc.verifyEqual(rid.defaultOptions().roi.probeMethod, "sobolCorners");
        end

        function theTraceRunsOverTheStochasticProbesOnly(tc)
            % THE DEFECT THIS GUARDS. The seeds are probed first and sit exactly
            % where the deepest holes are, so a trace over all probes is flat from
            % the first quarter onwards and `converged` is true for free -- the
            % convergence-trace safeguard retired by the feature.
            %
            % MEASURED ON THIS FIXTURE, at this seed: the all-probes trace reads
            % [49 98 196] -> climb [0 0], CONVERGED; the stochastic trace reads
            % [30 60 120] -> climb [8.0 %, 11.1 %], UNCONVERGED. The seeds carry
            % the maximum, so the old rule certified a probe count of 120 random
            % probes in 6-D as sufficient -- which it plainly is not.
            %
            % A tiny stochastic portion is NOT enough on its own to make this
            % visible: a running maximum is a step function, and over 20 or 40
            % stochastic probes it is flat by luck rather than by construction.
            % 120 is the smallest count in the scan at which the climb is well
            % clear of the tolerance in both increments.
            nDim    = 6;                                % 2^6 + 12 = 76 seeds
            nProbes = 196;                              % ... and 120 stochastic
            nrm     = rid.roi.Normalization([zeros(nDim, 1) ones(nDim, 1)], ...
                "d" + string(1:nDim));
            Q       = rid.roi.supportingPoints(nrm, 300, "method", "lhs", "seed", 0);

            [est, info] = rid.roi.estimateFillDistance(Q, nrm, ...
                "nProbes", nProbes, "seed", tc.Seed, "probeMethod", "sobolCorners");

            tc.verifyEqual(info.nSeededProbes, 76);
            tc.verifyEqual(info.nStochasticProbes, 120);
            tc.verifyEqual(info.trace.nProbes, [30 60 120], tc.seedDiagnostic());
            tc.verifyFalse(info.converged, tc.seedDiagnostic());
            tc.verifyGreaterThan(max(info.trace.relativeClimb), 0.02, ...
                tc.seedDiagnostic());

            % And the ESTIMATE still sees the seeds: only the trace narrowed. The
            % strict inequality is the point -- the quantity the trace is taken on
            % never reaches the value that is reported.
            tc.verifyLessThan(info.trace.value(end), est, tc.seedDiagnostic());
        end

        function rejectsAnUnknownProbeMethod(tc)
            tc.verifyError(@() rid.roi.estimateFillDistance([0.5 0.5], ...
                tc.unitSquare(), "probeMethod", "halton"), ...
                "rid:roi:unknownProbeMethod");
        end

        % --- the trace acts ------------------------------------------

        function tinyProbeCountIsReportedUnconverged(tc)
            % A 3-D grid probed with 40 points: the estimate lands 12 % below the
            % analytic supremum sqrt(3)/8 and the trace is still jumping. The
            % estimator must say so itself -- the convergence-trace safeguard's whole point is that no
            % reader is trusted to compare the three numbers. This fixture is also
            % the one where the last-increment rule reports CONVERGED, because the
            % running maximum happens to be flat over the final doubling.
            nrm = rid.roi.Normalization([0 1; 0 1; 0 1], ["a" "b" "c"]);
            Q   = tc.cellCentredGrid(4, [1 1 1]);
            [est, info] = rid.roi.estimateFillDistance(Q, nrm, ...
                "nProbes", 40, "seed", tc.Seed, "probeMethod", "sobol");
            tc.verifyFalse(info.converged, tc.seedDiagnostic());
            tc.verifyLessThan(est, 0.9 * sqrt(3) / 8, tc.seedDiagnostic());
        end

        function ampleProbeCountIsReportedConverged(tc)
            % The same point set, probed densely enough that doubling the probe
            % count from a half to all of the set moves the maximum by less than
            % the declared tolerance.
            nrm = rid.roi.Normalization([0 1; 0 1; 0 1], ["a" "b" "c"]);
            Q   = tc.cellCentredGrid(4, [1 1 1]);
            [~, info] = rid.roi.estimateFillDistance(Q, nrm, ...
                "nProbes", 2e5, "seed", tc.Seed, "probeMethod", "sobol");
            tc.verifyTrue(info.converged, tc.seedDiagnostic());
        end

        function verdictFollowsTheDeclaredTolerance(tc)
            % Same estimate, same probes: only the declared tolerance moves, so
            % the verdict is a comparison against it and not a property of the
            % probe count alone.
            nrm = tc.unitSquare();
            Q   = tc.cellCentredGrid(4, [1 1]);
            [~, strict] = rid.roi.estimateFillDistance(Q, nrm, ...
                "nProbes", 400, "seed", tc.Seed, "probeMethod", "sobol", ...
                "traceTolerance", 0);
            [~, loose]  = rid.roi.estimateFillDistance(Q, nrm, ...
                "nProbes", 400, "seed", tc.Seed, "probeMethod", "sobol", ...
                "traceTolerance", 1);
            tc.verifyFalse(strict.converged, tc.seedDiagnostic());
            tc.verifyTrue(loose.converged, tc.seedDiagnostic());
            tc.verifyEqual(strict.trace.value, loose.trace.value, ...
                tc.seedDiagnostic());
        end

        function relativeClimbIsTheTraceItself(tc)
            % The reported climb is the quantity the verdict is taken on, not a
            % separate calculation that could drift from the trace.
            nrm = tc.unitSquare();
            [~, info] = rid.roi.estimateFillDistance(tc.cellCentredGrid(4, [1 1]), ...
                nrm, "nProbes", 4000, "seed", tc.Seed, "probeMethod", "sobol");
            v = info.trace.value;
            tc.verifyEqual(info.trace.relativeClimb, diff(v) ./ v(2:end), ...
                "RelTol", 1e-12, tc.seedDiagnostic());
            tc.verifyEqual(info.converged, ...
                max(info.trace.relativeClimb) <= info.traceTolerance, ...
                tc.seedDiagnostic());
        end

        function aTraceThatCollapsedToOnePointIsNotConverged(tc)
            % nProbes = 1 makes the quarter, half and full counts the same number,
            % so nothing was compared. An estimate from a single probe must not be
            % able to certify its own convergence.
            [~, info] = rid.roi.estimateFillDistance([0.5 0.5], tc.unitSquare(), ...
                "nProbes", 1, "seed", tc.Seed, "probeMethod", "sobol", ...
                "traceTolerance", 1);
            tc.verifyEqual(numel(info.trace.value), 1);
            tc.verifyFalse(info.converged, tc.seedDiagnostic());
        end

        function rejectsNegativeTraceTolerance(tc)
            tc.verifyError(@() rid.roi.estimateFillDistance([0.5 0.5], ...
                tc.unitSquare(), "traceTolerance", -0.1), ...
                "rid:roi:invalidTraceTolerance");
        end

        % --- prefixes and the centre ---------------------------------

        function everyPrefixEqualsASeparateCall(tc)
            % The one-pass prefix path is a reorganization of the same minima and
            % maxima, so every per-prefix quantity must be BIT-identical to a
            % separate call on the leading block -- value, bracket, centre, trace
            % and verdict. Checked for each probe method, with a block size that
            % does not divide the probe count, on a weighted box and on points
            % partly outside it. Prefixes include 1 and steps of one sample, where a
            % carried distance that is never lowered would show first.
            nrm = tc.weightedSquare();
            rs  = RandStream("twister", "Seed", tc.Seed);
            P   = [2, 1] .* (1.2 * rand(rs, 300, 2) - 0.1);
            sizes = [1, 2, 3, 17, 18, 120, 299, 300];
            for method = ["sobolCorners", "sobol", "uniform"]
                args = {"nProbes", 2500, "seed", tc.Seed, "probeMethod", method, ...
                    "blockSize", 333};
                [eps, info] = rid.roi.estimateFillDistance(P, nrm, args{:}, ...
                    "prefixSizes", sizes);
                tc.verifySize(eps, [1, numel(sizes)]);
                tc.verifySize(info.centre, [numel(sizes), 2]);
                tc.verifyEqual(info.prefixSizes, sizes);
                for iPrefix = 1:numel(sizes)
                    [e1, i1] = rid.roi.estimateFillDistance(P(1:sizes(iPrefix), :), ...
                        nrm, args{:});
                    msg = sprintf("%s, prefix %d, %s", method, sizes(iPrefix), ...
                        tc.seedDiagnostic());
                    tc.verifyEqual(eps(iPrefix), e1, msg);
                    tc.verifyEqual(info.bracket(iPrefix, :), i1.bracket, msg);
                    tc.verifyEqual(info.centre(iPrefix, :), i1.centre, msg);
                    tc.verifyEqual(info.trace.value(iPrefix, :), i1.trace.value, msg);
                    tc.verifyEqual(info.trace.relativeClimb(iPrefix, :), ...
                        i1.trace.relativeClimb, msg);
                    tc.verifyEqual(info.converged(iPrefix), i1.converged, msg);
                    tc.verifyEqual(info.probeFillDistance, i1.probeFillDistance, msg);
                end
            end
        end

        function theFullPrefixIsTheDefaultCall(tc)
            % prefixSizes = nPoints is the default call, field for field.
            nrm = tc.unitSquare();
            P   = tc.cellCentredGrid(5, [1 1]);
            [a, infoA] = rid.roi.estimateFillDistance(P, nrm, "nProbes", 3000, "seed", tc.Seed);
            [b, infoB] = rid.roi.estimateFillDistance(P, nrm, "nProbes", 3000, ...
                "seed", tc.Seed, "prefixSizes", size(P, 1));
            tc.verifyEqual(a, b, tc.seedDiagnostic());
            tc.verifyEqual(infoA, infoB, tc.seedDiagnostic());
            tc.verifyEqual(infoA.prefixSizes, size(P, 1));
        end

        function theCentreIsTheProbeThatAttainsTheEstimate(tc)
            % Its distance to the nearest point IS the estimate, and it lies in the
            % normalized box. Anchored on the distance itself, so a centre taken
            % from the wrong probe, the wrong set or the wrong prefix fails here.
            nrm = tc.weightedSquare();
            rs  = RandStream("twister", "Seed", tc.Seed);
            P   = [2, 1] .* rand(rs, 40, 2);
            sizes = [5, 20, 40];
            [eps, info] = rid.roi.estimateFillDistance(P, nrm, "nProbes", 4000, ...
                "seed", tc.Seed, "prefixSizes", sizes);
            for iPrefix = 1:numel(sizes)
                c = info.centre(iPrefix, :);
                d = min(rid.ext.pairwiseDistance(c, P(1:sizes(iPrefix), :)));
                tc.verifyEqual(d, eps(iPrefix), tc.seedDiagnostic());
                tc.verifyGreaterThanOrEqual(c, [0 0]);
                tc.verifyLessThanOrEqual(c, [2 1]);
            end
        end

        function theCentreOfACornerPointIsTheOppositeVertex(tc)
            % Analytic: P = {(0,0)} on the unit square; the deepest point is the
            % vertex (1,1), which the default probe set seeds exactly.
            [eps, info] = rid.roi.estimateFillDistance([0 0], tc.unitSquare(), ...
                "nProbes", 1000, "seed", tc.Seed);
            tc.verifyEqual(info.centre, [1 1], tc.seedDiagnostic());
            tc.verifyEqual(eps, sqrt(2), tc.seedDiagnostic());
        end

        function aTieKeepsTheEarlierProbeWhateverTheBlocks(tc)
            % P = 0.5 on [0,1]: both ends are seeded probes at exactly 0.5. The
            % first in probe order, the vertex 0, is the centre -- also when the
            % two ends fall into different blocks.
            for blockSize = [1, 2048]
                [~, info] = rid.roi.estimateFillDistance(0.5, tc.unitInterval(), ...
                    "nProbes", 100, "seed", tc.Seed, "blockSize", blockSize);
                tc.verifyEqual(info.centre, 0, sprintf("blockSize %d", blockSize));
            end
        end

        function rejectsInvalidPrefixSizes(tc)
            nrm = tc.unitSquare();
            P   = tc.cellCentredGrid(3, [1 1]);
            for bad = {[0 3], [3 3], [5 4], [2 10], [1.5 3], [1 NaN]}
                tc.verifyError(@() rid.roi.estimateFillDistance(P, nrm, ...
                    "prefixSizes", bad{1}), "rid:roi:invalidPrefixSizes", ...
                    mat2str(bad{1}));
            end
        end

        % --- randomness contract ---------------------------------------------

        function deterministicUnderAFixedSeed(tc)
            nrm = tc.unitSquare();
            P   = tc.cellCentredGrid(4, [1 1]);
            a = rid.roi.estimateFillDistance(P, nrm, "nProbes", 3000, "seed", tc.Seed);
            b = rid.roi.estimateFillDistance(P, nrm, "nProbes", 3000, "seed", tc.Seed);
            tc.verifyEqual(a, b, tc.seedDiagnostic());
        end

        function differsAcrossSeeds(tc)
            % Pinned to a method with NO deterministic seeds, because that is what
            % this test is about. Under the seeded default the estimate can be
            % seed-INDEPENDENT and often is: the corners carry the maximum and the
            % corners are the same at every seed. That is the feature working --
            % cornerSeedingIsExactOnACellCentredGridAndThatIsTautological asserts
            % the same fact from the other side.
            nrm = tc.unitSquare();
            P   = tc.cellCentredGrid(4, [1 1]);
            a = rid.roi.estimateFillDistance(P, nrm, "nProbes", 3000, ...
                "seed", tc.Seed, "probeMethod", "sobol");
            b = rid.roi.estimateFillDistance(P, nrm, "nProbes", 3000, ...
                "seed", tc.OtherSeed, "probeMethod", "sobol");
            tc.verifyNotEqual(a, b, ...
                sprintf("seeds = %d, %d", tc.Seed, tc.OtherSeed));
        end

        function blockSizeChangesNothingButMemory(tc)
            % Probes are drawn in probe order, so blocking is not observable in the
            % result. Bit-identical, not merely close.
            nrm = tc.unitSquare();
            P   = tc.cellCentredGrid(4, [1 1]);
            [a, infoA] = rid.roi.estimateFillDistance(P, nrm, ...
                "nProbes", 3000, "seed", tc.Seed, "blockSize", 1);
            [b, infoB] = rid.roi.estimateFillDistance(P, nrm, ...
                "nProbes", 3000, "seed", tc.Seed, "blockSize", 997);
            [c, infoC] = rid.roi.estimateFillDistance(P, nrm, ...
                "nProbes", 3000, "seed", tc.Seed, "blockSize", 100000);
            tc.verifyEqual(a, b, tc.seedDiagnostic());
            tc.verifyEqual(a, c, tc.seedDiagnostic());
            tc.verifyEqual(infoA.trace.value, infoB.trace.value, tc.seedDiagnostic());
            tc.verifyEqual(infoA.trace.value, infoC.trace.value, tc.seedDiagnostic());
        end

        function doesNotTouchTheGlobalStream(tc)
            % Reproducible in isolation and inside parfor: the estimator builds its
            % own stream and neither reads nor advances the global one.
            RandStream.setGlobalStream(RandStream("philox", "Seed", 4242));
            before = RandStream.getGlobalStream();
            stateBefore = before.State;

            rid.roi.estimateFillDistance([0.5 0.5], tc.unitSquare(), ...
                "nProbes", 5000, "seed", tc.Seed);

            after = RandStream.getGlobalStream();
            tc.verifySameHandle(after, before, tc.seedDiagnostic());
            tc.verifyEqual(after.State, stateBefore, tc.seedDiagnostic());
        end

        % --- provenance ------------------------------------------------------

        function infoReportsMethodSeedAndProbeCount(tc)
            [~, info] = rid.roi.estimateFillDistance([0.5 0.5], tc.unitSquare(), ...
                "nProbes", 1234, "seed", tc.Seed, "blockSize", 256);
            tc.verifyEqual(info.method, "monteCarlo");
            tc.verifyEqual(info.nProbes, 1234);
            tc.verifyEqual(info.seed, tc.Seed);
            tc.verifyEqual(info.blockSize, 256);
            tc.verifyEqual(info.traceTolerance, ...
                rid.defaultOptions().roi.traceTolerance);
            tc.verifyEqual(info.probeMethod, ...
                rid.defaultOptions().roi.probeMethod);
        end

        % --- rejections ------------------------------------------------------

        function rejectsEmptyPointSet(tc)
            tc.verifyError(@() rid.roi.estimateFillDistance(zeros(0, 2), ...
                tc.unitSquare()), "rid:roi:emptyPointSet");
        end

        function rejectsDimensionMismatch(tc)
            tc.verifyError(@() rid.roi.estimateFillDistance([0.5 0.5 0.5], ...
                tc.unitSquare()), "rid:roi:dimensionMismatch");
        end

        function rejectsNonPositiveOrFractionalProbeCount(tc)
            nrm = tc.unitSquare();
            tc.verifyError(@() rid.roi.estimateFillDistance([0.5 0.5], nrm, ...
                "nProbes", 0), "rid:roi:invalidProbeCount");
            tc.verifyError(@() rid.roi.estimateFillDistance([0.5 0.5], nrm, ...
                "nProbes", 10.5), "rid:roi:invalidProbeCount");
            tc.verifyError(@() rid.roi.estimateFillDistance([0.5 0.5], nrm, ...
                "blockSize", 0), "rid:roi:invalidProbeCount");
        end

    end

end
