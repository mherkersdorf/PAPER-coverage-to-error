classdef EvaluateTest < matlab.unittest.TestCase
%EVALUATETEST  Analytic cases and properties P1-P7 for rid.crit.evaluate.
%
%   Layers 1 (analytic cases) and 2 (property tests) of the project's test plan, against the
%   permanent naive implementation. The properties are the test contract of
%   the method description; each test names the property it covers.
%
%   P1 permutation invariance      P5 saturation at psiMax, J in [0, psiMax]
%   P2 monotone under added samples P6 uniform demand, and the distance-based
%                                      limit where nothing is yet resolved
%   P3 epsReq0 scaling equivalence  P7 continuity, with kinks exactly at the
%   P4 J = 0 iff all dTilde <= 1       two documented places
%
%   Every test is seeded from a literal; randomness comes from a private
%   RandStream, never the global one.

    properties (Constant)
        Seed = 20260830
        Tol  = 1e-12
    end

    methods (Access = private)

        function s = stream(tc, offset)
            s = RandStream("twister", "Seed", tc.Seed + offset);
        end

        function d = seedDiagnostic(tc)
            d = sprintf("seed = %d", tc.Seed);
        end

        function Q = grid1d(~, nSup)
            Q = ((0:nSup - 1).' / (nSup - 1));
        end

    end

    methods (Test)

        % --- layer 1: analytic cases -----------------------------------------

        function singleSampleSingleSupportIsHandComputable(tc)
            % dRaw = 0.3, epsReq = 0.1 -> dTilde = 3, psi = min(2, 1) = 1, J = 1.
            J = rid.crit.evaluate(0.2, 0.5, 0.1, "psiMax", 1, "pNorm", 1);
            tc.verifyEqual(J, 1, "AbsTol", tc.Tol);

            % Unclipped, the same case gives psi = 2 exactly.
            J = rid.crit.evaluate(0.2, 0.5, 0.1, "psiMax", Inf, "pNorm", 1);
            tc.verifyEqual(J, 2, "RelTol", tc.Tol);
        end

        function uniformGridWithOneCentredSampleIsHandComputable(tc)
            % Q = {0, 0.25, 0.5, 0.75, 1}, one sample at 0.5, epsReq = 0.25.
            % dRaw   = [0.5 0.25 0 0.25 0.5]
            % dTilde = [2 1 0 1 2],  psi = [1 0 0 0 1] at psiMax = 1
            % J      = mean(psi) = 0.4
            Q      = (0:0.25:1).';
            epsReq = 0.25 * ones(5, 1);
            [J, info] = rid.crit.evaluate(0.5, Q, epsReq, "psiMax", 1, "pNorm", 1);
            tc.verifyEqual(info.dRaw,   [0.5; 0.25; 0; 0.25; 0.5], "AbsTol", tc.Tol);
            tc.verifyEqual(info.dTilde, [2; 1; 0; 1; 2],           "AbsTol", tc.Tol);
            tc.verifyEqual(info.psi,    [1; 0; 0; 0; 1],           "AbsTol", tc.Tol);
            tc.verifyEqual(J, 0.4, "AbsTol", tc.Tol);
        end

        function clippingIsActiveExactlyAtPsiMax(tc)
            % One remote sample: dTilde - 1 = 9, clipped to psiMax = 0.5.
            [J, info] = rid.crit.evaluate(10, 0, 1, "psiMax", 0.5, "pNorm", 1);
            tc.verifyEqual(info.dTilde, 10, "AbsTol", tc.Tol);
            tc.verifyEqual(info.psi, 0.5, "AbsTol", tc.Tol);
            tc.verifyEqual(J, 0.5, "AbsTol", tc.Tol);
        end

        function emptyDistributionScoresPsiMax(tc)
            % Nothing is covered before the first sample: dRaw = Inf, psi = psiMax
            % at every supporting point.
            Q = tc.grid1d(5);
            [J, info] = rid.crit.evaluate(zeros(0, 1), Q, 0.2 * ones(5, 1), ...
                "psiMax", 1, "pNorm", 1);
            tc.verifyEqual(info.dRaw, inf(5, 1));
            tc.verifyEqual(J, 1, "AbsTol", tc.Tol);
            tc.verifyEqual(info.fracResolved, 0);

            % Unclipped, the same state is legitimately infinite.
            tc.verifyEqual(rid.crit.evaluate(zeros(0, 1), Q, 0.2 * ones(5, 1), ...
                "psiMax", Inf), Inf);
        end

        function pNormInfinityIsTheWorstCase(tc)
            % The formula would give Inf^(1/Inf) = NaN, so p = Inf has its own
            % branch: J = max(psi).
            Q      = (0:0.25:1).';
            epsReq = 0.25 * ones(5, 1);
            [~, info] = rid.crit.evaluate(0.5, Q, epsReq, "psiMax", 1, "pNorm", 1);
            J = rid.crit.evaluate(0.5, Q, epsReq, "psiMax", 1, "pNorm", Inf);
            tc.verifyEqual(J, max(info.psi), "AbsTol", tc.Tol);
            tc.verifyFalse(isnan(J));
        end

        function largePNormApproachesTheWorstCase(tc)
            s      = tc.stream(1);
            Q      = rand(s, 40, 2);
            Znorm   = rand(s, 3, 2);
            epsReq = 0.05 * ones(40, 1);
            psiInf = rid.crit.evaluate(Znorm, Q, epsReq, "psiMax", 1, "pNorm", Inf);
            jBig   = rid.crit.evaluate(Znorm, Q, epsReq, "psiMax", 1, "pNorm", 256);
            tc.verifyLessThanOrEqual(jBig, psiInf + tc.Tol, tc.seedDiagnostic());
            tc.verifyGreaterThan(jBig, 0.9 * psiInf, tc.seedDiagnostic());
        end

        % --- layer 2: the properties of section 6 -----------------------------

        function p1_invariantUnderRowPermutationOfZnorm(tc)
            s      = tc.stream(2);
            Q      = rand(s, 30, 2);
            Znorm   = rand(s, 12, 2);
            epsReq = 0.05 + 0.1 * rand(s, 30, 1);
            perm   = randperm(s, 12);
            tc.verifyEqual( ...
                rid.crit.evaluate(Znorm(perm, :), Q, epsReq), ...
                rid.crit.evaluate(Znorm, Q, epsReq), ...
                "AbsTol", tc.Tol, tc.seedDiagnostic());
        end

        function p2_monotoneNonIncreasingWhenSamplesAreAdded(tc)
            % Holds with epsReq FIXED -- that is the whole content of the property.
            % Between resolution updates it is an invariant of the loop; across one
            % it is not, and nothing here claims otherwise.
            s      = tc.stream(3);
            Q      = rand(s, 40, 2);
            epsReq = 0.08 * ones(40, 1);
            Znorm   = zeros(0, 2);
            jPrev  = rid.crit.evaluate(Znorm, Q, epsReq);
            for kStep = 1:25
                Znorm = [Znorm; rand(s, 1, 2)]; %#ok<AGROW> the growth is the test
                jNow = rid.crit.evaluate(Znorm, Q, epsReq);
                tc.verifyLessThanOrEqual(jNow, jPrev + tc.Tol, ...
                    sprintf("%s, step %d", tc.seedDiagnostic(), kStep));
                jPrev = jNow;
            end
        end

        function p3_scalingEpsReqIsScalingAllDistances(tc)
            % Scaling epsReq0 by alpha is equivalent to scaling every distance by
            % 1/alpha: the criterion sees only the ratio.
            s      = tc.stream(4);
            Q      = rand(s, 30, 2);
            Znorm   = rand(s, 8, 2);
            epsReq = 0.05 + 0.1 * rand(s, 30, 1);
            for alpha = [0.25, 0.5, 2, 7]
                tc.verifyEqual( ...
                    rid.crit.evaluate(Znorm, Q, alpha * epsReq), ...
                    rid.crit.evaluate(Znorm / alpha, Q / alpha, epsReq), ...
                    "RelTol", tc.Tol, ...
                    sprintf("%s, alpha = %g", tc.seedDiagnostic(), alpha));
            end
        end

        function p4_isZeroExactlyWhenEveryPointIsResolved(tc)
            % Both directions, and the zero is exact -- J = 0 is the certificate,
            % so a value of 1e-17 would be a different statement.
            Q      = (0:0.25:1).';
            epsReq = 0.3 * ones(5, 1);
            [J, info] = rid.crit.evaluate(Q, Q, epsReq);       % a sample on every point
            tc.verifyTrue(all(info.dTilde <= 1));
            tc.verifyEqual(J, 0);

            % Exactly at the boundary dTilde == 1: still resolved.
            J = rid.crit.evaluate(0, 0.3, 0.3);
            tc.verifyEqual(J, 0);

            % One point just beyond it: strictly positive.
            J = rid.crit.evaluate(0, 0.3 * (1 + 1e-9), 0.3);
            tc.verifyGreaterThan(J, 0);

            % And the converse over a random configuration.
            s      = tc.stream(5);
            Q      = rand(s, 50, 2);
            Znorm   = rand(s, 5, 2);
            [J, info] = rid.crit.evaluate(Znorm, Q, 0.05 * ones(50, 1));
            tc.verifyEqual(J == 0, all(info.dTilde <= 1), tc.seedDiagnostic());
            tc.verifyGreaterThan(J, 0, tc.seedDiagnostic());
        end

        function p5_saturatesAtPsiMaxAndStaysInRange(tc)
            % An arbitrarily distant supporting point contributes exactly psiMax and
            % no more: that is what stops a dynamically unreachable region from
            % dominating the optimization.
            Q = [0; 1e6];
            for psiMax = [0.1, 1, 5]
                [J, info] = rid.crit.evaluate(0, Q, [1; 1], "psiMax", psiMax);
                tc.verifyEqual(info.psi, [0; psiMax], "AbsTol", tc.Tol);
                tc.verifyEqual(J, psiMax / 2, "AbsTol", tc.Tol);
            end

            % J stays inside [0, psiMax] over random configurations.
            s = tc.stream(6);
            for kTrial = 1:20
                Q      = rand(s, 25, 2);
                Znorm   = rand(s, 4, 2);
                psiMax = 0.1 + 2 * rand(s);
                J = rid.crit.evaluate(Znorm, Q, 0.02 + 0.1 * rand(s, 25, 1), ...
                    "psiMax", psiMax, "pNorm", 1 + 3 * rand(s));
                tc.verifyGreaterThanOrEqual(J, 0, ...
                    sprintf("%s, trial %d", tc.seedDiagnostic(), kTrial));
                tc.verifyLessThanOrEqual(J, psiMax + tc.Tol, ...
                    sprintf("%s, trial %d", tc.seedDiagnostic(), kTrial));
            end
        end

        function p6_nominalFieldsWithUnclippedPenaltyIsTheSpaceFillingLimit(tc)
            % With nominal fields epsReq is uniform, and with psiMax -> Inf the
            % criterion collapses to J = mean(dRaw)/epsReq0 - 1: the distance-based
            % nearest-neighbour criterion of the earlier work (section 10). That
            % collapse holds ONLY while every supporting point is unresolved --
            % psiMax removes the upper clip, never the positive part -- so the
            % regime is asserted below rather than assumed.
            s       = tc.stream(7);
            Q       = rand(s, 40, 2);
            Znorm    = rand(s, 6, 2);
            epsReq0 = 1e-3;
            epsReq  = rid.crit.demandedResolution(rid.fields.nominalFields(40), ...
                "epsReq0", epsReq0);
            tc.verifyEqual(epsReq, epsReq0 * ones(40, 1), "AbsTol", tc.Tol);

            [J, info] = rid.crit.evaluate(Znorm, Q, epsReq, "psiMax", Inf, "pNorm", 1);

            % The closed form holds where nothing is resolved, so that the positive
            % part is inactive; a baseline this fine guarantees it here.
            tc.assertTrue(all(info.dTilde >= 1), tc.seedDiagnostic());
            tc.verifyEqual(J, mean(info.dRaw) / epsReq0 - 1, ...
                "RelTol", tc.Tol, tc.seedDiagnostic());
        end

        function p7_isContinuousInTheSamplePositionsWithKinksOnlyAtTheTwoJoins(tc)
            % One supporting point at the origin, epsReq = 1, psiMax = 0.5. Moving a
            % single sample outward along a ray makes J piecewise linear in the
            % distance t: 0 for t <= 1, t - 1 for 1 <= t <= 1.5, 0.5 beyond.
            Q      = 0;
            epsReq = 1;
            psiMax = 0.5;
            jAt    = @(t) rid.crit.evaluate(t, Q, epsReq, "psiMax", psiMax, "pNorm", 1);

            tGrid = linspace(0, 3, 601);
            jGrid = arrayfun(jAt, tGrid);

            % Continuity: the increment over a step of the grid is bounded by the
            % largest slope, which is 1. No jump anywhere, including at the joins.
            tc.verifyLessThanOrEqual(max(abs(diff(jGrid))), ...
                (tGrid(2) - tGrid(1)) + tc.Tol);

            % Smooth away from the joins: the second difference vanishes.
            isInterior = tGrid > 1.6 & tGrid < 2.9;
            tc.verifyLessThanOrEqual(max(abs(diff(jGrid(isInterior), 2))), tc.Tol);

            % Non-smooth exactly at dTilde = 1 and at the clip: the one-sided slopes
            % differ there and nowhere else.
            h = 1e-6;
            for tKink = [1, 1.5]
                slopeLeft  = (jAt(tKink) - jAt(tKink - h)) / h;
                slopeRight = (jAt(tKink + h) - jAt(tKink)) / h;
                tc.verifyGreaterThan(abs(slopeRight - slopeLeft), 0.5, ...
                    sprintf("expected a kink at t = %g", tKink));
            end
            for tSmooth = [0.5, 1.25, 2]
                slopeLeft  = (jAt(tSmooth) - jAt(tSmooth - h)) / h;
                slopeRight = (jAt(tSmooth + h) - jAt(tSmooth)) / h;
                tc.verifyLessThan(abs(slopeRight - slopeLeft), 1e-4, ...
                    sprintf("expected smoothness at t = %g", tSmooth));
            end
        end

        % --- validation -------------------------------------------------------

        function rejectsNonPositiveOrNaNResolution(tc)
            Q = (0:0.25:1).';
            tc.verifyError(@() rid.crit.evaluate(0.5, Q, zeros(5, 1)), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.evaluate(0.5, Q, -0.1 * ones(5, 1)), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.evaluate(0.5, Q, [0.2; NaN; 0.2; 0.2; 0.2]), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.evaluate(0.5, Q, [0.2; Inf; 0.2; 0.2; 0.2]), ...
                "rid:crit:invalidResolution");
        end

        function rejectsShapeMismatch(tc)
            Q = (0:0.25:1).';
            tc.verifyError(@() rid.crit.evaluate(0.5, Q, 0.2 * ones(4, 1)), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.evaluate([0.5 0.5], Q, 0.2 * ones(5, 1)), ...
                "rid:crit:dimensionMismatch");
            tc.verifyError(@() rid.crit.evaluate(0.5, zeros(0, 1), zeros(0, 1)), ...
                "rid:crit:invalidResolution");
        end

        function rejectsInvalidClipOrNorm(tc)
            tc.verifyError(@() rid.crit.evaluate(0.5, 0, 1, "psiMax", 0), ...
                "rid:crit:invalidPenaltyClip");
            tc.verifyError(@() rid.crit.evaluate(0.5, 0, 1, "psiMax", -1), ...
                "rid:crit:invalidPenaltyClip");
            tc.verifyError(@() rid.crit.evaluate(0.5, 0, 1, "pNorm", 0.5), ...
                "rid:crit:invalidNorm");
        end

    end

end
