classdef DemandedResolutionTest < matlab.unittest.TestCase
%DEMANDEDRESOLUTIONTEST  Analytic cases and properties P8-P12 for the demanded resolution.
%
%   Layers 1 and 2 of the project's test plan against rid.crit.demandedResolution. The properties
%   are the test contract of the method description:
%
%   P8  max(epsReq) <= epsReq0, with equality iff some point has refinement == 1;
%       epsReq < epsReq0 exactly where refinement < 1
%   P9  LOCALITY -- epsReq(j) changes only when a field changes at q_j. This is the
%       property that fails if a normalization by a field statistic (a mean, median
%       or maximum over the field) is ever reintroduced, which is why it is tested
%       against every other supporting point bit-for-bit rather than to a tolerance
%   P10 noiseLevel nominal reproduces the bias-only design exactly, for ANY betaNoise
%   P11 all fields nominal give epsReq == epsReq0 at every supporting point
%   P12 certification mode floors at fillDistanceQ and reports the count; budget
%       mode applies no floor

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

        function f = fieldsFrom(~, relevance, complexity, noiseLevel)
            f = struct("relevance", relevance, "complexity", complexity, ...
                       "noiseLevel", noiseLevel);
        end

    end

    methods (Test)

        % --- layer 1: analytic cases -----------------------------------------

        function refinementIsHandComputable(tc)
            % beta = 0.5, betaNoise = 0.5, epsReq0 = 0.2.
            % (a c)^-0.5 * v^-0.5 with a = 1, c = 4, v = 1 gives 0.5 -> epsReq = 0.1.
            % With c = 1/4 the factor would be 2, clipped to 1.
            f = tc.fieldsFrom([1; 1; 1], [4; 1; 0.25], [1; 1; 1]);
            [epsReq, info] = rid.crit.demandedResolution(f, ...
                "epsReq0", 0.2, "beta", 0.5, "betaNoise", 0.5);
            tc.verifyEqual(info.refinement, [0.5; 1; 1], "AbsTol", tc.Tol);
            tc.verifyEqual(epsReq, [0.1; 0.2; 0.2], "AbsTol", tc.Tol);
            tc.verifyEqual(info.fracRefined, 1/3, "AbsTol", tc.Tol);
            tc.verifyEqual(info.minEpsReq, 0.1, "AbsTol", tc.Tol);
            tc.verifyEqual(info.minEpsReqRatio, 0.5, "AbsTol", tc.Tol);
        end

        function relevanceAndComplexityEnterAsOneProduct(tc)
            % The specification pairs them inside a single power (section 3), so
            % a = 2, c = 8 must give exactly what a = 16, c = 1 gives.
            fSplit  = tc.fieldsFrom([2; 1], [8; 1], [1; 1]);
            fJoined = tc.fieldsFrom([16; 1], [1; 1], [1; 1]);
            tc.verifyEqual( ...
                rid.crit.demandedResolution(fSplit,  "epsReq0", 0.3), ...
                rid.crit.demandedResolution(fJoined, "epsReq0", 0.3), ...
                "AbsTol", tc.Tol);
        end

        function betaZeroDisablesTargetShapingEntirely(tc)
            s = tc.stream(1);
            f = tc.fieldsFrom(1 + 9 * rand(s, 20, 1), 1 + 9 * rand(s, 20, 1), ...
                              1 + 9 * rand(s, 20, 1));
            epsReq = rid.crit.demandedResolution(f, ...
                "epsReq0", 0.15, "beta", 0, "betaNoise", 0);
            tc.verifyEqual(epsReq, 0.15 * ones(20, 1), "AbsTol", tc.Tol, ...
                tc.seedDiagnostic());
        end

        function absentMembersDefaultToNominal(tc)
            % A field that is not supplied is nominal; nSup is
            % inferred from the members that are present.
            explicit = rid.crit.demandedResolution( ...
                tc.fieldsFrom([1; 1; 1], [4; 1; 1], [1; 1; 1]), "epsReq0", 0.2);
            implicit = rid.crit.demandedResolution( ...
                struct("complexity", [4; 1; 1]), "epsReq0", 0.2);
            tc.verifyEqual(implicit, explicit, "AbsTol", tc.Tol);
        end

        % --- layer 2: the properties of section 6 -----------------------------

        function p8_maxIsBoundedByTheBaselineWithEqualityIffAPointIsUnrefined(tc)
            s = tc.stream(2);

            % A configuration with at least one unrefined point: equality holds.
            complexity = [1; 2; 5; 1; 3];
            [epsReq, info] = rid.crit.demandedResolution( ...
                struct("complexity", complexity), "epsReq0", 0.2);
            tc.verifyTrue(any(info.refinement == 1));
            tc.verifyEqual(max(epsReq), 0.2, "AbsTol", tc.Tol);

            % epsReq < epsReq0 exactly where refinement < 1, index by index.
            tc.verifyEqual(epsReq < 0.2, info.refinement < 1);

            % Every point refined: the maximum lies STRICTLY below the baseline and
            % the certification bound is merely conservative, not violated.
            [epsReq, info] = rid.crit.demandedResolution( ...
                struct("complexity", [2; 3; 7]), "epsReq0", 0.2);
            tc.verifyTrue(all(info.refinement < 1));
            tc.verifyLessThan(max(epsReq), 0.2);
            tc.verifyEqual(info.fracRefined, 1);

            % And the inequality itself over random fields.
            for kTrial = 1:20
                nSup = randi(s, [1 30]);
                f = tc.fieldsFrom(0.1 + 5 * rand(s, nSup, 1), ...
                                  0.1 + 5 * rand(s, nSup, 1), ...
                                  0.1 + 5 * rand(s, nSup, 1));
                [epsReq, info] = rid.crit.demandedResolution(f, "epsReq0", 0.2);
                tc.verifyLessThanOrEqual(max(epsReq), 0.2 + tc.Tol, ...
                    sprintf("%s, trial %d", tc.seedDiagnostic(), kTrial));
                tc.verifyEqual(epsReq < 0.2, info.refinement < 1, ...
                    sprintf("%s, trial %d", tc.seedDiagnostic(), kTrial));
            end
        end

        function p9_epsReqAtOnePointIsUntouchedByFieldsElsewhere(tc)
            % Locality. Perturb one field at one supporting point and assert that
            % every OTHER entry is bit-identical -- not merely close. A normalization
            % by any statistic of the field (mean, median, max) makes this fail,
            % which is exactly why it is not permitted.
            s     = tc.stream(3);
            nSup  = 25;
            f     = tc.fieldsFrom(1 + 3 * rand(s, nSup, 1), ...
                                  1 + 3 * rand(s, nSup, 1), ...
                                  1 + 3 * rand(s, nSup, 1));
            base  = rid.crit.demandedResolution(f, "epsReq0", 0.2);

            names = ["relevance", "complexity", "noiseLevel"];
            for name = names
                for jSup = [1, 7, nSup]
                    perturbed = f;
                    perturbed.(name)(jSup) = perturbed.(name)(jSup) * 37;
                    epsReq = rid.crit.demandedResolution(perturbed, "epsReq0", 0.2);

                    others = setdiff(1:nSup, jSup);
                    tc.verifyEqual(epsReq(others), base(others), ...
                        sprintf("%s: perturbing %s at %d moved another point", ...
                            tc.seedDiagnostic(), name, jSup));
                    tc.verifyLessThan(epsReq(jSup), base(jSup), ...
                        sprintf("%s: perturbing %s at %d did not refine it", ...
                            tc.seedDiagnostic(), name, jSup));
                end
            end

            % The same statement in its enlarge-the-ROI form: appending easy
            % supporting points must not change the demand at the existing ones.
            enlarged = tc.fieldsFrom([f.relevance;  ones(40, 1)], ...
                                     [f.complexity; ones(40, 1)], ...
                                     [f.noiseLevel; ones(40, 1)]);
            epsReq = rid.crit.demandedResolution(enlarged, "epsReq0", 0.2);
            tc.verifyEqual(epsReq(1:nSup), base, tc.seedDiagnostic());
        end

        function p10_nominalNoiseReproducesTheBiasOnlyDesignForAnyBetaNoise(tc)
            s      = tc.stream(4);
            nSup   = 20;
            f      = tc.fieldsFrom(1 + 3 * rand(s, nSup, 1), ...
                                   1 + 3 * rand(s, nSup, 1), ...
                                   ones(nSup, 1));
            biasOnly = rid.crit.demandedResolution( ...
                struct("relevance", f.relevance, "complexity", f.complexity), ...
                "epsReq0", 0.2, "beta", 0.5);

            for betaNoise = [0, 0.25, 0.5, 1, 4]
                epsReq = rid.crit.demandedResolution(f, ...
                    "epsReq0", 0.2, "beta", 0.5, "betaNoise", betaNoise);
                tc.verifyEqual(epsReq, biasOnly, ...
                    sprintf("%s, betaNoise = %g", tc.seedDiagnostic(), betaNoise));
            end
        end

        function p11_allFieldsNominalGiveTheBaselineEverywhere(tc)
            % The pure space-filling design, and the robust fallback whenever no
            % target information is available. Exact, not to a tolerance.
            for nSup = [1, 5, 50]
                [epsReq, info] = rid.crit.demandedResolution( ...
                    rid.fields.nominalFields(nSup), "epsReq0", 0.2);
                tc.verifyEqual(epsReq, 0.2 * ones(nSup, 1));
                tc.verifyEqual(info.refinement, ones(nSup, 1));
                tc.verifyEqual(info.fracRefined, 0);
                tc.verifyEqual(info.nFloored, 0);
            end
        end

        function p11_belowNominalFieldsAreClippedNotCoarsened(tc)
            % adaptation is strictly additive. Fields below nominal leave the
            % baseline untouched instead of coarsening it.
            f = tc.fieldsFrom([0.01; 1], [0.1; 1], [0.5; 1]);
            [epsReq, info] = rid.crit.demandedResolution(f, "epsReq0", 0.2);
            tc.verifyEqual(info.refinement, [1; 1]);
            tc.verifyEqual(epsReq, [0.2; 0.2]);
        end

        function p12_certificationModeFloorsAndReportsBudgetModeDoesNot(tc)
            % Complexity 100 with beta = 0.5 gives refinement 0.1 -> epsReq = 0.01,
            % below fillDistanceQ = 0.05, so certification mode projects it and says
            % so; budget mode leaves the shape alone.
            f = struct("complexity", [100; 4; 1]);

            [epsReqBudget, infoBudget] = rid.crit.demandedResolution(f, ...
                "epsReq0", 0.1, "beta", 0.5, "mode", "budget", "fillDistanceQ", 0.05);
            tc.verifyEqual(epsReqBudget, [0.01; 0.05; 0.1], "AbsTol", tc.Tol);
            tc.verifyEqual(infoBudget.nFloored, 0);
            tc.verifyEqual(infoBudget.fracFloored, 0);
            tc.verifyEqual(infoBudget.headroom, 2, "AbsTol", tc.Tol);

            [epsReqCert, infoCert] = tc.verifyWarning( ...
                @() rid.crit.demandedResolution(f, "epsReq0", 0.1, "beta", 0.5, ...
                    "mode", "certification", "fillDistanceQ", 0.05), ...
                "rid:crit:resolutionFloorActive");
            tc.verifyEqual(epsReqCert, [0.05; 0.05; 0.1], "AbsTol", tc.Tol);
            tc.verifyEqual(infoCert.nFloored, 1);
            tc.verifyEqual(infoCert.fracFloored, 1/3, "AbsTol", tc.Tol);
            tc.verifyGreaterThanOrEqual(min(epsReqCert), 0.05);

            % The refinement factor itself is reported unfloored, so the projection
            % is visible as a difference rather than absorbed.
            tc.verifyEqual(infoCert.refinement, [0.1; 0.5; 1], "AbsTol", tc.Tol);
        end

        function p12_certificationModeIsSilentWhenTheFloorDoesNotAct(tc)
            % No projection, no warning: the clamp is reported when it acts and
            % never otherwise.
            f = struct("complexity", [4; 1]);
            [epsReq, info] = rid.crit.demandedResolution(f, "epsReq0", 0.1, ...
                "beta", 0.5, "mode", "certification", "fillDistanceQ", 0.001);
            tc.verifyEqual(epsReq, [0.05; 0.1], "AbsTol", tc.Tol);
            tc.verifyEqual(info.nFloored, 0);
            tc.verifyEqual(info.headroom, 100, "AbsTol", tc.Tol);
        end

        % --- validation -------------------------------------------------------

        function requiresEpsReq0(tc)
            tc.verifyError(@() rid.crit.demandedResolution(rid.fields.nominalFields(3)), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.demandedResolution( ...
                rid.fields.nominalFields(3), "epsReq0", 0), "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.demandedResolution( ...
                rid.fields.nominalFields(3), "epsReq0", -1), "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.demandedResolution( ...
                rid.fields.nominalFields(3), "epsReq0", NaN), "rid:crit:invalidResolution");
        end

        function rejectsUnknownFieldMembers(tc)
            % A misspelt member would otherwise default to nominal and leave a
            % region unrefined with no symptom at all.
            tc.verifyError(@() rid.crit.demandedResolution( ...
                struct("noiselevel", [2; 1]), "epsReq0", 0.2), ...
                "rid:crit:unknownField");
            tc.verifyError(@() rid.crit.demandedResolution( ...
                struct("complexity", [2; 1], "relevanceWeight", [1; 1]), ...
                "epsReq0", 0.2), "rid:crit:unknownField");
        end

        function rejectsAnEmptyFieldStruct(tc)
            % nSup is inferred from the members present, so there must be one.
            tc.verifyError(@() rid.crit.demandedResolution(struct(), "epsReq0", 0.2), ...
                "rid:crit:invalidResolution");
        end

        function rejectsNonPositiveFields(tc)
            % Exclusion is expressed through Q, never through a zero field.
            tc.verifyError(@() rid.crit.demandedResolution( ...
                struct("complexity", [0; 1]), "epsReq0", 0.2), "rid:crit:invalidField");
            tc.verifyError(@() rid.crit.demandedResolution( ...
                struct("relevance", [-1; 1]), "epsReq0", 0.2), "rid:crit:invalidField");
            tc.verifyError(@() rid.crit.demandedResolution( ...
                struct("noiseLevel", [NaN; 1]), "epsReq0", 0.2), "rid:crit:invalidField");
            tc.verifyError(@() rid.crit.demandedResolution( ...
                struct("noiseLevel", [Inf; 1]), "epsReq0", 0.2), "rid:crit:invalidField");
        end

        function rejectsFieldShapeMismatch(tc)
            tc.verifyError(@() rid.crit.demandedResolution( ...
                struct("complexity", [1; 1; 1], "relevance", [1; 1]), "epsReq0", 0.2), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.demandedResolution( ...
                struct("complexity", [1 1 1]), "epsReq0", 0.2), ...
                "rid:crit:invalidResolution");
        end

        function rejectsAnUnknownModeOrAMissingFillDistance(tc)
            tc.verifyError(@() rid.crit.demandedResolution( ...
                rid.fields.nominalFields(3), "epsReq0", 0.2, "mode", "certified"), ...
                "rid:crit:unknownMode");

            % Certification without a measured ruler: the floor is undefined.
            tc.verifyError(@() rid.crit.demandedResolution( ...
                rid.fields.nominalFields(3), "epsReq0", 0.2, "mode", "certification"), ...
                "rid:crit:invalidResolution");
            tc.verifyError(@() rid.crit.demandedResolution( ...
                rid.fields.nominalFields(3), "epsReq0", 0.2, "mode", "certification", ...
                "fillDistanceQ", 0), "rid:crit:invalidResolution");
        end

        function rejectsNegativeExponents(tc)
            % A negative beta inverts the design: with the clip at 1 it would refine
            % the EASY regions and leave the demanding ones at the baseline.
            tc.verifyError(@() rid.crit.demandedResolution( ...
                rid.fields.nominalFields(3), "epsReq0", 0.2, "beta", -0.5), ...
                "MATLAB:validators:mustBeNonnegative");
            tc.verifyError(@() rid.crit.demandedResolution( ...
                rid.fields.nominalFields(3), "epsReq0", 0.2, "betaNoise", -0.5), ...
                "MATLAB:validators:mustBeNonnegative");
        end

        function reportsProvenanceInInfo(tc)
            [~, info] = rid.crit.demandedResolution(struct("complexity", [4; 1]), ...
                "epsReq0", 0.2, "beta", 0.5, "betaNoise", 0.25, "mode", "budget");
            tc.verifyEqual(info.epsReq0, 0.2);
            tc.verifyEqual(info.beta, 0.5);
            tc.verifyEqual(info.betaNoise, 0.25);
            tc.verifyEqual(info.mode, "budget");
            tc.verifyEqual(info.nSup, 2);
            tc.verifyTrue(isnan(info.headroom));   % fillDistanceQ not supplied
        end

    end

end
