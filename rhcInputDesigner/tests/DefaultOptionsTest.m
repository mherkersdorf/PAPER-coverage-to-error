classdef DefaultOptionsTest < matlab.unittest.TestCase
%DEFAULTOPTIONSTEST  rid.defaultOptions is the single source of default values.
%
%   Guards the rule that no function carries a literal default of its own:
%   every option a +roi or +crit function
%   defaults must be present here, so that a run record's effective options struct
%   is complete and a changed default is one edit.
%
%   The ABSENCE tests matter as much as the presence ones. `epsReq0`, the penalty
%   weight and the rate limit are each absent on purpose, and each for a different
%   reason; without a test, "absent because it is derived" and "absent because
%   someone forgot" look identical to the next reader.

    methods (Test)

        function returnsAStructWithOneSectionPerImplementedPackage(tc)
            opts = rid.defaultOptions();
            tc.verifyClass(opts, "struct");
            tc.verifyTrue(isfield(opts, "roi"));
            tc.verifyTrue(isfield(opts, "crit"));
            tc.verifyTrue(isfield(opts, "rep"));
            tc.verifyTrue(isfield(opts, "surrogate"));
            tc.verifyTrue(isfield(opts, "gen"));
            tc.verifyTrue(isfield(opts, "opt"));
        end

        function carriesEveryRoiOptionTheRoiFunctionsDefault(tc)
            expected = ["method" "supportingSet" "improveMaxTime" "progress" ...
                        "seed" "probeSeed" "probeMethod" "nProbes" ...
                        "traceTolerance" "blockSize" ...
                        "lhsCriterion" "lhsIterations" "lhsSmooth"];
            actual = string(fieldnames(rid.defaultOptions().roi)).';
            tc.verifyEqual(sort(actual), sort(expected));
        end

        function roiDefaultsHaveTheDocumentedValuesAndTypes(tc)
            roi = rid.defaultOptions().roi;
            tc.verifyEqual(roi.method, "lhs");
            tc.verifyEqual(roi.supportingSet, "");          % none by default
            tc.verifyEqual(roi.improveMaxTime, 1800);      % the Latin-hypercube method option
            tc.verifyEqual(roi.progress, "auto");          % the liveness-signal rule
            tc.verifyEqual(roi.lhsCriterion, "maximin");
            tc.verifyClass(roi.lhsSmooth, "logical");
            tc.verifyGreaterThanOrEqual(roi.seed, 0);
            tc.verifyGreaterThanOrEqual(roi.probeSeed, 0);
            tc.verifyGreaterThan(roi.nProbes, 0);
            tc.verifyGreaterThan(roi.blockSize, 0);
            tc.verifyGreaterThanOrEqual(roi.traceTolerance, 0);   % the convergence-trace safeguard
            tc.verifyEqual(roi.probeMethod, "sobolCorners");      % the probe bracket
        end

        function carriesEveryCritOptionTheCritFunctionsDefault(tc)
            expected = ["pNorm" "psiMax" "convOrder" "beta" "betaNoise" "mode"];
            actual = string(fieldnames(rid.defaultOptions().crit)).';
            tc.verifyEqual(sort(actual), sort(expected));
        end

        function critDefaultsHaveTheDocumentedValuesAndTypes(tc)
            crit = rid.defaultOptions().crit;
            tc.verifyEqual(crit.pNorm, 1);
            tc.verifyEqual(crit.psiMax, 1);             % open
            tc.verifyEqual(crit.convOrder, 2);          % open
            tc.verifyEqual(crit.mode, "budget");        % the mode strings
        end

        function carriesEveryRepOptionTheRepClassesDefault(tc)
            expected = ["pastSource" "nk" "inputDifference" "nInputs"];
            actual = string(fieldnames(rid.defaultOptions().rep)).';
            tc.verifyEqual(sort(actual), sort(expected));
        end

        function repDefaultsAreTheDefaultCoordinatesOfD043(tc)
            % [u(k), du(k), yHat(k-1)]: nk = 0 represents the current input, and
            % the difference basis makes a rate limit an axis-aligned box rather
            % than a diagonal band.
            rep = rid.defaultOptions().rep;
            tc.verifyEqual(rep.nk, 0);
            tc.verifyTrue(rep.inputDifference);
            tc.verifyEqual(rep.pastSource, "measured");     % the pastSource option
            tc.verifyEqual(rep.nInputs, 1);
        end

        function carriesEverySurrogateOptionTheSurrogateClassesDefault(tc)
            expected = ["updateInterval" "minUpdateSamples" "steadyStateTol" ...
                        "steadyStateMaxSamples" "steadyStateBlock"];
            actual = string(fieldnames(rid.defaultOptions().surrogate)).';
            tc.verifyEqual(sort(actual), sort(expected));
        end

        function theUpdatePolicyIsOffByDefault(tc)
            % the update MECHANISM exists from the start, the POLICY does
            % not. A finite default would presuppose one.
            tc.verifyEqual(rid.defaultOptions().surrogate.updateInterval, Inf);
            tc.verifyGreaterThan(rid.defaultOptions().surrogate.minUpdateSamples, 0);
        end

        function betaIsDerivedFromConvOrderInThisFile(tc)
            % beta = 1/convOrder is derived at the point of definition, and
            % betaNoise defaults to beta. The consequence, which is the
            % reason the test exists: overriding convOrder at a CALL SITE does not
            % change beta -- a caller wanting another convergence order must
            % override beta too.
            crit = rid.defaultOptions().crit;
            tc.verifyEqual(crit.beta, 1 / crit.convOrder, "AbsTol", 1e-15);
            tc.verifyEqual(crit.betaNoise, crit.beta);
        end

        function epsReq0HasNoDefault(tc)
            % the baseline commitment is required. A default here would hide
            % the single most consequential choice in the design.
            tc.verifyFalse(isfield(rid.defaultOptions().crit, "epsReq0"));
        end

        function carriesEveryGenOptionTheGenClassesDefault(tc)
            % affineProbeStep and affineCheckTol joined this section when they
            % moved out of rid.gen.affineMap's arguments block: as literals there
            % they were absent from a record's effective options, so a replay could
            % not reconstruct which tolerance produced the affine/notAffine verdict
            % it was comparing against.
            expected = ["violationScale" "tailWeight" "seed" ...
                        "affineProbeStep" "affineCheckTol"];
            actual = string(fieldnames(rid.defaultOptions().gen)).';
            tc.verifyEqual(sort(actual), sort(expected));
        end

        function carriesEveryOptOptionTheOptimizersDefault(tc)
            expected = ["nStarts" "algorithm" "maxIterations" ...
                        "maxEvaluationsPerVariable" "optimalityTolerance" ...
                        "stepTolerance" "constraintTolerance" ...
                        "gridLevels" "maxGridPoints"];
            actual = string(fieldnames(rid.defaultOptions().opt)).';
            tc.verifyEqual(sort(actual), sort(expected));
        end

        function theSolverSettingsAreTheSolversOwnDocumentedDefaults(tc)
            % Citing a source is not inventing a number. These are fmincon's own
            % documented defaults at R2025b, reproduced here so that the effective
            % options struct of a run says what the search actually did.
            %
            % Fails if: a tolerance is quietly tightened or loosened, which would
            % change every run's result with nothing in the record to show it.
            opt = rid.defaultOptions().opt;
            tc.verifyEqual(opt.maxIterations, 400);
            tc.verifyEqual(opt.maxEvaluationsPerVariable, 100);
            tc.verifyEqual(opt.optimalityTolerance, 1e-6);
            tc.verifyEqual(opt.stepTolerance, 1e-6);
            tc.verifyEqual(opt.constraintTolerance, 1e-6);
            % The one project choice: the penalty is C1 and not C2, and
            % the method description records that SQP tolerates it.
            tc.verifyEqual(opt.algorithm, "sqp");
        end

        function theStartCountIsDeclaredWhereTheDecisionLives(tc)
            % how many local searches to run is a search-strategy setting,
            % so it lives under opt and NOT under gen, and
            % rid.gen.Direct.initialGuesses takes it as an argument with no default
            % of its own. One declaration, one decider.
            %
            % Fails if: a second copy reappears under gen, where it could drift.
            tc.verifyEqual(rid.defaultOptions().opt.nStarts, 4);
            tc.verifyFalse(isfield(rid.defaultOptions().gen, "nStarts"));

            gen  = rid.gen.Direct(rid.gen.Constraints(rid.rep.Narx(2, 1), ...
                       struct("u1_k", [-2 3], "du1_k", [-0.4 0.4], "y_km1", [-2 3])), 4);
            tc.verifyError(@() gen.initialGuesses(wellFormedHistory(0), zeros(0, 1)), ...
                "rid:gen:missingStartCount");
        end

        function genDefaultsHaveTheDocumentedValuesAndTypes(tc)
            gen = rid.defaultOptions().gen;
            tc.verifyEqual(gen.violationScale, 0.05);
            tc.verifyEqual(gen.tailWeight, 0.1);
            tc.verifyGreaterThanOrEqual(gen.seed, 0);
            % nStarts is NOT here: it is a search-strategy setting and lives at
            % opts.opt.nStarts, so +opt decides the count and +gen answers.
            tc.verifyFalse(isfield(gen, "nStarts"));
            % violationScale is a FRACTION of the admissible width; a value at or
            % above 1 would price a violation of the whole box at less than psiMax.
            tc.verifyGreaterThan(gen.violationScale, 0);
            tc.verifyLessThan(gen.violationScale, 1);
            % tailWeight is a RATIO between horizon positions: at 1 the weighting
            % vanishes, above 1 a violation at the horizon end would outweigh one
            % at the applied sample, which inverts the hard/soft constraint split's whole argument.
            tc.verifyGreaterThan(gen.tailWeight, 0);
            tc.verifyLessThanOrEqual(gen.tailWeight, 1);
        end

        function thePenaltyWeightIsNotAnOption(tc)
            % the weight is DERIVED as w = psiMax / violationScale^2 at the
            % point of use. A stored weight would keep its old value when psiMax
            % changed, silently decalibrating every penalty -- the specific failure
            % the violationScale parametrization exists to prevent.
            gen = rid.defaultOptions().gen;
            for forbidden = ["w" "weight" "penaltyWeight" "violationWeight"]
                tc.verifyFalse(isfield(gen, forbidden), ...
                    "The penalty weight is derived from psiMax and violationScale, " + ...
                    "never stored.");
            end
        end

        function psiMaxIsNotDuplicatedUnderGen(tc)
            % The penalty reads psiMax from opts.crit. A second copy under opts.gen
            % is a second name for one number and would let
            % the two drift apart between the criterion and the objective.
            tc.verifyFalse(isfield(rid.defaultOptions().gen, "psiMax"));
            tc.verifyTrue(isfield(rid.defaultOptions().crit, "psiMax"));
        end

        function theRateLimitIsNotAnOption(tc)
            % the rate limit is part of the constraint
            % DECLARATION, passed to rid.gen.Constraints as a name-value argument.
            % A default here would make it a derived quantity with a value nobody
            % declared, and would compete with the default that comes from
            % the declared admissible region.
            gen = rid.defaultOptions().gen;
            tc.verifyFalse(isfield(gen, "rateLimit"));
            tc.verifyFalse(isfield(gen, "duMax"));
            tc.verifyFalse(isfield(gen, "amplitude"));
        end

        function isPureAndRepeatable(tc)
            tc.verifyEqual(rid.defaultOptions(), rid.defaultOptions());
        end

    end

end
