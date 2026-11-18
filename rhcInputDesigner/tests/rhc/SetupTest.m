classdef SetupTest < matlab.unittest.TestCase
%SETUPTEST  The problem object, and the coherence it enforces.
%
%   Two kinds of failure are checked, and only the second needs argument:
%
%     - a MISORDERED argument, which the concrete class validators of the single-setup rule turn
%       into an immediate typed error. Untestable in an interesting way -- it is
%       the language doing it -- but the point of the generator contract was to make it true of all
%       nine, so one case pins that the ninth is typed too;
%     - an INCOHERENT one, where every argument is of the right class and they
%       describe different problems. That is silent: the run proceeds, measures
%       coverage in the wrong units or against the wrong points, and reports a
%       resolution it never achieved. Each such case below would pass a
%       type-checking constructor.

    methods (Test)

        function aCoherentSetupIsAccepted(tc)
            % The other half of every refusal below: a check that always refuses
            % would pass those tests alone.
            f = SetupTest.fixture();
            setup = rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, f.fields, ...
                                  f.gen, f.sol, f.state, f.sch, f.fieldSch, f.qInfo);
            tc.verifyClass(setup, "rid.rhc.Setup");
            tc.verifyEqual(setup.Q, f.Q);
            tc.verifySameHandle(setup.state, f.state, ...
                "the coverage state is the one mutable object the loop advances");
            tc.verifySubstring(setup.describe(), "rid.rhc.Setup");
            tc.verifyEqual(setup.qInfo, f.qInfo);
        end

        function qInfoHasNoDefaultAndMustDeclareItsOrigin(tc)
            % An absent provenance is precisely the state the twelfth
            % argument removes: a Setup built without it does not construct.
            %
            % Fails if: qInfo is given a default, or an incomplete struct is
            % accepted silently.
            f = SetupTest.fixture();
            tc.verifyError(@() rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, ...
                f.fields, f.gen, f.sol, f.state, f.sch, f.fieldSch), ...
                "MATLAB:minrhs");
            tc.verifyError(@() rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, ...
                f.fields, f.gen, f.sol, f.state, f.sch, f.fieldSch, ...
                struct("origin", "supplied")), ...
                "rid:rhc:invalidQInfo");
            tc.verifyError(@() rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, ...
                f.fields, f.gen, f.sol, f.state, f.sch, f.fieldSch, ...
                struct("origin", "unknownValue", "method", "", "seed", NaN, ...
                       "hash", "")), ...
                "rid:rhc:invalidQInfo");
        end

        function aSuppliedOriginInventsNothing(tc)
            % "nothing is invented -- no method, no seed, NaN where a seed
            % would be." A method or a seed alongside origin = "supplied" would
            % assert a generator ran when none did.
            f = SetupTest.fixture();
            tc.verifyError(@() rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, ...
                f.fields, f.gen, f.sol, f.state, f.sch, f.fieldSch, ...
                struct("origin", "supplied", "method", "lhs", "seed", NaN, ...
                       "hash", "")), ...
                "rid:rhc:invalidQInfo");
            tc.verifyError(@() rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, ...
                f.fields, f.gen, f.sol, f.state, f.sch, f.fieldSch, ...
                struct("origin", "supplied", "method", "", "seed", 0, ...
                       "hash", "")), ...
                "rid:rhc:invalidQInfo");

            setup = rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, f.fields, ...
                f.gen, f.sol, f.state, f.sch, f.fieldSch, ...
                struct("origin", "supplied", "method", "", "seed", NaN, ...
                       "hash", ""));
            tc.verifyEqual(setup.qInfo.origin, "supplied");
        end

        function theGeneratorArgumentIsTypedLikeTheOtherEight(tc)
            % The generator contract exists for this: without rid.gen.Generator the ninth argument
            % could not be typed, and it sits between Q and state in a positional
            % signature -- the two neighbours most easily swapped with it.
            %
            % Fails if: the validator is relaxed to (1,1) so anything scalar passes.
            f = SetupTest.fixture();
            tc.verifyError(@() rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, ...
                f.fields, f.state, f.sol, f.state, f.sch, f.fieldSch, f.qInfo), ...
                "MATLAB:validation:UnableToConvert");
        end

        function aNormalizationOverADifferentBoxIsRefused(tc)
            % The normalization IS the ROI expressed as a metric. Built over
            % a different box it scales every distance while epsReq stays put, so
            % the run reports a resolution it never achieved -- and nothing raises.
            %
            % Fails if: only the DIMENSION COUNT is compared, which agrees here.
            f    = SetupTest.fixture();
            wide = rid.roi.Normalization([-2 2; -1 1; -4 4], f.nrm.dimName);
            tc.verifyError(@() rid.rhc.Setup(f.srg, f.rep, wide, f.roi, f.Q, ...
                f.fields, f.gen, f.sol, f.state, f.sch, f.fieldSch, f.qInfo), ...
                "rid:rhc:inconsistentSetup");
        end

        function supportingPointsOfTheWrongWidthAreRefused(tc)
            % Q is NORMALIZED coordinates with samples in rows. A transposed or
            % short Q would be measured against by pairwiseDistance without
            % complaint only if the widths happened to agree; they must be asserted.
            f = SetupTest.fixture();
            tc.verifyError(@() rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, ...
                f.Q(:, 1:2), f.fields, f.gen, f.sol, f.state, f.sch, f.fieldSch, ...
                f.qInfo), ...
                "rid:rhc:inconsistentSetup");
        end

        function aCoverageStateBuiltAgainstOtherPointsIsRefused(tc)
            % THE SILENT ONE. The state caches distances to the Q it was
            % constructed with; a different Q here leaves the loop scoring against
            % one set while the record, resolutionReport and every plot describe
            % another. Same shape, same class, different problem.
            %
            % Fails if: only size(Q) is compared rather than the points themselves.
            f      = SetupTest.fixture();
            otherQ = f.Q + 0.01;
            elsewhere = rid.crit.CoverageState(otherQ, 0.2 * ones(size(otherQ, 1), 1));
            tc.verifyError(@() rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, ...
                f.fields, f.gen, f.sol, elsewhere, f.sch, f.fieldSch, f.qInfo), ...
                "rid:rhc:inconsistentSetup");
        end

        function theOptimizerArgumentIsTypedAndSelectable(tc)
            % `sol` is inside `setup` because A REPLAY MUST HAVE THE SAME SOLVER OR
            % IT DOES NOT REPRODUCE -- the exact mirror of why `process` is outside
            %. Constructed inside `run` from `opts`, as it was, it
            % could not be swapped at all.
            %
            % Fails if: the validator is relaxed, or `sol` is ignored in favour of
            % a solver the loop builds for itself.
            f = SetupTest.fixture();
            tc.verifyError(@() rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, ...
                f.fields, f.gen, f.gen, f.state, f.sch, f.fieldSch, f.qInfo), ...
                "MATLAB:validation:UnableToConvert");

            grid  = rid.opt.GridSearch();
            setup = rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, f.fields, ...
                f.gen, grid, f.state, f.sch, f.fieldSch, f.qInfo);
            tc.verifyClass(setup.sol, "rid.opt.GridSearch");
        end

        function anIncompleteOrMisSizedFieldSetIsRefused(tc)
            % `fields` is (nSup x 1) DATA living beside Q, so it
            % is indexed by the supporting points and a mismatch is a different
            % problem silently substituted.
            %
            % Completeness comes from ONE definition, in the package that
            % constructs fields: an absent member would default to nominal
            % downstream and discard a declared field.
            f = SetupTest.fixture();

            partial = rmfield(f.fields, "relevance");
            tc.verifyError(@() rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, ...
                partial, f.gen, f.sol, f.state, f.sch, f.fieldSch, f.qInfo), ...
                "rid:fields:incompleteFields");

            short            = f.fields;
            short.complexity = short.complexity(1:end - 1);
            tc.verifyError(@() rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, ...
                short, f.gen, f.sol, f.state, f.sch, f.fieldSch, f.qInfo), ...
                "rid:rhc:inconsistentSetup");
        end

        function anIncompleteRoiDeclarationIsRefusedByTheRepresentation(tc)
            % every dimension is declared separately and nothing is
            % inherited. The diagnostic must come from the representation that owns
            % the coordinate names, not from a generic check here.
            f       = SetupTest.fixture();
            partial = rmfield(f.roi, "y_km1");
            tc.verifyError(@() rid.rhc.Setup(f.srg, f.rep, f.nrm, partial, f.Q, ...
                f.fields, f.gen, f.sol, f.state, f.sch, f.fieldSch, f.qInfo), ...
                "rid:rep:missingRange");
        end

        function aSurrogateWhoseStatesAreNotThisSpaceIsRefusedAtAssembly(tc)
            % The joint-space representation contract, the assembly clause. `rid.rep.StateSpace` declares
            % that its state coordinates ARE the surrogate's states, so a
            % surrogate of a different order cannot build this space at all.
            % NOTHING ELSE IN THIS CONSTRUCTOR NOTICES: nDim, the names, the box,
            % Q, the coverage state and the fields are all coherent here -- the
            % space is well declared and the model is simply not its. That is the
            % third silent failure of this file.
            %
            % REFUSED HERE RATHER THAN AT STEP 1 is the point of the separate
            % member: `build` asks the same question again, but by
            % then a run has started and, online, has already touched a process.
            %
            % Fails if: the `rep.assertCompatible(srg)` call is removed from the
            % constructor -- the incoherent pairing then constructs in silence.
            f = SetupTest.stateSpaceFixture(3);
            tc.verifyError(@() rid.rhc.Setup( ...
                SetupTest.stateSpaceSurrogate(2), f.rep, f.nrm, f.roi, f.Q, ...
                f.fields, f.gen, f.sol, f.state, f.sch, f.fieldSch, f.qInfo), ...
                "rid:rep:stateDimensionMismatch");

            % The other half, as everywhere here: the matched pairing constructs,
            % so what the check refuses is the MISMATCH and not the class.
            setup = SetupTest.setupFrom(f);
            tc.verifyClass(setup, "rid.rhc.Setup");
            tc.verifyEqual(setup.rep.nStates, 3);
        end

    end

    methods (Static)

        function f = fixture(opts)
            %FIXTURE  A coherent SISO problem: first-order surrogate, NARX space.
            %
            %   Shared with RunTest, which needs exactly the same objects. The ROI
            %   and the normalization are built from ONE box so that the coherence
            %   the constructor checks holds by construction rather than by
            %   coincidence.
            arguments
                opts.nHorizon   (1,1) double = 4
                opts.nSup       (1,1) double = 12
                opts.epsReq0    (1,1) double = 0.25
                opts.pastSource (1,1) string = "simulated"
                opts.mode       (1,1) string = "budget"
                opts.seed       (1,1) double = 0
                opts.nk         (1,1) double = rid.defaultOptions().rep.nk
            end

            % The ROI is keyed by the representation's OWN names, so that an input
            % delay moves the declaration with the coordinates instead of leaving a
            % stale one behind. At the default nk = 0 this is exactly
            % struct("u1_k", [-1 1], "du1_k", [-0.5 0.5], "y_km1", [-2 2]).
            f.rep = rid.rep.Narx(2, 1, "nk", opts.nk, "pastSource", opts.pastSource);
            f.roi = SetupTest.declareFor(f.rep);
            [box, dimName] = f.rep.bounds(f.roi);
            f.nrm = rid.roi.Normalization(box, dimName);
            f.srg = rid.surrogate.LinearFirstOrder(0.5, 0.1, 1);

            [f.Q, ~, f.qInfo] = rid.roi.supportingPoints(f.nrm, opts.nSup, ...
                "seed", opts.seed, "progress", "off");
            f.fields = rid.fields.nominalFields(opts.nSup);
            % The MEASURED fill distance, with the same probe settings rid.rhc.run
            % uses. It has to be the same number: in certification mode the floor
            % of the certifiability floor is defined by it, so an epsReq built from a different one
            % would be refused by the loop's consistency check -- correctly.
            f.fillDistanceQ = rid.roi.estimateFillDistance(f.Q, f.nrm);
            f.epsReq = rid.crit.demandedResolution(f.fields, ...
                "epsReq0",       opts.epsReq0, ...
                "mode",          opts.mode, ...
                "fillDistanceQ", f.fillDistanceQ);
            f.state  = rid.crit.CoverageState(f.Q, f.epsReq);

            f.gen      = rid.gen.Direct(rid.gen.Constraints(f.rep, f.roi), opts.nHorizon);
            f.sol      = rid.opt.MultiStart("seed", opts.seed);
            f.sch      = rid.surrogate.ImmediateScheduler();
            f.fieldSch = rid.fields.ImmediateScheduler();
        end

        function f = stateSpaceFixture(nStates, opts)
            %STATESPACEFIXTURE  The same coherent problem, in the joint input--state space.
            %
            %   The counterpart of `fixture` for rid.rep.StateSpace, built the
            %   same way: ONE box feeds the ROI and the normalization, so the
            %   coherence the constructor checks holds by construction, and the
            %   surrogate's order matches the declared state coordinates. What a
            %   test varies is then the one thing it is about.
            %
            %   `opts.srg` SUPPLIES A DIFFERENT MODEL OF THE SAME ORDER -- the
            %   nonlinear rid.surrogate.StateSpaceFunction, for the tests about
            %   which gradient tier the loop takes. The order still has to match
            %   the declared state coordinates, and `rid.rhc.Setup` says so if it
            %   does not.
            %
            %   `opts.stateRange` NARROWS THE DECLARED STATE BOX so that the
            %   model's steady-state locus leaves it, which is the case the
            %   guided levels of the joint-space representation contract must clip against. Empty means the
            %   shared default, so the one place that number lives stays
            %   `declareFor`.
            arguments
                nStates         (1,1) double {mustBeInteger, mustBePositive}
                opts.nHorizon   (1,1) double = 4
                opts.nSup       (1,1) double = 8
                opts.epsReq0    (1,1) double = 0.5
                opts.pastSource (1,1) string = "simulated"
                opts.seed       (1,1) double = 0
                opts.srg                     = []
                opts.stateRange (1,:) double = []
            end

            f.rep = rid.rep.StateSpace(nStates, "pastSource", opts.pastSource);
            if isempty(opts.stateRange)
                f.roi = SetupTest.declareFor(f.rep);
            else
                f.roi = SetupTest.declareFor(f.rep, opts.stateRange);
            end
            [box, dimName] = f.rep.bounds(f.roi);
            f.nrm = rid.roi.Normalization(box, dimName);
            f.srg = opts.srg;
            if isempty(f.srg)
                f.srg = SetupTest.stateSpaceSurrogate(nStates);
            end

            [f.Q, ~, f.qInfo] = rid.roi.supportingPoints(f.nrm, opts.nSup, ...
                "seed", opts.seed, "progress", "off");
            f.fields        = rid.fields.nominalFields(opts.nSup);
            f.fillDistanceQ = rid.roi.estimateFillDistance(f.Q, f.nrm);
            f.epsReq = rid.crit.demandedResolution(f.fields, ...
                "epsReq0",       opts.epsReq0, ...
                "mode",          "budget", ...
                "fillDistanceQ", f.fillDistanceQ);
            f.state  = rid.crit.CoverageState(f.Q, f.epsReq);

            f.gen      = rid.gen.Direct(rid.gen.Constraints(f.rep, f.roi), opts.nHorizon);
            f.sol      = rid.opt.MultiStart("seed", opts.seed);
            f.sch      = rid.surrogate.ImmediateScheduler();
            f.fieldSch = rid.fields.ImmediateScheduler();
        end

        function srg = stateSpaceSurrogate(nStates)
            %STATESPACESURROGATE  A stable diagonal model of the requested order.
            %
            %   Distinct poles well inside the unit circle and no zero entry in
            %   Bd, so every state is reached and the ORDER is the only thing that
            %   changes with the argument.
            arguments
                nStates (1,1) double {mustBeInRange(nStates, 1, 4)}
            end
            Ad  = diag(0.5 + 0.1 * (0:nStates - 1));
            Bd  = 0.1 * (1:nStates).';
            Cd  = [1, zeros(1, nStates - 1)];
            srg = rid.surrogate.LinearStateSpace(Ad, Bd, Cd, 0.1);
        end

        function roi = declareFor(rep, stateRange)
            %DECLAREFOR  One range per dimension, keyed by the coordinate's name.
            %
            %   Per-dimension and nothing derived; the range is chosen by
            %   the coordinate's KIND, never by parsing its name, so the
            %   same three widths follow the coordinates when nk moves them.
            %
            %   `stateRange` is an argument because the locus tests need a state
            %   box the model's rest locus does not fit inside; the DEFAULT is
            %   declared here and nowhere else.
            arguments
                rep
                stateRange (1,2) double = [-3 3]
            end
            roi   = struct();
            width = struct("inputLevel",      [-1 1], ...
                           "inputDifference", [-0.5 0.5], ...
                           "output",          [-2 2], ...
                           "state",           stateRange);
            for e = rep.coordinateInfo()
                roi.(e.name) = width.(e.kind);
            end
        end

        function setup = setupFrom(f)
            %SETUPFROM  The twelve arguments, in order.
            setup = rid.rhc.Setup(f.srg, f.rep, f.nrm, f.roi, f.Q, f.fields, ...
                                  f.gen, f.sol, f.state, f.sch, f.fieldSch, f.qInfo);
        end

    end

end
