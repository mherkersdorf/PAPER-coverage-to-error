classdef GoldenVectorsTest < matlab.unittest.TestCase
%GOLDENVECTORSTEST  The criterion against frozen values from an independent reading.
%
%   Layer 4 of the project's test plan. The cases in tests/data/golden/ were produced in Python, before this implementation existed, and they catch the
%   one class of error no property test can reach: an implementation that is
%   CONSISTENTLY wrong. Drop the `- 1` from the penalty and permutation invariance,
%   monotonicity, the scaling equivalence, boundedness and locality all still hold;
%   only the value is wrong, everywhere.
%
%   NEVER hand-edit a value to make a test pass, and never consult the generator to
%   make the MATLAB agree. If a golden value looks wrong, say so and stop: one of
%   the two readings of the specification is wrong, and moving the number hides
%   which (tests/data/golden/README.md).
%
%   Each case is compared through BOTH implementations -- the naive
%   rid.crit.evaluate and the cached rid.crit.CoverageState -- so a golden vector
%   that agrees with one and not the other is caught here rather than left to the
%   oracle test.
%
%   Format, from the README: `inputs`, `expected` and `source`. Non-finite values
%   appear only as the scalar string "Inf"; where an intermediate would be infinite
%   the intermediate arrays are omitted and only J is given.

    properties (Constant)
        RelTol = 1e-12
    end

    methods (Static)

        function folder = goldenFolder()
            here   = fileparts(mfilename("fullpath"));
            folder = fullfile(fileparts(here), "data", "golden");
        end

        function value = numericOrInf(raw)
            %NUMERICORINF  A scalar that JSON carries either as a number or as "Inf".
            if isstring(raw) || ischar(raw)
                value = str2double(raw);
            else
                value = double(raw);
            end
        end

    end

    methods (Test)

        function everyGoldenCaseReproduces(tc)
            folder = GoldenVectorsTest.goldenFolder();
            files  = dir(fullfile(folder, "*.json"));

            % A HARD FAILURE, not an assumption. This was `assumeNotEmpty`, which
            % filters the test out and leaves it Incomplete -- and `assertSuccess`
            % checks only `Failed`, so a tree without the vectors reported GREEN
            % with the entire golden layer absent. Combined with tests/data/ having
            % been git-ignored, that is what this layer did on every machine except
            % the one the vectors were generated on.
            %
            % Missing golden data is a broken checkout, not a configuration: layer 4
            % of the project's test plan is not optional, and a test layer that can silently
            % remove itself is worse than no test layer, because the green suite is
            % read as evidence.
            tc.assertNotEmpty(files, ...
                sprintf("no golden vectors in %s -- the golden layer cannot " + ...
                    "silently skip", folder));

            for iFile = 1:numel(files)
                tc.verifyOneCase(fullfile(folder, files(iFile).name));
            end
        end

    end

    methods (Access = private)

        function verifyOneCase(tc, path)
            golden = jsondecode(fileread(path));
            [~, name] = fileparts(path);
            what = @(part) sprintf("golden case ""%s"" (%s): %s", name, ...
                golden.source.generator, part);

            in     = golden.inputs;
            expect = golden.expected;

            Q      = in.Q;
            nSup   = size(Q, 1);
            nDim   = size(Q, 2);
            Znorm   = in.Znorm;
            if isempty(Znorm)
                Znorm = zeros(0, nDim);
            end
            psiMax = GoldenVectorsTest.numericOrInf(in.psiMax);
            pNorm  = GoldenVectorsTest.numericOrInf(in.pNorm);

            % --- the demanded resolution -------------------------------------
            fields = struct("relevance", in.relevance, "complexity", in.complexity, ...
                            "noiseLevel", in.noiseLevel);
            args = {"epsReq0", in.epsReq0, "beta", in.beta, "betaNoise", in.betaNoise, ...
                    "mode", string(in.mode), "fillDistanceQ", in.fillDistanceQ};

            % The floor warns whenever it acts; that is asserted in
            % DemandedResolutionTest, and here it would only add noise.
            warnState = warning("off", "rid:crit:resolutionFloorActive");
            restore   = onCleanup(@() warning(warnState));
            [epsReq, info] = rid.crit.demandedResolution(fields, args{:});
            clear restore

            tc.verifyEqual(info.refinement, expect.refinement, ...
                "RelTol", tc.RelTol, what("refinement"));
            tc.verifyEqual(epsReq, expect.epsReq, ...
                "RelTol", tc.RelTol, what("epsReq"));

            % --- the criterion, through the naive path -----------------------
            [J, critInfo] = rid.crit.evaluate(Znorm, Q, epsReq, ...
                "psiMax", psiMax, "pNorm", pNorm);

            if isfield(expect, "dRaw")
                tc.verifyEqual(critInfo.dRaw,   expect.dRaw,   "RelTol", tc.RelTol, what("dRaw"));
                tc.verifyEqual(critInfo.dTilde, expect.dTilde, "RelTol", tc.RelTol, what("dTilde"));
                tc.verifyEqual(critInfo.psi,    expect.psi,    "RelTol", tc.RelTol, ...
                    "AbsTol", 1e-14, what("psi"));
            end
            tc.verifyEqual(J, expect.J, "RelTol", tc.RelTol, "AbsTol", 1e-14, what("J"));

            % --- and through the cached path ---------------------------------
            state = rid.crit.CoverageState(Q, epsReq, "psiMax", psiMax, "pNorm", pNorm);
            state.add(Znorm);
            tc.verifyEqual(state.value(), expect.J, ...
                "RelTol", tc.RelTol, "AbsTol", 1e-14, what("J via CoverageState"));

            tc.verifyEqual(nSup, numel(expect.epsReq), what("nSup"));
        end

    end

end
