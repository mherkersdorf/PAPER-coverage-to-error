classdef LinearStateSpaceTest < matlab.unittest.TestCase
%LINEARSTATESPACETEST  The order-n linear surrogate in physical units.
%
%       x(k+1)  = Ad x(k) + Bd u(k).' + ed
%       yHat(k) = Cd x(k) + directTerm u(k).' + fd
%
%   The model this class exists for is INTERPRETABLE: its states are the
%   process's, which is what lets a measured state be an initial condition for it
%   and a simulated state be compared with a measured one. The code
%   checks the dimension; nothing can check the basis. That is why no surrogate
%   offers a state-space representation as its `defaultRepresentation`.
%
%   The analytic case is the one that matters: a continuous second-order system
%   discretized by `fromContinuous` must reproduce its continuous step response
%   EXACTLY at the sample instants, both states. Euler would pass a "roughly
%   right" test and fail this one by three orders of magnitude.

    properties (Constant)
        Tol   = 1e-12
        W0    = 2      % undamped natural frequency of the analytic fixture
        Zeta  = 0.3    % damping ratio: underdamped, so the modes are complex
        TS    = 0.05
    end

    methods (Test)

        % --- the analytic case  -------

        function zohDiscretizationReproducesTheContinuousStepResponse(tc)
            % ZOH-EXACT, not approximate. `fromContinuous` exponentiates the
            % augmented matrix, so the recursion reproduces the continuous
            % solution at t = (i-1)*tS to machine precision -- for BOTH states,
            % which is what makes the state coordinates of rid.rep.StateSpace
            % mean what they say.
            %
            % Fails if: the discretization is Euler, Ad = I + Ac*tS. That is
            % wrong by O(tS^2) per step -- about 1e-3 here, nine orders above the
            % tolerance.
            [Ac, Bc, Cc] = LinearStateSpaceTest.secondOrder();
            srg = rid.surrogate.LinearStateSpace.fromContinuous( ...
                Ac, Bc, Cc, tc.TS);

            nSamples = 60;
            [yHat, ~, xTraj] = srg.simulateStates(ones(nSamples, 1), [0 0]);

            t = (0:nSamples - 1).' * tc.TS;
            [yRef, dyRef] = LinearStateSpaceTest.stepResponse(t, tc.W0, tc.Zeta);

            % x1 = y and x2 = dy/dt, both at the sample instants.
            tc.verifyEqual(xTraj, [yRef, dyRef], "AbsTol", tc.Tol);
            % Cc = [1 0] and the model is strictly proper, so yHat IS x1.
            tc.verifyEqual(yHat, yRef, "AbsTol", tc.Tol);
        end

        function expmIsUsedRatherThanAFirstOrderHold(tc)
            % The other half of the test above, stated as a property rather than
            % against a closed form: over ONE sample from an arbitrary state the
            % transition must be exp(Ac*tS), which pins the discretization
            % independently of the step-response algebra.
            [Ac, Bc, Cc] = LinearStateSpaceTest.secondOrder();
            srg = rid.surrogate.LinearStateSpace.fromContinuous( ...
                Ac, Bc, Cc, tc.TS);
            tc.verifyEqual(srg.Ad, expm(Ac * tc.TS), "AbsTol", tc.Tol);
            % And the input map is the ZOH integral, not Bc*tS.
            tc.verifyEqual(srg.Bd, Ac \ (srg.Ad - eye(2)) * Bc, "AbsTol", tc.Tol);
            tc.verifyEqual(srg.ed, zeros(2, 1));   % none declared
        end

        % --- the row convention --------------------------------

        function rowIIsTheStateBeforeTheIthInputActs(tc)
            % Row 1 is xStart BIT-EXACTLY, and xEnd -- the state after the last
            % row -- is NOT a row. The second assertion is the sharp one: it ties
            % the last row to xEnd through the model's own recursion, so a
            % trajectory shifted by one sample cannot satisfy both.
            %
            % Fails if: the rows hold x(k+1) instead of x(k).
            srg  = LinearStateSpaceTest.randomStable(41);
            rs   = RandStream("twister", "Seed", 42);
            uSeq = 2 * rand(rs, 7, 1) - 1;
            x0   = [0.3 -0.7];

            [~, xEnd, xTraj] = srg.simulateStates(uSeq, x0);

            tc.verifyEqual(xTraj(1, :), x0);
            tc.verifyEqual(size(xTraj), [7 2]);
            tc.verifyEqual(xEnd.', ...
                srg.Ad * xTraj(end, :).' + srg.Bd * uSeq(end, :).' + srg.ed, ...
                "AbsTol", tc.Tol);
        end

        function simulateAndSimulateStatesRunOnOneCore(tc)
            % the first two outputs ARE those of `simulate`, so the
            % two members must not be two implementations. Bit-for-bit.
            srg  = LinearStateSpaceTest.randomStable(41);
            rs   = RandStream("twister", "Seed", 43);
            uSeq = 2 * rand(rs, 20, 1) - 1;

            [ySim, xSim]          = srg.simulate(uSeq, [0.3 -0.7]);
            [yTraj, xTrajEnd, ~]  = srg.simulateStates(uSeq, [0.3 -0.7]);
            tc.verifyEqual(yTraj,    ySim);
            tc.verifyEqual(xTrajEnd, xSim);
        end

        % --- steadyState ----------------------------------------------

        function theRestStateIsAnalyticAndASimulationStartedThereStays(tc)
            % For the fixture, x_ss(u) = [u, 0]: position u, zero velocity. The
            % SECOND half is what makes this more than an algebra check -- a rest
            % state the recursion does not actually hold at is not a rest state.
            %
            % Fails if: the sign of `ed` is flipped in the solve. It is zero for
            % this fixture, so the affine case is checked separately below.
            [Ac, Bc, Cc] = LinearStateSpaceTest.secondOrder();
            srg = rid.surrogate.LinearStateSpace.fromContinuous( ...
                Ac, Bc, Cc, tc.TS);

            tc.verifyEqual(srg.steadyState(0.8), [0.8 0], "AbsTol", tc.Tol);

            xRest = srg.steadyState(0.8);
            [~, ~, xTraj] = srg.simulateStates(0.8 * ones(40, 1), xRest);
            tc.verifyEqual(xTraj, repmat(xRest, 40, 1), "AbsTol", tc.Tol);
        end

        function theAffineTermEntersTheRestStateWithTheRightSign(tc)
            % x_ss = (I - Ad) \ (Bd u + ed): `ed` is an INPUT to the solve, on the
            % same side as Bd*u.
            %
            % Fails if: the sign of ed is flipped, which the zero-ed fixture above
            % cannot see.
            srg = rid.surrogate.LinearStateSpace([0.5 0; 0 0.25], [1; 1], [1 0], ...
                0.1, "ed", [0.2; -0.4]);
            xRest = srg.steadyState(2);
            tc.verifyEqual(xRest, [(2 + 0.2) / 0.5, (2 - 0.4) / 0.75], ...
                "AbsTol", tc.Tol);

            [~, ~, xTraj] = srg.simulateStates(2 * ones(30, 1), xRest);
            tc.verifyEqual(xTraj, repmat(xRest, 30, 1), "AbsTol", tc.Tol);
        end

        % --- stability, refused at construction ----------------

        function anUnstableModelIsRefusedAtConstruction(tc)
            % An unstable surrogate makes the horizon simulation diverge and the
            % criterion meaningless, so it is
            % refused rather than carried.
            %
            % Fails if: the test is `<=` where it should be `<` -- a pole exactly
            % on the unit circle is not stable, and that is the case the strict
            % inequality exists for.
            tc.verifyError(@() rid.surrogate.LinearStateSpace( ...
                [1.01 0; 0 0.5], [1; 1], [1 0], 0.1), ...
                "rid:surrogate:unstableModel");
            tc.verifyError(@() rid.surrogate.LinearStateSpace( ...
                [1 0; 0 0.5], [1; 1], [1 0], 0.1), ...
                "rid:surrogate:unstableModel");
            % A complex pair on the unit circle: marginal, and equally refused.
            tc.verifyError(@() rid.surrogate.LinearStateSpace( ...
                [0 1; -1 0], [0; 1], [1 0], 0.1), ...
                "rid:surrogate:unstableModel");

            % The other half: a stable one constructs and says so.
            srg = rid.surrogate.LinearStateSpace([0.9 0; 0 0.5], [1; 1], [1 0], 0.1);
            tc.verifyTrue(srg.isStable());
        end

        % --- timeConstant: guidance, and it must be the SLOWEST mode ----------

        function theTimeConstantIsTheDecayTimeOfTheSlowestMode(tc)
            % `timeConstant` is guidance for choosing nHorizon, never used to
            % compute it -- but wrong guidance is worse than none, and the
            % horizon must be long enough for the SLOWEST transient to be visible.
            %
            % Fails if: `min` is used in place of `max` (this gives 0.5 here), or
            % the real part in place of the modulus (which the complex fixture
            % below separates).
            tS  = 0.05;
            srg = rid.surrogate.LinearStateSpace( ...
                diag([exp(-tS / 2), exp(-tS / 0.5)]), [1; 1], [1 0], tS);
            tc.verifyEqual(srg.timeConstant, 2, "RelTol", 1e-12);

            % COMPLEX modes: |exp((-zeta*w0 +- i*wd)*tS)| = exp(-zeta*w0*tS), so
            % the decay time is 1/(zeta*w0). The REAL PART is
            % exp(-zeta*w0*tS)*cos(wd*tS), a different number.
            [Ac, Bc, Cc] = LinearStateSpaceTest.secondOrder();
            osc = rid.surrogate.LinearStateSpace.fromContinuous(Ac, Bc, Cc, tc.TS);
            tc.verifyEqual(osc.timeConstant, 1 / (tc.Zeta * tc.W0), "RelTol", 1e-12);
        end

        function aDeadbeatModelHasATimeConstantOfZero(tc)
            % Spectral radius 0: the model settles in n steps and -tS/log(0) is
            % the limit 0. Stated rather than left to IEEE arithmetic.
            srg = rid.surrogate.LinearStateSpace([0 1; 0 0], [0; 1], [1 0], 0.1);
            tc.verifyEqual(srg.timeConstant, 0);
            tc.verifyTrue(srg.isStable());
        end

        % --- shapes and the affine terms --------------------------------------

        function misshapedMatricesAreRefusedByName(tc)
            good = {[0.5 0; 0 0.25], [1; 1], [1 0], 0.1};
            tc.verifyError(@() rid.surrogate.LinearStateSpace( ...
                [0.5 0 0; 0 0.25 0], good{2:4}), ...
                "rid:surrogate:dimensionMismatch", "Ad must be square");
            tc.verifyError(@() rid.surrogate.LinearStateSpace( ...
                good{1}, [1; 1; 1], good{3:4}), ...
                "rid:surrogate:dimensionMismatch", "Bd rows are states");
            tc.verifyError(@() rid.surrogate.LinearStateSpace( ...
                good{1:2}, [1 0 0], good{4}), ...
                "rid:surrogate:dimensionMismatch", "Cd is (1 x nStates)");
            tc.verifyError(@() rid.surrogate.LinearStateSpace( ...
                good{:}, "ed", [1; 2; 3]), ...
                "rid:surrogate:dimensionMismatch", "ed is (nStates x 1)");
            tc.verifyError(@() rid.surrogate.LinearStateSpace( ...
                good{:}, "directTerm", [1 2 3]), ...
                "rid:surrogate:dimensionMismatch", "one entry per channel");
        end

        function aScalarDirectTermExpandsAcrossChannels(tc)
            % The LinearFirstOrder convention, kept rather than re-invented.
            srg = rid.surrogate.LinearStateSpace([0.5 0; 0 0.25], ...
                [1 0; 0 1], [1 1], 0.1, "directTerm", 0.25);
            tc.verifyEqual(srg.directTerm, [0.25 0.25]);
            tc.verifyEqual(srg.nInputs, 2);
            tc.verifyEqual(srg.nStates, 2);
        end

        function theDirectTermReachesTheOutputAndNotTheState(tc)
            % A direct term acts on the OUTPUT, never on the state -- which is
            % exactly why a process can report x(N+1) before u(N+1) is decided
            %.
            plain  = rid.surrogate.LinearStateSpace([0.5 0; 0 0.25], [1; 1], [1 0], 0.1);
            direct = rid.surrogate.LinearStateSpace([0.5 0; 0 0.25], [1; 1], [1 0], ...
                0.1, "directTerm", 0.75, "fd", -0.2);

            uSeq = [0.4; -1.1; 0.9];
            [yPlain, ~, xPlain]  = plain.simulateStates(uSeq, [0.1 0.2]);
            [yDirect, ~, xDirect] = direct.simulateStates(uSeq, [0.1 0.2]);

            tc.verifyEqual(xDirect, xPlain, "AbsTol", tc.Tol);
            tc.verifyEqual(yDirect, yPlain + 0.75 * uSeq - 0.2, "AbsTol", tc.Tol);
        end

        % --- fitFromStates -------------------------------------

        function fitFromStatesRecoversANoiseFreeSystem(tc)
            % S1 is this fit and must be described as exactly this fit (the joint-space representation contract). Noise-free, so recovery is exact up to conditioning.
            %
            % Fails if: the intercept column is dropped -- `ed` and `fd` then
            % leak into Ad, Bd and Cd, and a model fitted around a non-zero
            % operating point comes out wrong while still looking plausible.
            truth = LinearStateSpaceTest.randomStable(41, "ed", [0.2; -0.15], ...
                "fd", 0.4);
            [u, x, y] = LinearStateSpaceTest.dataFrom(truth, 200, 44);

            fitted = rid.surrogate.LinearStateSpace.fitFromStates(u, x, y, truth.tS);

            tc.verifyEqual(fitted.Ad, truth.Ad, "AbsTol", 1e-10);
            tc.verifyEqual(fitted.Bd, truth.Bd, "AbsTol", 1e-10);
            tc.verifyEqual(fitted.ed, truth.ed, "AbsTol", 1e-10);
            tc.verifyEqual(fitted.Cd, truth.Cd, "AbsTol", 1e-10);
            tc.verifyEqual(fitted.fd, truth.fd, "AbsTol", 1e-10);
            tc.verifyEqual(fitted.tS, truth.tS);
            % NOT fitted unless asked for, and exactly zero when it is not: a
            % feedthrough path the user did not declare must not appear.
            tc.verifyEqual(fitted.directTerm, 0);
        end

        function fitFromStatesFindsADeclaredDirectTerm(tc)
            truth = LinearStateSpaceTest.randomStable(45, "directTerm", 0.6, ...
                "fd", -0.3);
            [u, x, y] = LinearStateSpaceTest.dataFrom(truth, 200, 46);

            fitted = rid.surrogate.LinearStateSpace.fitFromStates(u, x, y, ...
                truth.tS, "directTerm", true);
            tc.verifyEqual(fitted.directTerm, truth.directTerm, "AbsTol", 1e-10);
            tc.verifyEqual(fitted.Cd, truth.Cd, "AbsTol", 1e-10);
            tc.verifyEqual(fitted.fd, truth.fd, "AbsTol", 1e-10);
        end

        function fitFromStatesIsDeterministic(tc)
            truth = LinearStateSpaceTest.randomStable(41);
            [u, x, y] = LinearStateSpaceTest.dataFrom(truth, 120, 47);
            a = rid.surrogate.LinearStateSpace.fitFromStates(u, x, y, truth.tS);
            b = rid.surrogate.LinearStateSpace.fitFromStates(u, x, y, truth.tS);
            tc.verifyEqual(a.Ad, b.Ad);   % bit-identical, no tolerance
            tc.verifyEqual(a.Cd, b.Cd);
        end

        function aConstantInputCannotIdentifyTheModel(tc)
            % The unidentifiable case is DETECTED FROM THE DATA, not assumed from
            % the structure (the discipline): a constant input is collinear
            % with the intercept column, and the states have settled besides.
            truth = LinearStateSpaceTest.randomStable(41);
            u = ones(60, 1);
            [~, ~, x] = truth.simulateStates(u, truth.steadyState(1));
            y = x * truth.Cd.' + truth.fd;
            tc.verifyError(@() rid.surrogate.LinearStateSpace.fitFromStates( ...
                u, x, y, truth.tS), "rid:surrogate:insufficientData");
        end

        function tooFewRowsIsRefusedBeforeTheSolve(tc)
            truth = LinearStateSpaceTest.randomStable(41);
            [u, x, y] = LinearStateSpaceTest.dataFrom(truth, 6, 48);
            tc.verifyError(@() rid.surrogate.LinearStateSpace.fitFromStates( ...
                u, x, y, truth.tS), "rid:surrogate:insufficientData");
            % The other half: the same data pass once the floor allows them.
            tc.verifyWarningFree(@() rid.surrogate.LinearStateSpace.fitFromStates( ...
                u, x, y, truth.tS, "minSamples", 5));
        end

        function aFitThatComesOutUnstableIsRefused(tc)
            % The refusal comes from the CONSTRUCTOR, so the fit cannot hand back
            % a model the class would not accept from anyone else.
            nSamples = 60;
            rs = RandStream("twister", "Seed", 49);
            u  = 2 * rand(rs, nSamples, 1) - 1;
            x  = zeros(nSamples, 1);
            for kStep = 1:nSamples - 1
                x(kStep + 1) = 1.05 * x(kStep) + u(kStep);
            end
            tc.verifyError(@() rid.surrogate.LinearStateSpace.fitFromStates( ...
                u, x, x, 0.1), "rid:surrogate:unstableModel");
        end

        function fitFromStatesChecksItsShapes(tc)
            truth = LinearStateSpaceTest.randomStable(41);
            [u, x, y] = LinearStateSpaceTest.dataFrom(truth, 60, 50);
            tc.verifyError(@() rid.surrogate.LinearStateSpace.fitFromStates( ...
                u(1:end-1, :), x, y, truth.tS), "rid:surrogate:dimensionMismatch");
            tc.verifyError(@() rid.surrogate.LinearStateSpace.fitFromStates( ...
                u, x, y(1:end-1, :), truth.tS), "rid:surrogate:dimensionMismatch");
            tc.verifyError(@() rid.surrogate.LinearStateSpace.fitFromStates( ...
                u, x, [y, y], truth.tS), "rid:surrogate:dimensionMismatch");
        end

        % --- what this class deliberately does NOT do -------------------------

        function itCannotBeRefittedInTheLoop(tc)
            % canUpdate = false: an input--output refit would return
            % an arbitrary state basis, and this model's states must stay in the
            % PROCESS's coordinates or a measured state is not an initial condition
            % for it. The fit is `fitFromStates`, a static factory used OUTSIDE the
            % loop, on states that were genuinely measured.
            srg = LinearStateSpaceTest.randomStable(41);
            tc.verifyFalse(srg.canUpdate());
            same = srg.update(struct("u", 1, "x0", [0 0], "tS", srg.tS), struct());
            tc.verifyEqual(same.Ad, srg.Ad);
        end

        function itKeepsTheInheritedNarxDefaultRepresentation(tc)
            % a class cannot know whether its states are
            % INTERPRETABLE, so no surrogate offers a state-space representation as
            % its default. The inherited NARX default stands, and a state-space
            % representation is always constructed explicitly.
            srg = LinearStateSpaceTest.randomStable(41);   % nStates = 2
            rep = srg.defaultRepresentation();
            tc.verifyClass(rep, "rid.rep.Narx");
            tc.verifyEqual(rep.nLagU, 3);
            tc.verifyEqual(rep.nLagY, 2);
        end

        function theBatchIsTheInheritedFallbackAndAgreesWithIt(tc)
            % simulateBatch is NOT overridden: no caller needs a
            % vectorized batch yet, and an override would have to match the base
            % loop bit for bit. This pins that the inherited one works, including
            % the third output.
            srg = LinearStateSpaceTest.randomStable(41);
            rs  = RandStream("twister", "Seed", 51);
            uSeqBatch = 2 * rand(rs, 10, 1, 4) - 1;

            [Y, Xend, XTraj] = srg.simulateBatch(uSeqBatch, [0.3 -0.7]);
            tc.verifyEqual(size(Y),     [10 4]);
            tc.verifyEqual(size(Xend),  [4 2]);
            tc.verifyEqual(size(XTraj), [10 2 4]);
            for iSeq = 1:4
                [y, xEnd, xTraj] = srg.simulateStates(uSeqBatch(:, :, iSeq), [0.3 -0.7]);
                tc.verifyEqual(Y(:, iSeq),       y,     "AbsTol", tc.Tol);
                tc.verifyEqual(Xend(iSeq, :),    xEnd,  "AbsTol", tc.Tol);
                tc.verifyEqual(XTraj(:, :, iSeq), xTraj, "AbsTol", tc.Tol);
            end
        end

        function describeNamesTheOrderAndTheParameters(tc)
            srg = LinearStateSpaceTest.randomStable(41);
            s   = srg.describe();
            tc.verifyClass(s, "string");
            tc.verifyTrue(contains(s, "LinearStateSpace"));
            tc.verifyTrue(contains(s, "nStates = 2"));
        end

        function isValueClassNotHandle(tc)
            tc.verifyFalse(isa(LinearStateSpaceTest.randomStable(41), "handle"));
        end

    end

    methods (Static, Access = private)

        function [Ac, Bc, Cc] = secondOrder()
            %SECONDORDER  ydotdot + 2*zeta*w0*ydot + w0^2*y = w0^2*u, unit gain.
            %
            %   States are POSITION and VELOCITY -- physically interpretable, which
            %   is the premise a state-space representation rests on.
            w0   = LinearStateSpaceTest.W0;
            zeta = LinearStateSpaceTest.Zeta;
            Ac = [0, 1; -w0^2, -2 * zeta * w0];
            Bc = [0; w0^2];
            Cc = [1, 0];
        end

        function [y, dy] = stepResponse(t, w0, zeta)
            %STEPRESPONSE  The closed-form underdamped unit-step response and its
            %derivative -- written from the ODE, not from the code under test.
            wd = w0 * sqrt(1 - zeta^2);
            y  = 1 - exp(-zeta * w0 * t) .* ...
                     (cos(wd * t) + (zeta * w0 / wd) * sin(wd * t));
            dy = (w0 / sqrt(1 - zeta^2)) * exp(-zeta * w0 * t) .* sin(wd * t);
        end

        function srg = randomStable(seed, opts)
            %RANDOMSTABLE  A seeded two-state model with spectral radius 0.8.
            arguments
                seed (1,1) double
                opts.ed         (:,1) double = []
                opts.fd         (1,1) double = 0
                opts.directTerm (1,:) double = 0
            end
            rs = RandStream("twister", "Seed", seed);
            Ad = randn(rs, 2);
            Ad = 0.8 * Ad / max(abs(eig(Ad)));
            Bd = randn(rs, 2, 1);
            Cd = randn(rs, 1, 2);
            srg = rid.surrogate.LinearStateSpace(Ad, Bd, Cd, 0.1, ...
                "ed", opts.ed, "fd", opts.fd, "directTerm", opts.directTerm);
        end

        function [u, x, y] = dataFrom(srg, nSamples, seed)
            %DATAFROM  Noise-free input, state and output records of `srg`.
            %
            %   `x` is exactly what `fitFromStates` documents: row k is the state
            %   BEFORE u(k) acts.
            rs = RandStream("twister", "Seed", seed);
            u  = 2 * rand(rs, nSamples, 1) - 1;
            [y, ~, x] = srg.simulateStates(u, srg.steadyState(0));
        end

    end

end
