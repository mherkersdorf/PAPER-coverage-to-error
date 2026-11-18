classdef StateSpaceFunctionTest < matlab.unittest.TestCase
%STATESPACEFUNCTIONTEST  Any discrete-time state-space model as a surrogate.
%
%       x(k+1) = stepFcn(x, u, params)
%       y(k)   = outputFcn(x, u, params)
%
%   THE CLASS EXISTS SO THAT S0 -- the paper's reference process, used as its own
%   surrogate -- is not a class in `src/` and not an ad-hoc wrapper in
%   `experiments/`. The first would put paper-specific models into the library;
%   the second would be untested, and S0 produces the paper's reference results.
%
%   THE ORACLE IS `LinearStateSpace`. Given linear step and output functions with
%   the same matrices, the two classes must agree exactly -- which is what makes
%   the evaluation ORDER inside the wrapper observable rather than a matter of
%   opinion.

    properties (Constant)
        Tol = 1e-12
    end

    methods (Test)

        function itReproducesLinearStateSpaceExactly(tc)
            % The whole class in one assertion: with a linear map in the functions
            % it IS LinearStateSpace, so any difference is the wrapper's doing.
            %
            % Fails if: the output is evaluated AFTER the step. y(k) reads x(k),
            % the state before u(k) acts; evaluating it after would report
            % y(k) = Cd x(k+1) + ... and shift the whole response by one sample.
            [wrapped, oracle] = StateSpaceFunctionTest.pair();
            rs   = RandStream("twister", "Seed", 61);
            uSeq = 2 * rand(rs, 25, 1) - 1;
            x0   = [0.3 -0.7];

            [yW, xW]      = wrapped.simulate(uSeq, x0);
            [yO, xO]      = oracle.simulate(uSeq, x0);
            tc.verifyEqual(yW, yO, "AbsTol", tc.Tol);
            tc.verifyEqual(xW, xO, "AbsTol", tc.Tol);

            [~, ~, xTrajW] = wrapped.simulateStates(uSeq, x0);
            [~, ~, xTrajO] = oracle.simulateStates(uSeq, x0);
            tc.verifyEqual(xTrajW, xTrajO, "AbsTol", tc.Tol);
            tc.verifyEqual(xTrajW(1, :), x0);   % the row convention, bit-exactly
        end

        function simulateAndSimulateStatesRunOnOneLoop(tc)
            % the first two outputs ARE those of `simulate`.
            [wrapped, ~] = StateSpaceFunctionTest.pair();
            rs   = RandStream("twister", "Seed", 62);
            uSeq = 2 * rand(rs, 15, 1) - 1;

            [ySim, xSim]         = wrapped.simulate(uSeq, [0.3 -0.7]);
            [yTraj, xTrajEnd, ~] = wrapped.simulateStates(uSeq, [0.3 -0.7]);
            tc.verifyEqual(yTraj,    ySim);     % bit-identical
            tc.verifyEqual(xTrajEnd, xSim);
        end

        function theSteadyStateFallbackAgreesWithTheAnalyticOne(tc)
            % `steadyState` is the BASE-CLASS fallback here -- simulate a constant
            % input to convergence -- and it must land where the analytic
            % solve of the equivalent model does. To the fallback's own tolerance,
            % not to 1e-12: it stops when the state stops moving.
            [wrapped, oracle] = StateSpaceFunctionTest.pair();
            tc.verifyEqual(wrapped.steadyState(0.8), oracle.steadyState(0.8), ...
                "RelTol", 1e-8);
        end

        function repeatedCallsAreBitIdentical(tc)
            % `simulate` is PURE and DETERMINISTIC. A wrapper around a user function is where a hidden
            % persistent variable would most easily get in.
            [wrapped, ~] = StateSpaceFunctionTest.pair();
            rs   = RandStream("twister", "Seed", 63);
            uSeq = 2 * rand(rs, 30, 1) - 1;
            [yA, xA] = wrapped.simulate(uSeq, [0.3 -0.7]);
            [yB, xB] = wrapped.simulate(uSeq, [0.3 -0.7]);
            tc.verifyEqual(yA, yB);
            tc.verifyEqual(xA, xB);
        end

        % --- what a run record must be able to rebuild ---------

        function anAnonymousHandleIsRefused(tc)
            % A record cannot see the values a closure captures, so an object built
            % from one could never be rebuilt on replay. Refused at construction,
            % where the user can still do something about it.
            %
            % Fails if: the check is removed -- the object then constructs happily
            % and the loss surfaces only when someone tries to replay the run.
            params = StateSpaceFunctionTest.linearParams();
            tc.verifyError(@() rid.surrogate.StateSpaceFunction( ...
                @(x, u, p) x, @linearOutputFcn, params, 2, 1, 0.1), ...
                "rid:surrogate:anonymousFunction");
            tc.verifyError(@() rid.surrogate.StateSpaceFunction( ...
                @linearStepFcn, @(x, u, p) 0, params, 2, 1, 0.1), ...
                "rid:surrogate:anonymousFunction");
        end

        function aHandleThatIsNotOnThePathIsRefused(tc)
            % The other half of the same requirement: a NAMED handle whose name
            % nothing resolves is no more rebuildable than a closure. The handle is
            % constructed by str2func, which does not require the function to
            % exist -- exactly how such a handle reaches a constructor in practice.
            params = StateSpaceFunctionTest.linearParams();
            tc.verifyError(@() rid.surrogate.StateSpaceFunction( ...
                str2func("noSuchStepFunctionAnywhere"), @linearOutputFcn, ...
                params, 2, 1, 0.1), "rid:surrogate:unresolvableFunction");

            % The other half: a resolvable named handle constructs.
            tc.verifyWarningFree(@() rid.surrogate.StateSpaceFunction( ...
                @linearStepFcn, @linearOutputFcn, params, 2, 1, 0.1));
        end

        function describeNamesBothFunctionsAndPrintsTheParameters(tc)
            % The record shows the parameters, and a reader can see WHICH functions
            % produced a run without opening the code.
            [wrapped, ~] = StateSpaceFunctionTest.pair();
            s = wrapped.describe();
            tc.verifyClass(s, "string");
            tc.verifyTrue(contains(s, "linearStepFcn"));
            tc.verifyTrue(contains(s, "linearOutputFcn"));
            tc.verifyTrue(contains(s, "Ad"), "the parameters must be printed");
        end

        % --- a misbehaving user function is caught, not propagated ------------

        function aStepFunctionReturningTheWrongShapeIsCaught(tc)
            % The user's function is the one thing this class cannot validate in
            % advance, so it is validated per sample. Without the check a (n x 1)
            % return would transpose silently into a different trajectory.
            params = StateSpaceFunctionTest.linearParams();
            srg = rid.surrogate.StateSpaceFunction(@transposingStepFcn, ...
                @linearOutputFcn, params, 2, 1, 0.1);
            tc.verifyError(@() srg.simulate([0.4; 0.2], [0.3 -0.7]), ...
                "rid:surrogate:dimensionMismatch");
        end

        function aNonFiniteStateIsCaughtWhereItHappens(tc)
            % A model that blows up mid-horizon must say so at the sample it
            % happened, not return a trajectory full of NaN for the criterion to
            % score. `rid.assertHistory` calls a non-finite value corruption for
            % the same reason.
            params = StateSpaceFunctionTest.linearParams();
            srg = rid.surrogate.StateSpaceFunction(@divergingStepFcn, ...
                @linearOutputFcn, params, 2, 1, 0.1);
            tc.verifyError(@() srg.simulate(ones(5, 1), [0.3 -0.7]), ...
                "rid:surrogate:nonFiniteSimulation");
        end

        % --- the inherited fallbacks ----------------------------------

        function theOptionalMembersAreTheDocumentedFallbacks(tc)
            % steadyState and isStable use the base-class fallbacks
            % and canUpdate is false. `timeConstant` is DECLARED by the user and
            % defaults to NaN -- never detected, because this class knows nothing
            % about the model it wraps. `directTerm` is declared for the same
            % reason.
            params = StateSpaceFunctionTest.linearParams();
            plain  = rid.surrogate.StateSpaceFunction(@linearStepFcn, ...
                @linearOutputFcn, params, 2, 1, 0.1);
            tc.verifyTrue(isnan(plain.timeConstant));
            tc.verifyEqual(plain.directTerm, 0);
            tc.verifyFalse(plain.canUpdate());
            tc.verifyTrue(plain.isStable());   % the unverified default

            declared = rid.surrogate.StateSpaceFunction(@linearStepFcn, ...
                @linearOutputFcn, params, 2, 1, 0.1, ...
                "timeConstant", 1.7, "directTerm", 0.25);
            tc.verifyEqual(declared.timeConstant, 1.7);
            tc.verifyEqual(declared.directTerm, 0.25);
        end

        function itKeepsTheInheritedNarxDefaultRepresentation(tc)
            % The joint-space representation contract, as for LinearStateSpace: a class cannot know whether
            % its states are interpretable, so it offers no state-space default.
            [wrapped, ~] = StateSpaceFunctionTest.pair();
            rep = wrapped.defaultRepresentation();
            tc.verifyClass(rep, "rid.rep.Narx");
            tc.verifyEqual(rep.nLagU, 3);   % nStates = 2, nk = 0
            tc.verifyEqual(rep.nLagY, 2);
        end

        function theBatchIsTheInheritedFallbackAndAgreesWithIt(tc)
            [wrapped, ~] = StateSpaceFunctionTest.pair();
            rs = RandStream("twister", "Seed", 64);
            uSeqBatch = 2 * rand(rs, 8, 1, 3) - 1;

            [Y, Xend, XTraj] = wrapped.simulateBatch(uSeqBatch, [0.3 -0.7]);
            tc.verifyEqual(size(Y),     [8 3]);
            tc.verifyEqual(size(Xend),  [3 2]);
            tc.verifyEqual(size(XTraj), [8 2 3]);
            for iSeq = 1:3
                [y, ~, xTraj] = wrapped.simulateStates(uSeqBatch(:, :, iSeq), [0.3 -0.7]);
                tc.verifyEqual(Y(:, iSeq),        y,     "AbsTol", tc.Tol);
                tc.verifyEqual(XTraj(:, :, iSeq), xTraj, "AbsTol", tc.Tol);
            end
        end

        function isValueClassNotHandle(tc)
            [wrapped, ~] = StateSpaceFunctionTest.pair();
            tc.verifyFalse(isa(wrapped, "handle"));
        end

    end

    methods (Static, Access = private)

        function params = linearParams()
            %LINEARPARAMS  One stable two-state model, as a parameter struct.
            params = struct("Ad", [0.6 0.2; -0.1 0.45], "Bd", [1; 0.5], ...
                            "ed", [0.1; -0.2], "Cd", [1 -0.3], ...
                            "directTerm", 0.15, "fd", 0.3);
        end

        function [wrapped, oracle] = pair()
            %PAIR  The same model as a wrapped function and as LinearStateSpace.
            p = StateSpaceFunctionTest.linearParams();
            wrapped = rid.surrogate.StateSpaceFunction(@linearStepFcn, ...
                @linearOutputFcn, p, 2, 1, 0.05, "directTerm", p.directTerm);
            oracle = rid.surrogate.LinearStateSpace(p.Ad, p.Bd, p.Cd, 0.05, ...
                "ed", p.ed, "fd", p.fd, "directTerm", p.directTerm);
        end

    end

end
