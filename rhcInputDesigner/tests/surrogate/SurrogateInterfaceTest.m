classdef SurrogateInterfaceTest < matlab.unittest.TestCase
%SURROGATEINTERFACETEST  The minimal interface and its documented fallbacks.
%
%   The claim under test is that a model which can ONLY simulate is a usable
%   surrogate. MinimalSurrogate implements the two abstract members and nothing
%   else; everything exercised here comes from rid.surrogate.Surrogate's fallbacks.
%   If this file needs a new member added to MinimalSurrogate, the interface has
%   grown and "any model architecture" has quietly become false.
%
%   MinimalSurrogate has TWO states throughout, so that anything assuming
%   nStates == 1 fails here rather than in the first experiment with a
%   higher-order model.

    properties (Constant)
        Tol = 1e-12
    end

    methods (Test)

        % --- simulateBatch: the loop fallback and its override --------

        function theBatchFallbackEqualsRepeatedSimulate(tc)
            srg  = MinimalSurrogate();
            rs   = RandStream("twister", "Seed", 21);
            uSeqBatch = rand(rs, 15, 1, 7);

            [Y, Xend] = srg.simulateBatch(uSeqBatch, [0.3 -0.2]);
            for iCand = 1:7
                [y, xEnd] = srg.simulate(uSeqBatch(:, :, iCand), [0.3 -0.2]);
                tc.verifyEqual(Y(:, iCand), y, "AbsTol", tc.Tol);
                tc.verifyEqual(Xend(iCand, :), xEnd, "AbsTol", tc.Tol);
            end
            tc.verifyEqual(size(Y), [15 7]);
            tc.verifyEqual(size(Xend), [7 2]);   % (nSeq x nStates), nStates = 2
        end

        function theVectorizedOverrideEqualsTheBaseLoopExactly(tc)
            % The batch-simulation contract lets a surrogate override the fallback when it can do better.
            % The override is only correct if it agrees with the loop it replaces,
            % so both are run on the same model and compared. ALL THREE outputs:
            % XTraj is part of the same contract, and an override that got
            % the first two right and stacked the third wrongly would otherwise
            % pass here.
            fast = rid.surrogate.LinearFirstOrder(1.3, 0.05, 2, "offset", -0.4, ...
                "directTerm", 0.25);
            slow = LinearFirstOrderLoopBatch(1.3, 0.05, 2, "offset", -0.4, ...
                "directTerm", 0.25);

            rs = RandStream("twister", "Seed", 22);
            uSeqBatch = 3 * rand(rs, 40, 1, 25) - 1;

            [Yfast, XfastEnd, XTrajFast] = fast.simulateBatch(uSeqBatch, 0.6);
            [Yslow, XslowEnd, XTrajSlow] = slow.simulateBatch(uSeqBatch, 0.6);

            tc.verifyEqual(Yfast, Yslow, "AbsTol", 1e-12, "RelTol", 1e-12);
            tc.verifyEqual(XfastEnd, XslowEnd, "AbsTol", 1e-12, "RelTol", 1e-12);
            tc.verifyEqual(XTrajFast, XTrajSlow, "AbsTol", 1e-12, "RelTol", 1e-12);
        end

        % --- simulateStates: the state trajectory --------------

        function theStateTrajectoryFallbackIsAnalyticAndRowOneIsXStart(tc)
            % THE ROW CONVENTION, which is the whole: row i is the
            % state BEFORE uSeq(i,:) acts, so row 1 is xStart and xEnd -- the state
            % after the last row -- is NOT a row. This pairs each state with the
            % input applied at the same instant, as yHat(k) = C x(k) + D u(k) and
            % x(k+1) = f(x(k), u(k)) do.
            %
            % MinimalSurrogate: x(k+1) = x(k).*[0.5 0.8] + u(k), y(k) = sum(x(k)).
            % From xStart = [1 2] under u = [0.5; -1; 2], by hand:
            %
            %   k=1  x = [1,    2    ]   y = 3
            %   k=2  x = [1,    2.1  ]   y = 3.1      [0.5, 1.6 ] + 0.5
            %   k=3  x = [-0.5, 0.68 ]   y = 0.18     [0.5, 1.68] - 1
            %   end  x = [1.75, 2.544]                [-0.25, 0.544] + 2
            %
            % Fails if: the loop records the state AFTER each step -- row 1 is then
            % [1, 2.1] and every row is shifted by one sample, an off-by-one with
            % no symptom of its own.
            srg = MinimalSurrogate();
            [yHat, xEnd, xTraj] = srg.simulateStates([0.5; -1; 2], [1 2]);

            tc.verifyEqual(xTraj, [1, 2; 1, 2.1; -0.5, 0.68], "AbsTol", tc.Tol);
            tc.verifyEqual(yHat,  [3; 3.1; 0.18],             "AbsTol", tc.Tol);
            tc.verifyEqual(xEnd,  [1.75, 2.544],              "AbsTol", tc.Tol);

            % BIT-EXACTLY, not to a tolerance: row 1 is xStart itself, never a
            % value that has been through the recursion.
            tc.verifyEqual(xTraj(1, :), [1 2]);
        end

        function simulateStatesAgreesWithSimulateBitForBit(tc)
            % The first two outputs of simulateStates ARE those of simulate
            %. A class that overrides one of them runs both on one
            % private core precisely so that this cannot drift; the fallback gets
            % it by construction. Bit-for-bit, since there is no reason for the two
            % to differ by even an ulp.
            %
            % Fails if: an override's yHat path is perturbed -- a direct term
            % applied in one and not the other, say.
            rs = RandStream("twister", "Seed", 31);
            uSeq = 2 * rand(rs, 12, 1) - 1;

            [models, starts] = SurrogateInterfaceTest.everyModel();
            for iModel = 1:numel(models)
                srg = models{iModel};
                [ySim, xSim]        = srg.simulate(uSeq, starts{iModel});
                [yTraj, xTrajEnd, ~] = srg.simulateStates(uSeq, starts{iModel});
                tc.verifyEqual(yTraj,    ySim, class(srg) + " yHat");
                tc.verifyEqual(xTrajEnd, xSim, class(srg) + " xEnd");
            end
        end

        function theStateTrajectoryFallbackNeedsNStates(tc)
            % The fallback must preallocate (nHorizon x nStates), so a model that
            % cannot say how many states it has gets the typed refusal rather than
            % a generic MATLAB error.
            srg = MinimalSurrogate(false);   % declares no state dimension
            tc.verifyError(@() srg.simulateStates(rand(4, 1), 0), ...
                "rid:surrogate:nStatesRequired");
        end

        function aZeroLengthInputIsLegalOnEveryPath(tc)
            % A representation that looks back zero samples has an
            % EMPTY step window, and the loop asks for exactly this once per step
            %. Narx(1, 0) already reached the case; the contract never said
            % it was allowed.
            %
            % Fails if: a class indexes u(1,:) before testing the length.
            [models, starts] = SurrogateInterfaceTest.everyModel();
            for iModel = 1:numel(models)
                srg = models{iModel};
                x0  = starts{iModel};
                n   = srg.nStates;
                who = class(srg);

                [yHat, xEnd] = srg.simulate(zeros(0, 1), x0);
                tc.verifyEqual(size(yHat), [0 1], who);
                tc.verifyEqual(xEnd, x0, who);

                [yHat, xEnd, xTraj] = srg.simulateStates(zeros(0, 1), x0);
                tc.verifyEqual(size(yHat),  [0 1], who);
                tc.verifyEqual(xEnd, x0, who);
                tc.verifyEqual(size(xTraj), [0 n], who);

                [Y, Xend, XTraj] = srg.simulateBatch(zeros(0, 1, 3), x0);
                tc.verifyEqual(size(Y),     [0 3],   who);
                tc.verifyEqual(size(Xend),  [3 n],   who);
                tc.verifyEqual(size(XTraj), [0 n 3], who);
            end
        end

        function theBatchThirdOutputStacksTrajectoriesAlongDimensionThree(tc)
            % The uSeqBatch dimension rule applied to XTraj: dimension 2 is the STATE,
            % dimension 3 the SEQUENCE. The two models differ in nStates -- 2 and 1
            % -- so a transposed or (nSeq x L x n) stacking cannot pass both.
            %
            % Fails if: the batch is stacked (nSeq x L x n), or the pages are
            % filled in the wrong order.
            rs = RandStream("twister", "Seed", 32);
            uSeqBatch = 2 * rand(rs, 9, 1, 3) - 1;

            [models, starts] = SurrogateInterfaceTest.everyModel();
            for iModel = 1:numel(models)
                srg = models{iModel};
                x0  = starts{iModel};
                n   = srg.nStates;
                who = class(srg);

                [~, ~, XTraj] = srg.simulateBatch(uSeqBatch, x0);
                tc.verifyEqual(size(XTraj), [9 n 3], who);
                for iSeq = 1:3
                    [~, ~, xTraj] = srg.simulateStates(uSeqBatch(:, :, iSeq), x0);
                    tc.verifyEqual(XTraj(:, :, iSeq), xTraj, ...
                        "AbsTol", tc.Tol, who + sprintf(" sequence %d", iSeq));
                end

                % A SINGLE sequence: MATLAB trims the trailing singleton and no
                % code may assert ndims == 3. What survives is (L x n).
                [~, ~, XOne] = srg.simulateBatch(uSeqBatch(:, :, 1), x0);
                tc.verifyEqual(size(XOne), [9 n], who);
            end
        end

        function theBatchLayoutIsAlwaysThreeDimensional(tc)
            % a 2-D (nHorizon x nSeq) convention would be ambiguous at
            % nInputs = 1 -- (20 x 3) could be three SISO candidates or one
            % three-input candidate -- and the ambiguity would stay invisible until
            % the first MISO run.
            srg = rid.surrogate.LinearFirstOrder(1, 0.1, 1);
            tc.verifyError(@() srg.simulateBatch(rand(10, 4), 0), ...
                "rid:surrogate:dimensionMismatch");
            tc.verifyWarningFree(@() srg.simulateBatch(rand(10, 1, 4), 0));
        end

        function anEmptyBatchIsHandled(tc)
            srg = rid.surrogate.LinearFirstOrder(1, 0.1, 1);
            [Y, Xend] = srg.simulateBatch(zeros(10, 1, 0), 0.5);
            tc.verifyEqual(size(Y), [10 0]);
            tc.verifyEqual(size(Xend), [0 1]);
        end

        % --- steadyState fallback ---------------------------------------------

        function theSteadyStateFallbackSimulatesToConvergence(tc)
            % x(k+1) = x(k).*[0.5 0.8] + u  =>  x_ss = u ./ [0.5 0.2].
            srg = MinimalSurrogate();
            x   = srg.steadyState(1);
            tc.verifyEqual(x, [1/0.5, 1/0.2], "RelTol", 1e-9);

            [~, xEnd] = srg.simulate(ones(500, 1), x);
            tc.verifyEqual(xEnd, x, "RelTol", 1e-9);
        end

        % --- the remaining fallbacks ------------------------------------------

        function isStableDefaultsToTrueUnverified(tc)
            % A model that cannot check its own stability cannot be assumed
            % unstable either; refusing to run every such model would exclude
            % exactly the vendor-supplied models the open interface exists for.
            tc.verifyTrue(MinimalSurrogate().isStable());
        end

        function canUpdateDefaultsToFalseAndUpdateIsANoOp(tc)
            % explicit rather than implied, so the event log records WHY the
            % surrogate never changes. Without the flag, "the model never updated"
            % is indistinguishable from a bug.
            srg  = MinimalSurrogate();
            tc.verifyFalse(srg.canUpdate());

            hist = wellFormedHistory(ones(50, 1), "y", ones(50, 1), "tS", 0.1);
            same = srg.update(hist, struct());
            tc.verifyEqual(class(same), class(srg));
            tc.verifyEqual(same.nStates, srg.nStates);
        end

        function timeConstantDefaultsToNaNAndIsOnlyEverGuidance(tc)
            % nHorizon is user-defined and fixed for
            % the run; the time constant is guidance for CHOOSING it, never used to
            % compute it. A NaN here must therefore never break anything.
            tc.verifyTrue(isnan(MinimalSurrogate().timeConstant));
        end

        function theDirectTermDefaultsToStrictlyProper(tc)
            tc.verifyEqual(MinimalSurrogate().directTerm, 0);
        end

        % --- defaultRepresentation is order-agnostic --------------------------

        function defaultRepresentationDerivesLagsFromNStatesNotFromOne(tc)
            % nStates is never assumed to be 1. With two states the default
            % representation must be nLagY = 2, nLagU = 3 -- if anything hard-codes
            % the first-order form, this is where it shows.
            srg = MinimalSurrogate();
            rep = srg.defaultRepresentation();
            tc.verifyEqual(rep.nLagY, 2);
            tc.verifyEqual(rep.nLagU, 3);
            tc.verifyEqual(rep.nDim, 1 * 3 + 2);
        end

        function theExtraInputCoordinateIsUOfKAndExistsOnlyWhenNkIsZero(tc)
            % The minimal surrogate interface. The "+ 1" in nLagU = nStates + 1 is u(k), the
            % coordinate that carries local PE outside the delay window -- so it
            % exists only when the representation contains u(k). At nk >= 1 the
            % input block is the delay window alone, u(k-nk) ... u(k-nk-m+1), and
            % nLagU = m: FOUR coordinates at nStates = 2, not five.
            %
            % Fails if: nStates + 1 is applied at every nk, which added a spurious
            % u(k-3) at nk = 1 and a spurious u(k-4) at nk = 2.
            srg = MinimalSurrogate();          % nStates = 2, nInputs = 1

            % The lag basis, so the names read as the input block itself rather
            % than as differences of it.
            rep0 = srg.defaultRepresentation("nk", 0, "inputDifference", false);
            tc.verifyEqual(rep0.nLagU, 3);
            tc.verifyEqual(rep0.requiredNames(), ...
                ["u1_k", "u1_km1", "u1_km2", "y_km1", "y_km2"]);

            rep1 = srg.defaultRepresentation("nk", 1, "inputDifference", false);
            tc.verifyEqual(rep1.nLagU, 2);
            tc.verifyEqual(rep1.requiredNames(), ...
                ["u1_km1", "u1_km2", "y_km1", "y_km2"]);

            rep2 = srg.defaultRepresentation("nk", 2, "inputDifference", false);
            tc.verifyEqual(rep2.nLagU, 2);
            tc.verifyEqual(rep2.requiredNames(), ...
                ["u1_km2", "u1_km3", "y_km1", "y_km2"]);

            % nLagY is untouched: the output block is the delay window and nk is
            % an INPUT delay.
            tc.verifyEqual([rep0.nLagY, rep1.nLagY, rep2.nLagY], [2 2 2]);
        end

        function aFirstOrderModelAtNkOneHasASingleInputCoordinate(tc)
            % The amendment's corner: nLagU = nStates = 1 at nk >= 1, so the
            % DEFAULT difference basis has no difference to form. It is refused
            % rather than switched off silently -- a setting that quietly does
            % nothing says something false. Without it the model gives
            % [u(k-1), y(k-1)], the first candidate space.
            srg = rid.surrogate.LinearFirstOrder(1, 0.1, 1);   % nStates = 1
            tc.verifyError(@() srg.defaultRepresentation("nk", 1), ...
                "rid:rep:invalidRepresentation");

            rep = srg.defaultRepresentation("nk", 1, "inputDifference", false);
            tc.verifyEqual(rep.nLagU, 1);
            tc.verifyEqual(rep.requiredNames(), ["u1_km1", "y_km1"]);
        end

        function defaultRepresentationNeedsNStatesAndSaysSo(tc)
            srg = MinimalSurrogate(false);   % declares no state dimension
            tc.verifyError(@() srg.defaultRepresentation(), ...
                "rid:surrogate:nStatesRequired");
            tc.verifyError(@() srg.steadyState(1), ...
                "rid:surrogate:nStatesRequired");
        end

        % --- value semantics ---------------------------------------------------

        function surrogatesAreValueObjects(tc)
            tc.verifyFalse(isa(MinimalSurrogate(), "handle"));
            tc.verifyFalse(isa(rid.surrogate.LinearFirstOrder(1, 0.1, 1), "handle"));
        end

    end

    methods (Static, Access = private)

        function [models, starts] = everyModel()
            %EVERYMODEL  Every concrete surrogate, with a start state for each.
            %
            %   The contract clauses of the joint-space representation contract hold for ALL of them, so the list
            %   lives in one place: a class added to the project and not to this
            %   list is a class the shared contract is never checked against.
            %
            %   The list deliberately spans the two implementation routes -- the
            %   base-class FALLBACK (MinimalSurrogate, which implements nothing but
            %   `simulate`) and every OVERRIDE -- and two state dimensions, so a
            %   transposed trajectory cannot pass on all of them.
            %
            %   The last two entries are the SAME MODEL by two routes: the matrices
            %   directly, and the matrices in `params` behind named step and output
            %   functions. Any clause they disagree on is the wrapper's doing.
            p = struct("Ad", [0.6 0.2; -0.1 0.45], "Bd", [1; 0.5], ...
                       "ed", [0.1; -0.2], "Cd", [1 -0.3], ...
                       "directTerm", 0.15, "fd", 0.3);
            models = { ...
                MinimalSurrogate(), ...                                  % n = 2
                rid.surrogate.LinearFirstOrder(1.3, 0.05, 2, ...         % n = 1
                    "offset", -0.4, "directTerm", 0.25), ...
                rid.surrogate.LinearStateSpace(p.Ad, p.Bd, p.Cd, 0.05, ...
                    "ed", p.ed, "fd", p.fd, "directTerm", p.directTerm), ...
                rid.surrogate.StateSpaceFunction(@linearStepFcn, ...      % n = 2
                    @linearOutputFcn, p, 2, 1, 0.05, ...
                    "directTerm", p.directTerm)};
            starts = {[0.3 -0.2], 0.6, [0.3 -0.7], [0.3 -0.7]};
        end

    end

end
