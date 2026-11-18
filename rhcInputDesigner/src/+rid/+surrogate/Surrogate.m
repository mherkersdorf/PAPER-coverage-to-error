classdef Surrogate
%SURROGATE  Abstract contract for a surrogate model. Deliberately minimal.
%
%   The surrogate maps a candidate input sequence to the FUTURE OUTPUT (and state)
%   TRAJECTORY from which the representation then builds the points in regressor
%   space that the criterion scores. It is
%   not assumed to be structurally correct and it is not a predictive model of the
%   process: it must place candidate points in roughly the right region, not
%   predict the output.
%
%   THE SURROGATE IS NOT THE REGRESSOR SPACE. Its internal structure -- order,
%   states, linearity, physical model or black box -- is its own concern and is
%   independent of the coordinates in which coverage is measured.
%   `defaultRepresentation` offers the representation matching a model's structure
%   as a CONVENIENCE; any representation may be paired with any surrogate.
%
%   REQUIRED VERSUS OPTIONAL MEMBERS. Only `simulate` and `describe` are
%   abstract, and `defaultRepresentation` is provided here in an order-agnostic
%   form that a subclass may override. A vendor-supplied model -- an FMU, a
%   compiled simulator, a lookup-table model -- can usually provide simulation and
%   little else, and a fat abstract class is precisely what would make "any model
%   architecture" false in practice. Everything else is a capability with a
%   documented fallback:
%
%       steadyState(uConst)  simulate a constant input to convergence
%       timeConstant         NaN -- guidance for choosing nHorizon, never used to
%                            compute it
%       isStable()           true, UNVERIFIED -- see the method
%       canUpdate()/update() false, and a no-op returning the model unchanged
%       nStates              required only when defaultRepresentation is not
%                            overridden, and NEVER assumed to be 1
%       directTerm           0, i.e. strictly proper
%       simulateStates       step `simulate` one sample at a time and record each
%                            state; needs a finite nStates
%
%   SIMULATION CONTRACT (section 3):
%     - `simulate` is PURE and DETERMINISTIC: state in, state out, no noise, no
%       persistent variables. Candidate scoring calls it thousands of times per
%       step; hidden state is a bug factory.
%     - The surrogate works in PHYSICAL units. Normalization happens once, at the
%       representation boundary, never inside a model.
%     - VALUE SEMANTICS: `update` returns a NEW object. Deliberately unlike
%       CoverageState, which is a handle. Surrogate changes are rare and must be
%       explicit and traceable.
%     - `simulateBatch` is part of the contract. The loop fallback below
%       makes it free for every subclass; a subclass overrides it only when it can
%       do better.
%
%   xStart IS REQUIRED, NEVER DEFAULTED. "State in, state out" is what
%   makes `simulate` pure, and a defaulted initial condition would silently decide
%   the first horizon's trajectory. The canonical way to produce one is
%   `steadyState(uConst)`: an experiment starting from rest starts at the rest
%   state of its first applied input. Zero would be a claim that the operating
%   point is at the origin, which is true only by accident of units.
%
%   A ZERO-LENGTH INPUT IS LEGAL on every path: `simulate` returns
%   (0 x 1) and `xStart`, `simulateStates` additionally (0 x nStates). A
%   representation that looks back zero samples has an empty step window, and the
%   loop asks for exactly this once per step.
%
%   SHAPES. Samples in rows throughout:
%       uSeq       (nHorizon x nInputs)
%       uSeqBatch  (nHorizon x nInputs x nSeq)   dimension 3 counts SEQUENCES
%       yHat       (nHorizon x 1)
%       Y          (nHorizon x nSeq)
%       xStart     (1 x nStates)
%       xEnd       (1 x nStates), or (nSeq x nStates) from simulateBatch
%       xTraj      (nHorizon x nStates)          row i is the state BEFORE
%                                                uSeq(i,:) acts, so row 1 is
%                                                xStart and xEnd is not a row
%       XTraj      (nHorizon x nStates x nSeq)   dimension 2 the STATE,
%                                                dimension 3 the SEQUENCE

    properties (SetAccess = protected)
        % Capabilities, with the fallbacks documented above. A subclass assigns
        % whichever of these it can honestly state.
        nStates      (1,1) double = NaN   % NEVER assumed to be 1
        nInputs      (1,1) double = NaN
        tS           (1,1) double = NaN   % sample time, PHYSICAL
        timeConstant (1,1) double = NaN   % dominant time constant, PHYSICAL
        directTerm   (1,:) double = 0     % 0 = strictly proper
    end

    methods (Abstract)
        [yHat, xEnd] = simulate(obj, uSeq, xStart)
        s              = describe(obj)
    end

    methods

        function [yHat, xEnd, xTraj] = simulateStates(obj, uSeq, xStart)
            %SIMULATESTATES  The output trajectory AND the state trajectory.
            %
            %   THE FIRST TWO OUTPUTS ARE EXACTLY THOSE OF `simulate`. This is not
            %   a second simulation path: the fallback below steps `simulate` one
            %   sample at a time, and a subclass that overrides this runs both on
            %   ONE private core, so the two cannot drift.
            %
            %   ROW CONVENTION: row `i` of `xTraj` is the state BEFORE `uSeq(i,:)`
            %   acts. So `xTraj(1,:) = xStart`, and `xEnd` -- the state after the
            %   last row, x(L) in the method description -- is NOT a row. This pairs each state
            %   with the input applied at the same instant, as
            %   yHat(k) = C x(k) + D u(k) and x(k+1) = f(x(k), u(k)) do. Rows
            %   holding x(k+1) were rejected: row k would then carry what a
            %   downstream model PREDICTS, not what it uses.
            %
            %   WHY A SEPARATE MEMBER rather than a third output on `simulate`: the
            %   base class cannot supply a fallback for an output of an ABSTRACT
            %   method, so every existing surrogate would have had to change, and
            %   an old one would fail with a generic MATLAB error instead of a typed
            %   one.
            %
            %   THE PER-SAMPLE LOOP IS THE FALLBACK, NOT THE PATH. It is correct for
            %   any surrogate but costs one call per sample; a model that can do
            %   better overrides this.
            %
            %   A ZERO-LENGTH INPUT IS LEGAL and returns (0 x 1), `xStart` and
            %   (0 x nStates). A representation that looks back zero samples has an
            %   empty step window and the loop asks for exactly this once per step
            %  .
            %
            %   Inputs:
            %       uSeq    (nHorizon x nInputs)  physical input
            %       xStart  (1 x nStates)         initial state; REQUIRED
            %   Outputs:
            %       yHat    (nHorizon x 1)        physical output
            %       xEnd    (1 x nStates)         state AFTER the last sample
            %       xTraj   (nHorizon x nStates)  row i = state before uSeq(i,:)
            arguments
                obj
                uSeq   (:,:) double {mustBeReal}
                xStart (1,:) double {mustBeReal}
            end
            % BEFORE assertStateShape, which would otherwise report a dimension
            % mismatch against NaN for a model whose real problem is that it
            % declares no state dimension at all.
            if ~isfinite(obj.nStates)
                error("rid:surrogate:nStatesRequired", ...
                    "simulateStates returns a (nHorizon x nStates) trajectory, " + ...
                    "and this model does not declare nStates. Either set " + ...
                    "nStates or override simulateStates.");
            end
            obj.assertStateShape(xStart);

            nHorizon = size(uSeq, 1);
            yHat     = zeros(nHorizon, 1);
            xTraj    = zeros(nHorizon, obj.nStates);
            x        = xStart;
            for kStep = 1:nHorizon
                % RECORDED BEFORE THE STEP, which is the row convention above.
                xTraj(kStep, :) = x;
                [y, x]          = obj.simulate(uSeq(kStep, :), x);
                yHat(kStep)     = y;
            end
            xEnd = x;
        end

        function [Y, Xend, XTraj] = simulateBatch(obj, uSeqBatch, xStart)
            %SIMULATEBATCH  Loop fallback over candidates.
            %
            %   Candidate scoring calls the surrogate thousands of times per step,
            %   and a surrogate whose per-call overhead dominates -- a Gaussian
            %   process, an ensemble -- is unusable without batching. This fallback
            %   makes the member free for every subclass; a subclass overrides it
            %   only when it can genuinely do better, and its override must agree
            %   with this loop exactly.
            %
            %   `uSeqBatch` is (nHorizon x nInputs x nSeq): the array is a batch of
            %   uSeq, so dimension 3 counts SEQUENCES.
            %
            %   IT IS NOT "ALWAYS 3-D", and it cannot be. MATLAB trims trailing
            %   singleton dimensions, so a single-sequence batch arrives here as a
            %   2-D (nHorizon x nInputs) array and no caller can make it otherwise.
            %   A middle singleton is never trimmed, so a SISO batch of several
            %   sequences is (nHorizon x 1 x nSeq) and genuinely 3-D.
            %
            %   What removes the ambiguity the uSeqBatch dimension rule is about is `size(uSeqBatch, 2) ==
            %   nInputs`, not the number of dimensions: it is that check which
            %   distinguishes three SISO sequences from one three-input sequence,
            %   and assertBatchShape has always tested exactly it.
            %
            %   THE OPTIONAL THIRD OUTPUT, `XTraj`, is a batch of `xTraj` stacked
            %   along dimension 3 -- the uSeqBatch dimension rule again, with the same caveat that
            %   MATLAB drops a trailing singleton. It is built ONLY
            %   when it is asked for: the loop below then steps `simulateStates`
            %   instead of `simulate`, which costs one call per sample for a model
            %   without an override, and no caller should pay that for a trajectory
            %   it does not want.
            %
            %   Inputs:
            %       uSeqBatch  (nHorizon x nInputs x nSeq)  physical inputs
            %       xStart     (1 x nStates)  shared initial state
            %   Outputs:
            %       Y      (nHorizon x nSeq)
            %       Xend   (nSeq x nStates)
            %       XTraj  (nHorizon x nStates x nSeq)  dimension 2 is the STATE,
            %                                           dimension 3 the SEQUENCE
            arguments
                obj
                uSeqBatch (:,:,:) double {mustBeReal}
                xStart    (1,:) double {mustBeReal}
            end
            obj.assertBatchShape(uSeqBatch);

            [nHorizon, ~, nSeq] = size(uSeqBatch);
            wantTraj = nargout >= 3;
            Y = zeros(nHorizon, nSeq);
            if wantTraj
                % simulateStates raises this too, but the preallocation below
                % needs the dimension first and would otherwise fail with a
                % generic error on a zero-sequence batch.
                if ~isfinite(obj.nStates)
                    error("rid:surrogate:nStatesRequired", ...
                        "The third output of simulateBatch is a " + ...
                        "(nHorizon x nStates x nSeq) array, and this model does " + ...
                        "not declare nStates.");
                end
                XTraj = zeros(nHorizon, obj.nStates, nSeq);
            end
            for iCand = 1:nSeq
                if wantTraj
                    [yHat, xEnd, xTraj] = obj.simulateStates( ...
                        uSeqBatch(:, :, iCand), xStart);
                    XTraj(:, :, iCand) = xTraj;
                else
                    [yHat, xEnd] = obj.simulate(uSeqBatch(:, :, iCand), xStart);
                end
                Y(:, iCand) = yHat;
                if iCand == 1
                    Xend = zeros(nSeq, numel(xEnd));
                end
                Xend(iCand, :) = xEnd;
            end
            if nSeq == 0
                Xend = zeros(0, numel(xStart));
            end
        end

        function rep = defaultRepresentation(obj, opts)
            %DEFAULTREPRESENTATION  The representation matching this model's order.
            %
            %   A CONVENIENCE, NOT A CONSTRAINT: any representation may be paired with any surrogate, and
            %   a fifty-state model behind a three-dimensional regressor space is a
            %   legitimate and expected pairing.
            %
            %   Order-agnostic by construction: the lag orders come from `nStates`,
            %   so a higher-order or cascade model is a pure addition later rather
            %   than a change to this interface. nStates is NEVER assumed to be 1
            %  .
            %
            %   THE INPUT BLOCK DECOMPOSES INTO TWO PARTS WITH DIFFERENT JOBS, and
            %   the rule follows from that rather than from an offset:
            %
            %       u(k)                        the CURRENT input -- carries local
            %                                   PE, and is OUTSIDE the delay window
            %       u(k-nk) ... u(k-nk-m+1),    the uniform delay window, m = nStates
            %       y(k-1)  ... y(k-m)          -- the state reconstruction
            %
            %   so nLagY = m always, and nLagU = m + 1 only where u(k) is actually a
            %   coordinate. At nk >= 1 it is not: the input block is the delay
            %   window alone and nLagU = m. Applying m + 1 at every nk added a
            %   spurious deepest input lag -- u(k-3) at nk = 1 for a two-state model
            %   -- and no test covered the case (the minimal surrogate interface).
            %
            %   A FIRST-ORDER MODEL AT nk >= 1 HAS A SINGLE INPUT COORDINATE, so the
            %   default difference basis has no difference to form and rid.rep.Narx
            %   refuses it (`rid:rep:invalidRepresentation`) rather than switching it
            %   off silently. With `inputDifference = false` it gives
            %   [u(k-1), y(k-1)], the first candidate space.
            %
            %   rid.rep.StateSpace is unaffected: it has no nk and always contains
            %   u(k).
            %
            %   Inputs:
            %       opts.nk / .inputDifference / .pastSource  as rid.rep.Narx
            %   Outputs:
            %       rep  (1x1) rid.rep.Narx
            arguments
                obj
                opts.nk              (1,1) double = rid.defaultOptions().rep.nk
                opts.inputDifference (1,1) logical ...
                    = rid.defaultOptions().rep.inputDifference
                opts.pastSource      (1,1) string ...
                    = rid.defaultOptions().rep.pastSource
            end
            if ~isfinite(obj.nStates)
                error("rid:surrogate:nStatesRequired", ...
                    "defaultRepresentation derives the lag orders from nStates, " + ...
                    "which this model does not declare. Either set nStates or " + ...
                    "override defaultRepresentation.");
            end
            if ~isfinite(obj.nInputs)
                error("rid:surrogate:nInputsRequired", ...
                    "defaultRepresentation needs nInputs, which this model does " + ...
                    "not declare.");
            end
            % The `+ 1` is u(k), so it exists only when the representation
            % contains it.
            nLagU = obj.nStates + (opts.nk == 0);
            rep = rid.rep.Narx(nLagU, obj.nStates, ...
                "nk", opts.nk, "nInputs", obj.nInputs, ...
                "inputDifference", opts.inputDifference, ...
                "pastSource", opts.pastSource);
        end

        function x = steadyState(obj, uConst, opts)
            %STEADYSTATE  Fallback: simulate a constant input to convergence.
            %
            %   The canonical producer of `xStart`. A model with an analytic
            %   rest state overrides this; a model that can only simulate gets it
            %   for free here.
            %
            %   Inputs:
            %       uConst  (1 x nInputs)  constant physical input
            %   Outputs:
            %       x       (1 x nStates)  state the model rests at under uConst
            arguments
                obj
                uConst (1,:) double {mustBeReal, mustBeFinite}
                opts.tol        (1,1) double = rid.defaultOptions().surrogate.steadyStateTol
                opts.maxSamples (1,1) double = rid.defaultOptions().surrogate.steadyStateMaxSamples
                opts.blockSize  (1,1) double = rid.defaultOptions().surrogate.steadyStateBlock
            end
            if ~isfinite(obj.nStates)
                error("rid:surrogate:nStatesRequired", ...
                    "The steadyState fallback starts from zeros(1, nStates), " + ...
                    "which this model does not declare. Either set nStates or " + ...
                    "override steadyState.");
            end
            x     = zeros(1, obj.nStates);
            uSeq  = repmat(uConst, opts.blockSize, 1);
            nDone = 0;
            while nDone < opts.maxSamples
                xPrev  = x;
                [~, x] = obj.simulate(uSeq, x);
                nDone  = nDone + opts.blockSize;
                if norm(x - xPrev) <= opts.tol * max(1, norm(x))
                    return
                end
            end
            error("rid:surrogate:steadyStateNotReached", ...
                "The constant-input simulation did not settle within %d " + ...
                "samples. The model may be unstable, marginally stable, or have " + ...
                "no steady state at this input.", opts.maxSamples);
        end

        function tf = isStable(~)
            %ISSTABLE  Fallback: true, UNVERIFIED.
            %
            %   A model that cannot check its own stability cannot be assumed
            %   unstable either, and refusing to run every model that lacks the
            %   check would exclude exactly the vendor-supplied models the open
            %   interface exists for. A caller cannot distinguish this default from
            %   a verified true; a model that CAN check must override.
            %
            %   Outputs:
            %       tf  (1x1) logical
            tf = true;
        end

        function tf = canUpdate(~)
            %CANUPDATE  Fallback: false -- this model cannot be refitted.
            %
            %   Explicit rather than implied, so that the scheduler does not
            %   schedule pointless work and the event log records WHY the surrogate
            %   never changes. Without the flag, "the model never updated" is
            %   indistinguishable from a bug.
            %
            %   Outputs:
            %       tf  (1x1) logical
            tf = false;
        end

        function obj = update(obj, ~, ~)
            %UPDATE  Fallback: a no-op returning the model unchanged.
            %
            %   `update` is a METHOD OF THE MODEL, never a generic fitter in the
            %   framework: refitting a linear ARX is least squares, refitting a
            %   physically derived model means deciding which of its parameters are
            %   free, and no framework can make that choice on a model it does not
            %   understand.
            %
            %   Outputs:
            %       obj  (1x1) the same model, unchanged
        end

    end

    methods (Access = protected)

        function assertBatchShape(obj, uSeqBatch)
            %ASSERTBATCHSHAPE  uSeqBatch is (nHorizon x nInputs x nSeq).
            %
            %   Dimension 2 is what carries the meaning; dimension 3 may be absent
            %   because MATLAB trims a trailing singleton, and that is legal.
            if isfinite(obj.nInputs) && size(uSeqBatch, 2) ~= obj.nInputs
                error("rid:surrogate:dimensionMismatch", ...
                    "uSeqBatch must be (nHorizon x nInputs x nSeq) with " + ...
                    "nInputs = %d (got %s). Dimension 2 is the input channel; " + ...
                    "dimension 3 counts sequences and may be absent for a single " + ...
                    "sequence, since MATLAB trims a trailing singleton.", ...
                    obj.nInputs, mat2str(size(uSeqBatch)));
            end
        end

        function assertInputShape(obj, uSeq)
            %ASSERTINPUTSHAPE  uSeq is (nHorizon x nInputs), samples in rows.
            if size(uSeq, 2) ~= obj.nInputs
                error("rid:surrogate:dimensionMismatch", ...
                    "uSeq must have nInputs = %d columns (got %d). " + ...
                    "Samples are in rows.", obj.nInputs, size(uSeq, 2));
            end
        end

        function assertStateShape(obj, xStart)
            %ASSERTSTATESHAPE  xStart is (1 x nStates).
            if ~isequal(size(xStart), [1 obj.nStates])
                error("rid:surrogate:dimensionMismatch", ...
                    "xStart must be (1 x nStates) = (1 x %d) (got %s). " + ...
                    "It is required, never defaulted; the canonical " + ...
                    "producer is steadyState(uConst).", ...
                    obj.nStates, mat2str(size(xStart)));
            end
        end

    end

end
