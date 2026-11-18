classdef LinearFirstOrder < rid.surrogate.Surrogate
%LINEARFIRSTORDER  ZOH-exact first-order lag; the entry-point surrogate.
%
%   The least prior knowledge that suffices for the loop to start: one time
%   constant. In PHYSICAL units throughout -- normalization happens once, at the
%   representation boundary, so a time constant and a gain mean what they say.
%
%       x(k+1) = a*x(k) + u(k,:)*b.' + e
%       y(k)   = x(k)   + u(k,:)*directTerm.'
%
%       a = exp(-tS/timeConstant)          ZOH-exact for a continuous first-order
%       b = (1-a) * (gain - directTerm)    so that y_ss = gain*u + offset
%       e = (1-a) * offset
%
%   The discrete recursion reproduces the continuous solution
%   y(t) = gain*u*(1 - exp(-t/T)) EXACTLY at the sample instants; that is what
%   "ZOH-exact" buys and it is asserted analytically in the tests.
%
%   THE GAIN IS DERIVED FROM THE ROI, not defaulted to one (the derived gain, superseding the
%   gain clause). Since the normalization derives from the declared ROI in
%   physical units, the gain decides WHERE IN THE ROI the simulated points
%   land: a gain-1 model against an ROI whose output range is ten times its input
%   range would place every predicted point in the bottom tenth of the output
%   dimension -- fatal offline, and a systematic distortion of candidate ranking
%   online, where the measured past coordinates stay correct and the mismatch is
%   therefore harder to notice. `fromRoi` derives
%
%       gain   = (yMax - yMin) / (nInputs * (uMax - uMin))
%       offset = yMid - gain * uMid.'
%
%   so the static map spans the ROI diagonal and the equilibrium passes through the
%   ROI centre by construction. The user still supplies only the time constant. The
%   1/nInputs factor generalizes "spans the diagonal" to several inputs; at
%   nInputs = 1 it is exactly the formula.
%
%   This does NOT make the gain meaningful as a process estimate. It is a scaling
%   choice that makes the surrogate's geometry consistent with the declared region,
%   nothing more.
%
%   STRICTLY PROPER BY DEFAULT: `directTerm = 0`, so the step response does
%   not jump. A non-zero direct term is available for a process with an algebraic
%   path; the static gain then splits as gain = b/(1-a) + directTerm, and that split
%   is a statement about the process, not something the framework can default.
%   Whether the REPRESENTATION contains u(k) is a separate and independently
%   adjustable choice.
%
%   Construction:
%       srg = rid.surrogate.LinearFirstOrder(timeConstant, tS, gain)
%       srg = rid.surrogate.LinearFirstOrder(timeConstant, tS, gain, ...
%                                            "offset", 0, "directTerm", 0)
%       srg = rid.surrogate.LinearFirstOrder.fromRoi(roi, timeConstant, tS)
%
%   The ROI format deliberately does not enter the main constructor: a model built
%   from process knowledge carries its own gain, and only the ROI-derived default
%   needs the declaration. `fromRoi` reads the LEVEL coordinate of each input
%   channel and the first output coordinate by name, because under the declared-bounds rule different
%   lags of the same channel may carry different ranges.
%
%   Inputs:
%       timeConstant       (1x1)  dominant time constant, physical, > 0
%       tS                 (1x1)  sample time, physical, > 0
%       gain               (1 x nInputs)  static gain per input channel
%       opts.offset        (1x1)  static-map offset, so y_ss = gain*u + offset
%       opts.directTerm    (1 x nInputs)  feedthrough; 0 = strictly proper
%
%   BOTH FORMS ARE REFITTABLE. `canUpdate` is true whether or not the model
%   declares a direct term: both are linear in their parameters, and the
%   feedthrough form is identifiable given excitation that separates u(k) from
%   u(k-1). Where the excitation does not separate them, the rank guard in `update`
%   rejects the fit -- the unidentifiable case is detected from the data, not
%   assumed from the structure.

    properties (SetAccess = immutable)
        gain   (1,:) double   % (1 x nInputs) static gain, PHYSICAL
        offset (1,1) double   % static-map offset: y_ss = gain*u + offset
        a      (1,1) double   % exp(-tS/timeConstant), the discrete pole
        b      (1,:) double   % (1 x nInputs) input matrix
        e      (1,1) double   % affine state term
    end

    methods

        function obj = LinearFirstOrder(timeConstant, tS, gain, opts)
            arguments
                timeConstant (1,1) double {mustBeReal, mustBePositive, mustBeFinite}
                tS           (1,1) double {mustBeReal, mustBePositive, mustBeFinite}
                gain         (1,:) double {mustBeReal, mustBeFinite}
                opts.offset     (1,1) double {mustBeReal, mustBeFinite} = 0
                opts.directTerm (1,:) double {mustBeReal, mustBeFinite} = 0
            end

            nInputs = numel(gain);
            d = opts.directTerm;
            if isscalar(d) && nInputs > 1
                d = repmat(d, 1, nInputs);
            end
            if numel(d) ~= nInputs
                error("rid:surrogate:dimensionMismatch", ...
                    "directTerm must have one entry per input channel " + ...
                    "(got %d for nInputs = %d).", numel(d), nInputs);
            end

            obj.timeConstant = timeConstant;
            obj.tS           = tS;
            obj.nInputs      = nInputs;
            obj.nStates      = 1;
            obj.directTerm   = d;

            obj.gain   = gain;
            obj.offset = opts.offset;
            obj.a      = exp(-tS / timeConstant);
            obj.b      = (1 - obj.a) * (gain - d);
            obj.e      = (1 - obj.a) * opts.offset;
        end

        function [yHat, xEnd] = simulate(obj, uSeq, xStart)
            %SIMULATE  Pure, deterministic forward simulation in physical units.
            %
            %   No noise, no persistent variables, no RNG: state in, state out.
            %   `filter` runs the linear recursion in compiled code, which matters
            %   because candidate scoring calls this thousands of times per step.
            %
            %   Inputs:
            %       uSeq    (nHorizon x nInputs)  physical input
            %       xStart  (1 x 1)               initial state; REQUIRED
            %   Outputs:
            %       yHat    (nHorizon x 1)  physical output
            %       xEnd    (1 x 1)         state after the last sample
            arguments
                obj
                uSeq   (:,:) double {mustBeReal}
                xStart (1,:) double {mustBeReal}
            end
            obj.assertInputShape(uSeq);
            obj.assertStateShape(xStart);
            [yHat, xEnd] = obj.recursion(uSeq, xStart);
        end

        function [yHat, xEnd, xTraj] = simulateStates(obj, uSeq, xStart)
            %SIMULATESTATES  The same simulation, with the trajectory it computes.
            %
            %   ONE PRIVATE CORE WITH `simulate`. The state
            %   trajectory is not extra work here -- the recursion produces it and
            %   `simulate` then discards all but the last value -- so the only thing
            %   an override must not do is compute it a SECOND way. Two paths that
            %   must agree sample for sample would not error when they diverged;
            %   they would return slightly wrong points and a slightly wrong design
            %   (the argument, one level down).
            %
            %   Row `i` is the state BEFORE `uSeq(i,:)` acts, so row 1 is `xStart`
            %   and `xEnd` is not a row.
            %
            %   Inputs:
            %       uSeq    (nHorizon x nInputs)  physical input
            %       xStart  (1 x 1)               initial state; REQUIRED
            %   Outputs:
            %       yHat    (nHorizon x 1)
            %       xEnd    (1 x 1)
            %       xTraj   (nHorizon x 1)        this model has one state
            arguments
                obj
                uSeq   (:,:) double {mustBeReal}
                xStart (1,:) double {mustBeReal}
            end
            obj.assertInputShape(uSeq);
            obj.assertStateShape(xStart);
            [yHat, xEnd, xTraj] = obj.recursion(uSeq, xStart);
        end

        function [Y, Xend, XTraj] = simulateBatch(obj, uSeqBatch, xStart)
            %SIMULATEBATCH  Vectorized override of the base loop fallback.
            %
            %   One columnwise `filter` over all sequences instead of nSeq calls.
            %   It must agree with the base-class loop exactly, and the tests assert
            %   that against rid.surrogate.Surrogate's fallback -- all three
            %   outputs, since XTraj is part of the same contract.
            %
            %   `X` IS ALREADY THE BATCH OF TRAJECTORIES, in the row convention of
            %   column s holds the states before each sample of sequence s.
            %   The third output is that array reshaped to put the single state on
            %   dimension 2 and the sequence on dimension 3 -- no second
            %   computation, so nothing can disagree with the first two outputs.
            %
            %   Inputs:
            %       uSeqBatch  (nHorizon x nInputs x nSeq)  dimension 3 counts
            %                                              sequences; it is absent
            %                                              for a single one
            %       xStart     (1 x 1)
            %   Outputs:
            %       Y      (nHorizon x nSeq)
            %       Xend   (nSeq x 1)
            %       XTraj  (nHorizon x 1 x nSeq)   nStates = 1 for this model
            arguments
                obj
                uSeqBatch (:,:,:) double {mustBeReal}
                xStart    (1,:) double {mustBeReal}
            end
            obj.assertBatchShape(uSeqBatch);
            obj.assertStateShape(xStart);

            [nHorizon, ~, nSeq] = size(uSeqBatch);
            if nHorizon == 0 || nSeq == 0
                Y     = zeros(nHorizon, nSeq);
                Xend  = repmat(xStart, nSeq, 1);
                XTraj = zeros(nHorizon, 1, nSeq);
                return
            end

            V     = reshape(pagemtimes(uSeqBatch, obj.b.'), nHorizon, nSeq) + obj.e;
            X     = filter(1, [1, -obj.a], ...
                        [repmat(xStart, 1, nSeq); V(1:end-1, :)], [], 1);
            Y     = X + reshape(pagemtimes(uSeqBatch, obj.directTerm.'), ...
                        nHorizon, nSeq);
            Xend  = (obj.a * X(end, :) + V(end, :)).';
            XTraj = reshape(X, nHorizon, 1, nSeq);
        end

        function x = steadyState(obj, uConst)
            %STEADYSTATE  Analytic rest state; the canonical producer of xStart.
            %
            %   x_ss = (u*b.' + e) / (1 - a), so that y_ss = gain*u + offset. An
            %   experiment starting from rest starts here.
            %
            %   Inputs:
            %       uConst  (1 x nInputs)  constant physical input
            %   Outputs:
            %       x       (1 x 1)        rest state
            arguments
                obj
                uConst (1,:) double {mustBeReal, mustBeFinite}
            end
            obj.assertInputShape(uConst);
            x = (uConst * obj.b.' + obj.e) / (1 - obj.a);
        end

        function tf = isStable(obj)
            %ISSTABLE  A genuine check, not the base class's unverified default.
            %
            %   0 < a < 1 is what a real, positive, finite time constant means. The
            %   constructor enforces timeConstant > 0, so a CONSTRUCTED object is
            %   always stable and this only ever bites on an update -- which is
            %   exactly where it must.
            %
            %   Outputs:
            %       tf  (1x1) logical
            tf = obj.a > 0 && obj.a < 1;
        end

        function tf = canUpdate(~)
            %CANUPDATE  This model can always be refitted.
            %
            %   Both the strictly proper and the feedthrough forms are linear in
            %   their parameters, so both are ordinary least squares. A feedthrough
            %   model is identifiable given excitation that separates u(k) from
            %   u(k-1); where the excitation does not, the rank guard in `update`
            %   detects the collinearity and rejects the fit. Declaring
            %   canUpdate = false for every feedthrough model would refuse the
            %   identifiable case along with the unidentifiable one.
            %
            %   Outputs:
            %       tf  (1x1) logical
            tf = true;
        end

        function obj = update(obj, hist, opts)
            %UPDATE  Refit (timeConstant, gain, offset) by least squares. Value semantics.
            %
            %   The fit is the equation-error (one-step) form. Strictly proper:
            %
            %       y(k) = a*y(k-1) + u(k-1,:)*b.' + e
            %
            %   With a direct term, substituting x(k) = y(k) - u(k,:)*d.' gives the
            %   same kind of regression with one block more:
            %
            %       y(k) = a*y(k-1) + u(k-1,:)*(b - a*d).' + u(k,:)*d.' + e
            %
            %   The FORM FOLLOWS THE MODEL'S DECLARED STRUCTURE: a strictly proper
            %   model is refitted strictly proper, so an update can never introduce
            %   a feedthrough path the user did not declare. Both forms are linear
            %   in their parameters and both are plain least squares; the second is
            %   identifiable exactly when the excitation separates u(k) from u(k-1),
            %   and where it does not, the rank guard below rejects the fit rather
            %   than the model refusing in advance.
            %
            %   Either form is fitted over ALL data acquired up to the snapshot,
            %   with NO forgetting factor:
            %   the process is assumed time-invariant and the surrogate is
            %   a geometric guide, so a forgetting factor would introduce drift
            %   indistinguishable from exploration. It fits INPUT-OUTPUT data only;
            %   this model's state IS its output, so no estimated state is involved
            %   -- an observer depends on the current model, and fitting to its
            %   output would make the update self-referential.
            %
            %   The one-step fit is used for multi-step horizon simulation; that
            %   tension is deliberate and open.
            %
            %   AN UPDATE THAT LEAVES THE MODEL UNUSABLE IS REJECTED, NEVER
            %   SILENTLY ACCEPTED, and the returned object is then the original:
            %
            %       abs(a) >= 1   rid:surrogate:unstableUpdate    diverges over the
            %                                                     horizon; J becomes
            %                                                     meaningless
            %       a <= 0        rid:surrogate:invalidUpdate     STABLE but has no
            %                                                     real first-order
            %                                                     time constant --
            %                                                     outside this model
            %                                                     class
            %       too little /  rid:surrogate:insufficientData
            %       rank-deficient
            %
            %   Inputs:
            %       hist              (1x1) struct  .u (nSteps x nInputs), .y (nSteps x 1)
            %       opts.minSamples   (1x1)  regression rows required
            %   Outputs:
            %       obj  (1x1) a NEW model; the original is untouched
            arguments
                obj
                hist (1,1) struct
                opts.minSamples (1,1) double {mustBeInteger, mustBePositive} ...
                    = rid.defaultOptions().surrogate.minUpdateSamples
            end

            rid.assertHistory(hist);
            % SEMANTIC: a well-formed offline history is fine, there is simply
            % nothing to refit on.
            if ~isfield(hist, "y") || isempty(hist.y)
                error("rid:surrogate:insufficientData", ...
                    "update needs measured hist.y. Offline there is nothing to " + ...
                    "refit on, which is why the offline row of the rebuild table " + ...
                    "says the surrogate never changes.");
            end
            u = double(hist.u);
            y = double(hist.y);
            obj.assertInputShape(u);
            if ~isequal(size(y), [size(u, 1), 1])
                error("rid:surrogate:dimensionMismatch", ...
                    "hist.y must be (nSteps x 1) = (%d x 1) (got %s).", ...
                    size(u, 1), mat2str(size(y)));
            end

            nRows = size(u, 1) - 1;
            if nRows < opts.minSamples
                error("rid:surrogate:insufficientData", ...
                    "The fit needs at least %d regression rows; the history " + ...
                    "gives %d. A first-order model fitted from a handful of " + ...
                    "samples is not worth swapping in.", opts.minSamples, ...
                    max(nRows, 0));
            end

            % Equation-error regression. The u(k) block is present only for a
            % model that declares a direct term, so the refit preserves the
            % model's structure.
            hasDirect = any(obj.directTerm ~= 0);
            if hasDirect
                Phi = [y(1:end-1), u(1:end-1, :), u(2:end, :), ones(nRows, 1)];
            else
                Phi = [y(1:end-1), u(1:end-1, :), ones(nRows, 1)];
            end
            target = y(2:end);
            if rank(Phi) < size(Phi, 2)
                error("rid:surrogate:insufficientData", ...
                    "The regression matrix is rank deficient (rank %d of %d). " + ...
                    "The history does not excite the model's parameters -- a " + ...
                    "constant input is the usual cause, and for a model with a " + ...
                    "direct term so is any input held constant across a step, " + ...
                    "which makes u(k) and u(k-1) collinear.", ...
                    rank(Phi), size(Phi, 2));
            end
            theta = Phi \ target;

            nIn  = obj.nInputs;
            aNew = theta(1);
            eNew = theta(end);
            if hasDirect
                dNew = theta(2 + nIn : 1 + 2 * nIn).';
                bNew = theta(2 : 1 + nIn).' + aNew * dNew;   % undo b - a*d
            else
                dNew = obj.directTerm;                        % stays zero
                bNew = theta(2 : 1 + nIn).';
            end

            if abs(aNew) >= 1
                error("rid:surrogate:unstableUpdate", ...
                    "The fit gives a discrete pole a = %.6g with |a| >= 1. An " + ...
                    "unstable surrogate makes the horizon simulation diverge and " + ...
                    "the criterion meaningless, so the update is rejected and the " + ...
                    "previous model kept.", aNew);
            end
            if aNew <= 0
                error("rid:surrogate:invalidUpdate", ...
                    "The fit gives a discrete pole a = %.6g. This is stable but " + ...
                    "has no real positive time constant (T = -tS/log(a) is not " + ...
                    "real), so it lies outside this model class. Rejected; the " + ...
                    "previous model is kept.", aNew);
            end

            % Back to physical parameters: the static gain splits as
            % K = b/(1-a) + d.
            obj = rid.surrogate.LinearFirstOrder( ...
                -obj.tS / log(aNew), obj.tS, bNew / (1 - aNew) + dNew, ...
                "offset", eNew / (1 - aNew), "directTerm", dNew);
        end

        function s = describe(obj)
            %DESCRIBE  Human-readable summary, for logs and run records.
            %
            %   Outputs:
            %       s  (1x1) string
            s = sprintf(...
                "rid.surrogate.LinearFirstOrder  T = %g, tS = %g, a = %.6g\n" + ...
                "  gain       : %s\n  offset     : %g\n  directTerm : %s\n" + ...
                "  nStates = %d, nInputs = %d, isStable = %d, canUpdate = %d", ...
                obj.timeConstant, obj.tS, obj.a, mat2str(obj.gain, 6), ...
                obj.offset, mat2str(obj.directTerm, 6), ...
                obj.nStates, obj.nInputs, obj.isStable(), obj.canUpdate());
        end

    end

    methods (Static)

        function obj = fromRoi(roi, timeConstant, tS, opts)
            %FROMROI  Entry-point model with the gain derived from the ROI.
            %
            %   THE GAIN COMES FROM THE LEVEL COORDINATE'S DECLARED RANGE. Under
            %   the declared-bounds rule every dimension carries its own bounds, so different lags of
            %   the same channel may legitimately carry different ranges and "the
            %   input range" is no longer unambiguous. This function therefore reads
            %   named coordinates rather than signals: the level coordinate of each
            %   input channel and the first output coordinate. The defaults are the
            %   default representation's (`u1_k`, `y_km1`); a representation with
            %   another input delay passes its own names.
            %
            %       gain   = (yMax - yMin) / (nInputs * (uMax - uMin))
            %       offset = yMid - gain * uMid.'
            %
            %   so the static map spans the ROI diagonal and the equilibrium passes
            %   through the ROI centre by construction. The 1/nInputs factor
            %   generalizes "spans the diagonal" to several inputs; at nInputs = 1 it
            %   is exactly the formula.
            %
            %   Names it is not given are ignored: this is a MODEL constructor, not a
            %   representation, and it consumes only the ranges its static map spans.
            %
            %   Inputs:
            %       roi              (1x1) struct  the declaration, keyed by
            %                                      coordinate name
            %       timeConstant     (1x1)  the only modelling input the user supplies
            %       tS               (1x1)  sample time
            %       opts.inputName   (1 x nInputs) string, level coordinate per channel
            %       opts.outputName  (1x1) string, the output coordinate
            %       opts.directTerm  (1 x nInputs)  feedthrough; 0 = strictly proper
            %   Outputs:
            %       obj  (1x1) rid.surrogate.LinearFirstOrder
            arguments
                roi          (1,1) struct
                timeConstant (1,1) double {mustBeReal, mustBePositive, mustBeFinite}
                tS           (1,1) double {mustBeReal, mustBePositive, mustBeFinite}
                opts.inputName  (1,:) string = rid.rep.Narx.coordinateName( ...
                    "u1", rid.defaultOptions().rep.nk)
                opts.outputName (1,1) string = rid.rep.Narx.coordinateName("y", 1)
                opts.directTerm (1,:) double {mustBeReal, mustBeFinite} = 0
            end

            nInputs = numel(opts.inputName);
            rY      = rid.surrogate.LinearFirstOrder.rangeOf(roi, opts.outputName);
            gain    = zeros(1, nInputs);
            uMid    = zeros(1, nInputs);
            for cIn = 1:nInputs
                rU = rid.surrogate.LinearFirstOrder.rangeOf(roi, opts.inputName(cIn));
                gain(cIn) = (rY(2) - rY(1)) / (nInputs * (rU(2) - rU(1)));
                uMid(cIn) = mean(rU);
            end
            offset = mean(rY) - gain * uMid.';

            obj = rid.surrogate.LinearFirstOrder(timeConstant, tS, gain, ...
                "offset", offset, "directTerm", opts.directTerm);
        end

    end

    methods (Access = private)

        function [yHat, xEnd, xTraj] = recursion(obj, uSeq, xStart)
            %RECURSION  The ZOH-exact recursion; the ONE core of both simulations.
            %
            %   `filter` runs it in compiled code, which matters because candidate
            %   scoring calls this thousands of times per step. Shapes are the
            %   callers' contract; nothing is validated here.
            nHorizon = size(uSeq, 1);
            if nHorizon == 0
                % Legal, and the loop asks for it once per step when the
                % representation looks back zero samples.
                yHat  = zeros(0, 1);
                xEnd  = xStart;
                xTraj = zeros(0, 1);
                return
            end

            v     = uSeq * obj.b.' + obj.e;             % (nHorizon x 1) state drive
            xTraj = filter(1, [1, -obj.a], [xStart; v(1:end-1)]);
            yHat  = xTraj + uSeq * obj.directTerm.';
            xEnd  = obj.a * xTraj(end) + v(end);
        end

    end

    methods (Static, Access = private)
        function r = rangeOf(roi, name)
            %RANGEOF  One declared range, by coordinate name.
            %
            %   The identifiers name THIS package, not +rep: an error identifier
            %   names the package that throws it. This is the surrogate
            %   reading a declaration, not a representation validating one.
            %
            %   The value check is kept here even though rid.roi.Normalization owns
            %   box validity in general, because this function DIVIDES by the
            %   span: an unordered or zero-width range would silently produce an
            %   infinite or negative gain, which is arithmetic of this class's own,
            %   not a statement about the box.
            if ~isfield(roi, name)
                error("rid:surrogate:missingRange", ...
                    "The declaration has no range for ""%s"", which is the " + ...
                    "coordinate this model's static map spans. Pass the name " + ...
                    "explicitly via ""inputName"" / ""outputName"" if the " + ...
                    "representation uses another input delay.", name);
            end
            r = roi.(name);
            if ~isnumeric(r) || ~isreal(r) || ~isequal(size(r), [1 2]) ...
                    || ~all(isfinite(r)) || r(2) <= r(1)
                error("rid:surrogate:invalidRange", ...
                    "The range declared for ""%s"" must be a finite (1 x 2) " + ...
                    "row [lower upper] with upper > lower; the ROI-derived gain " + ...
                    "divides by its width.", name);
            end
            r = double(r);
        end
    end

end
