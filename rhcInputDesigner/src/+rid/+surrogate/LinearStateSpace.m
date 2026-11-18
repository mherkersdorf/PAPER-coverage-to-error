classdef LinearStateSpace < rid.surrogate.Surrogate
%LINEARSTATESPACE  Order-n linear surrogate in physical units.
%
%       x(k+1)  = Ad x(k) + Bd u(k).' + ed
%       yHat(k) = Cd x(k) + directTerm u(k).' + fd
%
%   The surrogate that pairs with the joint input--state-space representation
%   `rid.rep.StateSpace`. What makes that pairing legitimate is a premise the code
%   cannot check: THE STATES ARE THE PROCESS'S OWN. Only then is a measured state
%   an initial condition for this model, and only then is a simulated state
%   comparable with a measured one. `assertCompatible` checks the
%   dimension; nothing can check the basis. That is why a state-space
%   representation needs an INTERPRETABLE model and why this class does
%   not offer one as its `defaultRepresentation` (part 2): a class cannot know
%   whether its own states are interpretable.
%
%   PARAMETER NAMES ARE THIS MODEL'S OWN. `Ad`, `Bd`, `ed`,
%   `Cd`, `fd` are local shorthand and do not refer to the symbols of the notation
%   tables -- `d` there is the dimension of the regressor space. The `d` suffix
%   follows the method description and keeps these clear of the affine
%   Jacobian. The direct term is the BASE CLASS's `directTerm`, not a second name
%   for it.
%
%   `ed` AND `fd` EXIST BECAUSE A FIT PRODUCES THEM. A least-squares fit with an
%   intercept around an operating point yields an affine state term and an output
%   offset; they correspond to `LinearFirstOrder`'s `e` and `offset`.
%
%   STABILITY IS REFUSED, NOT REPORTED. The constructor rejects
%   max(abs(eig(Ad))) >= 1 with `rid:surrogate:unstableModel`: an unstable
%   surrogate makes the horizon simulation diverge and the criterion meaningless.
%   A pole exactly on the unit circle is
%   refused too -- marginal is not stable over a long run.
%
%   NOT REFITTABLE IN THE LOOP: `canUpdate` is false and `update` is the inherited
%   no-op. An input--output refit would return an arbitrary state basis, which
%   breaks the premise above. The fit that IS available is the static
%   `fitFromStates`, used OUTSIDE the loop on states that were genuinely measured.
%
%   `simulateBatch` IS NOT OVERRIDDEN. No caller needs a vectorized batch yet, and
%   an override would have to match the base loop bit for bit. The
%   inherited fallback is correct and is asserted to be.
%
%   Construction:
%       srg = rid.surrogate.LinearStateSpace(Ad, Bd, Cd, tS)
%       srg = rid.surrogate.LinearStateSpace(Ad, Bd, Cd, tS, "ed", e, "fd", f, ...
%                                            "directTerm", d)
%       srg = rid.surrogate.LinearStateSpace.fromContinuous(Ac, Bc, Cc, tS)
%       srg = rid.surrogate.LinearStateSpace.fitFromStates(u, x, y, tS)
%
%   Inputs:
%       Ad               (nStates x nStates)  state transition, DISCRETE
%       Bd               (nStates x nInputs)  input matrix, DISCRETE
%       Cd               (1 x nStates)        output map; q = 1 project-wide
%       tS               (1x1)                sample time, physical, > 0
%       opts.ed          (nStates x 1)        affine state term; [] means zeros
%       opts.fd          (1x1)                output offset
%       opts.directTerm  (1 x nInputs)        feedthrough; scalar expands, 0 = strictly proper

    properties (SetAccess = immutable)
        Ad (:,:) double   % (nStates x nStates)
        Bd (:,:) double   % (nStates x nInputs)
        ed (:,1) double   % (nStates x 1)  affine state term
        Cd (1,:) double   % (1 x nStates)
        fd (1,1) double   % output offset
    end

    properties (SetAccess = immutable, GetAccess = private)
        % ONE source for the three things that read it: the constructor's refusal,
        % `isStable` and `timeConstant`. Computed once because Ad is immutable, and
        % private because it is arithmetic of this class rather than part of the
        % surrogate contract -- callers ask `isStable`.
        spectralRadius (1,1) double
    end

    methods

        function obj = LinearStateSpace(Ad, Bd, Cd, tS, opts)
            arguments
                Ad (:,:) double {mustBeReal, mustBeFinite}
                Bd (:,:) double {mustBeReal, mustBeFinite}
                Cd (:,:) double {mustBeReal, mustBeFinite}
                tS (1,1) double {mustBeReal, mustBePositive, mustBeFinite}
                opts.ed         (:,:) double {mustBeReal, mustBeFinite} = []
                opts.fd         (1,1) double {mustBeReal, mustBeFinite} = 0
                opts.directTerm (1,:) double {mustBeReal, mustBeFinite} = 0
            end

            nStates = size(Ad, 1);
            if size(Ad, 2) ~= nStates
                error("rid:surrogate:dimensionMismatch", ...
                    "Ad must be square (nStates x nStates); got %s.", ...
                    mat2str(size(Ad)));
            end
            nInputs = size(Bd, 2);
            if size(Bd, 1) ~= nStates
                error("rid:surrogate:dimensionMismatch", ...
                    "Bd must be (nStates x nInputs) = (%d x nInputs); got %s.", ...
                    nStates, mat2str(size(Bd)));
            end
            if ~isequal(size(Cd), [1 nStates])
                error("rid:surrogate:dimensionMismatch", ...
                    "Cd must be (1 x nStates) = (1 x %d); got %s. The output is " + ...
                    "(nSteps x 1) project-wide, so q = 1.", ...
                    nStates, mat2str(size(Cd)));
            end

            e = opts.ed;
            if isempty(e)
                e = zeros(nStates, 1);
            end
            if ~isequal(size(e), [nStates 1])
                error("rid:surrogate:dimensionMismatch", ...
                    "ed must be (nStates x 1) = (%d x 1); got %s.", ...
                    nStates, mat2str(size(e)));
            end

            % The LinearFirstOrder convention, kept rather than re-invented: one
            % scalar means the same feedthrough on every channel.
            d = opts.directTerm;
            if isscalar(d) && nInputs > 1
                d = repmat(d, 1, nInputs);
            end
            if numel(d) ~= nInputs
                error("rid:surrogate:dimensionMismatch", ...
                    "directTerm must have one entry per input channel " + ...
                    "(got %d for nInputs = %d).", numel(d), nInputs);
            end

            radius = max(abs(eig(Ad)));
            if isempty(radius)
                radius = 0;
            end
            % STRICTLY less than one. A pole ON the unit circle does not decay, so
            % the horizon simulation never settles and J measures a trajectory the
            % process would not follow.
            if radius >= 1
                error("rid:surrogate:unstableModel", ...
                    "Ad has spectral radius %.6g >= 1, so this model is not " + ...
                    "asymptotically stable. An unstable surrogate makes the " + ...
                    "horizon simulation diverge and the criterion meaningless " + ...
                    ".", radius);
            end

            obj.Ad = Ad;
            obj.Bd = Bd;
            obj.ed = e;
            obj.Cd = Cd;
            obj.fd = opts.fd;
            obj.spectralRadius = radius;

            obj.nStates    = nStates;
            obj.nInputs    = nInputs;
            obj.tS         = tS;
            obj.directTerm = d;
            % The decay time of the SLOWEST mode: |lambda|^k = exp(k*tS/T) gives
            % T = -tS/log(max|lambda|). Guidance for choosing nHorizon, never used
            % to compute it. A deadbeat model settles in n steps, and the
            % limit as the radius goes to zero is zero.
            if radius == 0
                obj.timeConstant = 0;
            else
                obj.timeConstant = -tS / log(radius);
            end
        end

        function [yHat, xEnd] = simulate(obj, uSeq, xStart)
            %SIMULATE  Pure, deterministic forward simulation in physical units.
            %
            %   Inputs:
            %       uSeq    (nHorizon x nInputs)  physical input
            %       xStart  (1 x nStates)         initial state; REQUIRED
            %   Outputs:
            %       yHat    (nHorizon x 1)
            %       xEnd    (1 x nStates)         state after the last sample
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
            %SIMULATESTATES  The same simulation, with the state trajectory.
            %
            %   ONE PRIVATE CORE WITH `simulate`. Two implementations that must
            %   agree sample for sample would not error when they diverged; they
            %   would return slightly wrong points and a slightly wrong design.
            %
            %   Row `i` is the state BEFORE `uSeq(i,:)` acts, so row 1 is `xStart`
            %   and `xEnd` is not a row.
            %
            %   Inputs:
            %       uSeq    (nHorizon x nInputs)  physical input
            %       xStart  (1 x nStates)         initial state; REQUIRED
            %   Outputs:
            %       yHat    (nHorizon x 1)
            %       xEnd    (1 x nStates)
            %       xTraj   (nHorizon x nStates)
            arguments
                obj
                uSeq   (:,:) double {mustBeReal}
                xStart (1,:) double {mustBeReal}
            end
            obj.assertInputShape(uSeq);
            obj.assertStateShape(xStart);
            [yHat, xEnd, xTraj] = obj.recursion(uSeq, xStart);
        end

        function x = steadyState(obj, uConst)
            %STEADYSTATE  Analytic rest state; the canonical producer of xStart.
            %
            %   x_ss = (I - Ad) \ (Bd u.' + ed), the fixed point of the recursion.
            %   `ed` enters on the same side as the input drive, which is what it
            %   is: a constant contribution to the state update.
            %
            %   (I - Ad) is invertible because the constructor refused a unit
            %   eigenvalue, so this cannot be reached with a singular solve.
            %
            %   Inputs:
            %       uConst  (1 x nInputs)  constant physical input
            %   Outputs:
            %       x       (1 x nStates)  state the model rests at under uConst
            arguments
                obj
                uConst (1,:) double {mustBeReal, mustBeFinite}
            end
            obj.assertInputShape(uConst);
            x = ((eye(obj.nStates) - obj.Ad) \ (obj.Bd * uConst.' + obj.ed)).';
        end

        function tf = isStable(obj)
            %ISSTABLE  A genuine check, not the base class's unverified default.
            %
            %   Always true for a CONSTRUCTED object -- the constructor refuses
            %   anything else -- which is deliberate: the refusal is where the
            %   information is, and this member exists so a caller can ask without
            %   knowing that.
            %
            %   Outputs:
            %       tf  (1x1) logical
            tf = obj.spectralRadius < 1;
        end

        function s = describe(obj)
            %DESCRIBE  Human-readable summary, for logs and run records.
            %
            %   Outputs:
            %       s  (1x1) string
            s = sprintf( ...
                "rid.surrogate.LinearStateSpace  nStates = %d, nInputs = %d, " + ...
                "tS = %g\n  Ad         : %s\n  Bd         : %s\n" + ...
                "  ed         : %s\n  Cd         : %s\n  fd         : %g\n" + ...
                "  directTerm : %s\n" + ...
                "  spectral radius = %.6g, timeConstant = %g, " + ...
                "isStable = %d, canUpdate = %d", ...
                obj.nStates, obj.nInputs, obj.tS, ...
                mat2str(obj.Ad, 6), mat2str(obj.Bd, 6), mat2str(obj.ed, 6), ...
                mat2str(obj.Cd, 6), obj.fd, mat2str(obj.directTerm, 6), ...
                obj.spectralRadius, obj.timeConstant, ...
                obj.isStable(), obj.canUpdate());
        end

    end

    methods (Access = private)

        function [yHat, xEnd, xTraj] = recursion(obj, uSeq, xStart)
            %RECURSION  The ONE core of both simulations. Shapes are the callers'.
            %
            %   A plain loop: n is small and the horizon short, and the vectorized
            %   alternative would be a second path to keep in step for no caller
            %   that has asked.
            nHorizon = size(uSeq, 1);
            yHat     = zeros(nHorizon, 1);
            xTraj    = zeros(nHorizon, obj.nStates);
            x        = xStart.';
            for kStep = 1:nHorizon
                % RECORDED BEFORE THE STEP: row i is the state before uSeq(i,:)
                % acts, which pairs each state with the input at the same instant.
                xTraj(kStep, :) = x.';
                u               = uSeq(kStep, :).';
                yHat(kStep)     = obj.Cd * x + obj.directTerm * u + obj.fd;
                x               = obj.Ad * x + obj.Bd * u + obj.ed;
            end
            xEnd = x.';
        end

    end

    methods (Static)

        function obj = fromContinuous(Ac, Bc, Cc, tS, opts)
            %FROMCONTINUOUS  Exact zero-order-hold discretization, via `expm`.
            %
            %   For xdot = Ac x + Bc u + ec with u held constant over a sample,
            %
            %       x(t+T) = e^{Ac T} x(t) + int_0^T e^{Ac s} ds (Bc u + ec)
            %
            %   and exponentiating the AUGMENTED matrix produces all three blocks
            %   at once:
            %
            %       M = expm([Ac, Bc, ec; 0, 0, 0] * tS)
            %       Ad = M(1:n, 1:n)   Bd = M(1:n, n+1:n+p)   ed = M(1:n, end)
            %
            %   EXACT AT THE SAMPLE INSTANTS, not an approximation -- which is what
            %   makes the state coordinates of rid.rep.StateSpace mean what they
            %   say, and is asserted analytically against a closed-form step
            %   response in the tests.
            %
            %   `expm` IS BASE MATLAB. `c2d` would need the Control System Toolbox,
            %   which this project does not have.
            %
            %   Inputs:
            %       Ac               (nStates x nStates)  continuous state matrix
            %       Bc               (nStates x nInputs)  continuous input matrix
            %       Cc               (1 x nStates)        output map; carried over
            %                                             unchanged, since a ZOH
            %                                             discretization does not
            %                                             touch the output equation
            %       tS               (1x1)                sample time, > 0
            %       opts.ec          (nStates x 1)        continuous affine term
            %       opts.fd          (1x1)                output offset
            %       opts.directTerm  (1 x nInputs)        feedthrough
            %   Outputs:
            %       obj  (1x1) rid.surrogate.LinearStateSpace
            arguments
                Ac (:,:) double {mustBeReal, mustBeFinite}
                Bc (:,:) double {mustBeReal, mustBeFinite}
                Cc (:,:) double {mustBeReal, mustBeFinite}
                tS (1,1) double {mustBeReal, mustBePositive, mustBeFinite}
                opts.ec         (:,:) double {mustBeReal, mustBeFinite} = []
                opts.fd         (1,1) double {mustBeReal, mustBeFinite} = 0
                opts.directTerm (1,:) double {mustBeReal, mustBeFinite} = 0
            end

            nStates = size(Ac, 1);
            if size(Ac, 2) ~= nStates
                error("rid:surrogate:dimensionMismatch", ...
                    "Ac must be square (nStates x nStates); got %s.", ...
                    mat2str(size(Ac)));
            end
            nInputs = size(Bc, 2);
            if size(Bc, 1) ~= nStates
                error("rid:surrogate:dimensionMismatch", ...
                    "Bc must be (nStates x nInputs) = (%d x nInputs); got %s.", ...
                    nStates, mat2str(size(Bc)));
            end

            ec = opts.ec;
            if isempty(ec)
                ec = zeros(nStates, 1);
            end
            if ~isequal(size(ec), [nStates 1])
                error("rid:surrogate:dimensionMismatch", ...
                    "ec must be (nStates x 1) = (%d x 1); got %s.", ...
                    nStates, mat2str(size(ec)));
            end

            M = expm([Ac, Bc, ec; zeros(nInputs + 1, nStates + nInputs + 1)] * tS);
            Ad = M(1:nStates, 1:nStates);
            Bd = M(1:nStates, nStates + 1 : nStates + nInputs);
            ed = M(1:nStates, end);

            obj = rid.surrogate.LinearStateSpace(Ad, Bd, Cc, tS, ...
                "ed", ed, "fd", opts.fd, "directTerm", opts.directTerm);
        end

        function obj = fitFromStates(u, x, y, tS, opts)
            %FITFROMSTATES  Least squares on MEASURED states.
            %
            %   Two independent linear regressions, both with an intercept:
            %
            %       x(k+1) ~ [x(k), u(k), 1]        k = 1..N-1   -> Ad, Bd, ed
            %       y(k)   ~ [x(k), (u(k)), 1]      k = 1..N     -> Cd, (D), fd
            %
            %   A STATIC FACTORY, NOT `update`, and used OUTSIDE the loop. The update/canUpdate contract
            %   objects to a GENERIC fitter in the framework; this one is specific
            %   to the model it builds, which is what that entry permits. Putting it
            %   here rather than in `experiments/` also puts it under unit test: S1
            %   is evidence for the paper, and `experiments/` has no test layer. The
            %   pilot signal and its data stay there.
            %
            %   THE STATES MUST BE GENUINELY MEASURED. Estimated states must not be
            %   used: an observer depends on the current model, so fitting to its
            %   output would make the fit self-referential.
            %
            %   STATED ONCE, FOR THE PAPER: an equation-error fit on NOISY measured
            %   states is BIASED -- the errors-in-variables problem, since x(k)
            %   appears on both sides. That is acceptable for a deliberately coarse
            %   S1, but S1 must be described as exactly this fit and the pilot noise
            %   level must be recorded.
            %
            %   THE OUTPUT FIT USES ALL N ROWS while the state fit uses N-1: y(k)
            %   needs only x(k), so discarding the last row there would throw away a
            %   measurement for symmetry's sake.
            %
            %   Inputs:
            %       u                (N x nInputs)  applied input, physical
            %       x                (N x nStates)  MEASURED state; row k is the
            %                                       state before u(k) acts
            %       y                (N x 1)        measured output, physical
            %       tS               (1x1)          sample time, > 0
            %       opts.directTerm  (1x1) logical  fit a feedthrough term; false
            %                                       fixes it at exactly 0, so a
            %                                       path the user did not declare
            %                                       cannot appear
            %       opts.minSamples  (1x1)          regression rows required
            %   Outputs:
            %       obj  (1x1) rid.surrogate.LinearStateSpace
            arguments
                u  (:,:) double {mustBeReal, mustBeFinite}
                x  (:,:) double {mustBeReal, mustBeFinite}
                y  (:,:) double {mustBeReal, mustBeFinite}
                tS (1,1) double {mustBeReal, mustBePositive, mustBeFinite}
                opts.directTerm (1,1) logical = false
                opts.minSamples (1,1) double {mustBeInteger, mustBePositive} ...
                    = rid.defaultOptions().surrogate.minUpdateSamples
            end

            nSamples = size(x, 1);
            nStates  = size(x, 2);
            nInputs  = size(u, 2);
            if size(u, 1) ~= nSamples
                error("rid:surrogate:dimensionMismatch", ...
                    "u has %d row(s) and x has %d. The input, the state and the " + ...
                    "output are sampled together.", size(u, 1), nSamples);
            end
            if ~isequal(size(y), [nSamples 1])
                error("rid:surrogate:dimensionMismatch", ...
                    "y must be (N x 1) = (%d x 1); got %s.", ...
                    nSamples, mat2str(size(y)));
            end

            nRows = nSamples - 1;
            if nRows < opts.minSamples
                error("rid:surrogate:insufficientData", ...
                    "The state fit needs at least %d regression rows; the data " + ...
                    "give %d. A model fitted from a handful of samples is not " + ...
                    "worth building.", opts.minSamples, max(nRows, 0));
            end

            % --- the state equation ------------------------------------------
            % The intercept column is what identifies `ed`. Without it an affine
            % term fitted around a non-zero operating point leaks into Ad and Bd,
            % and the model comes out wrong while still looking plausible.
            Phi = [x(1:end-1, :), u(1:end-1, :), ones(nRows, 1)];
            rid.surrogate.LinearStateSpace.assertFullRank(Phi, "state");
            Theta = Phi \ x(2:end, :);                 % (nStates+nInputs+1 x nStates)
            Ad = Theta(1:nStates, :).';
            Bd = Theta(nStates + 1 : nStates + nInputs, :).';
            ed = Theta(end, :).';

            % --- the output equation -----------------------------------------
            % THE FORM FOLLOWS THE DECLARED STRUCTURE, as it does for
            % LinearFirstOrder.update: the u(k) block is present only when a direct
            % term was asked for, so a fit can never introduce a feedthrough path
            % the user did not declare.
            if opts.directTerm
                Psi = [x, u, ones(nSamples, 1)];
            else
                Psi = [x, ones(nSamples, 1)];
            end
            rid.surrogate.LinearStateSpace.assertFullRank(Psi, "output");
            theta = Psi \ y;
            Cd = theta(1:nStates).';
            fd = theta(end);
            if opts.directTerm
                d = theta(nStates + 1 : nStates + nInputs).';
            else
                d = 0;
            end

            % The constructor refuses an unstable Ad, so the fit cannot hand back a
            % model the class would not accept from anyone else.
            obj = rid.surrogate.LinearStateSpace(Ad, Bd, Cd, tS, ...
                "ed", ed, "fd", fd, "directTerm", d);
        end

    end

    methods (Static, Access = private)

        function assertFullRank(M, which)
            %ASSERTFULLRANK  The unidentifiable case, DETECTED from the data.
            %
            %   Not assumed from the structure: a model is refused because THESE
            %   data cannot identify it, which is the discipline. A constant
            %   input is the usual cause -- it is collinear with the intercept
            %   column, and the states have settled besides, so neither Bd nor the
            %   affine term can be separated.
            if rank(M) < size(M, 2)
                error("rid:surrogate:insufficientData", ...
                    "The %s regression matrix is rank deficient (rank %d of " + ...
                    "%d). The data do not excite the model's parameters -- a " + ...
                    "constant input is the usual cause, being collinear with " + ...
                    "the intercept column that identifies the affine term.", ...
                    which, rank(M), size(M, 2));
            end
        end

    end

end
