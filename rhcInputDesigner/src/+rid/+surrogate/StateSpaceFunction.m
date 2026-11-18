classdef StateSpaceFunction < rid.surrogate.Surrogate
%STATESPACEFUNCTION  Any discrete-time state-space model as a surrogate.
%
%       x(k+1) = stepFcn(x, u, params)
%       y(k)   = outputFcn(x, u, params)
%
%   The wrapper that lets a process model be used as its own surrogate -- S0 in
%   the paper's terms. The model itself, a continuous model plus its
%   discretization, is written in `experiments/` as NAMED functions; this class
%   holds them and nothing else.
%
%   Two alternatives were rejected. A CLASS PER PROCESS in `src/` would put
%   paper-specific models into the library. WRAPPING THE PROCESS AD HOC in
%   `experiments/` would be untested, and S0 produces the paper's reference
%   results.
%
%   THE HANDLES MUST BE NAMED AND RESOLVABLE, and both refusals exist for one
%   reason: a run record must be able to rebuild this object. An ANONYMOUS handle
%   captures values the record cannot see (`rid:surrogate:anonymousFunction`); a
%   named handle whose name nothing on the path resolves is no more rebuildable
%   (`rid:surrogate:unresolvableFunction`). `params` is a plain struct, which the
%   record stores and `describe` prints, so the parameters of a run are visible
%   and a replay can reconstruct the model.
%
%   THE MODEL MUST BE A FIXED-STEP DISCRETE MAP. The loop
%   re-simulates the whole history every step, so every sample is a MATLAB
%   call: an adaptive ODE solve per sample would make a long offline run
%   impractical. The remedy -- carrying the trajectory forward -- is
%   adopted only on profiling evidence.
%
%   WHAT IS DECLARED AND NEVER DETECTED. `timeConstant` defaults to NaN and
%   `directTerm` to 0, and both are the user's statements about a model this class
%   knows nothing about. NaN is the honest answer for a black box, and a
%   direct term is a claim about the model's structure, not something a wrapper
%   can infer by probing.
%
%   The optional members are the base-class fallbacks: `steadyState` simulates a
%   constant input to convergence, `isStable` returns the unverified `true`, and
%   `canUpdate` is false with the inherited no-op `update` -- a framework cannot
%   refit a model it does not understand. `defaultRepresentation` stays
%   the inherited NARX one: a class cannot know whether its states are
%   interpretable.
%
%   Construction:
%       srg = rid.surrogate.StateSpaceFunction(stepFcn, outputFcn, params, ...
%                                              nStates, nInputs, tS)
%       srg = rid.surrogate.StateSpaceFunction(..., "timeConstant", T, ...
%                                              "directTerm", d)
%
%   Inputs:
%       stepFcn            (1x1) function_handle  xNext = stepFcn(x, u, params),
%                                                 (1 x nStates), (1 x nInputs),
%                                                 struct -> (1 x nStates)
%       outputFcn          (1x1) function_handle  y = outputFcn(x, u, params) -> (1x1)
%       params             (1x1) struct           passed to both, unchanged
%       nStates            (1x1)                  state dimension, >= 1
%       nInputs            (1x1)                  input channels, >= 1
%       tS                 (1x1)                  sample time, physical, > 0
%       opts.timeConstant  (1x1)                  DECLARED; NaN when unknown
%       opts.directTerm    (1 x nInputs)          DECLARED; 0 = strictly proper

    properties (SetAccess = immutable)
        % The two handles are declared WITHOUT validators, and the typing lives in
        % the `arguments` block below -- the precedent of rid.rhc.Setup and
        % rid.rep.Representation, and for the same reason: a property validated as
        % (1,1) function_handle needs a default value, and there is no meaningful
        % empty function handle to give it. Nothing is lost. The constructor is the
        % only writer, the properties are immutable, and both handles are validated
        % on the way in. Shapes:
        %
        %   stepFcn    (1x1) function_handle  xNext = stepFcn(x, u, params)
        %   outputFcn  (1x1) function_handle  y     = outputFcn(x, u, params)
        stepFcn
        outputFcn
        params    (1,1) struct
    end

    methods

        function obj = StateSpaceFunction(stepFcn, outputFcn, params, ...
                                          nStates, nInputs, tS, opts)
            arguments
                stepFcn   (1,1) function_handle
                outputFcn (1,1) function_handle
                params    (1,1) struct
                nStates   (1,1) double {mustBeInteger, mustBePositive}
                nInputs   (1,1) double {mustBeInteger, mustBePositive}
                tS        (1,1) double {mustBeReal, mustBePositive, mustBeFinite}
                opts.timeConstant (1,1) double {mustBeReal} = NaN
                opts.directTerm   (1,:) double {mustBeReal, mustBeFinite} = 0
            end

            rid.surrogate.StateSpaceFunction.assertRebuildable(stepFcn, "stepFcn");
            rid.surrogate.StateSpaceFunction.assertRebuildable(outputFcn, "outputFcn");

            % The LinearFirstOrder convention: one scalar means the same
            % feedthrough on every channel.
            d = opts.directTerm;
            if isscalar(d) && nInputs > 1
                d = repmat(d, 1, nInputs);
            end
            if numel(d) ~= nInputs
                error("rid:surrogate:dimensionMismatch", ...
                    "directTerm must have one entry per input channel " + ...
                    "(got %d for nInputs = %d).", numel(d), nInputs);
            end

            obj.stepFcn   = stepFcn;
            obj.outputFcn = outputFcn;
            obj.params    = params;

            obj.nStates      = nStates;
            obj.nInputs      = nInputs;
            obj.tS           = tS;
            obj.timeConstant = opts.timeConstant;
            obj.directTerm   = d;
        end

        function [yHat, xEnd] = simulate(obj, uSeq, xStart)
            %SIMULATE  Pure, deterministic forward simulation in physical units.
            %
            %   Purity is the user function's to keep: this class calls it once per
            %   sample with the state and the input, and holds no state of its own
            %   between calls. A persistent variable inside `stepFcn` would break
            %   the contract here as surely as one inside this file.
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
            [yHat, xEnd] = obj.march(uSeq, xStart);
        end

        function [yHat, xEnd, xTraj] = simulateStates(obj, uSeq, xStart)
            %SIMULATESTATES  The same simulation, with the state trajectory.
            %
            %   ONE PRIVATE LOOP WITH `simulate`, so the two cannot drift.
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
            [yHat, xEnd, xTraj] = obj.march(uSeq, xStart);
        end

        function s = describe(obj)
            %DESCRIBE  Both function names and the parameters.
            %
            %   The parameters are printed, not summarized: a run record showing
            %   "some parameters" would leave a reader unable to say which model
            %   produced the result, and a replay unable to rebuild it.
            %
            %   Outputs:
            %       s  (1x1) string
            try
                shown = string(jsonencode(obj.params));
            catch
                % A struct holding something jsonencode cannot render is still
                % worth describing by its field names; saying so beats an error
                % from a diagnostic.
                shown = "<not JSON-encodable; fields: " + ...
                    strjoin(string(fieldnames(obj.params)).', ", ") + ">";
            end
            s = sprintf( ...
                "rid.surrogate.StateSpaceFunction  nStates = %d, nInputs = %d, " + ...
                "tS = %g\n  stepFcn      : %s\n  outputFcn    : %s\n" + ...
                "  params       : %s\n  timeConstant : %g (declared)\n" + ...
                "  directTerm   : %s (declared)\n" + ...
                "  isStable = %d (unverified), canUpdate = %d", ...
                obj.nStates, obj.nInputs, obj.tS, ...
                func2str(obj.stepFcn), func2str(obj.outputFcn), shown, ...
                obj.timeConstant, mat2str(obj.directTerm, 6), ...
                obj.isStable(), obj.canUpdate());
        end

    end

    methods (Access = private)

        function [yHat, xEnd, xTraj] = march(obj, uSeq, xStart)
            %MARCH  The ONE loop of both simulations. Shapes are the callers'.
            %
            %   THE OUTPUT IS EVALUATED BEFORE THE STEP. y(k) is a function of
            %   x(k) -- the state before u(k) acts -- so evaluating it after would
            %   report y(k) from x(k+1) and shift the whole response by one
            %   sample: plausible numbers, a different model, and no error.
            nHorizon = size(uSeq, 1);
            yHat     = zeros(nHorizon, 1);
            xTraj    = zeros(nHorizon, obj.nStates);
            x        = xStart;
            for kStep = 1:nHorizon
                u = uSeq(kStep, :);
                xTraj(kStep, :) = x;

                y = obj.outputFcn(x, u, obj.params);
                obj.assertSample(y, [1 1], "outputFcn", kStep);
                yHat(kStep) = y;

                x = obj.stepFcn(x, u, obj.params);
                obj.assertSample(x, [1 obj.nStates], "stepFcn", kStep);
            end
            xEnd = x;
        end

        function assertSample(obj, value, expected, which, kStep) %#ok<INUSD>
            %ASSERTSAMPLE  The user's return, checked where it is produced.
            %
            %   The wrapped model is the one thing this class cannot validate in
            %   advance, so it is validated per sample. The cost is real -- this is
            %   the inner loop -- and it is paid deliberately: an (nStates x 1)
            %   return would otherwise transpose the trajectory silently, and a
            %   NaN would reach the criterion as a coordinate rather than as an
            %   error.
            if ~isnumeric(value) || ~isreal(value) || ~isequal(size(value), expected)
                error("rid:surrogate:dimensionMismatch", ...
                    "%s returned %s of class %s at sample %d; the contract is a " + ...
                    "real %s. Samples are in ROWS throughout.", ...
                    which, mat2str(size(value)), class(value), kStep, ...
                    mat2str(expected));
            end
            if ~all(isfinite(value))
                error("rid:surrogate:nonFiniteSimulation", ...
                    "%s returned a non-finite value at sample %d. A model that " + ...
                    "diverges mid-horizon must say so here rather than hand the " + ...
                    "criterion a trajectory of NaN to score.", which, kStep);
            end
        end

    end

    methods (Static, Access = private)

        function assertRebuildable(h, role)
            %ASSERTREBUILDABLE  A handle a run record could reconstruct.
            %
            %   Two refusals, one reason. An ANONYMOUS handle captures values that
            %   exist only in the workspace that built it, and a record storing
            %   `@(x,u,p) ...` stores a string nobody can evaluate back into the
            %   same function. A NAMED handle that nothing on the path resolves has
            %   the same defect one step later: the record names a function the
            %   replay cannot find.
            %
            %   Checked at CONSTRUCTION, where the user can still do something
            %   about it -- not at save time, by which point the run has happened.
            name = string(func2str(h));
            if startsWith(name, "@")
                error("rid:surrogate:anonymousFunction", ...
                    "%s is an anonymous function. A run record cannot see the " + ...
                    "values a closure captures, so the surrogate could never be " + ...
                    "rebuilt on replay. Write it as a named function and " + ...
                    "pass its parameters in `params`, which the record stores.", ...
                    role);
            end
            if isempty(which(name))
                error("rid:surrogate:unresolvableFunction", ...
                    "%s names ""%s"", which does not resolve on the path. A " + ...
                    "replay rebuilds the surrogate from this name, so a handle " + ...
                    "nothing resolves is as unrebuildable as a closure. Local " + ...
                    "and nested functions are not reachable by name either.", ...
                    role, name);
            end
        end

    end

end
