classdef MinimalSurrogate < rid.surrogate.Surrogate
%MINIMALSURROGATE  Test double: a surrogate that can ONLY simulate.
%
%   This class is the point. It implements the two abstract members --
%   `simulate` and `describe` -- and nothing else, standing in for the
%   vendor-supplied models the open interface exists for: an FMU, a compiled
%   simulator, a lookup-table model. Everything the tests exercise through it comes
%   from the base class's documented fallbacks: `simulateBatch` (loop),
%   `steadyState` (simulate to convergence), `isStable` (true, unverified),
%   `canUpdate` (false), `update` (no-op) and `defaultRepresentation`.
%
%   It has TWO states on purpose: nStates is never to be assumed 1, and
%   defaultRepresentation must derive the lag orders from it. Its order also has
%   nothing to do with the dimension of the representation it is paired with, which
%   is the separation the method description insists on.
%
%       x(k+1) = x(k) .* pole + u(k)
%       y(k)   = sum(x(k))

    properties (SetAccess = immutable)
        pole (1,2) double = [0.5, 0.8]
    end

    methods
        function obj = MinimalSurrogate(declareStates)
            % `declareStates = false` stands in for a model that cannot even say
            % how many states it has -- a black-box simulator. It must still be
            % constructible; only the two fallbacks that need the state dimension
            % may then refuse.
            arguments
                declareStates (1,1) logical = true
            end
            if declareStates
                obj.nStates = 2;
            end
            obj.nInputs = 1;
            obj.tS      = 0.1;
        end

        function [yHat, xEnd] = simulate(obj, uSeq, xStart)
            nHorizon = size(uSeq, 1);
            yHat     = zeros(nHorizon, 1);
            x        = xStart;
            for kStep = 1:nHorizon
                yHat(kStep) = sum(x);
                x = x .* obj.pole + uSeq(kStep, :);
            end
            xEnd = x;
        end

        function s = describe(obj)
            s = sprintf("MinimalSurrogate (test double), nStates = %d", obj.nStates);
        end
    end

end
