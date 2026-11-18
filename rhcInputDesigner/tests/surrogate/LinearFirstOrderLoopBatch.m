classdef LinearFirstOrderLoopBatch < rid.surrogate.LinearFirstOrder
%LINEARFIRSTORDERLOOPBATCH  Test double exposing the BASE simulateBatch fallback.
%
%   rid.surrogate.LinearFirstOrder overrides simulateBatch with a vectorized
%   columnwise `filter`. The batch-simulation contract requires an override to agree with the loop it
%   replaces, and the only way to compare the two on the same model is to reach
%   past the override to rid.surrogate.Surrogate's implementation.

    methods
        function obj = LinearFirstOrderLoopBatch(varargin)
            obj@rid.surrogate.LinearFirstOrder(varargin{:});
        end

        function [Y, Xend, XTraj] = simulateBatch(obj, uSeqBatch, xStart)
            % ALL THREE outputs are forwarded, and the third only when it is
            % asked for: the base fallback builds XTraj on `nargout >= 3`, so
            % requesting it unconditionally here would change the very branch
            % this double exists to expose.
            if nargout >= 3
                [Y, Xend, XTraj] = ...
                    simulateBatch@rid.surrogate.Surrogate(obj, uSeqBatch, xStart);
            else
                [Y, Xend] = ...
                    simulateBatch@rid.surrogate.Surrogate(obj, uSeqBatch, xStart);
            end
        end
    end

end
