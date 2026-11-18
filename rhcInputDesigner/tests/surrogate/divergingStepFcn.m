function xNext = divergingStepFcn(x, u, ~)
%DIVERGINGSTEPFCN  Test fixture: a step function that reaches a non-finite state.
%
%   Stands in for a model that blows up mid-horizon -- an unstable map, a division
%   by a state that reached zero, an overflow. The surrogate must say so at the
%   sample it happened rather than return a trajectory full of NaN for the
%   criterion to score: a NaN coordinate is corruption, not missing data, which is
%   the same judgement `rid.assertHistory` makes about a signal.
%
%   Inputs / Outputs: as `linearStepFcn`.

xNext = x ./ (u - u) + 0 * x;   % 0/0 -> NaN on the first sample

end
