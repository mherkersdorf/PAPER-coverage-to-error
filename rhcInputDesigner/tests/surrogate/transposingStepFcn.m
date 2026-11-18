function xNext = transposingStepFcn(x, u, params)
%TRANSPOSINGSTEPFCN  Test fixture: a step function that returns a COLUMN.
%
%   The user's function is the one thing rid.surrogate.StateSpaceFunction cannot
%   validate in advance, so it validates the return per sample. Without that check
%   an (nStates x 1) return would be accepted and silently transpose the
%   trajectory -- plausible numbers, a different model, no error.
%
%   Inputs / Outputs: as `linearStepFcn`, except that `xNext` is (nStates x 1).

xNext = params.Ad * x.' + params.Bd * u.' + params.ed;

end
