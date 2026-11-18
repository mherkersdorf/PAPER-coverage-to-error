function xNext = linearStepFcn(x, u, params)
%LINEARSTEPFCN  Test fixture: x(k+1) = Ad x(k) + Bd u(k).' + ed, as a NAMED function.
%
%   Named, not anonymous, because rid.surrogate.StateSpaceFunction refuses a
%   closure: a run record cannot see the values an anonymous handle captures, so a
%   replay could not rebuild the object. The parameters travel in
%   `params`, which the record stores and `describe` prints.
%
%   It reproduces rid.surrogate.LinearStateSpace exactly, so the two classes can be
%   compared on one model -- which is what pins the evaluation ORDER inside the
%   wrapper's loop.
%
%   Inputs:
%       x       (1 x nStates)  state before this sample's input acts
%       u       (1 x nInputs)  input applied at this sample
%       params  (1x1) struct   .Ad .Bd .ed
%   Outputs:
%       xNext   (1 x nStates)

xNext = (params.Ad * x.' + params.Bd * u.' + params.ed).';

end
