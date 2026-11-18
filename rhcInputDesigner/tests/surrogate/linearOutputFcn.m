function y = linearOutputFcn(x, u, params)
%LINEAROUTPUTFCN  Test fixture: y(k) = Cd x(k) + directTerm u(k).' + fd.
%
%   The partner of `linearStepFcn`. It reads x(k) -- the state BEFORE this
%   sample's input acts -- which is what makes the evaluation order inside
%   rid.surrogate.StateSpaceFunction observable: evaluating the output after the
%   step would silence one sample of the response.
%
%   Inputs:
%       x       (1 x nStates)  state before this sample's input acts
%       u       (1 x nInputs)  input applied at this sample
%       params  (1x1) struct   .Cd .directTerm .fd
%   Outputs:
%       y       (1x1)

y = params.Cd * x.' + params.directTerm * u.' + params.fd;

end
