function [values, nOut] = probeArity(scoreFcn, theta, nWanted)
%PROBEARITY  Call scoreFcn asking for nWanted outputs, and report how many it gave.
%
%   WHY THIS IS NOT `nargout(scoreFcn)`. An anonymous handle reports -1 there,
%   and it passes its own nargout through to what it wraps -- `@(t) obj(t)`
%   returns two outputs whenever `obj` does. So the declared arity classifies
%   neither the handle nor the function behind it, and the only reliable question
%   is the one asked by calling.
%
%   The catch is deliberately narrow. `MATLAB:maxlhs` (an anonymous handle) and
%   `MATLAB:TooManyOutputs` (a named function) are the two ways MATLAB says "that
%   many outputs do not exist"; anything else is a genuine failure inside
%   scoreFcn and is rethrown untouched, so a bug in the objective cannot be
%   mistaken for a missing gradient.
%
%   The call is a REAL evaluation and the caller is expected to use its value
%   rather than throw it away -- it goes through the counted handle, so it
%   appears in `nEvaluations` like any other.
%
%   Inputs:
%       scoreFcn  function_handle
%       theta     (nTheta x 1)
%       nWanted   (1x1)  outputs to ask for
%   Outputs:
%       values    (1 x nOut) cell  what came back
%       nOut      (1x1)            nWanted, or fewer if scoreFcn cannot supply it

for nTry = nWanted:-1:1
    try
        values = cell(1, nTry);
        [values{:}] = scoreFcn(theta);
        nOut = nTry;
        return
    catch err
        if err.identifier ~= "MATLAB:maxlhs" && ...
                err.identifier ~= "MATLAB:TooManyOutputs"
            rethrow(err);
        end
    end
end

error("rid:opt:scoreFcnReturnsNothing", ...
    "scoreFcn returned no output at all; it must return at least the objective.");

end
