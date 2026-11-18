function [counted, readCount] = countedScore(scoreFcn)
%COUNTEDSCORE  Wrap scoreFcn so that every call through it is counted.
%
%   The count is the diagnostic `nEvaluations` of the method description,
%   and it is also what makes contract item 1 testable: a test scoreFcn that
%   counts its own calls must agree with this number, so an optimizer that
%   reimplemented or short-circuited the criterion would disagree.
%
%   Nested rather than anonymous on purpose: an anonymous handle captures its
%   workspace BY VALUE, so `@() n` would return the count as it stood when the
%   handle was made -- always zero. Handles to nested functions share the live
%   workspace, so `readCount` sees what `counted` has done.
%
%   The wrapper is arity-transparent: it requests from `scoreFcn` exactly as many
%   outputs as its own caller asked for, so the [f], [f, g] and [f, g, info]
%   forms all pass through unchanged.
%
%   Inputs:
%       scoreFcn   function_handle
%   Outputs:
%       counted    function_handle  same signature, counted
%       readCount  function_handle  () -> (1x1) double

n         = 0;
counted   = @callAndCount;
readCount = @readN;

    function varargout = callAndCount(theta)
        varargout = cell(1, max(nargout, 1));
        % Counted only AFTER the call returns. An arity probe that asks for
        % one output too many raises at the call boundary without the
        % objective being entered, and counting it would make this number
        % disagree with a scoreFcn that counts itself -- which is the one
        % assertion contract item 1 rests on.
        [varargout{:}] = scoreFcn(theta);
        n = n + 1;
    end

    function c = readN()
        c = n;
    end

end
