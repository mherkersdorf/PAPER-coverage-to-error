classdef NonlinearChainDouble
%NONLINEARCHAINDOUBLE  A theta -> Znorm closure that is deliberately NOT affine.
%
%   Exists for one assertion: that the affine-consistency check of
%   the method description can FAIL. Without it the check could pass on
%   everything and would establish nothing. The nonlinearity is a squared term in
%   the output coordinate, which is what a Hammerstein-like static nonlinearity
%   would do to a tier-1 composition.

    methods (Static)
        function Znorm = evaluate(theta, nHorizon, nInputs)
            uSeq  = reshape(theta, nHorizon, nInputs);
            y     = cumsum(sum(uSeq, 2)) / nHorizon;
            Znorm = [uSeq, y.^2];        % the square is the whole point
        end
    end
end
