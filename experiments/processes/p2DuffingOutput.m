function y = p2DuffingOutput(x, u, p) %#ok<INUSD>
%P2DUFFINGOUTPUT  Output of the P2 Duffing process, S0's outputFcn.
%
%   y = x1, no direct term.
%
%   Inputs:
%       x  (1x2) double  [x1, x2]
%       u  (1x1) double  unused; the output has no direct term
%       p  (1x1) struct  p2DuffingParams()
%   Outputs:
%       y  (1x1) double
arguments
    x (1,2) double {mustBeReal, mustBeFinite}
    u (1,1) double {mustBeReal, mustBeFinite}
    p (1,1) struct
end

y = x(1);
end
