function srg = p2DuffingS0(p)
%P2DUFFINGS0  S0 -- the Duffing process itself, as a StateSpaceFunction surrogate.
%
%   Builds rid.surrogate.StateSpaceFunction the same way every script that needs
%   S0 for this process does, so no two callers can assemble it differently
%  .
%
%   Inputs:
%       p    (1x1) struct  p2DuffingParams()
%   Outputs:
%       srg  (1x1) rid.surrogate.StateSpaceFunction
arguments
    p (1,1) struct
end

srg = rid.surrogate.StateSpaceFunction(@p2DuffingStep, @p2DuffingOutput, p, ...
    2, 1, p.tS, "timeConstant", p.timeConstant, "directTerm", 0);
end
