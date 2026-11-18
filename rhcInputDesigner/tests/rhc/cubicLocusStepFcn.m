function xNext = cubicLocusStepFcn(x, u, params)
%CUBICLOCUSSTEPFCN  A one-state map whose steady-state locus is x1 = u^3.
%
%       x(k+1) = a x(k) + (1 - a) u(k)^3,    0 < a < 1
%
%   The fixed point under a constant input is x = u^3, so the steady-state locus
%   of the joint input--state space is a cubic through the ROI rather than a
%   straight line. THAT IS WHAT MAKES THE TWO CANDIDATE TILINGS DISTINGUISHABLE:
%   levels spread equally along the locus are clustered in `u`,
%   and levels spread equally in `u` are clustered along the locus. For a linear
%   model the two sets coincide, so a linear fixture cannot tell them apart and
%   would pass either implementation.
%
%   The map is a contraction in `x` at rate `a`, so the base-class `steadyState`
%   fallback settles from any start -- which is what the locus construction needs
%   from a black-box model.
%
%   NAMED, NOT ANONYMOUS, because rid.surrogate.StateSpaceFunction refuses a
%   closure: a run record cannot see the values an anonymous handle captures, so
%   a replay could not rebuild the object.
%
%   Inputs:
%       x       (1 x 1)      state before this sample's input acts
%       u       (1 x 1)      input applied at this sample
%       params  (1x1) struct .a, plus the .Cd .directTerm .fd that
%                            `linearOutputFcn` reads
%   Outputs:
%       xNext   (1 x 1)

xNext = params.a * x + (1 - params.a) * u^3;

end
