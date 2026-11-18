function xNext = gapLocusStepFcn(x, u, params)
%GAPLOCUSSTEPFCN  A two-state map whose rest locus LEAVES the state box and re-enters elsewhere.
%
%       x(k+1) = a x(k) + (1 - a) r(u(k)),    r(u) = [g1 u,  c - b u^2]
%
%   The fixed point under a constant input is r(u), so the steady-state locus of
%   the joint input--state space is a parabola: `x2` peaks near u = 0 and falls
%   away on both sides, while `x1` sweeps the input range. With the state box
%   declared narrower than that peak, the locus lies INSIDE the box over an
%   interval of negative u, leaves it, and comes back inside over an interval of
%   positive u -- re-entering at a DIFFERENT `x1`, because `x1` went on moving
%   while the locus was outside.
%
%   THAT IS WHAT MAKES THE GAP CLAUSE OF THE JOINT-SPACE REPRESENTATION CONTRACT PART 9 TESTABLE. Arc length is
%   accumulated over CONSECUTIVE probes only: the chord across the gap spans a
%   piece of locus that was never probed inside the box, so it is not arc length
%   and is not added. A construction that bridged it would measure a longer locus
%   and choose different levels. Were the exit and the re-entry at the same
%   point -- which is what a one-state locus is forced to do, since it cannot
%   leave an interval and return without crossing back through it -- the chord
%   would be zero and either rule would give the same answer. Two states are the
%   smallest fixture in which the two rules can disagree.
%
%   The map is a contraction in `x` at rate `a`, so the base-class `steadyState`
%   fallback settles from any start -- which is all the locus construction asks
%   of a black-box model.
%
%   NAMED, NOT ANONYMOUS, because rid.surrogate.StateSpaceFunction refuses a
%   closure: a run record cannot see the values an anonymous handle captures, so
%   a replay could not rebuild the object.
%
%   Inputs:
%       x       (1 x 2)      state before this sample's input acts
%       u       (1 x 1)      input applied at this sample
%       params  (1x1) struct .a .g1 .c .b, plus the .Cd .directTerm .fd that
%                            `linearOutputFcn` reads
%   Outputs:
%       xNext   (1 x 2)

xNext = params.a * x + (1 - params.a) * gapLocusRest(u, params);

end

% ---------------------------------------------------------------------------

function rest = gapLocusRest(u, params)
%GAPLOCUSREST  The rest state under a constant input: the locus itself.
rest = [params.g1 * u, params.c - params.b * u^2];

end
