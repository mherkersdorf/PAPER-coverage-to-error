function [f, g, info] = stepObjective(ctx, theta)
%STEPOBJECTIVE  The function the optimizer minimizes: J plus the penalties.
%
%       objective(theta) = J(theta) + sum_i w_i * p_i(theta)
%
%   `J` is NEVER a constraint term and the two are logged separately: `J` measures
%   distribution quality and must stay comparable across runs, while the penalty is
%   a property of one declaration.
%
%   THREE ARITIES, and which is in use is discovered by CALLING:
%
%       f            = objective(theta)   % value only; the hot path when tier 3
%       [f, g]       = objective(theta)   % the analytic path
%       [f, g, info] = objective(theta)   % .J .penalty .nShell, at the optimum
%
%   `nargout` is read HERE, inside a real function, which is reliable; the ban is
%   on `nargout(handle)` from outside, where an anonymous wrapper reports -1 and
%   classifies nothing.
%
%   THE VALUE AND THE GRADIENT COME FROM ONE CALL. `probeGradient` returns
%   `info.J`, so the analytic path never asks the criterion twice -- calling
%   `probe` separately for the value would double the dominant cost of every
%   iteration for nothing.
%
%   IT GOES THROUGH THE CACHE, NEVER THROUGH `rid.crit.evaluate`. The reference
%   path walks every sample, `O(nSup * nSamples)`, and `nSamples` grows for the
%   whole run; wiring it in here would violate algorithmic invariant 1, whose
%   symptom is a loop that grows steadily slower for no visible reason.
%
%   THE TWO GRADIENTS LIVE IN DIFFERENT SPACES and are converted before they are
%   added. The criterion's is already normalized; the penalty's is physical and
%   passes through `diag(scale ./ dimWeights)` on its way (section 4). The
%   correction is arithmetically trivial and structurally easy to omit.
%
%   THE HORIZON ROWS ARE SELECTED FROM `kIdx`, not from an assumed position, and
%   the count is asserted against `size(A,1)/nDim`. Here the selection happens to
%   be the identity -- the window yields horizon rows only -- but the assertion is
%   what makes relying on that safe, and it is the check section 3 requires.
%
%   Inputs:
%       ctx    (1x1) struct  .gen .rep .nrm .srg .state .win .nHorizon .A .penaltyOpts
%       theta  (nParams x 1) the decision vector
%   Outputs:
%       f      (1x1) double         the objective
%       g      (nParams x 1) double d(objective)/d(theta)
%       info   (1x1) struct         .J .penalty .nShell, SEPARATELY

wantGradient = nargout >= 2;

[Zphys, Znorm, kIdx] = candidatePoints(ctx, theta);

nShell = NaN;
if wantGradient
    [gradZnorm, critInfo] = ctx.state.probeGradient(Znorm);
    J      = critInfo.J;
    nShell = critInfo.nShell;
else
    J = ctx.state.probe(Znorm);
end

[pen, gradPenPhys] = ctx.gen.penalty(Zphys, ...
    "violationScale", ctx.penaltyOpts.violationScale, ...
    "tailWeight",     ctx.penaltyOpts.tailWeight, ...
    "psiMax",         ctx.penaltyOpts.psiMax);

f = J + pen;

if wantGradient
    gJ = rid.gen.horizonGradient(gradZnorm, kIdx, kIdx, ctx.A);
    gP = rid.gen.gradientToTheta(gradPenPhys, ctx.A, ...
                                 ctx.nrm.span, ctx.nrm.dimWeights);
    g = gJ + gP;
end

if nargout >= 3
    info = struct("J", J, "penalty", pen, "nShell", nShell);
end

end
