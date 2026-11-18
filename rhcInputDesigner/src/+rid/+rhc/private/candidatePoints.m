function [Zphys, Znorm, kIdx] = candidatePoints(ctx, theta)
%CANDIDATEPOINTS  Where a candidate would put samples in regressor space.
%
%   The composition `theta -> uSeq -> yHat -> Zphys -> Znorm` of
%   the method description, evaluated through the SAME `rep` and
%   `srg` the loop uses everywhere else. There is deliberately no second
%   construction path: `rep.buildFuture` was not added, because two paths per
%   representation would have to agree sample for sample and would not error when
%   they diverged -- they would return slightly wrong points and a slightly wrong
%   design.
%
%   THE ROW COUNT IS ASSERTED, not assumed. `build` drops the warm-up rows;
%   the window was sized so that exactly the horizon survives. A history
%   shorter than the warm-up depth returns fewer rows, and every downstream
%   consumer -- the penalty's `gamma`, the affine map's row count -- would then be
%   working on a horizon of a different length.
%
%   BOTH SPACES ARE RETURNED because the two consumers live in different ones: the
%   criterion works in NORMALIZED coordinates and the penalty in PHYSICAL ones.
%   Returning one and converting at the
%   call site is how the two arrays get confused, which is why they no longer share
%   a name.
%
%   THE STATES ARE SPLICED LIKE THE OUTPUTS when the representation reads measured
%   ones: the window's past states are the measured tail and the
%   horizon's are simulated forward from `xNow`, which `stepWindow` anchored on the
%   last measurement. Both halves come from ONE call -- `simulateStates` returns
%   the trajectory beside the outputs it already computed -- so the horizon is
%   never simulated twice for the two quantities.
%
%   Inputs:
%       ctx    (1x1) struct  .gen .rep .nrm .srg .win .nHorizon
%       theta  (nParams x 1) the decision vector
%   Outputs:
%       Zphys  (nHorizon x nDim)  PHYSICAL horizon coordinates
%       Znorm  (nHorizon x nDim)  NORMALIZED horizon coordinates
%       kIdx   (nHorizon x 1)     the window's time index for each row

uSeq = ctx.gen.toTrajectory(theta);

histWin = struct("u", [ctx.win.uTail; uSeq], "x0", ctx.win.xWin, "tS", ctx.win.tS);
if ctx.rep.pastSource == "measured"
    % Measured past, predicted future -- and the future is simulated from xNow,
    % the state the horizon starts from, not from xWin.
    if ctx.rep.needsMeasuredStates()
        [yFuture, ~, xFuture] = ctx.srg.simulateStates(uSeq, ctx.win.xNow);
        histWin.x             = [ctx.win.xTail; xFuture];
    else
        yFuture = ctx.srg.simulate(uSeq, ctx.win.xNow);
    end
    histWin.y = [ctx.win.yTail; yFuture];
end

[Zphys, kIdx] = ctx.rep.build(histWin, ctx.srg);

if size(Zphys, 1) ~= ctx.nHorizon
    error("rid:rhc:historyTooShort", ...
        "The candidate window produced %d regressor row(s) where the horizon " + ...
        "is %d. The window carries %d historical sample(s) and the " + ...
        "representation discards %d as warm-up, so the run needs a longer " + ...
        "priming history before the first horizon is fully representable.", ...
        size(Zphys, 1), ctx.nHorizon, ctx.win.w, ...
        size(histWin.u, 1) - size(Zphys, 1));
end

Znorm = ctx.nrm.toNormalized(Zphys);

end
