function win = stepWindow(hist, srg, rep, w)
%STEPWINDOW  The tail of the history a candidate horizon must be built on.
%
%   Candidate points are built by appending the horizon to the last `w` applied
%   samples and calling the ordinary `rep.build` -- `w` being the warm-up depth
%   `build` already discards, so the window returns exactly `nHorizon`
%   rows. This function assembles everything that does not depend on `theta`, ONCE
%   PER STEP, so the per-evaluation cost is independent of the history length.
%
%   THE STATE IS THE TRAP. The window needs the state at the START OF THE WINDOW,
%   not `hist.x0`. Using the run's initial state simulates as though the window
%   began at time zero: plausible numbers, wrong trajectory, and no error anywhere.
%   Both states are derived here, from `hist.x0` with the CURRENT surrogate, and
%   never carried between steps -- a carried state is not merely stale after a
%   surrogate update, it is type-incorrect.
%
%   ONE PASS OVER THE HISTORY, not two: the prefix is simulated to reach `xWin`,
%   and the tail is simulated onward from there to reach `xNow`. Together they
%   traverse `hist.u` exactly once.
%
%   `yTail` is read only when `pastSource == "measured"`, where the window's `y` is
%   a SPLICE of measured past and predicted future -- the seam.
%   With `"simulated"` the whole window is re-derived by `build` from `xWin` and
%   the measured tail is not consulted at all.
%
%   WITH MEASURED STATES THE HORIZON IS ANCHORED ON A MEASUREMENT.
%   A state-space representation looks back zero samples, so no measured value
%   would reach the horizon through the splice above: measured states would change
%   the realized points and never the design, and a coarse surrogate that had
%   drifted would start every horizon in the wrong place. So `xNow` is ONE
%   SURROGATE STEP FROM THE LAST MEASURED STATE,
%
%       [~, xNow] = srg.simulate(hist.u(end,:), hist.x(end,:))
%
%   and the states of the window's past are the measured tail `xTail`, spliced
%   exactly as `yTail` is. It is still derived from the history at every step and
%   never carried forward, and it is valid only under the assumption of
%   the joint-space representation contract, which no code can check: the surrogate's states are the
%   process's. Nothing simulates the prefix here -- the measurements already say
%   where the window begins -- so this branch is O(w) rather than O(N).
%
%   WHY THE FIRST HORIZON STATE IS PREDICTED AND NOT MEASURED. x(N+1) does not
%   depend on u(N+1) -- a direct term acts on the output, never on the state -- so
%   a process that senses its states could report it before the decision. The
%   optimizer needs time between samples, so the state at the next application
%   instant has to be predicted from the last one available; that is how a
%   receding-horizon loop absorbs its own computation time, and it keeps the process contract's
%   order, apply then measure, a single call. Measuring before deciding is a possible alternative.
%
%   Inputs:
%       hist  (1x1) struct  the signal history
%       srg   (1x1) object  the CURRENT surrogate
%       rep   (1x1) object  the representation; pastSource and
%                           needsMeasuredStates are read
%       w     (1x1) double  window depth in samples, from kIndex (see rid.rhc.run)
%   Outputs:
%       win   (1x1) struct
%           .uTail  (w x nInputs)   the last w applied samples
%           .yTail  (w x 1)         their measurements; empty unless "measured"
%           .xTail  (w x nStates)   their measured states; empty unless the
%                                   representation reads measured states
%           .xWin   (1 x nStates)   the state BEFORE the window
%           .xNow   (1 x nStates)   the state the horizon starts from
%           .tS     (1x1)
%           .w      (1x1)

nHist = size(hist.u, 1);
w     = min(w, nHist);

xTail = zeros(0, 1);
if rep.needsMeasuredStates()
    uTail     = hist.u(nHist - w + 1:end, :);
    xTail     = hist.x(nHist - w + 1:end, :);
    [~, xNow] = srg.simulate(hist.u(end, :), hist.x(end, :));
    if w > 0
        % The window begins at a sample whose state was MEASURED, so there is
        % nothing to simulate to reach it.
        xWin = hist.x(nHist - w + 1, :);
    else
        % An empty window has no state of its own. `xNow` is what `build` will be
        % handed as `x0`, and with no rows before the horizon it is also where the
        % horizon starts.
        xWin = xNow;
    end
else
    if w >= nHist
        % The whole history is shorter than the warm-up depth: the window IS the
        % history, and the state before it is the run's initial state.
        uTail = hist.u;
        xWin  = hist.x0;
    else
        [~, xWin] = srg.simulate(hist.u(1:nHist - w, :), hist.x0);
        uTail     = hist.u(nHist - w + 1:end, :);
    end
    [~, xNow] = srg.simulate(uTail, xWin);
end

yTail = zeros(0, 1);
if rep.pastSource == "measured"
    yTail = hist.y(nHist - w + 1:end, :);
end

win = struct("uTail", uTail, "yTail", yTail, "xTail", xTail, ...
             "xWin", xWin, "xNow", xNow, "tS", hist.tS, "w", w);

end
