function report = resolutionReport(epsReq0, fillDistanceQ, mode, opts)
%RESOLUTIONREPORT  The run-start check on the baseline against the discretization.
%
%   Called by +rhc once at run start. It states what this discretization
%   can resolve and refuses a certification run that claims more.
%
%   The finest resolution that can be certified is fillDistanceQ: a resolved
%   supporting point certifies coverage only within epsReq(j) of itself, so no
%   claim finer than the ruler can be made, whatever the data.
%   Since refinement only ever makes epsReq finer,
%   epsReq0 is the COARSEST demand in the field, so in certification mode
%   `epsReq0 < fillDistanceQ` means every epsReq(j) would be floored: the
%   refinement structure is erased and the run becomes a uniform design at the
%   discretization limit while claiming otherwise. That is an error
%   (`rid:crit:uncertifiableBaseline`), and it is the only check available before
%   the fields have been estimated. Partial flooring is NOT an error -- it is the
%   intended behaviour of the certifiability floor and is reported by rid.crit.demandedResolution as
%   the floored fraction. An error threshold on that fraction would mean inventing
%   a number. The remedy is a coarser epsReq0 or a larger nSup.
%
%   In budget mode the same configuration is a REPORT, not an error: nothing is
%   certified there, so nothing is over-claimed.
%
%   `fillDistanceQ` is a Monte-Carlo estimate and is biased low, so the
%   headroom and the bound are accurate to that estimate and are not a guarantee.
%   This is why the message says "up to the accuracy of that estimate": the
%   inequality eps <= epsReq0 + fillDistanceQ is the theoretical justification of
%   the method, not an operational certificate.
%
%   `estimateConverged` IS THAT ACCURACY, MADE TO ACT. The estimator draws
%   its own conclusion from the quarter/half/full probe trace, and this is where it
%   lands: a headroom computed from an estimate that was still climbing is a
%   headroom computed from the wrong number, so the run-start statement says so
%   rather than reporting the ratio as if it were settled. It does NOT change the
%   error condition -- an under-resolved probe set makes fillDistanceQ too small
%   and the floor look finer than it is, which is worth saying loudly, not worth
%   refusing to run over.
%
%   `relativeClimb` IS QUOTED BESIDE THE FLAG, and that is not decoration (the convergence-trace safeguard). `traceTolerance = 0.02` is a chosen value, not a sourced one, and at
%   defaults it fires on most runs -- which is the estimator telling the truth, not
%   noise. A boolean that fires on most runs becomes wallpaper; a magnitude stays
%   readable, and a 2.1 %% climb and a 40 %% climb are not the same news.
%
%   It defaults to `true` because the analytic grid path has no probe set
%   and nothing to converge. The Monte-Carlo path must pass the estimator's own
%   `info.converged`, and +rhc does: the default is not load-bearing anywhere in
%   `src/`.
%
%   `probeFillDistance` is the same discipline applied to the OTHER half of the
%   estimate's accuracy. The probe set has a fill distance of its own, and
%   because the distance-to-nearest function is 1-Lipschitz the true value cannot
%   exceed fillDistanceQ + hProbe. The run-start statement therefore reports a
%   BRACKET rather than a point value, and the coverage bound with it: the
%   pessimistic end is what a reader needs in order to judge whether a headroom is
%   real. Neither end changes the error condition -- the floor and the certifiability floor
%   projection continue to use the measured fillDistanceQ, so nothing downstream
%   moves; this is a statement, not a new threshold.
%
%   NOT A RIGOROUS BOUND: hProbe is itself a Monte-Carlo estimate and is biased low
%   in the same way, so the upper end is the best available statement and not a
%   certificate. Pass NaN -- the default -- when there is no probe set, as
%   on the analytic grid path, and nothing is claimed.
%
%   The flag is deliberately NOT threaded into rid.crit.demandedResolution, whose
%   `info.headroom` is the same run constant computed per step. One fact, one
%   source: the run-start statement is where a reader learns the estimate's
%   standing, and a copy per step would be the duplication this project keeps
%   paying for.
%
%   The function does not print. +crit may not write to the console,
%   so the required run-start statement travels as `report.message` for +rhc and
%   +io to place in the log and the run record.
%
%   Inputs:
%       epsReq0        (1x1)  the baseline commitment, finite and > 0
%       fillDistanceQ  (1x1)  measured fill distance of Q, finite and > 0
%       mode           (1x1) string  "certification" | "budget"
%       opts.estimateConverged (1x1) logical  the verdict of
%                             rid.roi.estimateFillDistance on its probe trace;
%                             true for the analytic grid value
%       opts.relativeClimb (1xN) the trace's consecutive relative increments, the
%                             quantity that verdict was taken on; [] when there is
%                             no trace
%       opts.probeFillDistance (1x1) hProbe, the fill distance of the probe set
%                             that produced fillDistanceQ; NaN when there is no
%                             probe set (the analytic grid value)
%
%   Outputs:
%       report  (1x1) struct with fields
%           .mode               (1x1) string
%           .epsReq0            (1x1)
%           .fillDistanceQ      (1x1)
%           .finestCertifiable  (1x1)  = fillDistanceQ
%           .headroom           (1x1)  epsReq0 / fillDistanceQ
%           .bound              (1x1)  epsReq0 + fillDistanceQ, reported alongside
%                                      every J = 0 certification
%           .baselineBelowFloor (1x1) logical  budget mode only; in certification
%                                      mode this configuration has already errored
%           .fillDistanceConverged (1x1) logical  as passed in
%           .relativeClimb      (1xN)  as passed in
%           .probeFillDistance  (1x1)  as passed in, or NaN
%           .fillDistanceBracket (1x2) [fillDistanceQ, fillDistanceQ + hProbe],
%                                      or [NaN NaN] when hProbe was not supplied
%           .message            (1x1) string   the run-start statement

arguments
    epsReq0       (1,1) double {mustBeReal}
    fillDistanceQ (1,1) double {mustBeReal}
    mode          (1,1) string
    opts.estimateConverged (1,1) logical = true
    opts.probeFillDistance (1,1) double {mustBeReal} = NaN
    opts.relativeClimb (1,:) double {mustBeReal} = []
end

if ~isfinite(epsReq0) || epsReq0 <= 0
    error("rid:crit:invalidResolution", ...
        "epsReq0 must be finite and strictly positive (got %g).", epsReq0);
end
if ~isfinite(fillDistanceQ) || fillDistanceQ <= 0
    error("rid:crit:invalidResolution", ...
        "fillDistanceQ must be finite and strictly positive (got %g). " + ...
        "It is measured from Q, never assumed.", fillDistanceQ);
end
if ~isnan(opts.probeFillDistance) && ...
        (~isfinite(opts.probeFillDistance) || opts.probeFillDistance < 0)
    error("rid:crit:invalidResolution", ...
        "probeFillDistance must be NaN (no probe set), or finite and " + ...
        "nonnegative (got %g).", opts.probeFillDistance);
end
if ~ismember(mode, ["certification", "budget"])
    error("rid:crit:unknownMode", ...
        "mode must be ""certification"" or ""budget"" (got ""%s"").", mode);
end

baselineBelowFloor = epsReq0 < fillDistanceQ;

if mode == "certification" && baselineBelowFloor
    error("rid:crit:uncertifiableBaseline", ...
        "epsReq0 = %g is finer than fillDistanceQ = %g: this discretization " + ...
        "cannot certify the baseline, so every epsReq would be floored and the " + ...
        "refinement structure erased. Coarsen epsReq0 or grow nSup.", ...
        epsReq0, fillDistanceQ);
end

headroom = epsReq0 / fillDistanceQ;
bound    = epsReq0 + fillDistanceQ;

if mode == "certification"
    message = sprintf( ...
        "certification-driven run: the finest resolution this discretization " + ...
        "can resolve is fillDistanceQ = %g, so the achieved resolution is " + ...
        "reported up to the accuracy of that Monte-Carlo estimate, never as a " + ...
        "guarantee. Refinement headroom epsReq0/fillDistanceQ = %.3g. " + ...
        "On J = 0 the coverage bound is eps <= epsReq0 + fillDistanceQ = %g.", ...
        fillDistanceQ, headroom, bound);
else
    message = sprintf( ...
        "budget-driven run: nothing is certified, the floor is not " + ...
        "applied, and epsReq only shapes the distribution. fillDistanceQ = %g, " + ...
        "refinement headroom epsReq0/fillDistanceQ = %.3g.", ...
        fillDistanceQ, headroom);
    if baselineBelowFloor
        message = message + sprintf( ...
            " Note: epsReq0 = %g is finer than fillDistanceQ, so Q cannot " + ...
            "resolve the baseline and the criterion is partially blind -- " + ...
            "samples may fall between supporting points and register no gain. " + ...
            "This is a report, not an error, because nothing is certified here.", ...
            epsReq0);
    end
end

% the shortfall of a probe maximum is bounded by the probe set's own fill
% distance, so the floor is stated as a bracket rather than as a point value.
hasBracket = ~isnan(opts.probeFillDistance);
if hasBracket
    fillDistanceBracket = fillDistanceQ + [0, opts.probeFillDistance];
    message = message + sprintf( ...
        " The probe set's own fill distance is %g, and the distance-to-nearest " + ...
        "function is 1-Lipschitz, so fillDistanceQ lies in [%g, %g] and the " + ...
        "coverage bound in [%g, %g]. hProbe is itself estimated and " + ...
        "biased low, so the upper end is the best available statement, not a " + ...
        "certificate.", ...
        opts.probeFillDistance, fillDistanceBracket(1), fillDistanceBracket(2), ...
        epsReq0 + fillDistanceBracket(1), epsReq0 + fillDistanceBracket(2));
else
    fillDistanceBracket = [NaN NaN];
end

% the estimate's standing travels with the headroom, in the same sentence,
% and the MAGNITUDE travels with the flag -- a boolean that fires on most runs is
% wallpaper, a percentage is not.
climbText = "";
if ~isempty(opts.relativeClimb)
    climbText = " Trace climb " + ...
        join(compose("%.2f %%", 100 * opts.relativeClimb(:).'), ", ") + ...
        ", taken over the stochastic probes; the verdict is on the largest.";
end

if ~opts.estimateConverged
    message = message + sprintf( ...
        " WARNING: fillDistanceQ = %g is an UNCONVERGED estimate -- its probe " + ...
        "maximum was still climbing across the quarter/half/full trace by more " + ...
        "than the declared tolerance. The true fill distance is larger, " + ...
        "so the certifiable floor is reported finer than it is and the headroom " + ...
        "above is optimistic. Raise opts.roi.nProbes and re-measure.", ...
        fillDistanceQ) + climbText;
else
    message = message + climbText;
end

report = struct( ...
    "mode",               mode, ...
    "epsReq0",            epsReq0, ...
    "fillDistanceQ",      fillDistanceQ, ...
    "finestCertifiable",  fillDistanceQ, ...
    "headroom",           headroom, ...
    "bound",              bound, ...
    "baselineBelowFloor", baselineBelowFloor, ...
    "fillDistanceConverged", opts.estimateConverged, ...
    "relativeClimb",      opts.relativeClimb, ...
    "probeFillDistance",  opts.probeFillDistance, ...
    "fillDistanceBracket", fillDistanceBracket, ...
    "message",            string(message));

end
