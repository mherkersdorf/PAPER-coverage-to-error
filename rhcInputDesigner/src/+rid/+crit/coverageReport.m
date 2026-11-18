function report = coverageReport(state)
%COVERAGEREPORT  The three-state classification of a COMMITTED coverage state.
%
%       report = rid.crit.coverageReport(state)
%
%   THE OWNER OF THE THREE-STATE LABEL. `rid.crit.evaluate` returns
%   `dTilde` and `psi` but needs `Znorm`, which is derived and never stored;
%   `CoverageState.probe` returns them for a CANDIDATE. Nothing
%   returned them for the state the loop has already committed, so every
%   consumer that wanted to say "resolved" formed the boundary itself -- which
%   is criterion maths outside `+crit`, and a second implementation of the
%   definition that explains a run's whole behaviour. The classification is not
%   a display choice: it is what "resolved" MEANS, and a consumer that
%   re-derives it is one that can disagree with the criterion about the state of
%   the experiment.
%
%   THE THREE STATES, and why they are three and not two:
%
%       resolved   dTilde <= 1               psi = 0; contributes nothing
%       active     1 < dTilde < 1 + psiMax   THE SHELL: the only place a
%                                            gradient exists
%       saturated  dTilde >= 1 + psiMax      psi is at its clip; the gradient is
%                                            zero here too, for the opposite
%                                            reason
%
%   Merging the last two into "unresolved" hides the distinction that explains
%   why a run stalls: an empty active set with everything resolved and an empty
%   active set with everything saturated are opposite diagnoses.
%
%   BOTH BOUNDARIES ARE NON-STRICT, on the side the criterion uses. A point at
%   `dTilde = 1` has `psi = 0` and is not in the shell, so it is RESOLVED --
%   the same comparison `CoverageState.refreshActive` makes. A point at
%   `dTilde = 1 + psiMax` has reached the clip and has no gradient left, so it
%   is SATURATED.
%
%   AN UNCOVERED POINT IS SATURATED. `dRaw` is `Inf` where nothing has been
%   sampled, so `dTilde` is `Inf` and `psi` sits at its clip. This stays correct
%   under `psiMax = Inf`, where the test is `Inf >= Inf`.
%
%   IT READS THE STATE, IT DOES NOT RECOMPUTE IT. The distances come from the
%   cache the run actually used, so this reports the run's own view rather than
%   a second opinion about it. `rid.crit.evaluate` remains the independent
%   recomputation, and `tests/crit/CoverageReportTest` holds the two together.
%
%   Inputs:
%       state  (1x1) rid.crit.CoverageState  a committed state; for a run's
%                                            FINAL coverage this is the state
%                                            after `rid.rhc.run`'s end-of-run
%                                            rebuild
%   Outputs:
%       report  (1x1) struct
%           .nSup         (1x1)
%           .psiMax       (1x1)
%           .dTilde       (nSup x 1)  dRaw ./ epsReq
%           .psi          (nSup x 1)  min(max(dTilde - 1, 0), psiMax)
%           .label        (nSup x 1) string  "resolved" | "active" | "saturated"
%           .isResolved   (nSup x 1) logical
%           .isActive     (nSup x 1) logical
%           .isSaturated  (nSup x 1) logical
%           .stateName    (1 x 3) string  the three, in report order
%           .count        (1 x 3) double  how many in each; sums to nSup

arguments
    state (1,1) rid.crit.CoverageState
end

dTilde = state.dRaw ./ state.epsReq;
psi    = min(max(dTilde - 1, 0), state.psiMax);

% ONE SOURCE FOR THE PARTITION. The masks are derived from the label rather than
% computed alongside it, so the string and the logicals cannot disagree -- which
% would be this function's own reason for existing, reintroduced one level down.
label = repmat("active", state.nSup, 1);
label(dTilde <= 1)                  = "resolved";
label(dTilde >= 1 + state.psiMax)   = "saturated";

report = struct( ...
    "nSup",        state.nSup, ...
    "psiMax",      state.psiMax, ...
    "dTilde",      dTilde, ...
    "psi",         psi, ...
    "label",       label, ...
    "isResolved",  label == "resolved", ...
    "isActive",    label == "active", ...
    "isSaturated", label == "saturated", ...
    "stateName",   ["resolved", "active", "saturated"], ...
    "count",       [sum(label == "resolved"), sum(label == "active"), ...
                    sum(label == "saturated")]);

end
