function assertHistory(hist)
%ASSERTHISTORY  The one structural validator of the signal history.
%
%   The history is a plain struct and gains no methods, so the one thing a class
%   would have bought -- a single place to validate -- comes from here instead.
%   It lives at the package root beside `defaultOptions` because `+rep`, `+gen`
%   and `+surrogate` all consume a history, and having them depend on `+rhc`
%   would invert the dependency direction.
%
%   IT ANSWERS ONE QUESTION: *is this a well-formed history?* -- and answers it
%   the same way everywhere. Whether a well-formed history has what a particular
%   caller needs is a different question with a package-specific answer, and those
%   guards stay where they are:
%
%       rid.gen.Direct              needs a last row of `u`, and says why the
%                                   actuator's starting position is UNRECORDED
%                                   rather than undefined (rid:gen:missingInitialInput)
%       rid.rep.Representation      needs `.y` when pastSource is "measured" and a
%                                   non-empty `.x0` when it is "simulated" -- the
%                                   condition is pastSource, and only the
%                                   representation knows it
%       rid.surrogate.LinearFirstOrder  needs measurements to refit on at all
%
%   CALL IT FIRST. A malformed history that reaches a semantic check produces a
%   confusing symptom instead of a clear one -- "pastSource = measured needs
%   hist.y" is the wrong thing to say about a struct whose `u` is a string.
%
%   THERE IS NO "WHICH FIELDS DO YOU REQUIRE" ARGUMENT. What `pastSource` implies
%   is one fact, and an argument would put it at every call site instead of at the
%   representation that owns it -- the same objection that made the starts
%   provider a closure rather than forwarded name-value pairs.
%
%   REQUIRED, ALWAYS:
%       .u   (nSteps x nInputs)  applied input, physical. `u(1,:)` is the pre-run
%                                actuator position by declared convention.
%       .x0  (1 x nStates)       the RUN's initial surrogate state. Required in
%                                both pastSource modes, because the xStart re-derivation re-derives
%                                the current state from it at every step rather
%                                than carrying `xEnd` forward -- a carried state is
%                                not merely stale after a surrogate update, it is
%                                type-incorrect.
%       .tS  (1x1)               sample time, physical, finite and > 0. The history-as-single-source rule lists
%                                it among the stored signals and `+rep` already
%                                demands it.
%
%   VALIDATED WHEN PRESENT:
%       .y   (nSteps x nOutputs) measured output. ABSENT means offline, and that is
%                                the only way to say it: an empty `y` would be a
%                                sentinel, and sentinels invite `isempty` checks
%                                that conflate "no measurement" with "not yet".
%       .x   (nSteps x nStates)  measured states, physical. Row `k`
%                                is the state measured at the instant `u(k)` was
%                                applied, so it counts the same steps as `u` for
%                                the same reason `y` does. ABSENT means no states
%                                were measured, and is read the same way as an
%                                absent `y`.
%
%   THE WIDTH OF .x IS NOT CHECKED HERE. How many states there are to measure is
%   the representation's question, not a history's: `rid.rep.StateSpace` compares
%   it with its own `nStates` (rid:rep:dimensionMismatch) and decides whether the
%   field has to be there at all (rid:rep:missingStateMeasurement) -- the same
%   split that leaves "pastSource = measured needs y" with `+rep`.
%
%   WHAT IT DELIBERATELY DOES NOT CHECK. Emptiness is not malformation. A history
%   with zero rows is structurally sound; whether a caller can work with one is
%   that caller's question, and `rid.gen.Direct` answers it with a better message
%   than a generic validator could. The single-output convention `size(y,2) == 1`
%   is likewise a property of the current model class rather than of a history, so
%   it stays with the packages that assume it.
%
%   Inputs:
%       hist  (1x1) struct  the signal history
%   Outputs:
%       none; throws on a malformed history

if ~isstruct(hist) || ~isscalar(hist)
    error("rid:invalidHistory", ...
        "hist must be a scalar struct (got %s of size %s).", ...
        class(hist), mat2str(size(hist)));
end

% The hints are not decoration. Absorbing a semantic guard into a structural one
% must not throw away what that guard told the user to DO: the "absent x0" branch
% of rid:rep:missingInitialState is unreachable once x0 is required everywhere, so
% its pointer to the canonical producer moves here.
hints = struct( ...
    "u",  "It carries the applied input, and u(1,:) is the actuator's position " + ...
          "before the run began -- unrecorded, not undefined.", ...
    "x0", "It is the run's initial surrogate state, from which the current state " + ...
          "is re-derived at every step. Obtain it with " + ...
          "srg.steadyState(hist.u(1,:)), the canonical producer.", ...
    "tS", "It is the sample time, in physical units.");

for name = ["u", "x0", "tS"]
    if ~isfield(hist, name)
        error("rid:missingHistoryField", ...
            "hist is missing the required field %s. %s A history carries " + ...
            "u, x0 and tS in every run; y and x are present only when a " + ...
            "process supplied those measurements.", ...
            """" + name + """", hints.(name));
    end
end

assertRealMatrix(hist.u, "u");

if ~isnumeric(hist.tS) || ~isreal(hist.tS) || ~isscalar(hist.tS) || ...
        ~isfinite(hist.tS) || hist.tS <= 0
    error("rid:invalidHistory", ...
        "hist.tS must be a finite, strictly positive scalar sample time.");
end

assertRealMatrix(hist.x0, "x0");
if size(hist.x0, 1) ~= 1
    % The state is a ROW, matching simulate(uSeq, xStart) and the samples-in-rows
    % convention. A column here transposes silently into a multi-sample state.
    error("rid:invalidHistory", ...
        "hist.x0 must be (1 x nStates); got %s. The initial state is a row, " + ...
        "as srg.simulate takes it.", mat2str(size(hist.x0)));
end

if isfield(hist, "y")
    assertRealMatrix(hist.y, "y");
    if size(hist.y, 1) ~= size(hist.u, 1)
        error("rid:invalidHistory", ...
            "hist.y has %d row(s) and hist.u has %d. The two signals are " + ...
            "sampled together, so they carry the same number of steps.", ...
            size(hist.y, 1), size(hist.u, 1));
    end
end

if isfield(hist, "x")
    assertRealMatrix(hist.x, "x");
    if size(hist.x, 1) ~= size(hist.u, 1)
        error("rid:invalidHistory", ...
            "hist.x has %d row(s) and hist.u has %d. Row k of x is the state " + ...
            "measured at the instant u(k) was applied, so the two " + ...
            "carry the same number of steps.", ...
            size(hist.x, 1), size(hist.u, 1));
    end
end

end

% ---------------------------------------------------------------------------

function assertRealMatrix(value, name)
%ASSERTREALMATRIX  A signal is a real, finite, two-dimensional numeric array.
%
%   NaN and Inf are malformation, not missing data: an absent measurement is an
%   absent field (see the header), so a non-finite entry is a corrupted one.
if ~isnumeric(value) || ~isreal(value) || ~ismatrix(value)
    error("rid:invalidHistory", ...
        "hist.%s must be a real numeric matrix (got %s of size %s).", ...
        name, class(value), mat2str(size(value)));
end
if ~all(isfinite(value), "all")
    error("rid:invalidHistory", ...
        "hist.%s contains a non-finite value. An absent measurement is an " + ...
        "absent field, never a NaN.", name);
end

end
