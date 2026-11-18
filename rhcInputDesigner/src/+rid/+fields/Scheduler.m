classdef Scheduler < handle
%SCHEDULER  Abstract seam for refinement-field update policies.
%
%   The refinement fields decide how fine a resolution each supporting point
%   demands, and they improve with DATA. Re-estimating them runs behind a
%   scheduler, deliberately the same seam shape as `rid.surrogate.Scheduler`,
%   so there is one pattern to learn rather than two:
%
%       fieldSch             = rid.fields.ImmediateScheduler(...)   % the default
%       [fieldsNew, info]    = fieldSch.poll(hist, srg, rep, nrm, Q, fields)
%
%   `fieldsNew` is `[]` unless an update is ready. When one is, it is a COMPLETE
%   set of three, ready for `rid.crit.demandedResolution`.
%
%   IT IS NOT TIED TO THE SURROGATE'S SCHEDULE, and that is the point of a second
%   seam. The two rates are independent in principle and in practice: complexity
%   and noise estimates improve with data, which arrives every step, while the
%   surrogate is refit rarely by design. Coupling them would let the
%   resolution sharpen only when the model is refit -- backwards, since "where is
%   this process hard?" is answerable from data long before a new model is
%   warranted, and it makes the cheap update wait on the expensive one.
%
%   `poll` RECEIVES THE CURRENT `fields` AND OVERRIDES ONLY WHAT IT ESTIMATES.
%   Only `complexity` and `noiseLevel` are scheduled; `relevance` is declared a
%   priori and never estimated. Letting a
%   scheduler HOLD a declared relevance would put one value in two places -- here
%   and wherever the run's initial `epsReq` was built -- and a divergence would be
%   invisible: a wrong relevance does not raise, it silently weights regions
%   incorrectly and yields a worse experiment forever. Passing the current fields
%   in removes the duplication structurally. That is the coordinate-meaning rule's disease and the coordinate-meaning rule's
%   remedy.
%
%   ORDERING WHEN BOTH SCHEDULERS FIRE IN ONE STEP: surrogate first, then fields.
%   Estimators take `srg` and must see the new one. The cache follows the same
%   order -- `rebuild` (O(nSup * nAll), invariant 3) before `setResolution`
%   (O(nSup)), each recomputing the active set. Two events are logged, in
%   that order. The ordering is the LOOP's to honour; it is written here because
%   this is where a reader looking at the seam will ask.
%
%   `info` carries the fields the event log needs from ANY field scheduler, so the
%   log shape does not change when a second one arrives:
%
%       .status     "none" | "applied" | "rejected"
%       .updated    (1x1) logical
%       .k0 .k1     (1x1)  the step the estimate was made from, and applied at
%       .estimated  (1 x n) string, which members were re-estimated
%       .message    (1x1) string, empty unless something is worth recording
%
%   THERE IS NO `modelAge`. The surrogate seam reports it because a fit can take
%   longer than an iteration and the design then works from a stale model, which
%   must be visible. Field estimation has no asynchronous path, so k0 == k1
%   always; carrying a field that is structurally zero would suggest a latency
%   that does not exist.
%
%   A FIELD UPDATE IS NOT FREE, and the cost is not computational. `epsReq`
%   changes, so `J` is not comparable across it and the monotonicity assert holds
%   only between such changes. A scheduler that
%   fires often makes the `J` trace unreadable -- a reason to schedule
%   deliberately, not a reason to couple.
%
%   THE POLICY IS NOT THE MECHANISM. This seam exists from the start; the
%   trigger that fires it is switched on later and must not be presupposed here.
%   The default interval is Inf: nothing is re-estimated unless a run asks for it,
%   so v1 behaviour is the pure space-filling design of `rid.fields.Nominal`.
%
%   Schedulers are HANDLE objects. They and `rid.crit.CoverageState` are the only
%   places in this design where mutable state is intended.

    properties (SetAccess = protected)
        nUpdates       (1,1) double = 0    % accepted re-estimations so far
        nRejected      (1,1) double = 0    % attempted and refused
        lastUpdateStep (1,1) double = 0    % step of the last accepted update
    end

    methods

        function reset(obj)
            %RESET  Return the scheduler to its pre-run state.
            %
            %   A scheduler is one of the two intended homes of mutable state, so
            %   it accumulates across a run: `nUpdates`, `nRejected` and above all
            %   `lastUpdateStep`, which the interval policy measures from. Reusing
            %   one in a second run would carry that clock forward and silently
            %   change WHEN the second run updates.
            %
            %   It exists for `rid.io.replay`, which must re-run from a record's
            %   own objects and needs them at the state the original run started
            %   from -- not the state it finished in. Restoring by hand at the call
            %   site would put the scheduler's internals outside the scheduler.
            obj.nUpdates       = 0;
            obj.nRejected      = 0;
            obj.lastUpdateStep = 0;
        end

    end

    methods (Abstract)
        [fieldsNew, info] = poll(obj, hist, srg, rep, nrm, Q, fields)
    end

    methods (Static, Access = protected)

        function info = emptyInfo()
            %EMPTYINFO  The event-log record shape, shared by every field scheduler.
            info = struct("status", "none", "updated", false, ...
                "k0", NaN, "k1", NaN, "estimated", strings(1, 0), "message", "");
        end

    end

    methods (Static)

        function assertCompleteFields(fields)
            %ASSERTCOMPLETEFIELDS  All three members present, before any override.
            %
            %   An absent member is NOT an error to `rid.crit.demandedResolution`:
            %   it defaults to nominal. That is right there and wrong here.
            %   A scheduler handed a `fields` without `relevance` would return a set
            %   without one, the declared relevance would silently become ones, and
            %   the run would keep going with a target it was never given. The
            %   member the scheduler must NOT estimate is exactly the member whose
            %   loss cannot be detected downstream.
            %
            %   PUBLIC, because `rid.rhc.Setup` needs the same rule: `fields` is
            %   now one of the eleven arguments, and a set arriving there
            %   incomplete would reach a scheduler that could only pass through
            %   what it was given. One definition, in the package that constructs
            %   fields.
            %
            %   IT CHECKS COMPLETENESS AND NOT VALUES. Shape, finiteness and strict
            %   positivity belong to `rid.crit.demandedResolution`, which owns them
            %   and reports them with the message about excluding a region
            %   through `Q`; a second copy here would be a second definition of a
            %   valid field.
            names = ["relevance", "complexity", "noiseLevel"];
            if ~isstruct(fields) || ~isscalar(fields)
                error("rid:fields:incompleteFields", ...
                    "fields must be a scalar struct carrying %s (got %s).", ...
                    strjoin("""" + names + """", ", "), class(fields));
            end
            missing = names(~isfield(fields, names));
            if ~isempty(missing)
                error("rid:fields:incompleteFields", ...
                    "fields is missing %s. A field scheduler overrides only the " + ...
                    "members it estimates and passes the rest through, so it must " + ...
                    "be handed a COMPLETE set of %s -- an absent member would " + ...
                    "default to nominal downstream and silently discard a " + ...
                    "declared field. Use rid.fields.nominalFields(nSup) " + ...
                    "for the nominal case.", ...
                    strjoin("""" + missing + """", ", "), ...
                    strjoin("""" + names + """", ", "));
            end
        end

    end

end
