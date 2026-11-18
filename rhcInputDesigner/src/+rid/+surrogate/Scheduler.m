classdef Scheduler < handle
%SCHEDULER  Abstract seam for surrogate update policies.
%
%   Fitting a surrogate can take longer than one iteration of the loop, and the
%   process does not stop meanwhile. Updating therefore runs behind a scheduler:
%
%       sch            = rid.surrogate.ImmediateScheduler(opts)   % blocking, default
%       sch            = rid.surrogate.AsyncScheduler(opts)       % parfeval
%       [srgNew, info] = sch.poll(hist, srg)   % srgNew = [] unless one is ready
%
%   The asynchronous cycle is SNAPSHOT -- FIT -- SWAP: the data are frozen at step
%   k0, the fit runs in the background while the loop keeps designing with the
%   current surrogate, the new surrogate is swapped in at step k1 when it is ready,
%   and the next fit starts from the data acquired in between. Value semantics are
%   what make the swap safe -- nothing mutates under the loop.
%
%   `info` carries the fields the event log needs from EITHER scheduler, so the log
%   shape does not change when AsyncScheduler arrives:
%
%       .status    "none" | "applied" | "rejected" | "unsupported"
%       .updated   (1x1) logical
%       .k0 .k1    (1x1)  snapshot and swap step
%       .modelAge  k1 - k0, in steps; 0 for a blocking scheduler
%       .message   (1x1) string, empty unless something is worth recording
%
%   MODEL AGE is reported on every swap, never inferred: a fit that is slow
%   relative to the process dynamics means the design has been working from a stale
%   model for a long time, and that must be visible.
%
%   THE POLICY IS NOT THE MECHANISM. This seam exists from the start; the
%   trigger that fires it -- fixed interval, residual-based, information-based -- is
%   switched on later and must not be presupposed here. The interval default is
%   therefore Inf: nothing refits unless a run asks for it.
%
%   Schedulers are HANDLE objects. They and rid.crit.CoverageState are the only two
%   places in this design where mutable state is intended.

    properties (SetAccess = protected)
        nUpdates       (1,1) double = 0    % accepted swaps so far
        nRejected      (1,1) double = 0    % fits attempted and refused
        lastUpdateStep (1,1) double = 0    % step of the last accepted swap
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
        [srgNew, info] = poll(obj, hist, srg)
    end

    methods (Static, Access = protected)
        function info = emptyInfo()
            %EMPTYINFO  The event-log record shape, shared by every scheduler.
            info = struct("status", "none", "updated", false, ...
                "k0", NaN, "k1", NaN, "modelAge", NaN, "message", "");
        end
    end

end
