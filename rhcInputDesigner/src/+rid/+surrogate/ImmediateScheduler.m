classdef ImmediateScheduler < rid.surrogate.Scheduler
%IMMEDIATESCHEDULER  Blocking surrogate updates; the default.
%
%   The fit runs inside `poll`, so the snapshot and the swap are the same step:
%   k0 == k1 and the model age is always 0. Those fields are reported anyway,
%   because AsyncScheduler fills them with something else and the event log
%   must not change shape when it arrives.
%
%   POLICY. A fit is attempted when `updateInterval` steps have passed since the
%   last accepted swap. The default is Inf -- nothing refits unless a run asks for
%   it -- because the mechanism exists from the start but the policy does not
%  .
%
%   A MODEL THAT CANNOT REFIT IS NOT POLLED. `canUpdate = false` returns
%   status "unsupported" rather than "none", so the event log records WHY the
%   surrogate never changes; without that, "the model never updated" is
%   indistinguishable from a bug.
%
%   A REJECTED FIT IS REPORTED, NOT RETHROWN. An update that would leave
%   the model unstable, unrepresentable or under-determined is refused by the model
%   itself; the scheduler catches that, keeps the previous surrogate, and returns
%   status "rejected" with the model's message for the event log. One bad fit --
%   a stretch of history that happens not to excite the parameters -- must not end
%   a run, and the rejection is visible in the log rather than on the console,
%   which +surrogate may not write to.
%
%   Construction:
%       sch = rid.surrogate.ImmediateScheduler()
%       sch = rid.surrogate.ImmediateScheduler("updateInterval", 50)
%
%   Inputs:
%       opts.updateInterval  (1x1)  steps between fit attempts; Inf = never
%       opts.minSamples      (1x1)  passed through to the model's update

    properties (SetAccess = immutable)
        updateInterval (1,1) double
        minSamples     (1,1) double
    end

    methods

        function obj = ImmediateScheduler(opts)
            arguments
                opts.updateInterval (1,1) double {mustBeReal, mustBePositive} ...
                    = rid.defaultOptions().surrogate.updateInterval
                opts.minSamples (1,1) double {mustBeInteger, mustBePositive} ...
                    = rid.defaultOptions().surrogate.minUpdateSamples
            end
            obj.updateInterval = opts.updateInterval;
            obj.minSamples     = opts.minSamples;
        end

        function [srgNew, info] = poll(obj, hist, srg)
            %POLL  Attempt an update if the policy says so; [] otherwise.
            %
            %   Inputs:
            %       hist    (1x1) struct  the signal history, single source of truth
            %       srg     (1x1) object  the current surrogate
            %   Outputs:
            %       srgNew  (1x1) object or []  a NEW surrogate, or [] if none
            %       info    (1x1) struct        the event-log record; see Scheduler
            arguments
                obj
                hist (1,1) struct
                srg  (1,1)
            end

            srgNew = [];
            info   = rid.surrogate.Scheduler.emptyInfo();
            % This class had no history validation at all: it reads hist.u and
            % trusts it. The history struct puts one definition of well-formed everywhere.
            rid.assertHistory(hist);
            kStep  = size(hist.u, 1);
            info.k0 = kStep;
            info.k1 = kStep;

            if ~srg.canUpdate()
                info.status  = "unsupported";
                info.message = "the surrogate declares canUpdate = false";
                return
            end
            if kStep - obj.lastUpdateStep < obj.updateInterval
                return
            end

            try
                srgNew = srg.update(hist, "minSamples", obj.minSamples);
            catch err
                % Rejected by the model, not by the scheduler: keep the previous
                % surrogate and record why.
                obj.nRejected = obj.nRejected + 1;
                srgNew        = [];
                info.status   = "rejected";
                info.message  = string(err.identifier) + ": " + string(err.message);
                return
            end

            obj.nUpdates       = obj.nUpdates + 1;
            obj.lastUpdateStep = kStep;
            info.status        = "applied";
            info.updated       = true;
            info.modelAge      = 0;   % blocking: snapshot and swap are the same step
        end

        function s = describe(obj)
            %DESCRIBE  One line for the run record.
            %
            %   The record stores what produced a run, and the update POLICY is
            %   part of that: "the surrogate never changed" means something
            %   different at updateInterval = Inf than at 50.
            %
            %   Outputs:
            %       s  (1x1) string
            s = sprintf("rid.surrogate.ImmediateScheduler (updateInterval = %g, " + ...
                "minSamples = %d)", obj.updateInterval, obj.minSamples);
        end

    end

end
