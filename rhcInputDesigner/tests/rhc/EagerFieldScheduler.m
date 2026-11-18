classdef EagerFieldScheduler < rid.fields.Scheduler
%EAGERFIELDSCHEDULER  A field scheduler that fires on EVERY step.
%
%   The default never fires, which is right for v1 and useless for testing
%   what happens when one does. This one re-estimates every step, returning the
%   fields it was given unchanged, so the RESOLUTION-CHANGE EVENT is produced
%   without the values moving.
%
%   That separation is the point: it makes "the monotonicity assert must not span a
%   logged discontinuity" testable on its own, with nothing else different. It also
%   pins the amended the refinement-field scheduler contract from the caller's side -- it can only pass the
%   fields through because it was given them.

    methods

        function [fieldsNew, info] = poll(obj, hist, ~, ~, ~, ~, fields)
            rid.assertHistory(hist);
            rid.fields.Scheduler.assertCompleteFields(fields);
            kStep              = size(hist.u, 1);
            fieldsNew          = fields;
            obj.nUpdates       = obj.nUpdates + 1;
            obj.lastUpdateStep = kStep;
            info               = rid.fields.Scheduler.emptyInfo();
            info.status        = "applied";
            info.updated       = true;
            info.k0            = kStep;
            info.k1            = kStep;
            info.estimated     = ["complexity", "noiseLevel"];
        end

        function s = describe(~)
            s = "EagerFieldScheduler (fires every step)";
        end

    end

end
