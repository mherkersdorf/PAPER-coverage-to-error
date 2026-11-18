classdef ImmediateScheduler < rid.fields.Scheduler
%IMMEDIATESCHEDULER  Blocking refinement-field updates; the default.
%
%   The estimation runs inside `poll`, so the step it was made from and the step
%   it takes effect at are the same: k0 == k1. Those fields are reported anyway,
%   for the same reason `rid.surrogate.ImmediateScheduler` reports them -- the
%   event log must not change shape if an asynchronous variant ever arrives.
%
%   The name matches the surrogate's on purpose. Two seams with one shape is the
%   whole argument; two seams with two names would make a reader look for
%   a difference that is not there.
%
%   POLICY. An estimate is attempted when `updateInterval` steps have passed since
%   the last accepted update. The default is Inf -- nothing is re-estimated unless
%   a run asks for it -- because the mechanism exists from the start but the policy
%   does not. With the default the run is the pure space-filling design:
%   `rid.fields.Nominal` everywhere, `epsReq = epsReq0`, and this class never
%   fires.
%
%   IT ESTIMATES TWO MEMBERS AND PASSES THE THIRD THROUGH. `complexity` and
%   `noiseLevel` are estimated; `relevance` is declared a priori and arrives in the
%   `fields` argument, so this class never holds a copy of it. Holding one would
%   put a declared value in two places, and a divergence would be silent -- see the
%   base class.
%
%   TWO ESTIMATORS, NOT ONE, and they are kept apart deliberately: `complexity` is
%   a BIAS quantity and `noiseLevel` a VARIANCE one, they enter `epsReq` through
%   different exponents (`beta` and `betaNoise`), and an estimator built on
%   residuals must declare which of the two it measures. One object
%   serving both would make that declaration impossible to state.
%
%   A REJECTED ESTIMATE IS REPORTED, NOT RETHROWN, on the precedent. An
%   estimator that cannot honour its contract on the data so far -- too few
%   samples, a region with no residuals in it -- raises; this class catches that,
%   keeps the previous fields, and returns status "rejected" with the message for
%   the event log. One bad estimate must not end a run, and the rejection is
%   visible in the log rather than on the console, which +fields may not write to.
%
%   THE VALUES THEMSELVES ARE NOT VALIDATED HERE. `rid.crit.demandedResolution`
%   owns that check -- shape, finiteness and strict positivity, with the exclusion-through-Q rule
%   message about excluding a region through `Q` rather than through a field -- and
%   a second copy here would be a second definition of a valid field. What IS
%   checked here is completeness of the incoming set, which is the one failure
%   downstream cannot see (base class, `assertCompleteFields`).
%
%   Construction:
%       fieldSch = rid.fields.ImmediateScheduler()
%       fieldSch = rid.fields.ImmediateScheduler("updateInterval", 25, ...
%                      "complexity", myComplexityEstimator)
%
%   Inputs:
%       opts.complexity     (1x1) object  field estimator; default rid.fields.Nominal
%       opts.noiseLevel     (1x1) object  field estimator; default rid.fields.Nominal
%       opts.updateInterval (1x1)         steps between attempts; Inf = never

    properties (SetAccess = immutable)
        complexity                      % (1x1) object, the complexity estimator
        noiseLevel                      % (1x1) object, the noise-level estimator
        updateInterval (1,1) double
    end

    methods

        function obj = ImmediateScheduler(opts)
            arguments
                opts.complexity     (1,1) = rid.fields.Nominal()
                opts.noiseLevel     (1,1) = rid.fields.Nominal()
                opts.updateInterval (1,1) double {mustBeReal, mustBePositive} ...
                    = rid.defaultOptions().fields.updateInterval
            end
            obj.complexity     = opts.complexity;
            obj.noiseLevel     = opts.noiseLevel;
            obj.updateInterval = opts.updateInterval;
        end

        function [fieldsNew, info] = poll(obj, hist, srg, rep, nrm, Q, fields)
            %POLL  Re-estimate the scheduled fields if the policy says so; [] otherwise.
            %
            %   Inputs:
            %       hist       (1x1) struct  the signal history, single source of truth
            %       srg        (1x1) object  the CURRENT surrogate -- after a
            %                                surrogate update in the same step, the
            %                                new one (the refinement-field scheduler ordering)
            %       rep        (1x1) object  the representation
            %       nrm        (1x1) object  the normalization
            %       Q          (nSup x nDim) supporting points, normalized
            %       fields     (1x1) struct  the CURRENT complete set of three
            %   Outputs:
            %       fieldsNew  (1x1) struct or []  a complete set of three, or []
            %       info       (1x1) struct        the event-log record; see Scheduler
            arguments
                obj
                hist   (1,1) struct
                srg    (1,1)
                rep    (1,1)
                nrm    (1,1)
                Q      (:,:) double {mustBeReal, mustBeFinite}
                fields (1,1) struct
            end

            fieldsNew = [];
            info      = rid.fields.Scheduler.emptyInfo();
            % Structural first: a malformed history must not reach an
            % estimator, which would then say something specific and confusing
            % about a struct that is simply not a history.
            rid.assertHistory(hist);
            rid.fields.Scheduler.assertCompleteFields(fields);

            kStep   = size(hist.u, 1);
            info.k0 = kStep;
            info.k1 = kStep;

            if kStep - obj.lastUpdateStep < obj.updateInterval
                return
            end

            try
                complexityValue = obj.complexity.estimate(hist, srg, rep, nrm, Q);
                noiseValue      = obj.noiseLevel.estimate(hist, srg, rep, nrm, Q);
            catch err
                % Refused by the estimator, not by the scheduler: keep the previous
                % fields and record why.
                obj.nRejected = obj.nRejected + 1;
                fieldsNew     = [];
                info.status   = "rejected";
                info.message  = string(err.identifier) + ": " + string(err.message);
                return
            end

            % `relevance` is carried through untouched: it is declared, not
            % estimated, and this class holds no copy of it.
            fieldsNew            = fields;
            fieldsNew.complexity = complexityValue;
            fieldsNew.noiseLevel = noiseValue;

            obj.nUpdates       = obj.nUpdates + 1;
            obj.lastUpdateStep = kStep;
            info.status        = "applied";
            info.updated       = true;
            info.estimated     = ["complexity", "noiseLevel"];
        end

        function s = describe(obj)
            %DESCRIBE  One line for the run record.
            %
            %   Outputs:
            %       s  (1x1) string
            s = sprintf("rid.fields.ImmediateScheduler (updateInterval = %g, " + ...
                "complexity = %s, noiseLevel = %s)", obj.updateInterval, ...
                obj.complexity.describe(), obj.noiseLevel.describe());
        end

    end

end
