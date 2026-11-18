classdef ConstraintlessGeneratorDouble < rid.gen.Generator
%CONSTRAINTLESSGENERATORDOUBLE  A parametrization that omits `constraints`.
%
%   The mirror of `IncompleteGeneratorDouble`, and deliberately a separate class:
%   that one omits `applied` and shows the abstract METHOD block is enforced,
%   this one omits `constraints` and shows the abstract PROPERTY block is. A
%   single double missing both could not tell which block did the refusing.
%
%   `constraints` is the member left out because its absence is the one that
%   USED to be silent. Before the constraints contract the box lived on `rid.gen.Direct` alone, so
%   a generator without it did not raise -- `rid.viz.signal` simply drew no
%   bounds, and a run pinned against its actuator limits became
%   indistinguishable from one that never approached them. The limits belong to
%   the actuator, not to the parametrization; only how they are ENFORCED in a
%   given parameter space is parametrization-specific, and that is what
%   `feasibleSet` returns.

    properties (SetAccess = immutable)
        nHorizon
    end

    methods

        function obj = ConstraintlessGeneratorDouble()
            obj.nHorizon = 1;
        end

        function p = nInputs(~),                      p = 1;                end
        function n = nParams(obj),                    n = obj.nHorizon;     end
        function uSeq = toTrajectory(~, theta),       uSeq = theta;         end
        function uApp = applied(~, uSeq),             uApp = uSeq(1, :);    end
        function theta = shift(~, thetaPrev),         theta = thetaPrev;    end
        function [lb, ub, A, b] = feasibleSet(obj, ~)
            lb = -ones(obj.nParams(), 1);
            ub =  ones(obj.nParams(), 1);
            A  = zeros(0, obj.nParams());
            b  = zeros(0, 1);
        end
        function [Theta0, info] = initialGuesses(obj, ~, ~, opts)
            arguments
                obj
                ~
                ~
                opts.nStarts (1,1) double = 1
            end
            Theta0 = zeros(opts.nStarts, obj.nParams());
            info   = struct("seed", 0, "nStarts", opts.nStarts, "warmStart", "none");
        end
        function [pen, gradPen, info] = penalty(~, Zphys, ~)
            pen     = 0;
            gradPen = zeros(size(Zphys));
            info    = struct();
        end
        function s = describe(~),                     s = "constraintless"; end

    end

end
