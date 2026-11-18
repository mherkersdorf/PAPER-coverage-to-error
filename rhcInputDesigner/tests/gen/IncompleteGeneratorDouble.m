classdef IncompleteGeneratorDouble < rid.gen.Generator
%INCOMPLETEGENERATORDOUBLE  A parametrization that omits `applied`.
%
%   It implements every member of the contract except the one that declares the
%   applied unit, so MATLAB must refuse to instantiate it. That refusal is the
%   only evidence that `rid.gen.Generator`'s abstract block is doing work: with an
%   empty block this class would construct, and a parametrization that never says
%   what it commits would reach the loop.
%
%   `applied` is the member left out deliberately. It is the one whose absence
%   would otherwise be filled in by the loop assuming "the first sample",
%   which is exactly the step-parametrization assumption the interface exists to
%   remove.

    properties (SetAccess = immutable)
        nHorizon
    end

    methods

        function obj = IncompleteGeneratorDouble()
            obj.nHorizon = 1;
        end

        function p = nInputs(~),                      p = 1;                end
        function n = nParams(obj),                    n = obj.nHorizon;     end
        function uSeq = toTrajectory(~, theta),       uSeq = theta;         end
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
        function s = describe(~),                     s = "incomplete";     end

    end

end
