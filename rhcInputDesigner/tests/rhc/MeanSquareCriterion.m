classdef MeanSquareCriterion < rid.crit.Criterion
%MEANSQUARECRITERION  A criterion that keeps every default of rid.crit.Criterion.
%
%   J = mean over committed and candidate points of ||z - target||^2. Not a
%   design criterion anyone would use; it exists so the loop can be driven by a
%   criterion that has NONE of the coverage criterion's features -- no
%   certification, no monotonicity, no resolution, no active set, no supporting
%   points. It is deliberately NOT monotone: committing a point farther
%   from the target than the current mean raises J, so the loop's invariant-3
%   assert would fire if it were not gated by `isMonotone()`.
%
%   With nothing committed, J is 0: the mean over an empty set is taken as the
%   empty sum, which also guarantees that the first commit raises J.

    properties (SetAccess = immutable)
        target (1,:) double
    end

    properties (SetAccess = private)
        sumSq    (1,1) double = 0   % sum of ||z - target||^2 over committed points
        nSamples (1,1) double = 0
    end

    methods

        function obj = MeanSquareCriterion(target)
            obj.target = target;
        end

        function add(obj, Znew)
            obj.sumSq    = obj.sumSq + sum(sum((Znew - obj.target).^2, 2));
            obj.nSamples = obj.nSamples + size(Znew, 1);
        end

        function J = probe(obj, Zcand)
            n = obj.nSamples + size(Zcand, 1);
            if n == 0
                J = 0;
                return
            end
            J = (obj.sumSq + sum(sum((Zcand - obj.target).^2, 2))) / n;
        end

        function [gradCand, info] = probeGradient(obj, Zcand)
            n        = obj.nSamples + size(Zcand, 1);
            gradCand = zeros(size(Zcand));
            if n > 0
                gradCand = 2 * (Zcand - obj.target) / n;
            end
            info = struct("J", obj.probe(Zcand), "nShell", size(Zcand, 1));
        end

        function J = value(obj)
            J = obj.probe(zeros(0, numel(obj.target)));
        end

        function rebuild(obj, Zall)
            obj.sumSq    = 0;
            obj.nSamples = 0;
            obj.add(Zall);
        end

        function fresh = freshCopy(obj)
            fresh = MeanSquareCriterion(obj.target);
        end

        function s = describe(obj)
            s = "MeanSquareCriterion: target " + mat2str(obj.target);
        end

    end

end
