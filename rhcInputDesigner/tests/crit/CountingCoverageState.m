classdef CountingCoverageState < rid.crit.CoverageState
%COUNTINGCOVERAGESTATE  A CoverageState that counts the distances it computes.
%
%   The O(nSup * nHorizon) contract of the method description -- algorithmic
%   invariant 1 -- cannot be asserted from outside. A violation is INVISIBLE
%   behaviourally: an implementation that scanned the whole history and threw the
%   result away would return exactly the right number, just slowly, and the run
%   would grow steadily slower for no reason a plot would show.
%
%   `rid.crit.CoverageState` routes every distance it computes through one
%   protected seam, `distanceMatrix`, precisely so the count can be taken here.
%   The override changes nothing about the values: it forwards to the base
%   implementation and adds `numel(D)` to a running total.

    properties
        nDistances (1,1) double = 0
    end

    methods (Access = protected)

        function D = distanceMatrix(obj, Znorm, rows)
            D = distanceMatrix@rid.crit.CoverageState(obj, Znorm, rows);
            obj.nDistances = obj.nDistances + numel(D);
        end

    end

    methods

        function reset(obj)
            %RESET  Start counting again, so a setup phase is not charged to a probe.
            obj.nDistances = 0;
        end

    end

end
