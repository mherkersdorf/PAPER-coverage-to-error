classdef RisingCoverageState < rid.crit.CoverageState
%RISINGCOVERAGESTATE  A coverage state whose J RISES as samples are committed.
%
%   The loop's monotonicity assert (invariant 3) can only be exercised by an input
%   that violates it, and the loop is correct: adding samples shrinks
%   nearest-neighbour distances, so no real run makes J rise within a segment.
%   This double produces the violation the assert exists to catch -- the symptom of
%   a stale cache, or of rows committed that were never applied -- and drives the
%   REAL assert through the REAL loop.
%
%   Reproducing the comparison inside the test instead would assert against the
%   test's own copy of it, which is worth nothing.
%
%   Only `value` is altered. `probe` and `probeGradient` are untouched, so the
%   optimizer still sees the true landscape and the run is otherwise ordinary --
%   which is what makes the failure attributable to the assert rather than to a
%   broken search.

    properties (SetAccess = private)
        bump (1,1) double = 0
    end

    properties (SetAccess = immutable)
        increment (1,1) double
    end

    methods

        function obj = RisingCoverageState(Q, epsReq, increment)
            arguments
                Q         (:,:) double
                epsReq    (:,1) double
                increment (1,1) double = 0.1
            end
            obj@rid.crit.CoverageState(Q, epsReq);
            obj.increment = increment;
        end

        function add(obj, Znew)
            add@rid.crit.CoverageState(obj, Znew);
            obj.bump = obj.bump + obj.increment;
        end

        function J = value(obj)
            J = value@rid.crit.CoverageState(obj) + obj.bump;
        end

    end

end
