classdef SurrogateDependentRep < rid.rep.Representation
%SURROGATEDEPENDENTREP  Test double: coordinates that DO depend on the surrogate.
%
%   Stands in for the fourth row of the rebuild table:
%   a representation whose COORDINATES are
%   reconstructed using surrogate parameters -- an estimated state, an OBF basis
%   whose poles come from the surrogate, a transformed axis -- so that refitting
%   moves the historical points even though the output was measured.
%
%   No such representation exists in this layer, and that is exactly why the double is
%   needed: without it, `requiresRebuildOnSurrogateUpdate` could return
%   `pastSource == "simulated"` and every test would still pass.

    properties (SetAccess = immutable)
        nDim               = 1
        dependsOnSurrogate = true
    end

    methods
        function obj = SurrogateDependentRep(pastSource)
            arguments
                pastSource (1,1) string = "measured"
            end
            obj@rid.rep.Representation(pastSource);
        end

        function info = coordinateInfo(~)
            % requiredNames is inherited and derived from this.
            info = SurrogateDependentRep.makeEntry("u1_k", "inputLevel", 1, 0);
        end
        function [box, dimName] = bounds(obj, roi)
            dimName = obj.requiredNames();
            box     = obj.lookupRanges(roi, dimName);
        end
        function [Zphys, kIndex] = build(~, hist, ~)
            Zphys = hist.u;  kIndex = (1:size(hist.u, 1)).';
        end
        function s = describe(~),                 s = "SurrogateDependentRep";  end
    end

    methods (Static)
        function e = makeEntry(name, kind, channel, lag)
            %MAKEENTRY  Public reach-through to the inherited protected factory.
            %
            %   Two callers, one spelling. This double builds its own entry with
            %   it, and the coordinateInfo contract suite uses it to show that an
            %   unknown kind cannot be CONSTRUCTED at all -- the closed set is
            %   enforced where entries are made, not merely where they are read
            %  .
            %
            %   It moved here when the shape double was retired. Writing the
            %   four fields out by hand instead would make this file a second
            %   producer of the entry format, which is the drift the coordinate-meaning rule exists to
            %   prevent one level up.
            e = rid.rep.Representation.coordinateEntry(name, kind, channel, lag);
        end
    end

end
