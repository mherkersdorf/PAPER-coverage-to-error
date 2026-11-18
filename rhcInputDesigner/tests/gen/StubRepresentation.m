classdef StubRepresentation < rid.rep.Representation
%STUBREPRESENTATION  Test double: coordinates declared entry by entry.
%
%   Exists so that +gen can be tested against the Representation INTERFACE rather
%   than against rid.rep.Narx. Two things need it:
%
%   - configurations Narx cannot produce -- a channel with a rate but no level, a
%     representation with no input coordinate at all, a "state" coordinate (which
%     arrives with the StateSpace);
%   - names that follow no convention, which is what proves the classification is a
%     switch on `.kind` and not a name match.
%
%   Constructed from {name, kind, channel, lag} cells, one per coordinate:
%
%       rep = StubRepresentation({"u1_k", "inputLevel", 1, 0}, ...
%                                {"y_km1", "output",  NaN, 1})

    properties (SetAccess = immutable)
        nDim
        dependsOnSurrogate = false
        coordInfo (1,:) struct
    end

    methods
        function obj = StubRepresentation(varargin)
            obj@rid.rep.Representation("measured");
            info = repmat(struct("name", "", "kind", "", ...
                                 "channel", NaN, "lag", NaN), 1, nargin);
            for iDim = 1:nargin
                entry = varargin{iDim};
                info(iDim) = struct("name", string(entry{1}), ...
                                    "kind", string(entry{2}), ...
                                    "channel", entry{3}, "lag", entry{4});
            end
            obj.coordInfo = info;
            obj.nDim      = nargin;
        end

        function info = coordinateInfo(obj),  info = obj.coordInfo;  end

        function [box, dimName] = bounds(obj, roi)
            dimName = obj.requiredNames();
            box     = obj.lookupRanges(roi, dimName);
        end

        function [Zphys, kIndex] = build(~, hist, ~)
            Zphys = hist.u;  kIndex = (1:size(hist.u, 1)).';
        end

        function s = describe(~),  s = "StubRepresentation";  end
    end

end
