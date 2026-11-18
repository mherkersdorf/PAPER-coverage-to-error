classdef MisshapedProcess < handle
%MISSHAPEDPROCESS  A process that returns the wrong thing, for the guard.
%
%   `mode` selects which way it is wrong: a whole horizon instead of the applied
%   unit, or a non-finite reading. Both must raise `rid:rhc:invalidMeasurement`,
%   and the loop must NOT catch it -- what to do about a failing process is a
%   policy, and policies stay out until there is hardware to have opinions about.

    properties (SetAccess = immutable)
        mode (1,1) string
    end

    methods

        function obj = MisshapedProcess(mode)
            arguments
                mode (1,1) string {mustBeMember(mode, ["tooManyRows", "notFinite"])}
            end
            obj.mode = mode;
        end

        function yMeasured = apply(obj, uApplied)
            switch obj.mode
                case "tooManyRows"
                    yMeasured = zeros(size(uApplied, 1) + 3, 1);
                case "notFinite"
                    yMeasured = NaN(size(uApplied, 1), 1);
            end
        end

    end

end
