classdef ScriptedStream < handle
%SCRIPTEDSTREAM  Test double: `rand` returns a scripted sequence, cyclically.
%
%   Lets a test pin the two draws `targetSequence` makes -- the target and the
%   speed fraction `g` -- without `initialGuesses` gaining a knob it would not
%   otherwise have. `g = 1 - rand`, so a scripted 0 forces g = 1.
%
%   The draw order per target is: target, then g. A test that wants g = 1 with
%   varying targets therefore supplies [t1 0 t2 0 ...].

    properties (SetAccess = immutable)
        values (1,:) double
    end
    properties (Access = private)
        nDrawn (1,1) double = 0
    end

    methods
        function obj = ScriptedStream(values)
            arguments
                values (1,:) double {mustBeGreaterThanOrEqual(values, 0), ...
                                     mustBeLessThan(values, 1)}
            end
            obj.values = values;
        end

        function v = rand(obj, varargin)
            obj.nDrawn = obj.nDrawn + 1;
            v = obj.values(1 + mod(obj.nDrawn - 1, numel(obj.values)));
        end
    end
end
