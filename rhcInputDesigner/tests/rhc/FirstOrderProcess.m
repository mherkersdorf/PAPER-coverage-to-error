classdef FirstOrderProcess < handle
%FIRSTORDERPROCESS  A synthetic plant satisfying the contract.
%
%   The applied unit in, its measurement out. It is a first-order lag with a
%   deliberately DIFFERENT gain and time constant from the surrogate the design
%   uses, so an online test exercises the real situation -- the model is a
%   geometric guide, not a predictive model -- rather than
%   the degenerate one where the plant and the surrogate agree exactly.
%
%   It is a handle because it carries the plant's state between calls, which is
%   what a real process does. Nothing in `src/` may depend on that: the loop calls
%   `apply` and nothing else.
%
%   Deterministic and seedless: the noise channel is out of scope until a
%   noise-level estimator to feed, and a random plant here would make every run
%   test irreproducible for no gain.

    properties (SetAccess = immutable)
        a    (1,1) double   % discrete pole
        b    (1,1) double   % input gain
    end

    properties (SetAccess = private)
        x    (1,1) double = 0
        nApplications (1,1) double = 0
    end

    methods

        function obj = FirstOrderProcess(timeConstant, tS, gain)
            arguments
                timeConstant (1,1) double {mustBePositive}
                tS           (1,1) double {mustBePositive}
                gain         (1,1) double
            end
            obj.a = exp(-tS / timeConstant);
            obj.b = (1 - obj.a) * gain;
        end

        function yMeasured = apply(obj, uApplied)
            %APPLY  (nApplied x 1) in, (nApplied x 1) out.
            nApplied  = size(uApplied, 1);
            yMeasured = zeros(nApplied, 1);
            for iRow = 1:nApplied
                obj.x            = obj.a * obj.x + obj.b * uApplied(iRow, 1);
                yMeasured(iRow)  = obj.x;
            end
            obj.nApplications = obj.nApplications + 1;
        end

    end

end
