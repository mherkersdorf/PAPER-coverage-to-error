classdef ConstantFieldEstimator
%CONSTANTFIELDESTIMATOR  A field estimator returning one fixed value everywhere.
%
%   Distinguishable from rid.fields.Nominal, which is the whole point: a test that
%   scheduled Nominal estimators could not tell "the scheduler overrode the member"
%   from "the scheduler left it alone", because both give ones.
%
%   Inputs:
%       value  (1x1)  the constant returned for every supporting point; > 0

    properties (SetAccess = immutable)
        value        (1,1) double
        calibration  (1,1) string
        isOptimistic (1,1) logical
    end

    methods

        function obj = ConstantFieldEstimator(value)
            arguments
                value (1,1) double {mustBePositive}
            end
            obj.value        = value;
            obj.calibration  = "constant test estimator; 1 means nominal";
            obj.isOptimistic = false;
        end

        function [out, info] = estimate(obj, ~, ~, ~, ~, Q)
            out  = obj.value * ones(size(Q, 1), 1);
            info = struct("calibration", obj.calibration, ...
                          "isOptimistic", obj.isOptimistic, "nSup", size(Q, 1));
        end

        function s = describe(obj)
            s = sprintf("ConstantFieldEstimator(%g)", obj.value);
        end

    end

end
