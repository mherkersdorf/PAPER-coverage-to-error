classdef FailingFieldEstimator
%FAILINGFIELDESTIMATOR  A field estimator that refuses, for the path.
%
%   An estimator that cannot honour its contract on the data so far must fail
%   loudly rather than clamp. The scheduler catches
%   that, keeps the previous fields and logs the refusal; this double is what makes
%   that path reachable in a test.

    properties (SetAccess = immutable)
        calibration  (1,1) string = "never produces a value"
        isOptimistic (1,1) logical = false
    end

    methods

        function [out, info] = estimate(~, ~, ~, ~, ~, ~) %#ok<STOUT>
            error("rid:fields:estimatorRefused", ...
                "not enough residuals in this region to estimate anything.");
        end

        function s = describe(~)
            s = "FailingFieldEstimator";
        end

    end

end
