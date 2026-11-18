function lme = calculate_lme(y, yLinear, opts)
%CALCULATE_LME  Normalized linear model error.
%
%   DEFAULT (opts.offset = false), UNCHANGED from the original formula, so
%   existing numbers stay reproducible:
%
%       LME = 100 * sqrt(MSE(y, yLinear) / var(y))
%
%   opts.offset = true replaces yLinear with an AFFINE reference: yLinear fitted
%   to y by least squares WITH AN INTERCEPT, y ~ a + b*yLinear. This is never a
%   blind removal of means (a plain y - mean(y), yLinear - mean(yLinear)) -- in an
%   output comparison a constant offset can be the static nonlinearity itself
%   (e.g. a Hammerstein nonlinearity that shifts the mean), and subtracting it
%   unconditionally would remove exactly what the metric is meant to detect. An
%   explicit least-squares affine fit is the alternative: it is chosen to
%   minimize the residual the metric then reports, not to erase a nuisance mean.
%   var(y) is unaffected by opts.offset; only the reference is refitted.
%
%   Inputs:
%       y           (N x 1) or (N x d)  process output (or one-step increment)
%       yLinear     (N x 1) or (N x d)  linear reference, same shape as y
%       opts.offset (1x1) logical  default false; see above
%   Output:
%       lme  (1x1) or (1 x d)  LME [%], one value per column of y
arguments
    y       (:,:) double {mustBeReal, mustBeFinite}
    yLinear (:,:) double {mustBeReal, mustBeFinite}
    opts.offset (1,1) logical = false
end

if ~isequal(size(y), size(yLinear))
    error("experiments:calculateLme:dimensionMismatch", ...
        "y and yLinear must have the same shape; got %s and %s.", ...
        mat2str(size(y)), mat2str(size(yLinear)));
end

if opts.offset
    N = size(y, 1);
    for jCol = 1:size(y, 2)
        Phi = [yLinear(:, jCol), ones(N, 1)];
        theta = Phi \ y(:, jCol);
        yLinear(:, jCol) = Phi * theta;
    end
end

mse      = mean((y - yLinear).^2, 1);
variance = var(y, 1, 1);

if any(variance == 0)
    error("experiments:calculateLme:zeroVariance", ...
        "Variance of y is zero in at least one column. LME is undefined.");
end

lme = 100 * sqrt(mse ./ variance);
end
