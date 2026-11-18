function model = gpMaternFit(Z, T, lengthScale, jitterRel)
%GPMATERNFIT  Fit the reference GP: affine mean by least squares, Matern 5/2 on the residuals.
%
%   The downstream reference model: a Gaussian process with FIXED,
%   deliberately local hyperparameters. NOTHING IS OPTIMIZED -- the length
%   scale is an input, fixed before any data and identical across data sets,
%   because a model whose hyperparameters adapt to the data spacing would adapt
%   to exactly the quantity under test.
%
%       t(z) = [z, 1] * beta  +  GP(0, k)        k(r) = (1 + s + s^2/3) exp(-s),
%                                                s = sqrt(5) r / lengthScale
%
%   beta is the least-squares affine fit to T; the GP interpolates its
%   residuals with an isotropic Matern 5/2 kernel of unit signal variance. With
%   noise-free data the predicted mean does not depend on the signal variance.
%
%   AN INPUT WITHOUT SPREAD in Z (e.g. an APRBS holding its first level over the
%   whole training set) gets slope 0 in beta, exactly: its column is left out of
%   the least-squares fit, which would otherwise be rank-deficient and split the
%   offset between it and the constant arbitrarily. Its distance still enters the
%   kernel. With every input varying, the fit is unchanged. The local model
%   network treats such an input the same way.
%
%   JITTER. K + jitter*I (added on the diagonal in place) with jitter = jitterRel (relative to the unit signal
%   variance). If the Cholesky factorization fails, jitter is raised by a
%   factor of 10 and the factorization retried; every raise is returned in
%   model.jitterRaises. Above a relative jitter of 1 the model would no longer
%   interpolate in any sense, so that is an error, not a further raise.
%
%   Z is expected in the caller's normalized coordinates (the gate uses the
%   ROI mapped to [-1, 1]^d); this function does not normalize.
%
%   Inputs:
%       Z            (N x d) double  training inputs, samples in rows
%       T            (N x m) double  targets, one column per output; the columns
%                                    share the kernel and its factorization
%       lengthScale  (1x1) positive  Matern length scale, in Z's units
%       jitterRel    (1x1) positive  initial relative jitter, e.g. 1e-8
%   Outputs:
%       model  (1x1) struct
%           .Z             (N x d)      training inputs
%           .beta          ((d+1) x m)  affine coefficients, [z, 1] * beta; the
%                                       rows of inputs without spread are 0
%           .alpha         (N x m)      (K + jitter I) \ residuals
%           .lengthScale   (1x1)
%           .jitter        (1x1)        the jitter finally used
%           .jitterRaises  (1x1)        number of factor-10 raises (0 = none)
arguments
    Z           (:,:) double {mustBeReal, mustBeFinite}
    T           (:,:) double {mustBeReal, mustBeFinite}
    lengthScale (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
    jitterRel   (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
end

[N, d] = size(Z);
if size(T, 1) ~= N
    error("experiments:gpMaternFit:dimensionMismatch", ...
        "Z has %d rows but T has %d.", N, size(T, 1));
end
if N < d + 1
    error("experiments:gpMaternFit:tooFewPoints", ...
        "The affine mean needs at least d + 1 = %d points; got %d.", d + 1, N);
end

varies    = max(Z, [], 1) > min(Z, [], 1);         % inputs with spread get a slope
fitted    = [varies, true];                         % ... and the constant always
Phi       = [Z, ones(N, 1)];
beta      = zeros(d + 1, size(T, 2));
beta(fitted, :) = Phi(:, fitted) \ T;
residuals = T - Phi * beta;

K = gpMaternKernel(Z, Z, lengthScale);

% The jitter is added to K's diagonal in place, from the unjittered diagonal at
% every attempt: the same matrix, bit for bit, as K + jitter*eye(N), without the
% two further N x N arrays that expression builds (0.5 GB each at N = 8000).
diagIdx  = 1:(N + 1):N^2;
diagBare = K(diagIdx);
jitter       = jitterRel;
jitterRaises = 0;
while true
    K(diagIdx) = diagBare + jitter;
    [L, flag] = chol(K, "lower");
    if flag == 0
        break
    end
    jitter       = 10 * jitter;
    jitterRaises = jitterRaises + 1;
    if jitter > 1
        error("experiments:gpMaternFit:notPositiveDefinite", ...
            "Cholesky failed up to a relative jitter of 1 (N = %d, lengthScale = %g).", ...
            N, lengthScale);
    end
end

alpha = L.' \ (L \ residuals);

model = struct("Z", Z, "beta", beta, "alpha", alpha, "lengthScale", lengthScale, ...
    "jitter", jitter, "jitterRaises", jitterRaises);
end
