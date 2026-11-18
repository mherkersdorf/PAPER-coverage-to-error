function K = gpMaternKernel(A, B, lengthScale)
%GPMATERNKERNEL  Isotropic Matern 5/2 kernel, unit variance: (nA x d), (nB x d) -> (nA x nB).
%   Squared distances by expansion (base MATLAB), clipped at 0 against
%   round-off; the kernel is flat at r = 0, so that round-off is harmless.
D2 = sum(A.^2, 2) + sum(B.^2, 2).' - 2 * (A * B.');
s  = sqrt(5 * max(D2, 0)) / lengthScale;
K  = (1 + s + s.^2 / 3) .* exp(-s);
end
