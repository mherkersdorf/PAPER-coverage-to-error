function T = gpMaternPredict(model, Zq, chunkSize)
%GPMATERNPREDICT  Posterior mean of the reference GP fitted by gpMaternFit, in chunks.
%
%   t(zq) = [zq, 1] * beta + k(zq, Z) * alpha. The query points are processed
%   chunkSize rows at a time, so memory stays at chunkSize x N kernel entries
%   whatever the number of query points.
%
%   Inputs:
%       model      (1x1) struct         from gpMaternFit
%       Zq         (nQ x d) double      query points, in the model's coordinates
%       chunkSize  (1x1) positive int   query rows per chunk
%   Outputs:
%       T  (nQ x m) double  predicted mean, one column per fitted output
arguments
    model     (1,1) struct
    Zq        (:,:) double {mustBeReal, mustBeFinite}
    chunkSize (1,1) double {mustBeInteger, mustBePositive}
end

if size(Zq, 2) ~= size(model.Z, 2)
    error("experiments:gpMaternPredict:dimensionMismatch", ...
        "Zq has %d columns; the model was fitted on %d.", size(Zq, 2), size(model.Z, 2));
end

nQ = size(Zq, 1);
T  = zeros(nQ, size(model.alpha, 2));
for iStart = 1:chunkSize:nQ
    rows = iStart:min(iStart + chunkSize - 1, nQ);
    Kq   = gpMaternKernel(Zq(rows, :), model.Z, model.lengthScale);
    T(rows, :) = [Zq(rows, :), ones(numel(rows), 1)] * model.beta + Kq * model.alpha;
end
end
