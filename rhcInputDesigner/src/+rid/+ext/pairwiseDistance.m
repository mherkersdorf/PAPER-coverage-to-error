function D = pairwiseDistance(A, B)
%PAIRWISEDISTANCE  Euclidean distances between every pair of points of two sets.
%
%   Toolbox: Statistics and Machine Learning Toolbox (pdist2).
%
%   Isolates the project's pairwise-distance functionality. The criterion's
%   nearest-neighbour distance is this matrix reduced along its second dimension,
%   exactly as written:
%
%       dRaw = min(rid.ext.pairwiseDistance(Q, Znorm), [], 2);
%
%   The reduction stays with the caller: it is the criterion's logic, not the
%   toolbox's. Distances are Euclidean in the normalized coordinates of
%   the method description. The anisotropic extension
%   is deliberately not a parameter here -- a per-point metric changes which
%   neighbour is nearest, and therefore changes the caching contract of
%   the method description, not merely this call.
%
%   Empty point sets are legal and return the correspondingly empty matrix:
%   the acquired set is empty before the first sample is committed.
%
%   Inputs:
%       A  (nA x nDim)  point set, samples in rows
%       B  (nB x nDim)  point set, samples in rows
%
%   Outputs:
%       D  (nA x nB)    D(ii,jj) = norm(A(ii,:) - B(jj,:), 2)

arguments
    A (:,:) double {mustBeReal, mustBeFinite}
    B (:,:) double {mustBeReal, mustBeFinite}
end

if size(A, 2) ~= size(B, 2)
    error("rid:ext:dimensionMismatch", ...
        "A and B must have the same number of columns (got %d and %d).", ...
        size(A, 2), size(B, 2));
end

D = pdist2(A, B);

end
