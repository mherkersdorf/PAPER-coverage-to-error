function h = contentHash(Q)
%CONTENTHASH  A SHA-256 hash of Q's bytes, for identity by CONTENT, not path.
%
%   Paths move; a run record names a stored supporting-point set by this hash
%   instead, for the same reason it stores a git SHA rather than a filename.
%   Shared by rid.roi.saveSupportingSet and the "stored" branch of
%   rid.roi.supportingPoints, so the identity a file is saved under and the one
%   a caller reads back cannot drift apart (the one-producer precedent).
%
%   Inputs:
%       Q  (nSup x nDim) double
%   Outputs:
%       h  (1x1) string  lowercase hex SHA-256 of Q's raw IEEE-754 bytes,
%                        column-major, or "unknown" if the digest is unavailable

arguments
    Q (:,:) double {mustBeReal, mustBeFinite}
end

try
    digest = java.security.MessageDigest.getInstance("SHA-256");
    h = string(lower(reshape(dec2hex(typecast( ...
        digest.digest(typecast(Q(:), "uint8")), "uint8")).', 1, [])));
catch
    h = "unknown";
end

end
