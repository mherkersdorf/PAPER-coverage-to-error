function [srg, provenance] = p2DuffingS1()
%P2DUFFINGS1  S1 -- the linear state-space surrogate of the Duffing process, from S1.json.
%
%   Rebuilds rid.surrogate.LinearStateSpace(Ad, Bd, Cd, tS, "ed", ed, "fd", fd,
%   "directTerm", directTerm) from experiments/p2/surrogates/S1.json, the
%   committed file. The numbers carry
%   17 significant digits, so the object is bit-identical to the fitted one.
%   Every script that needs S1 builds it here, so no two callers can assemble it
%   differently -- the counterpart of p2DuffingS0.
%
%   S1 is the equation-error least-squares fit on the measured states of the
%   noise-free pilot APRBS; its states are the process's coordinates [x1, x2].
%
%   Usage (with experiments/processes and
%   rhcInputDesigner/setupPath on the path):
%
%       S1 = p2DuffingS1();
%       xSs = S1.steadyState(0.5);                       % rest state at u = 0.5
%       [y, xEnd, x] = S1.simulateStates(u, [0 0]);      % u (N x 1), from rest
%
%   Inputs:  none
%   Outputs:
%       srg         (1x1) rid.surrogate.LinearStateSpace
%       provenance  (1x1) struct  the file's provenance block: pilotRecordId,
%                                 fitRecordId, seed, pilotNoiseLevel, fitMethod,
%                                 script, commit, date

jsonPath = fullfile(fileparts(fileparts(mfilename("fullpath"))), "p2", "surrogates", "S1.json");
if ~isfile(jsonPath)
    error("experiments:p2DuffingS1:missingFile", ...
        "%s does not exist; it is part of the repository.", jsonPath);
end
s = jsondecode(fileread(jsonPath));

srg = rid.surrogate.LinearStateSpace(s.Ad, s.Bd, s.Cd, s.tS, ...
    "ed", s.ed, "fd", s.fd, "directTerm", s.directTerm);

% timeConstant is derived from Ad, so the stored value is a consistency check
% of the file, not an input.
if ~isequal(srg.timeConstant, s.timeConstant)
    error("experiments:p2DuffingS1:inconsistentFile", ...
        "S1.json states timeConstant = %.17g, but its Ad gives %.17g.", ...
        s.timeConstant, srg.timeConstant);
end
provenance = s.provenance;
end
