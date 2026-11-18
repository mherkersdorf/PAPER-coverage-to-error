function [srg, info] = p3BoucWenBla(jsonPath)
%P3BOUCWENBLA  The P3 surrogate -- the order-3 BLA of the Bouc-Wen process (750 Hz interpolation path), from bla.json.
%
%   Rebuilds rid.surrogate.LinearStateSpace(Ad, Bd, Cd, tS, "ed", ed, "fd", fd,
%   "directTerm", directTerm) from experiments/p3/surrogates/bla.json, the
%   committed file.
%   The numbers carry 17 significant digits, so the object is bit-identical to
%   the fitted one. Every P3 script that needs the BLA, or the APRBS hold
%   derived from it, takes it from here -- the counterpart of p2DuffingS1.
%
%   The BLA is the surrogate of the 750 Hz path of p3BoucWenSimulate (750 Hz
%   samples interpolated band-limited to 15 kHz, interp(u, 20); the test data's
%   path), u [N] -> y [m], estimated by the robust method on own
%   random-phase multisines, never on the test data.
%   Its states are the subspace estimate's coordinates, NOT the process's: the
%   surrogate suits a representation built from u and y (the NARX spaces of
%   the study), not rid.rep.StateSpace.
%
%   Usage (with experiments/processes and
%   rhcInputDesigner/setupPath on the path):
%
%       [bla, info] = p3BoucWenBla();
%       y = bla.simulate(u, zeros(1, bla.nStates));   % u (N x 1) [N], from rest
%       T = info.aprbsT;                               % APRBS.m's T [samples]
%
%   Inputs:
%       jsonPath  (1x1) string  optional; the file to read. Default (and the
%                 only use outside the tests): experiments/p3/surrogates/bla.json
%   Outputs:
%       srg   (1x1) rid.surrogate.LinearStateSpace
%       info  (1x1) struct
%           .tauDominantSamples  (1x1) -1/ln|z| of the pole of largest |z| [samples at 750 Hz]
%           .aprbsT              (1x1) tauDominantSamples / 2, APRBS.m's T
%           .provenance          (1x1) struct  the file's provenance block

arguments
    jsonPath (1,1) string = fullfile(fileparts(fileparts(mfilename("fullpath"))), "p3", "surrogates", "bla.json")
end
if ~isfile(jsonPath)
    error("experiments:p3BoucWenBla:missingFile", ...
        "%s does not exist; it is part of the repository.", jsonPath);
end
s = jsondecode(fileread(jsonPath));

srg = rid.surrogate.LinearStateSpace(s.Ad, s.Bd, s.Cd, s.tS, ...
    "ed", s.ed, "fd", s.fd, "directTerm", s.directTerm);

% The time constants are derived from Ad, so recomputing them checks the file.
% eig differs in the last bits across platforms and releases (cluster job
% 12330325, x86 R2024b, refused the file written on Apple Silicon R2025b), so the
% check is relative to 1e-12 -- a last-bit change of Ad moves tau by ~1e-14, a
% change of Ad itself by far more -- and the STORED values are returned: every
% machine takes the same APRBS hold from the same file. srg.timeConstant is the
% object's own, recomputed here; it is guidance only and enters no result.
tauDominant = -1 / log(max(abs(eig(srg.Ad))));
checkDerived("timeConstant", srg.timeConstant, s.timeConstant);
checkDerived("tauDominantSamples", tauDominant, s.tauDominantSamples);
checkDerived("aprbsT", tauDominant / 2, s.aprbsT);
info = struct("tauDominantSamples", s.tauDominantSamples, "aprbsT", s.aprbsT, ...
    "provenance", s.provenance);
end

% ---------------------------------------------------------------------------

function checkDerived(name, recomputed, stored)
%CHECKDERIVED  The value derived from Ad on this machine agrees with bla.json's to a relative 1e-12.
relTol = 1e-12;
if ~(abs(recomputed - stored) <= relTol * abs(stored))
    error("experiments:p3BoucWenBla:inconsistentFile", ...
        "bla.json's %s does not follow from its Ad: stored %.17g, recomputed %.17g " + ...
        "(relative difference %.3g, tolerance %g).", name, stored, recomputed, ...
        abs(recomputed - stored) / abs(stored), relTol);
end
end
