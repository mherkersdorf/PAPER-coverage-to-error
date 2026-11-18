function Znorm = realizedPoints(v)
%REALIZEDPOINTS  The run's realized regressor points, rebuilt from the record.
%
%   `Z` IS DERIVED AND NEVER STORED, so recovering it means calling
%   `rep.build(hist, srg)`. The question was always WHICH surrogate: a `Setup`
%   is immutable, so `setup.srg` is the model the run BEGAN with, and after a
%   refit the loop remaps every past coordinate under the new one (invariant 3).
%   Rebuilding with the initial model would then show a distribution the run
%   never visited -- the worse failure, because it is a plausible one.
%
%   THE RECORD NOW CARRIES THE FINAL SURROGATE, in `outputs`, which is
%   exactly sufficient: invariant 3 makes the remap total, so at the end of a
%   run every realized coordinate is a function of the history and the final
%   model. `viewOf` reads it from there, so this function has no choice left to
%   make and no case left to withhold.
%
%   Inputs:
%       v  (1x1) struct  a viewOf() view; `v.srg` is the run's FINAL surrogate
%   Outputs:
%       Znorm  (nRealized x nDim)  normalized coordinates

[Zphys, ~] = v.rep.build(v.hist, v.srg);
Znorm      = v.nrm.toNormalized(Zphys);

end
