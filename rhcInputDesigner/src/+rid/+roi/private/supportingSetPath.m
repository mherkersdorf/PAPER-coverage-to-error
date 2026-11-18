function path = supportingSetPath(name)
%SUPPORTINGSETPATH  Where a stored supporting-point set lives.
%
%   ONE PRODUCER, shared by rid.roi.saveSupportingSet and the "stored" branch of
%   rid.roi.supportingPoints: a convention with two producers drifts (the lag-suffix naming's
%   precedent). `name` is the basename of the file, not a path -- the location is
%   fixed at `supportingSets/` under the project root, resolved through
%   rid.io.projectRoot (the project-root rule, no second path definition).
%
%   Inputs:
%       name  (1x1) string  the file's basename, e.g. "narx2d_lhsFull_5000.mat"
%   Outputs:
%       path  (1x1) string  absolute path under <projectRoot>/supportingSets/

arguments
    name (1,1) string
end

path = string(fullfile(rid.io.projectRoot(), "supportingSets", name));

end
