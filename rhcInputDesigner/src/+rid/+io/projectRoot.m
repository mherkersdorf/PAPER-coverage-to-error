function root = projectRoot()
%PROJECTROOT  The folder holding src/, tests/, docs/ and runs/.
%
%   ONE DEFINITION, because it was three: `saveRun` and `loadRun` each carried a
%   local copy and `rid.rhc.run` computed it inline for the git capture. Three
%   copies of one path is three chances for them to disagree about where the
%   project is -- and the disagreement shows up as a run that writes a record
%   somewhere `loadRun` does not look, which reads as a missing record rather than
%   as a path bug.
%
%   It is derived from THIS FILE'S location and never from `pwd`: the working
%   directory is global state, and an experiment must not depend on where it was
%   launched from.
%
%   THE INVARIANT IS CHECKED, NOT MERELY DOCUMENTED. The derivation is a fixed
%   number of `fileparts` calls, so it is silently wrong the moment this file
%   moves -- one package deeper and it returns `src/`, which exists, is readable,
%   and is not the project root. Nothing downstream would notice: `runsRoot` would
%   happily resolve `src/runs`, `saveRun` would create it, and the failure would
%   surface much later as a record reported missing immediately after it was
%   written. An error at the point of computation costs one `isfolder` call per
%   invocation and names the actual defect.
%
%   `src` AND `tests` ARE THE MARKERS, not `docs` or `runs`. Both are required for
%   the project to be usable at all and neither is optional or generated: `runs/`
%   does not exist until the first record is written, and a checkout without
%   `docs/` is imaginable in a way that one without `src/` is not. Two markers
%   rather than one, because a single common name could be matched by an
%   unrelated ancestor directory.
%
%   Inputs:  none
%   Outputs:
%       root  (1x1) string  absolute path to the project root

% src/+rid/+io/projectRoot.m -> +io -> +rid -> src -> the project root.
root = string(fileparts(fileparts(fileparts(fileparts(mfilename("fullpath"))))));

markers = ["src", "tests"];
missing = markers(~arrayfun(@(name) isfolder(fullfile(root, name)), markers));
if ~isempty(missing)
    error("rid:io:invalidProjectRoot", ...
        "Derived the project root as ""%s"", which has no %s. The root is " + ...
        "obtained by walking up from this file (src/+rid/+io/projectRoot.m), " + ...
        "so a fixed number of steps is wrong the moment the file moves -- and a " + ...
        "plausible wrong root surfaces later as a record reported missing " + ...
        "rather than as a path error. Check where this file lives.", ...
        root, strjoin(missing + "/", " or "));
end

end
