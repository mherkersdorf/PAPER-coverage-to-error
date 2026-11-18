function rec = loadRun(runId, opts)
%LOADRUN  Read a run record back from disk.
%
%   The `.mat` is the complete record and the one this reads; `run.json` is the
%   searchable summary and is deliberately NOT the source of truth, because it
%   carries only the headline fields.
%
%   IT MERGES NOTHING. The record stores the EFFECTIVE options struct precisely so
%   that a later change of default cannot alter what an old run means; loading
%   defaults and overlaying the stored overrides would reintroduce exactly that
%   dependency. What is stored is what is returned.
%
%   THE ROOT IS RESOLVED BY `rid.io.runsRoot`, the same function the writer used.
%   It was computed here independently, and with an absolute `opts.io.runsDir` the
%   two disagreed: the run wrote to the right directory and this looked under
%   `<projectRoot>/<abs/path>`, reporting a record it had just written as missing.
%
%   Inputs:
%       runId      (1x1) string  the record's directory name, or a full path
%       opts.root  (1x1) string  the runs directory; default from rid.defaultOptions
%   Outputs:
%       rec        (1x1) struct  the record as it was written

arguments
    runId     (1,1) string
    opts.root (1,1) string = ""
end

path = runId;
if ~isfolder(path)
    root = opts.root;
    if strlength(root) == 0
        root = rid.io.runsRoot();
    end
    path = string(fullfile(root, runId));
end

file = fullfile(path, "run.mat");
if ~isfile(file)
    error("rid:io:recordNotFound", ...
        "No run.mat at %s. A record is a DIRECTORY holding run.mat, run.json " + ...
        "and, when the tree was dirty, dirty.patch (written by " + ...
        "rid.io.saveRun).", path);
end

loaded = load(file, "rec");
rec    = loaded.rec;

end
