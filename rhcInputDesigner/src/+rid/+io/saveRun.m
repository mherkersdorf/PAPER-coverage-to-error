function path = saveRun(rec, note, opts)
%SAVERUN  Write a run record to disk, immutably.
%
%       runs/<runId>/run.mat        the record, for MATLAB
%       runs/<runId>/run.json       options, seeds, git state, status, headline
%                                   metrics -- searchable without starting MATLAB
%       runs/<runId>/dirty.patch    only when the tree was dirty
%
%   RECORDS ARE IMMUTABLE. An existing
%   `runs/<runId>/` is never overwritten; re-running produces a new record. A
%   collision therefore raises rather than silently replacing a result that a
%   figure somewhere already cites.
%
%   TWO FORMATS, and the JSON is not a convenience. `run.mat` is git-ignored
%   repository-wide, so a record that existed only as `.mat` would exist only on
%   the machine that made it; the JSON keeps the run searchable and diffable. Both
%   describe the same run, and the `.mat` is the complete one.
%
%   `run.mat` IS `-v7`, AND `-v7.3` ONLY WHEN IT MUST BE. `-v7.3` (HDF5)
%   pays per struct field, and a record holds a struct per step: a 500-step record
%   was 12.8 MB in `-v7.3` and 0.35 MB in `-v7`. `-v7` cannot hold a variable of
%   2 GB or more, and `save` does not fail on one: it warns and writes a file
%   WITHOUT the variable. So the file is checked for `rec` after writing, and a
%   record that is not in it is written again as `-v7.3`. The check reads the file,
%   not the warning, which a caller may have switched off.
%
%   THE DIFF IS WRITTEN OUT SEPARATELY AS WELL AS EMBEDDED. It is in the record so
%   the record stays self-contained; it is on disk as `dirty.patch` so it can be
%   applied with `git apply` without a MATLAB session.
%
%   Inputs:
%       rec        (1x1) struct  the record built by rid.rhc.run
%       note       (1x1) string  optional free text stored in the record
%       opts.root  (1x1) string  the runs directory; default from rid.defaultOptions
%   Outputs:
%       path       (1x1) string  the directory written

arguments
    rec       (1,1) struct
    note      (1,1) string = ""
    opts.root (1,1) string = ""
end

if ~isfield(rec, "identity") || ~isfield(rec.identity, "runId")
    error("rid:io:invalidRecord", ...
        "rec has no identity.runId. A record is built by rid.rhc.run and " + ...
        "carries its own identity.");
end

rec.note = note;

root = opts.root;
if strlength(root) == 0
    % No root given means this was called on its own, with no run to take an
    % effective `opts.io.runsDir` from. `rid.rhc.run` passes the effective value,
    %   and the resolution lives in one place so
    % that write and read cannot disagree.
    root = rid.io.runsRoot();
end
path = string(fullfile(root, rec.identity.runId));

if isfolder(path)
    error("rid:io:recordExists", ...
        "%s already exists. Records are immutable: re-running produces a new " + ...
        "record and never overwrites one a figure may already cite " + ...
        ".", path);
end
mkdir(path);

writeMat(fullfile(path, "run.mat"), rec);

summary = struct( ...
    "runId",        rec.identity.runId, ...
    "createdAt",    rec.identity.createdAt, ...
    "host",         rec.identity.host, ...
    "matlabRelease", rec.identity.matlabRelease, ...
    "commit",       rec.code.commit, ...
    "branch",       rec.code.branch, ...
    "isDirty",      rec.code.isDirty, ...
    "srcHash",      rec.code.srcHash, ...
    "options",      rec.configuration.effective, ...
    "overrides",    rec.configuration.overrides, ...
    "seeds",        rec.seeds, ...
    "status",       rec.outputs.status, ...
    "nSteps",       rec.outputs.nSteps, ...
    "J",            rec.outputs.J, ...
    "fillDistanceQ", rec.problem.fillDistanceQ, ...
    "nSup",         rec.problem.nSup, ...
    "storedQ",      rec.problem.storedQ, ...
    "note",         rec.note);
writeText(fullfile(path, "run.json"), jsonencode(summary, "PrettyPrint", true));

if rec.code.isDirty && strlength(rec.code.diff) > 0
    writeText(fullfile(path, "dirty.patch"), rec.code.diff);
end

end

% ---------------------------------------------------------------------------

function writeMat(file, rec) %#ok<INUSD> saved by name
%WRITEMAT  `rec` as `-v7`, or as `-v7.3` when `-v7` could not hold it.
%
%   `save -v7` of a variable of 2 GB or more leaves a file with no variable in it
%   and only a warning (`MATLAB:save:sizeTooBigForMATFile`). The fallback is
%   decided on what the file holds, so it acts whether or not the warning is shown.
save(file, "rec", "-v7");
if ~holdsRecord(file)
    delete(file);
    save(file, "rec", "-v7.3");
    if ~holdsRecord(file)
        error("rid:io:writeFailed", ...
            "%s does not hold the record after saving as -v7 and as -v7.3.", file);
    end
end

end

% ---------------------------------------------------------------------------

function tf = holdsRecord(file)
%HOLDSRECORD  Does the MAT-file contain the variable `rec`?
tf = isfile(file) && ismember("rec", string(who("-file", file)));

end

% ---------------------------------------------------------------------------

function writeText(file, text)
%WRITETEXT  Write a string to a file, or say which file could not be written.
fid = fopen(file, "w");
if fid < 0
    error("rid:io:writeFailed", "Could not open %s for writing.", file);
end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, char(text));

end
