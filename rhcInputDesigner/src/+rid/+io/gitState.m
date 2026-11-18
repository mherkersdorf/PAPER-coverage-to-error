function code = gitState(projectRoot)
%GITSTATE  The commit a run was made from, and what differed from it.
%
%   A plot is reproducible if and only if its record says which code produced it.
%   This reads that from git and **degrades
%   gracefully**: a missing git installation never aborts an experiment, and a
%   dirty tree is allowed but recorded, with the diff embedded.
%
%   THE DIFF IS EMBEDDED, NOT REFERENCED. A path to a patch file elsewhere is not
%   provenance -- the file moves, the branch advances, and the record then names a
%   commit whose working tree nobody can reconstruct. The whole point of section 1
%   is that the record is self-contained.
%
%   THE FALLBACK IS A CONTENT HASH OF `src/`, not silence. "unknown" alone would
%   make two runs from different code indistinguishable; a hash of the sources at
%   least separates them, even when it cannot name them. It is computed from the
%   `.m` files' bytes, so it is stable across machines and ignores mtimes.
%
%   CAPTURED AT RUN START, so that `requireCleanTree` can refuse before an
%   experiment rather than after it, and so that the recorded state is the one the
%   run actually used.
%
%   Inputs:
%       projectRoot  (1x1) string  the folder to interrogate; the project root
%   Outputs:
%       code  (1x1) struct
%           .commit      (1x1) string   full hash, or "unknown"
%           .shortCommit (1x1) string   first 7 characters, or "unknown"
%           .branch      (1x1) string   or "unknown"
%           .isDirty     (1x1) logical  false when git is unavailable, with
%                                       .gitAvailable false to say so
%           .diff        (1x1) string   the full diff when dirty, "" otherwise
%           .srcHash     (1x1) string   content hash of src/**/*.m
%           .gitAvailable (1x1) logical

arguments
    projectRoot (1,1) string
end

code = struct("commit", "unknown", "shortCommit", "unknown", ...
              "branch", "unknown", "isDirty", false, "diff", "", ...
              "srcHash", contentHash(fullfile(projectRoot, "src")), ...
              "gitAvailable", false);

[failed, hash] = runGit(projectRoot, "rev-parse HEAD");
if failed
    return
end
code.gitAvailable = true;
code.commit       = strtrim(hash);
code.shortCommit  = extractBefore(code.commit + "        ", 8);

[failed, branch] = runGit(projectRoot, "rev-parse --abbrev-ref HEAD");
if ~failed
    code.branch = strtrim(branch);
end

% `git status --porcelain` covers staged, unstaged and untracked alike; `git diff`
% alone would call a tree with only staged changes clean.
[failed, status] = runGit(projectRoot, "status --porcelain");
if failed
    return
end
code.isDirty = strlength(strtrim(status)) > 0;
if ~code.isDirty
    return
end

% HEAD, not the index: the record must describe what differed from the COMMIT it
% names, which is what a reader will check out.
[failed, diffText] = runGit(projectRoot, "diff HEAD");
if ~failed
    code.diff = string(diffText);
end

end

% ---------------------------------------------------------------------------

function [failed, out] = runGit(projectRoot, args)
%RUNGIT  One git invocation, confined to the project folder, that cannot BLOCK.
%
%   `-C` rather than `cd`: the working directory is global state and an experiment
%   must not depend on where it was launched from.
%
%   `--no-pager` IS LOAD-BEARING, not hygiene. MATLAB's `system` hands the child a
%   terminal, so `git diff` starts its pager and waits for a keypress that will
%   never come -- and the symptom is a run that hangs with no output at all,
%   before the first step, on any dirty tree. `GIT_TERMINAL_PROMPT=0` closes the
%   same door on credential prompts.
cmd = "GIT_TERMINAL_PROMPT=0 git --no-pager " + ...
      sprintf("-C ""%s"" %s", projectRoot, args);
try
    [code, out] = system(cmd);
    failed = code ~= 0;
    out    = string(out);
catch
    failed = true;
    out    = "";
end

end

% ---------------------------------------------------------------------------

function h = contentHash(srcRoot)
%CONTENTHASH  A hash of every .m file under srcRoot, by CONTENT.
%
%   Sorted by relative path and built from the bytes, so it is identical on two
%   machines with the same sources and unaffected by checkout times.
h = "unknown";
try
    files = dir(fullfile(srcRoot, "**", "*.m"));
    if isempty(files)
        return
    end
    names  = string({files.folder}) + filesep + string({files.name});
    names  = sort(names);
    digest = java.security.MessageDigest.getInstance("SHA-256");
    for iFile = 1:numel(names)
        fid = fopen(names(iFile), "r");
        if fid < 0
            continue
        end
        bytes = fread(fid, Inf, "*uint8");
        fclose(fid);
        digest.update(uint8(char(extractAfter(names(iFile), strlength(srcRoot)))));
        digest.update(bytes);
    end
    h = string(lower(reshape(dec2hex(typecast(digest.digest(), "uint8")).', 1, [])));
catch
    h = "unknown";
end

end
