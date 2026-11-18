function root = runsRoot(runsDir)
%RUNSROOT  Where run records live: the one resolution of `opts.io.runsDir`.
%
%   The method description documents `runsDir` as **relative to the
%   project root**. That says how a relative value is RESOLVED; it is not a
%   prohibition on an absolute one, and an absolute path is resolved to itself.
%
%   WHY IT IS ONE FUNCTION. The resolution existed in three places -- `rid.rhc.run`
%   handled the absolute case, `rid.io.saveRun` and `rid.io.loadRun` each joined
%   the project root on unconditionally -- and the three disagreed. With an
%   absolute `runsDir` the run wrote to the right directory while `loadRun` looked
%   under `<projectRoot>/<abs/path>` and found nothing: **write and read
%   disagreed**, and the symptom was a record that had just been written reported
%   as missing. A path that three functions compute is a path three functions can
%   get differently right.
%
%   IT IS NOT A DEFAULT LOOKUP. The caller passes the value it is acting under --
%   `rid.rhc.run` passes the EFFECTIVE `opts.io.runsDir`, so an override reaches
%   the directory the run writes to and the record names where it actually went.
%   Omitting the argument falls back to `rid.defaultOptions`, which is right only
%   for `saveRun` and `loadRun` called on their own, with no run to take an
%   effective value from.
%
%   Inputs:
%       runsDir  (1x1) string  relative to the project root, or absolute;
%                              omitted means rid.defaultOptions().io.runsDir
%   Outputs:
%       root     (1x1) string  absolute path to the runs directory

arguments
    runsDir (1,1) string = rid.defaultOptions().io.runsDir
end

if isAbsolute(runsDir)
    root = runsDir;
else
    root = string(fullfile(rid.io.projectRoot(), runsDir));
end

end

% ---------------------------------------------------------------------------

function tf = isAbsolute(path)
%ISABSOLUTE  A leading separator, or a Windows drive letter.
%
%   Joining a project root onto an absolute path does not fail: it produces a
%   third directory that is neither, and writes the run there. That is why this is
%   a test and not a `fullfile` call.
tf = startsWith(path, filesep) || ~isempty(regexp(path, "^[A-Za-z]:", "once"));

end
