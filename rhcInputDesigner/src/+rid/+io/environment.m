function env = environment()
%ENVIRONMENT  What this run was executed on, for the run record.
%
%   A run that behaves differently after a toolbox upgrade must be identifiable as
%   such, which is only possible if the toolbox
%   versions were recorded at the time. `ver` is the source; nothing here is
%   inferred.
%
%   IT IS CAPTURED AT RUN START, NOT AT SAVE TIME. The record describes the
%   run, not the moment it finished, and the two differ whenever anything changes
%   under a long experiment.
%
%   IT NEVER THROWS. An experiment must not die because a version string could not
%   be read, so every field degrades to "unknown" rather than raising.
%
%   Inputs:  none
%   Outputs:
%       env  (1x1) struct
%           .matlabRelease  (1x1) string   e.g. "R2025b"
%           .matlabVersion  (1x1) string
%           .host           (1x1) string
%           .platform       (1x1) string   computer("arch")
%           .toolboxes      (1xn) struct   .Name .Version .Release, from `ver`
%           .capturedAt     (1x1) string   ISO-8601, local time with offset

env = struct( ...
    "matlabRelease", "unknown", ...
    "matlabVersion", "unknown", ...
    "host",          "unknown", ...
    "platform",      "unknown", ...
    "toolboxes",     struct("Name", {}, "Version", {}, "Release", {}), ...
    "capturedAt",    string(datetime("now", "TimeZone", "local", ...
                            "Format", "uuuu-MM-dd'T'HH:mm:ssZZZZZ")));

try
    env.matlabRelease = string(version("-release"));
    env.matlabVersion = string(version());
catch
    % leave the placeholders
end

try
    env.platform = string(computer("arch"));
catch
end

try
    [failed, name] = system("hostname");
    if failed == 0
        env.host = strtrim(string(name));
    end
catch
end

try
    installed = ver();
    keep      = struct("Name", {}, "Version", {}, "Release", {});
    for iTool = 1:numel(installed)
        keep(end + 1) = struct("Name",    string(installed(iTool).Name), ...
                               "Version", string(installed(iTool).Version), ...
                               "Release", string(installed(iTool).Release)); %#ok<AGROW>
    end
    env.toolboxes = keep;
catch
end

end
