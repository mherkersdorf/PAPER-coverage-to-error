function [opts, overrides] = resolveOptions(overrides)
%RESOLVEOPTIONS  Defaults plus the caller's overrides, with typos refused.
%
%   The run record stores the EFFECTIVE options struct and the overrides
%   separately: the resolved struct makes the record
%   self-contained while defaults keep changing, and the overrides preserve the
%   INTENT, which is what a reader needs later. This function produces both.
%
%   AN UNRECOGNIZED NAME IS AN ERROR, never a silently ignored field. A dropped
%   override does not raise anywhere else: `"psiMax"` written under `opts.gen`
%   instead of `opts.crit` simply does nothing, and the only symptom is a run that
%   is slightly worse than the one the user thinks they configured. That is the
%   objection the starts-provider rule raised against forwarding name-value pairs, arriving through
%   the options struct instead.
%
%   TWO NAMES MAY BE ADDED THAT NO DEFAULT DECLARES, and they are enumerated below
%   rather than allowed by a wildcard: `epsReq0` and `nSteps` are
%   commitments with no defensible default, so they exist only when a caller
%   supplies them. Any other new name is a typo by construction, because a real
%   option would have a default in `rid.defaultOptions`.
%
%   THREE USED TO BE HERE AND ARE NOT ANY MORE, all for one reason: they carried
%   DATA through a channel meant for values. The refinement
%   `fields` are `(nSup x 1)` vectors indexed by `Q` and belong in
%   `rid.rhc.Setup`; `u0` and `y0` were the experiment's initial condition and
%   belong in the `hist0` argument, where `rid.assertHistory` validates them.
%   `nSteps` stays, because a budget length genuinely is a setting.
%
%   Inputs:
%       overrides  (1x1) struct  nested like rid.defaultOptions, partial
%   Outputs:
%       opts       (1x1) struct  the effective options
%       overrides  (1x1) struct  the input, returned unchanged for the record

arguments
    overrides (1,1) struct
end

opts = rid.defaultOptions();

% Dotted paths that may be ADDED to the effective struct. Everything else must
% already exist, which is what turns a misspelling into an error.
addable = ["crit.epsReq0", "rhc.nSteps"];

sections = string(fieldnames(overrides)).';
for section = sections
    if ~isfield(opts, section)
        error("rid:rhc:unknownOption", ...
            "opts has no section ""%s"". The sections are %s.", ...
            section, strjoin("""" + string(fieldnames(opts)).' + """", ", "));
    end
    if ~isstruct(overrides.(section)) || ~isscalar(overrides.(section))
        error("rid:rhc:unknownOption", ...
            "opts.%s must be a scalar struct of settings; options are nested " + ...
            "one level, by package.", section);
    end

    names = string(fieldnames(overrides.(section))).';
    for name = names
        path = section + "." + name;
        if ~isfield(opts.(section), name) && ~ismember(path, addable)
            error("rid:rhc:unknownOption", ...
                "opts.%s is not an option. A real option has a default in " + ...
                "rid.defaultOptions; the only names that may be added are %s, " + ...
                "which are commitments with no defensible default. An " + ...
                "unrecognized name is refused here because a dropped override " + ...
                "raises nowhere else -- it merely produces a worse run.", ...
                path, strjoin("""" + addable + """", ", "));
        end
        opts.(section).(name) = overrides.(section).(name);
    end
end

end
