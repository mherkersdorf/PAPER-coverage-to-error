function report = reachabilityReport(box, dimName, constraints, opts)
%REACHABILITYREPORT  What the declarations imply, reported and never corrected.
%
%   Because nothing in this design is derived from anything else, a user can
%   declare a set of bounds that is internally inconsistent, and the consequences
%   are not cosmetic: supporting points placed where the signal may never go are
%   never covered, saturate at psiMax, and dilute J permanently -- and in
%   CERTIFICATION mode J = 0 becomes unattainable, so the run cannot terminate.
%   Run assembly therefore reads the declarations and reports what they imply.
%
%   IT CORRECTS NOTHING AND CLIPS NOTHING. A clip only ever narrows, which
%   looks safe, but it changes the normalization and therefore the meaning of every
%   coordinate in every stored run. The inputs come back untouched.
%
%   Two checks:
%
%       ROI inside admissible   the ROI box against the admissible box, dimension
%                               by dimension, matched by name. Reports the fraction
%                               of ROI box VOLUME lying outside. WARNING in budget
%                               mode, ERROR in certification mode.
%       rate consistency        a `du<X>_k` admissible bound against channel X's
%                               declared rate limit, where both exist. WARNING
%                               ONLY, in BOTH modes -- they are separate
%                               declarations by design, and escalating this
%                               in certification mode would be wrong.
%
%   Coordinate names come from rid.rep.Narx.coordinateName, never from a
%   string assembled here.
%
%   WHAT THIS DOES NOT CLAIM. It compares DECLARED BOXES. It cannot see the
%   REACHABLE region: a part of the ROI inside the admissible box may still be
%   dynamically unreachable, and no box comparison will say so -- an asymmetric rate
%   limit shrinks the horizon-reachable set far below the amplitude box. That is
%   deliberately not implemented here. The report detects a declared
%   inconsistency, not unreachability, and must not be read as more.
%
%   Inputs:
%       box              (nDim x 2)   the ROI box, physical
%       dimName          (1 x nDim)   its coordinate names
%       constraints      (1x1)        rid.gen.Constraints
%       opts.mode        (1x1) string "budget" (default) or "certification"
%   Outputs:
%       report  (1x1) struct
%           .fractionOutside   (1x1)      ROI box volume outside the admissible box
%           .perDimOverlap     (1 x nDim) overlap fraction per dimension
%           .roiOutside        (1x1) logical
%           .rateDisagreement  (1 x nChecked) struct .channel .declared .admissible
%           .mode              (1x1) string
%           .messages          (1 x nMsg) string, everything reported

arguments
    box         (:,2) double {mustBeReal}
    dimName     (1,:) string
    constraints (1,1) rid.gen.Constraints
    opts.mode   (1,1) string {mustBeMember(opts.mode, ...
                    ["budget", "certification"])} = rid.defaultOptions().crit.mode
end

nDim = size(box, 1);
if numel(dimName) ~= nDim
    error("rid:gen:invalidRoiBox", ...
        "dimName has %d entries; box has %d row(s).", numel(dimName), nDim);
end

messages = strings(1, 0);

% --- check 1: the ROI box against the admissible box -------------------------
%
% A NAME PRESENT IN ONE BOX AND NOT THE OTHER IS AN ERROR, NOT A SKIP. The declaration-validation rule puts
% the missing-name and unknown-name checks on `bounds`, which both boxes already
% passed; two boxes that disagree about which coordinates exist therefore did not
% come from one representation, and comparing the dimensions they happen to share
% would report a fraction computed over a different space than the one declared.
% Silently skipping is the acceptance the project's coding standard forbids, and the fraction would be
% quietly wrong rather than absent.
missing = setdiff(dimName, constraints.dimName, "stable");
extra   = setdiff(constraints.dimName, dimName, "stable");
if ~isempty(missing) || ~isempty(extra)
    error("rid:gen:dimensionNameMismatch", ...
        "The ROI box and the admissible declaration describe different " + ...
        "coordinates: %s only in the ROI, %s only in the admissible region. " + ...
        "Both are declared over one representation, so this is a " + ...
        "declaration error, not a subset to compare.", ...
        listOrNone(missing), listOrNone(extra));
end

perDimOverlap = zeros(1, nDim);
for iDim = 1:nDim
    iAdm = find(constraints.dimName == dimName(iDim), 1);
    adm  = constraints.admissibleBox(iAdm, :);
    roi  = box(iDim, :);
    width = roi(2) - roi(1);
    overlap = min(roi(2), adm(2)) - max(roi(1), adm(1));
    if width <= 0
        perDimOverlap(iDim) = 1;          % degenerate ROI extent: nothing outside
    else
        perDimOverlap(iDim) = max(0, min(1, overlap / width));
    end
end

fractionOutside = 1 - prod(perDimOverlap);
roiOutside      = fractionOutside > 0;

if roiOutside
    worst = find(perDimOverlap < 1);
    detail = strjoin(arrayfun(@(k) sprintf("%s %.3g%% inside", ...
        dimName(k), 100 * perDimOverlap(k)), worst), ", ");
    msg = sprintf("%.4g%% of the ROI box volume lies outside the admissible " + ...
        "region (%s). Supporting points there can never be covered: they " + ...
        "saturate at psiMax and dilute J permanently.", ...
        100 * fractionOutside, detail);
    messages(end+1) = msg;                                         %#ok<AGROW>
    if opts.mode == "certification"
        % In certification mode J = 0 is the termination condition, and an
        % uncoverable supporting point makes it unattainable: the run cannot
        % terminate. That is a stop, not a warning.
        error("rid:gen:roiOutsideAdmissible", ...
            "%s In certification mode J = 0 is then unattainable and the run " + ...
            "cannot terminate.", msg);
    end
    warning("rid:gen:roiOutsideAdmissible", "%s", msg);
end

% --- check 2: the declared rate limit against the du<X>_k admissible bound ----
%
% WARNING ONLY, IN BOTH MODES. The two are separate declarations by design
%; a disagreement is worth surfacing but is not an inconsistency the run
% cannot survive, and escalating it in certification mode would be wrong.
rateDisagreement = struct("channel", {}, "declared", {}, "admissible", {});
for cIn = 1:numel(constraints.channelName)
    duName = rid.rep.Narx.coordinateName("d" + constraints.channelName(cIn), 0);
    iAdm   = find(constraints.dimName == duName, 1);
    if isempty(iAdm) || ~constraints.rateIsDeclared(cIn)
        continue                       % the check needs BOTH to exist
    end
    admRate = constraints.admissibleBox(iAdm, :);
    decRate = constraints.rate(cIn, :);
    if ~isequal(admRate, decRate)
        rateDisagreement(end+1) = struct( ...                        %#ok<AGROW>
            "channel", constraints.channelName(cIn), ...
            "declared", decRate, "admissible", admRate);
        msg = sprintf("Channel ""%s"": the declared rate limit [%g, %g] " + ...
            "disagrees with the admissible bound on ""%s"", [%g, %g]. These " + ...
            "are separate declarations by design; neither is " + ...
            "corrected.", constraints.channelName(cIn), decRate(1), decRate(2), ...
            duName, admRate(1), admRate(2));
        messages(end+1) = msg;                                       %#ok<AGROW>
        warning("rid:gen:rateDisagreement", "%s", msg);
    end
end

report = struct("fractionOutside", fractionOutside, ...
                "perDimOverlap", perDimOverlap, ...
                "roiOutside", roiOutside, ...
                "rateDisagreement", rateDisagreement, ...
                "mode", opts.mode, "messages", messages);

end

function text = listOrNone(names)
%LISTORNONE  A readable name list, or "(none)" when there is nothing to list.
if isempty(names)
    text = "(none)";
else
    text = strjoin("""" + names + """", ", ");
end
end
