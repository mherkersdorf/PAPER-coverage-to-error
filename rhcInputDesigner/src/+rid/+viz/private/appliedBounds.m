function b = appliedBounds(gen, nInputs)
%APPLIEDBOUNDS  The hard box on u and on du, from the generator's contract.
%
%   ON THE CONTRACT, SO IT IS ALWAYS THERE. `rid.gen.Generator` declares
%   `constraints` alongside `nHorizon`, because the amplitude and rate limits
%   belong to the ACTUATOR: every parametrization that designs an input signal
%   is subject to the same physical box, and what varies between them is only
%   how that box is enforced in their own parameter space, which is what
%   `feasibleSet` returns.
%
%   THIS FUNCTION USED TO ASK WHETHER THE BOX EXISTED, and that was the defect.
%   A generator without it did not raise; the plot simply drew no bounds, and a
%   run pinned against its actuator limits became indistinguishable from one
%   that never approached them. The refusal now happens at construction, where a
%   missing box is a broken generator rather than a figure with one annotation
%   fewer, so there is nothing left here to ask.
%
%   WHAT IS STILL CHECKED is coherence, not presence: a box declaring a
%   different number of channels than the history carries cannot be drawn
%   against it, and drawing it against the wrong channel would put a false
%   statement on the axes.
%
%   Inputs:
%       gen      (1x1) rid.gen.Generator
%       nInputs  (1x1) double  channels in the signal history
%   Outputs:
%       b  (1x1) struct
%           .channelName  (1 x nInputs) string
%           .amplitude    (nInputs x 2) double  [uMin uMax], physical
%           .rate         (nInputs x 2) double  [duMin duMax], physical
%           .rateSource   (1 x nInputs) string  "declared" | "duCoordinate"
%                                               | "amplitudeWidth"

con = gen.constraints;

if numel(con.channelName) ~= nInputs
    error("rid:viz:channelMismatch", ...
        "The generator's constraints declare %d channel(s) and the history " + ...
        "carries %d. Drawing one against the other would put a bound on the " + ...
        "wrong axes.", numel(con.channelName), nInputs);
end

b = struct( ...
    "channelName", con.channelName, ...
    "amplitude",   con.amplitude, ...
    "rate",        con.rate, ...
    "rateSource",  con.rateSource);

end
