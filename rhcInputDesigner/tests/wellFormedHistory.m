function hist = wellFormedHistory(u, opts)
%WELLFORMEDHISTORY  A history that satisfies rid.assertHistory, for tests.
%
%   Most tests care about one field. `rid.assertHistory` requires `.u`, `.x0` and
%   `.tS` in every run, so a test that only cares about the input still
%   has to supply the other two -- and writing them out at eighty call sites would
%   bury the one value each test is actually about.
%
%   The placeholders are DELIBERATELY BORING: `x0 = 0` and `tS = 0.1`. A test that
%   depends on either passes it explicitly, which then reads as the point of the
%   test rather than as scenery.
%
%   `.y` is set only when asked for. Absent means offline, which is the only way
%   to say it -- an empty `y` would be a sentinel.
%
%   Inputs:
%       u          (nSteps x nInputs)  the applied input
%       opts.y     (nSteps x nOutputs) measured output; omitted means offline
%       opts.x0    (1 x nStates)       initial surrogate state
%       opts.tS    (1x1)               sample time
%   Outputs:
%       hist       (1x1) struct

arguments
    u        (:,:) double
    opts.y   (:,:) double = zeros(0, 0)
    opts.x0  (1,:) double = 0
    opts.tS  (1,1) double = 0.1
end

hist = struct("u", u, "x0", opts.x0, "tS", opts.tS);
if ~isempty(opts.y)
    hist.y = opts.y;
end

end
