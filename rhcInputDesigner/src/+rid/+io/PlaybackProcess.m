classdef PlaybackProcess < handle
%PLAYBACKPROCESS  A record's stored measurements, handed back one unit at a time.
%
%   `rid.io.replay` must reproduce an ONLINE run without a device. The measurements
%   are in the record, so this stands where the process stood and returns them in
%   order. It is not hardware and cannot become hardware: it has no connection to
%   anything, and it runs out.
%
%   IT IGNORES THE INPUT IT IS GIVEN, and that is correct rather than sloppy. What
%   a replay checks is whether the DESIGN reproduces; the plant's response is a
%   fixed record of what happened. While the reproduced input matches the stored
%   one the playback is exactly faithful, and the step at which it stops matching
%   is the divergence `replay` reports -- so a mismatch is caught by the comparison
%   rather than papered over by a re-simulation here.
%
%   RUNNING OUT IS AN ERROR, not a wrap-around or a zero. A replay that asked for
%   more measurements than the run took has already diverged, and quietly
%   supplying something would turn that into a subtle numerical difference many
%   steps later.
%
%   It is a handle because it has a cursor. That makes it a third piece of mutable
%   state, where the project's coding standard names two -- and the exception is narrow and stated: it
%   holds a read position into an immutable array, not a derived quantity that can
%   go stale, and it is created and destroyed inside one `replay` call. It is
%   public only because MATLAB does not permit a class definition in a `private`
%   folder; nothing outside `rid.io.replay` needs it, though replaying a stored
%   experiment against its stored measurements is a legitimate use of it.
%
%   IT PLAYS BACK THE MEASURED STATES TOO, through the two-output
%   form of the process contract. They are a second measurement of the same
%   applied samples, so they are stored and handed back alongside `y` rather than
%   reconstructed: a state re-simulated here would be the surrogate's answer, and
%   what the run was anchored on was the plant's. A playback built from a record
%   that measured none REFUSES the second output rather than returning zeros --
%   an invented state would anchor every horizon at the origin, and the
%   divergence would surface many steps later as a puzzling `u`.
%
%   Construction:
%       process = rid.io.PlaybackProcess(y)
%       process = rid.io.PlaybackProcess(y, x)
%
%   Inputs:
%       y  (nSteps x nOutputs)  the stored measurements, in order
%       x  (nSteps x nStates)   the stored measured states, in order; zeros(0,0)
%                               when the run measured none

    properties (SetAccess = immutable)
        y  (:,:) double   % (nSteps x nOutputs) the stored measurements
        x  (:,:) double   % (nSteps x nStates)  the stored states, or empty
    end

    properties (SetAccess = private)
        cursor (1,1) double = 0   % rows already handed back
    end

    methods

        function obj = PlaybackProcess(y, x)
            arguments
                y (:,:) double {mustBeReal}
                x (:,:) double {mustBeReal} = zeros(0, 0)
            end
            % THE TWO ARE MEASUREMENTS OF THE SAME APPLIED SAMPLES, so a state
            % block of a different length cannot be a record of the same run.
            % Refused here, where the mismatch is still legible, rather than at
            % the step where one of the two runs out first.
            if ~isempty(x) && size(x, 1) ~= size(y, 1)
                error("rid:io:dimensionMismatch", ...
                    "The playback holds %d output row(s) and %d state row(s). " + ...
                    "Both are measurements of the same applied samples " + ...
                    ", so they have the same length.", ...
                    size(y, 1), size(x, 1));
            end
            obj.y = y;
            obj.x = x;
        end

        function [yMeasured, xMeasured] = apply(obj, uApplied)
            %APPLY  The next nApplied stored measurements (the shape).
            %
            %   The second output is the joint-space representation contract's, and exists only when the
            %   recorded run measured states.
            %
            %   Inputs:
            %       uApplied   (nApplied x nInputs)  ignored; see the class header
            %   Outputs:
            %       yMeasured  (nApplied x nOutputs)
            %       xMeasured  (nApplied x nStates)
            if nargout >= 2 && isempty(obj.x)
                error("rid:io:noStatePlayback", ...
                    "This playback was built from a record that measured no " + ...
                    "states, so it has none to hand back. A zero would be a " + ...
                    "state the run never visited: every horizon is anchored " + ...
                    "on the last measured one, so the replay " + ...
                    "would diverge silently. Replay a run whose " + ...
                    "representation read measured states, or ask for one " + ...
                    "output.");
            end
            nApplied = size(uApplied, 1);
            first    = obj.cursor + 1;
            last     = obj.cursor + nApplied;
            if last > size(obj.y, 1)
                error("rid:io:playbackExhausted", ...
                    "The replay asked for measurement rows %d-%d and the " + ...
                    "record holds %d. It has designed more steps than the " + ...
                    "original run, so it had already diverged before this " + ...
                    "point.", first, last, size(obj.y, 1));
            end
            yMeasured  = obj.y(first:last, :);
            xMeasured  = zeros(0, 0);
            if ~isempty(obj.x)
                xMeasured = obj.x(first:last, :);
            end
            obj.cursor = last;
        end

    end

end
