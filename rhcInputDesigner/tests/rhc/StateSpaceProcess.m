classdef StateSpaceProcess < handle
%STATESPACEPROCESS  A plant that measures its states, for the joint-space representation contract.
%
%   The two-output form of the process contract:
%
%       [yMeasured, xMeasured] = process.apply(uApplied)
%
%   with `xMeasured` `(nApplied x nStates)` and row `k` the state at the instant
%   `u(k)` was applied -- the state BEFORE that input acts, which is the row
%   convention of the joint-space representation contract and the one `simulateStates` returns.
%
%   ITS DYNAMICS DIFFER FROM THE SURROGATE'S ON PURPOSE, in the same way
%   `FirstOrderProcess` does: the surrogate is a geometric guide, not a predictive
%   model. Its states are in the SAME basis, because
%   that is the scientific assumption a measured-state run rests on and cannot
%   check -- only the pole and the gain move. A plant identical to
%   the surrogate would make the anchor of the joint-space representation contract invisible: the anchored and
%   the re-simulated state would agree, and a test could not tell them apart.
%
%   `lastNargout` RECORDS HOW MANY OUTPUTS THE LOOP ASKED FOR. That is the sharp
%   assertion for the other half of the contract: a representation that does not
%   read measured states must leave every existing one-output adapter untouched,
%   and the only way to see that the loop did not ask is to ask the process.
%
%   It is a handle because it carries the plant's state between calls, which is
%   what a real process does. Nothing in `src/` may depend on that: the loop calls
%   `apply` and nothing else.

    properties (SetAccess = immutable)
        Ad      (:,:) double    % (nStates x nStates)
        Bd      (:,:) double    % (nStates x nInputs)
        Cd      (1,:) double    % (1 x nStates)
        defect  (1,1) string    % how the measured states are wrong, if at all
    end

    properties (SetAccess = private)
        x             (1,:) double     % the state BEFORE the next input acts
        nApplications (1,1) double = 0
        lastNargout   (1,1) double = NaN
    end

    methods

        function obj = StateSpaceProcess(Ad, Bd, Cd, x0, opts)
            %STATESPACEPROCESS  A linear plant in the surrogate's state basis.
            %
            %   Inputs:
            %       Ad   (nStates x nStates)
            %       Bd   (nStates x nInputs)
            %       Cd   (1 x nStates)
            %       x0   (1 x nStates)   the state measured before the run begins
            %       opts.defect  (1x1) string  "none", "wrongWidth" or "notFinite"
            arguments
                Ad (:,:) double
                Bd (:,:) double
                Cd (1,:) double
                x0 (1,:) double
                opts.defect (1,1) string ...
                    {mustBeMember(opts.defect, ["none", "wrongWidth", "notFinite"])} = "none"
            end
            obj.Ad     = Ad;
            obj.Bd     = Bd;
            obj.Cd     = Cd;
            obj.x      = x0;
            obj.defect = opts.defect;
        end

        function [yMeasured, xMeasured] = apply(obj, uApplied)
            %APPLY  (nApplied x nInputs) in; its measurement and its states out.
            obj.lastNargout   = nargout;
            obj.nApplications = obj.nApplications + 1;

            nApplied  = size(uApplied, 1);
            nStates   = numel(obj.x);
            yMeasured = zeros(nApplied, 1);
            xMeasured = zeros(nApplied, nStates);

            for iRow = 1:nApplied
                % Recorded BEFORE the step: row k is the state at the instant
                % u(k) was applied, the quantity `y(k) = C x(k)` also reads.
                xMeasured(iRow, :) = obj.x;
                yMeasured(iRow, 1) = obj.Cd * obj.x.';
                obj.x = (obj.Ad * obj.x.' + obj.Bd * uApplied(iRow, :).').';
            end

            switch obj.defect
                case "wrongWidth"
                    % One state short: the shape the loop must refuse rather than
                    % append to a history whose width then disagrees with the
                    % representation's (needs nStates >= 2 to be a real defect).
                    xMeasured = xMeasured(:, 1:end - 1);
                case "notFinite"
                    xMeasured(1, 1) = NaN;
            end
        end

    end

end
