function y = duffingOutputFcn(x, ~, ~)
%DUFFINGOUTPUTFCN  The measured output of the Duffing process: its position.
%
%   The partner of `duffingStepFcn`. It reads x(k) -- the state BEFORE this
%   sample's input acts -- and has no direct term, which is what a position
%   measurement of a second-order mechanical system looks like.
%
%   Inputs:
%       x  (1 x 2)  state before this sample's input acts
%   Outputs:
%       y  (1x1)

y = x(1);

end
