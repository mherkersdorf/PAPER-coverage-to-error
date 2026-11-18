function xNext = duffingStepFcn(x, u, params)
%DUFFINGSTEPFCN  One fixed-step RK4 step of a damped, driven Duffing oscillator.
%
%   The smoothly nonlinear test process, in the state basis
%   `x = [position, velocity]`:
%
%       dx1/dt = x2
%       dx2/dt = -delta x2 - alpha x1 - beta x1^3 + gain u
%
%   SMOOTHLY NONLINEAR IS THE POINT (the joint-space representation contract, and the third known risk of
%   the method description). A piecewise-linear process whose probe and
%   check points share one linear region PASSES the affineness check, and the loop
%   would then take the analytic path and reuse one Jacobian for the whole run --
%   a test built on such a process would report the value-only tier only by luck.
%   The cubic term makes the map non-affine everywhere except at the origin.
%
%   A FIXED-STEP DISCRETE MAP, NOT AN ADAPTIVE SOLVE. The loop
%   re-simulates the whole history three times per step and every sample here is a
%   MATLAB call, so an adaptive integration per sample would put the run's cost
%   somewhere nobody can predict. RK4 at the sample time is deterministic and
%   costs four field evaluations, always.
%
%   NAMED, NOT ANONYMOUS, because rid.surrogate.StateSpaceFunction refuses a
%   closure: a run record cannot see the values an anonymous handle captures, so a
%   replay could not rebuild the object.
%
%   Inputs:
%       x       (1 x 2)      state before this sample's input acts
%       u       (1 x 1)      input applied at this sample
%       params  (1x1) struct .tS .delta .alpha .beta .gain
%   Outputs:
%       xNext   (1 x 2)

h  = params.tS;
k1 = duffingField(x,               u, params);
k2 = duffingField(x + (h / 2) * k1, u, params);
k3 = duffingField(x + (h / 2) * k2, u, params);
k4 = duffingField(x + h * k3,       u, params);

xNext = x + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

end

% ---------------------------------------------------------------------------

function dx = duffingField(x, u, params)
%DUFFINGFIELD  The continuous vector field, as a (1 x 2) row.
dx = [x(2), ...
      -params.delta * x(2) - params.alpha * x(1) - params.beta * x(1)^3 + ...
       params.gain * u];

end
