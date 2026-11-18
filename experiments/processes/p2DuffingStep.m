function xNext = p2DuffingStep(x, u, p)
%P2DUFFINGSTEP  One sample of the P2 Duffing process, S0's stepFcn.
%
%   Classical RK4, u held constant over the sample, p.nSub substeps -- a
%   fixed-step discrete map, no adaptive ODE solve.
%
%   Inputs:
%       x  (1x2) double  [x1, x2]
%       u  (1x1) double
%       p  (1x1) struct  p2DuffingParams()
%   Outputs:
%       xNext  (1x2) double
arguments
    x (1,2) double {mustBeReal, mustBeFinite}
    u (1,1) double {mustBeReal, mustBeFinite}
    p (1,1) struct
end

h  = p.tS / p.nSub;
xc = x(:);
for iSub = 1:p.nSub
    k1 = duffingDerivative(xc,          u, p);
    k2 = duffingDerivative(xc + h/2*k1, u, p);
    k3 = duffingDerivative(xc + h/2*k2, u, p);
    k4 = duffingDerivative(xc + h*k3,   u, p);
    xc = xc + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end
xNext = xc.';
end

function xDot = duffingDerivative(xc, u, p)
xDot = [ ...
    xc(2); ...
    -p.omega0^2 * (xc(1) + p.alpha*xc(1)^3) - 2*p.zeta*p.omega0*xc(2) + p.omega0^2*u ...
    ];
end
