function p = p2DuffingParams()
%P2DUFFINGPARAMS  Parameters of the P2 Duffing test process.
%
%   x1' = x2
%   x2' = -omega0^2*(x1 + alpha*x1^3) - 2*zeta*omega0*x2 + omega0^2*u
%   y   = x1
%
%   The unique, globally asymptotically stable equilibrium exists for every
%   constant u (energy argument). timeConstant is the envelope
%   decay 1/(zeta*omega0); it does not depend on the operating point because the
%   damping is linear, so it is declared once here rather than derived per call.
%
%   Outputs:
%       p  (1x1) struct with fields:
%           omega0        (1x1)  natural frequency [rad/s]
%           zeta          (1x1)  damping ratio
%           alpha         (1x1)  hardening-spring coefficient (starting value 1;
%                                settled by the LME procedure)
%           tS            (1x1)  sample time [s]; T/tS ~ 17
%           nSub          (1x1)  RK4 substeps per sample, u held constant
%           timeConstant  (1x1)  T = 1/(zeta*omega0) [s]

    p.omega0 = 1.0;
    p.zeta   = 0.3;
    p.alpha  = 1.0;
    p.tS     = 0.2;
    p.nSub   = 1;

    p.timeConstant = 1 / (p.zeta * p.omega0);
end
