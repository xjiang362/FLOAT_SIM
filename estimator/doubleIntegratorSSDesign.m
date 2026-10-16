function ss = doubleIntegratorSSDesign(b, Tisr, lambdar, lambdae)
%DOUBLEINTEGRATORSSDESIGN Discrete-time state-space regulator + observer
%design for one SISO double-integrator channel (heave, roll, or pitch
%each have this same structure: position/rate driven by a generalized
%force/torque through 1/mass or 1/inertia).
%
%   Continuous plant:  xdot = [0 1;0 0]*x + [0;b]*u,  y = [1 0]*x
%   (x = [position; rate], u = generalized force/torque, b = 1/m or 1/I)
%
%   Discretized exactly via zero-order hold at the ISR sample period Tisr
%   (matches the true actuation model: u is genuinely held constant
%   between ISR ticks -- see simulateTrapDevice.m):
%     Ad = [1 Tisr; 0 1],  Bd = [b*Tisr^2/2; b*Tisr],  Cd = [1 0]
%
%   Regulator gain K: state feedback u = -K*xhat placing the closed-loop
%   poles at the discrete image of a repeated continuous pole at
%   -lambdar (critically damped, bandwidth lambdar rad/s):
%     p_r = exp(-lambdar*Tisr) (double root)
%
%   Observer gain L: predictor-form Luenberger observer
%     xhat[k+1] = Ad*xhat[k] + Bd*u[k] + L*(y[k] - Cd*xhat[k])
%   placing eig(Ad - L*Cd) at the discrete image of a repeated continuous
%   pole at -lambdae (bandwidth lambdae rad/s; should be faster than
%   lambdar, typically 3-6x, for the separation principle to hold up
%   well in practice despite this only being approximately decoupled).

Ad = [1, Tisr; 0, 1];
Bd = [b * Tisr^2 / 2; b * Tisr];
Cd = [1, 0];

pr = exp(-lambdar * Tisr);
pe = exp(-lambdae * Tisr);

K = ackermannSISO(Ad, Bd, [pr, pr]);
L = ackermannSISO(Ad', Cd', [pe, pe])';

ss.Ad = Ad;
ss.Bd = Bd;
ss.Cd = Cd;
ss.K  = K;
ss.L  = L;

end
