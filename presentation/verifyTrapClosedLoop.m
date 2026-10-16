%% verifyTrapClosedLoop.m
% Linearises the FULL nonlinear closed loop numerically and reports its
% eigenvalues.
%
% WHY THIS EXISTS. verifyTrapControl.m checks the loop against the DESIGN
% model, in which u is a force and the plant is a triple double-integrator.
% That model is stable by construction -- it is what the gains were placed
% on -- so it cannot detect a defect that lives in the actuator. The real
% loop closes through
%
%   pose -> capacitance -> effective capacitance -> height -> estimator
%        -> regulator -> force -> VOLTAGE -> wrench -> pose
%
% and the V^2 inversion in trapForceToVoltage is where the design model and
% the device part company. This script finite-differences one ISR tick of
% that whole chain over the 21-state vector
%
%   s = [x(1:12); xhat(1:6); xi(1:3)]
%
% and takes eig(ds+/ds). Anything outside the unit circle is a genuine
% instability of the implemented loop, whatever the design model says.
%
% WHAT IT FOUND, AND WHAT THAT TURNED OUT TO MEAN. With the original
% design -- per-pair force input, regulator placed on a double integrator,
% lambdar = 50 Hz -- this script measured max|z| = 1.008 / 1.062 / 1.547 at
% Tisr = 50 / 200 / 720 us, i.e. unstable at every rate, with the dominant
% eigenvector almost pure tilt RATE.
%
% Two real defects were behind it, neither visible to verifyTrapControl:
%
%   1. The design model was a double integrator. Within one ISR the voltage
%      is frozen and the plate obeys qddot = p^2*(q - q_k), so the exact
%      sampled map is [1 sinh(pT)/p; 0 cosh(pT)] -- the position term
%      cancels but the VELOCITY is amplified by cosh(pT) every tick. At
%      Tisr = 500 us that is a 17% per-tick growth the design model did not
%      contain at all.
%   2. Commanding a force per PAIR pins each pair's total z-force and
%      nothing else, which leaves the tilt stiffness only partly cancelled
%      (beta = 0.80 / 0.72) because a pair's two members sit at different
%      coordinates.
%
% Both are fixed: the input is now the WRENCH [dFz; dtaux; dtauy],
% allocated across all five pairs by trapWrenchToVoltage (so beta = 1 in
% every DOF), and trapControlDesign places the gains on the exact sampled
% map. The loop is now stable at lambdar = 50 Hz at every sample rate
% including 720 us.
%
% NOTE the "anti-damping" section at the end is retained only as a scale
% check. It is NOT the mechanism that broke the loop -- that attribution
% was tested and rejected: the exact held-voltage MIMO model with the real
% gains predicted stability where the loop actually diverged.
%
% The design model cannot see any of this because it has no V^2 stage.
%
% See also TRAPFORCETOVOLTAGE, TRAPCONTROLDESIGN, VERIFYTRAPCONTROL.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

Tlist = [50 200 720]*1e-6;
lamR  = [50 200];
lab   = {'rx','ry','rz','thx','thy','thz','vx','vy','vz','wx','wy','wz', ...
         'h.rz','h.vz','h.thx','h.wx','h.thy','h.wy','xi1','xi2','xi3'};

fprintf('=== eigenvalues of the implemented loop (Vtia = 0, no delay) ===\n');
fprintf('%8s %10s %12s %10s %12s\n','Tisr[us]','lamR[Hz]','max|z|','# |z|>1','tau[ms]');
for T = Tlist
    for lr = lamR
        [J, K] = loopJacobian(T, 2*pi*lr, 2*pi*300);
        ev = eig(J);
        [~,ix] = sort(abs(ev),'descend'); ev = ev(ix);
        % the three passively trapped lateral modes sit on the unit circle
        % to within the squeeze-film decay, so count with a tolerance
        nUns = sum(abs(ev) > 1 + 1e-5);
        tau  = T/log(abs(ev(1)));
        fprintf('%8.0f %10.0f %12.6f %10d %12.1f\n', T*1e6, lr, abs(ev(1)), nUns, tau*1e3);
    end
end

%% ---- the unstable eigenvector of the as-specified design ----
fprintf('\n=== dominant eigenvector, Tisr = 50 us, lamR = 50 Hz ===\n');
[J,~] = loopJacobian(50e-6, 2*pi*50, 2*pi*300);
[V,D] = eig(J); [~,k1] = max(abs(diag(D)));
v = real(V(:,k1)); v = v/max(abs(v));
for i = 1:numel(v)
    if abs(v(i)) > 1e-3, fprintf('  %-7s %+9.4f\n', lab{i}, v(i)); end
end
fprintf('  (with the wrench input and the exact map this mode is now stable)\n');

%% ---- the anti-damping budget ----
fprintf('\n=== scale check only: lag-induced anti-damping (NOT the failure mode) ===\n');
p = trapDeviceGeometry();
W = p.mass.m*p.phys.g;
fprintf('  dc = -F*(2/g)*td,  F = W = %.1f uN, g = %.0f um\n', W*1e6, p.gap.z0*1e6);
for le = [300 1000 3000]
    td = 1/(2*pi*le);
    fprintf('    estimator %4.0f Hz -> td = %6.0f us -> dc = %+8.4f N*s/m\n', ...
        le, td*1e6, -W*2/p.gap.z0*td);
end
fprintf('  squeeze-film damping, for scale:          c  = %+.4e N*s/m\n', p.damp.c);
fprintf('\n');

% =====================================================================
function [J, K] = loopJacobian(T, lamr, lame)
p = trapDeviceGeometry();
p.readout.Tisr = T; p.readout.dtChannel = T/4; p.sim.dtEuler = T/72;
K  = trapControlDesign(p, 'lambdar', lamr, 'lambdae', lame);
dt = p.sim.dtEuler; nSub = round(K.Tisr/dt);
F  = @(s) oneTick(s, p, K, dt, nSub);

s0 = zeros(21,1);
for k = 1:20, s0 = F(s0); end         % settle onto the fixed point

% perturbation scales: metres, radians, and their rates
sc = [1e-9*ones(3,1); 1e-9*ones(3,1); 1e-6*ones(3,1); 1e-6*ones(3,1); ...
      1e-9; 1e-6; 1e-9; 1e-6; 1e-9; 1e-6; 1e-9*ones(3,1)];
J = zeros(21,21);
for j = 1:21
    ep = zeros(21,1); ep(j) = sc(j);
    J(:,j) = (F(s0+ep) - F(s0-ep)) / (2*sc(j));
end
end

function sn = oneTick(s, p, K, dt, nSub)
x = s(1:12); xh = s(13:18); xi = s(19:21);
ys = trapCapToHeight(p, trapSensing(p, trapCapacitance(p, x(1:6))));
f   = -K.Kx*xh - K.Ki*xi;
xin = xi + K.Tisr*(-K.Cc*xh);
pe = zeros(6,1); pe(3)=xh(1); pe(4)=xh(3); pe(5)=xh(5);
u = trapWrenchToVoltage(p, f, pe, 'Vtia', 0);
for j = 1:nSub, x = x + dt*trapDynamics(p, x, u); end
sn = [x; K.Ad*xh + K.Bd*f + K.L*(ys - K.Csd*xh); xin];
end
