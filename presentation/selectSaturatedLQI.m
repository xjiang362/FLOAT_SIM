%% selectSaturatedLQI.m
% Full-state LQI redesigned to live inside a +/-1g acceleration envelope,
% with anti-windup on the integrator.
%
% WHY +/-1g. The electrode array is ABOVE the plate only, so every pair
% pulls upward and the sole downward force is gravity. The achievable
% acceleration is therefore [-g, +22g] at Vmax = 120 V -- wildly asymmetric.
% Braking is the scarce direction, so the usable SYMMETRIC envelope is
% +/-1g, and a controller that commands more than that is writing cheques
% the actuator cannot cash.
%
% THE HARD LIMIT THIS IMPOSES. The uncancelled electrostatic stiffness
% contributes an acceleration kappa_z * r_z away from the setpoint. Merely
% holding station at displacement r_z therefore costs
%       |a_z| = kappa_z |r_z|
% so with |a_z| <= g the plate can only be held inside
%       |r_z| <= g/kappa_z = %.2f um
% Beyond that the spring beats the actuator and no control law recovers.
% Delay shrinks it further: the plate grows by cosh(p_z t_d) = 1.96 before
% the oldest channel is even read.
%
% ANTI-WINDUP. Back-calculation:
%       x_I[k+1] = x_I[k] + T(r - H y[k]) + K_aw (u_sat[k] - u[k])
% so the integrator accumulates only what the actuator actually delivered.
% Without it the integrator charges through every saturated interval and
% the loop fails even inside the envelope above.
%
% See also TRAPFULLSTATEMODEL, SIMULATEFULLSTATELQI.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p  = trapDeviceGeometry();
g  = p.phys.g;  m = p.mass.m;
Jd = [m; p.mass.I(1,1); p.mass.I(2,2)];
Mref = trapFullStateModel(p);
kz = Mref.kappa(1);  pz = Mref.omega(1);  td = max(Mref.age);

rHold = g/kz;
fprintf('kappa_z = %.4e 1/s^2  ->  static hold envelope |rz| <= g/kappa_z = %.2f um\n', kz, rHold*1e6);
fprintf('growth before the oldest channel is read: cosh(p_z t_d) = %.3f\n', cosh(pz*td));
fprintf('  => releasable |rz0| <~ %.2f um before the spring already outruns 1g\n\n', rHold/cosh(pz*td)*1e6);

aMax = g;             % +/-1g symmetric envelope
Kaw  = 1.0;           % back-calculation gain

%% ---- 1. control weight vs peak demand, WITH saturation + anti-windup ----
fprintf('=== peak |a_z| in the linear loop, rz0 = -3 um, sat +/-1g, AW on ===\n');
fprintf('%8s %10s %12s %12s %12s %10s\n','rhoU','max|z|','min az[g]','max az[g]','peak rz[um]','sat ticks');
best = [];
for ru = [1e2 1e3 1e4 1e5 1e6 1e7]
    M = trapFullStateModel(p,'rhoU',ru,'rhoI',1e2);
    S = linSat(M, -3e-6, aMax, Kaw);
    fprintf('%8.0e %10.4f %12.2f %12.2f %12.2f %9.0f%%\n', ru, max(abs(M.clR)), ...
        min(S.a)/g, max(S.a)/g, S.pkZ*1e6, 100*S.fsat);
    if isempty(best) && max(abs(S.a)) <= aMax*1.001, best = ru; end
end

%% ---- 2. capture envelope at a chosen weight ----
ru = 1e5;
fprintf('\n=== release envelope at rhoU = %.0e (sat +/-1g, AW on) ===\n', ru);
M = trapFullStateModel(p,'rhoU',ru,'rhoI',1e2);
fprintf('%10s %12s %12s %12s %10s %16s\n','rz0[um]','min az[g]','max az[g]','peak rz[um]','sat ticks','nonlinear plant');
for z0 = [-1 -2 -3 -4 -5 -6]*1e-6
    S = linSat(M, z0, aMax, Kaw);
    fate = nlSat(p, M, Jd, z0, aMax, Kaw, 0.3);
    fprintf('%10.1f %12.2f %12.2f %12.2f %9.0f%% %16s\n', z0*1e6, ...
        min(S.a)/g, max(S.a)/g, S.pkZ*1e6, 100*S.fsat, fate);
end

%% ---- 3. does anti-windup matter? ----
fprintf('\n=== anti-windup on vs off, rhoU = %.0e ===\n', ru);
fprintf('%10s %18s %18s\n','rz0[um]','AW on','AW off');
for z0 = [-1 -2 -3]*1e-6
    fprintf('%10.1f %18s %18s\n', z0*1e6, ...
        nlSat(p, M, Jd, z0, aMax, Kaw, 0.3), ...
        nlSat(p, M, Jd, z0, aMax, 0,   0.3));
end

% =====================================================================
function S = linSat(M, rz0, aMax, Kaw)
T = M.T;
x = [rz0;0;0;0;0;0]; z = M.Phi*x;
xh = zeros(9,1); uPrev = zeros(3,1); xI = zeros(3,1);
S.a = []; S.pkZ = abs(rz0); nsat = 0; N = round(0.3/T);
for k = 1:N
    y = M.Cbx*x + M.Cbz*z + M.Dbar*uPrev;
    u = -M.Kx*xh(1:6) - M.Kz*xh(7:9) - M.Ku*uPrev - M.KI*xI;
    us = max(min(u, aMax), -aMax);
    if any(abs(us-u) > 1e-12), nsat = nsat+1; end
    xI = xI + T*(0 - M.H*y) + Kaw*(us - u);
    S.a(end+1) = us(1);
    xn = M.Ad*x + M.Bd*us;  zn = M.Phi*x + M.Gam*us;
    xh = M.Abar*xh + M.Bbar*us + M.L*(y - M.Cbar*xh - M.Dbar*uPrev);
    uPrev = us; x = xn; z = zn;
    S.pkZ = max(S.pkZ, abs(x(1)));
end
S.fsat = nsat/N;
end

function fate = nlSat(p, M, Jd, rz0, aMax, Kaw, tEnd)
T = M.T; dt = p.sim.dtEuler; nSub = round(T/dt); nT = floor(tEnd/T);
x = zeros(12,1); x(3) = rz0;
xh = zeros(9,1); uPrev = zeros(3,1); xI = zeros(3,1);
s = trapSamplerInit(p, T); fate = '';
for k = 1:nT
    t = (k-1)*T;
    if ~checkPoseValid(p, x(1:6)), fate = sprintf('contact %.1f ms', t*1e3); break; end
    try
        [s, CeD, si] = trapSamplerStep(s, x(1:6), t);
        if si.valid, Ce = CeD; else, Ce = si.CeffNow; end
        y = trapCapToHeight(p, Ce);
    catch
        fate = sprintf('contact %.1f ms', t*1e3); break
    end
    u  = -M.Kx*xh(1:6) - M.Kz*xh(7:9) - M.Ku*uPrev - M.KI*xI;
    us = max(min(u, aMax), -aMax);
    xI = xI + T*(0 - M.H*y) + Kaw*(us - u);
    V  = trapWrenchToVoltage(p, Jd.*us, zeros(6,1), 'Vtia', 0, 'Vmax', 120);
    hit = false;
    for j = 1:nSub
        if ~checkPoseValid(p, x(1:6)), hit=true; fate=sprintf('contact %.1f ms',(t+(j-1)*dt)*1e3); break; end
        x = x + dt*trapDynamics(p, x, V);
    end
    if hit, break; end
    xh = M.Abar*xh + M.Bbar*us + M.L*(y - M.Cbar*xh - M.Dbar*uPrev);
    uPrev = us;
end
if isempty(fate)
    if abs(x(3)) < 2e-7, fate = 'settles';
    else, fate = sprintf('drift %.2f um', x(3)*1e6); end
end
end
