%% selectFullStateWeights.m
% Weight selection for the full-state LQI design. Sweeps the control weight
% rhoU, the integrator weight rhoI and the observer measurement weight
% rMeas, and records for each: the closed-loop spectral radius, the
% equivalent bandwidth, the peak commanded acceleration and peak excursion
% from a standard release, and whether the REAL nonlinear plant survives.
%
% The binding constraints are physical, not spectral:
%   (i)  the actuator is UNILATERAL -- it can only pull toward the
%        electrodes, so the most negative achievable a_z is -g = -9.81
%        m/s^2, while upward it is limited by Vmax;
%   (ii) the plate has +/-15 um of travel, so any design whose transient
%        overshoots that is useless however well-damped it is on paper.
% A spectrally excellent design (rho_u small, max|z| ~ 0.23) violates both.
%
% Writes presentation/fullstate_selection_table.tex.
%
% See also TRAPFULLSTATEMODEL, SIMULATEFULLSTATELQI.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p  = trapDeviceGeometry();
m  = p.mass.m; Ix = p.mass.I(1,1); Iy = p.mass.I(2,2);
Jd = [m; Ix; Iy];
W  = m*p.phys.g;

% available acceleration envelope
Mref = trapFullStateModel(p);
[~,Jc] = trapCapacitance(p, zeros(6,1));
[~,~,Bq] = trapWrench(p, Jc, zeros(6,1));
aUp = (sum(Bq(3,1:5))*120^2 - W)/m;
fprintf('actuator envelope: a_z in [%.2f, %+.1f] m/s^2 at Vmax = 120 V\n', -p.phys.g, aUp);
fprintf('travel: r_z in [%.0f, %+.0f] um\n\n', -p.gap.z0*1e6, p.gap.z0*1e6);

rz0 = -5e-6;   tEnd = 0.30;  Vmax = 120;

%% ---- sweep ----
rhoUs = [1e-2 1e-1 1e0 1e1 1e2 1e3 1e4];
fprintf('=== control weight sweep (rhoI = 1e2, rMeas = 1e-2) ===\n');
fprintf('%8s %10s %10s %12s %12s %16s\n', ...
    'rhoU','max|z|','BW[Hz]','peak az','peak rz[um]','nonlinear plant');
Srow = {};
for ru = rhoUs
    M = trapFullStateModel(p, 'rhoU', ru);
    [pkA, pkZ] = linPeak(M, rz0);
    fate = nlFate(p, M, Jd, rz0, tEnd, Vmax);
    bw = -log(max(abs(M.clR)))/M.T/2/pi;
    fprintf('%8.0e %10.4f %10.1f %12.2f %12.2f %16s\n', ru, max(abs(M.clR)), bw, pkA, pkZ*1e6, fate);
    Srow(end+1,:) = {ru, max(abs(M.clR)), bw, pkA, pkZ*1e6, fate}; %#ok<SAGROW>
end

%% ---- pick, then refine the integrator and observer weights ----
fprintf('\n=== integrator weight, at the chosen rhoU = 1e2 ===\n');
fprintf('%8s %10s %12s %12s %16s\n','rhoI','max|z|','peak az','peak rz[um]','nonlinear plant');
for ri = [1e0 1e1 1e2 1e3 1e4]
    M = trapFullStateModel(p, 'rhoU', 1e2, 'rhoI', ri);
    [pkA, pkZ] = linPeak(M, rz0);
    fprintf('%8.0e %10.4f %12.2f %12.2f %16s\n', ri, max(abs(M.clR)), pkA, pkZ*1e6, ...
        nlFate(p, M, Jd, rz0, tEnd, Vmax));
end

fprintf('\n=== observer measurement weight, at rhoU = 1e2, rhoI = 1e2 ===\n');
fprintf('%8s %12s %12s %16s\n','rMeas','obs max|z|','peak rz[um]','nonlinear plant');
for rm = [1e-4 1e-3 1e-2 1e-1 1e0]
    M = trapFullStateModel(p, 'rhoU', 1e2, 'rhoI', 1e2, 'rMeas', rm);
    [~, pkZ] = linPeak(M, rz0);
    fprintf('%8.0e %12.4f %12.2f %16s\n', rm, max(abs(M.clE)), pkZ*1e6, ...
        nlFate(p, M, Jd, rz0, tEnd, Vmax));
end

% =====================================================================
function [pkA, pkZ] = linPeak(M, rz0)
% peak commanded acceleration and peak heave on M's OWN linear model
T = M.T;
x = [rz0;0;0;0;0;0]; z = M.Phi*x;
xh = zeros(9,1); uPrev = zeros(3,1); xI = zeros(3,1);
pkA = 0; pkZ = abs(rz0);
for k = 1:round(0.2/T)
    y = M.Cbx*x + M.Cbz*z + M.Dbar*uPrev;
    u = -M.Kx*xh(1:6) - M.Kz*xh(7:9) - M.Ku*uPrev - M.KI*xI;
    xI = xI + T*(0 - M.H*y);
    pkA = max(pkA, abs(u(1)));
    xn = M.Ad*x + M.Bd*u;  zn = M.Phi*x + M.Gam*u;
    xh = M.Abar*xh + M.Bbar*u + M.L*(y - M.Cbar*xh - M.Dbar*uPrev);
    uPrev = u; x = xn; z = zn;
    pkZ = max(pkZ, abs(x(1)));
end
end

function fate = nlFate(p, M, Jd, rz0, tEnd, Vmax)
T = M.T; dt = p.sim.dtEuler; nSub = round(T/dt); nT = floor(tEnd/T);
x = zeros(12,1); x(3) = rz0;
xh = zeros(9,1); uPrev = zeros(3,1); xI = zeros(3,1);
s = trapSamplerInit(p, T);
fate = '';
for k = 1:nT
    t = (k-1)*T;
    if ~checkPoseValid(p, x(1:6)), fate = sprintf('contact %.1f ms', t*1e3); break; end
    try
        [s, CeD, sInfo] = trapSamplerStep(s, x(1:6), t);
        if sInfo.valid, Ce = CeD; else, Ce = sInfo.CeffNow; end
        y = trapCapToHeight(p, Ce);
    catch
        fate = sprintf('contact %.1f ms', t*1e3); break
    end
    u  = -M.Kx*xh(1:6) - M.Kz*xh(7:9) - M.Ku*uPrev - M.KI*xI;
    xI = xI + T*(0 - M.H*y);
    V  = trapWrenchToVoltage(p, Jd.*u, zeros(6,1), 'Vtia', 0, 'Vmax', Vmax);
    hit = false;
    for j = 1:nSub
        if ~checkPoseValid(p, x(1:6)), hit=true; fate=sprintf('contact %.1f ms',(t+(j-1)*dt)*1e3); break; end
        x = x + dt*trapDynamics(p, x, V);
    end
    if hit, break; end
    xh = M.Abar*xh + M.Bbar*u + M.L*(y - M.Cbar*xh - M.Dbar*uPrev);
    uPrev = u;
end
if isempty(fate)
    if abs(x(3)) < 1e-7, fate = 'settles';
    else, fate = sprintf('drift %.2f um', x(3)*1e6); end
end
end
