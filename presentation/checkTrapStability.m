%% checkTrapStability.m
% Open-loop stiffness of the levitated plate at the hover pose, split into
% the three UNSTABLE degrees of freedom (heave and the two tilts) and the
% three that should be passively trapped (the two lateral translations and
% yaw).
%
% Method: the generalised force Q(q) = B(q)*u.^2 is differentiated by
% central difference in each pose coordinate, giving the 6-by-6 stiffness
%       K = -dQ/dq
% with the sign chosen so POSITIVE K is a restoring spring. Eigenvalues of
% the mass-normalised blocks then give either an oscillation frequency
% (stable) or a divergence rate (unstable).
%
% NOTE the lateral block rests entirely on the plate-electrode overlap
% derivation, since all lateral force data in this project is corrupted by
% meshing and cannot be used to check it.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'));
p = trapDeviceGeometry();

%% ---- hover voltage, four quadrants equal ----
pose0 = zeros(6,1);
[~,J0]   = trapCapacitance(p, pose0);
[~,~,B0] = trapWrench(p, J0, zeros(6,1));
sel = [1 1 1 1 0 0]';
W   = p.mass.m*p.phys.g;
vHover = W/(B0(3,:)*sel);
u = sqrt(vHover)*sel;
fprintf('hover: %.2f V on each of the four quadrant pairs (weight %.1f uN)\n', u(1), W*1e6);

%% ---- stiffness by central difference ----
h = [1e-9 1e-9 1e-9 1e-7 1e-7 1e-7];     % m, m, m, rad, rad, rad
K = zeros(6,6);
for j = 1:6
    qp = pose0; qp(j) = qp(j)+h(j);
    qm = pose0; qm(j) = qm(j)-h(j);
    [~,Jp] = trapCapacitance(p,qp); [Fp,Tp] = trapWrench(p,Jp,u);
    [~,Jm] = trapCapacitance(p,qm); [Fm,Tm] = trapWrench(p,Jm,u);
    K(:,j) = -([Fp;Tp]-[Fm;Tm])/(2*h(j));
end

lbl = {'rx','ry','rz','thx','thy','thz'};
fprintf('\n=== stiffness K = -dQ/dq  (positive = restoring) ===\n');
fprintf('%8s','');
fprintf('%12s', lbl{:}); fprintf('\n');
un = {'N/m','N/m','N/m','Nm/rad','Nm/rad','Nm/rad'};
for i=1:6
    fprintf('%8s', lbl{i});
    fprintf('%12.4g', K(i,:));
    fprintf('   [%s]\n', un{i});
end

%% ---- diagonal summary ----
M = [p.mass.m p.mass.m p.mass.m diag(p.mass.I)'];
fprintf('\n=== per-axis behaviour ===\n');
fprintf('%6s %14s %10s %16s\n','dof','K_ii','sign','rate / freq');
for i=1:6
    k = K(i,i); om = sqrt(abs(k)/M(i));
    if k > 0
        s = 'restoring'; r = sprintf('%8.1f rad/s osc (%.1f Hz)', om, om/2/pi);
    else
        s = 'DESTABIL.'; r = sprintf('%8.1f rad/s growth (%.1f Hz)', om, om/2/pi);
    end
    fprintf('%6s %14.4g %10s   %s\n', lbl{i}, k, s, r);
end

%% ---- lateral block, properly ----
lat = [1 2 6];
Klat = K(lat,lat);
Mlat = diag(M(lat));
fprintf('\n=== lateral block (rx, ry, thz) ===\n');
disp(Klat);
ev = eig(Klat, Mlat);
fprintf('  generalised eigenvalues (K,M): %s\n', sprintf('%.4g  ', ev));
if all(real(ev) > 0)
    fprintf('  ALL POSITIVE -> the lateral 3 DOF are PASSIVELY STABLE\n');
    fprintf('  trap frequencies: %s Hz\n', sprintf('%.1f  ', sqrt(ev)/2/pi));
else
    fprintf('  NOT all positive -> lateral motion is not passively trapped\n');
end

%% ---- unstable block ----
vert = [3 4 5];
Kv = K(vert,vert); Mv = diag(M(vert));
evv = eig(Kv, Mv);
fprintf('\n=== heave / tilt block (rz, thx, thy) ===\n');
fprintf('  generalised eigenvalues: %s\n', sprintf('%.4g  ', evv));
fprintf('  divergence rates: %s rad/s\n', sprintf('%.1f  ', sqrt(abs(evv))));
fprintf('  -> %d of 3 unstable\n', sum(real(evv)<0));
