%% verifyTrapPipeline.m
% End-to-end checks of the levitation pipeline:
%   trapDeviceGeometry -> trapCapacitance -> trapWrench -> trapDynamics
%   with trapSensing on the readout branch and checkPoseValid on the
%   constraint branch.
%
% These are SELF-CONSISTENCY and PHYSICS checks, not fits. In particular
% none of the lateral results are compared against simulation: all Fx, Fy
% and tauz data in this project is corrupted by meshing, so the lateral
% model rests on its geometric derivation only.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'estimator'), fullfile(root,'presentation'));

p = trapDeviceGeometry();
nFail = 0; verdict = {'FAIL','PASS'};
chk = @(ok,msg) fprintf('  [%s] %s\n', verdict{logical(ok)+1}, msg);

%% ---- 1. hover: what voltage balances gravity? ----
fprintf('=== 1. hover equilibrium ===\n');
W = p.mass.m*p.phys.g;
pose0 = zeros(6,1);
[~,dCdq] = trapCapacitance(p, pose0);
[~,~,B]  = trapWrench(p, dCdq, zeros(6,1));
% all four quadrants at the same voltage, centrals off
sel = [1 1 1 1 0 0]';
Fz_per_v = B(3,:)*sel;                      % N per V^2
vHover = W/Fz_per_v;  uHover = sqrt(vHover);
fprintf('  weight %.2f uN;  dFz/d(u^2) = %.4e N/V^2\n', W*1e6, Fz_per_v);
fprintf('  hover voltage (4 quadrants) = %.2f V\n', uHover);
u = uHover*sel;
[F,tau] = trapWrench(p, dCdq, u);
chk(abs(F(3)-W)/W < 1e-9, sprintf('Fz balances weight to %.1e relative', abs(F(3)-W)/W));
chk(norm(F(1:2)) < 1e-12*W, 'no lateral force at the centred pose, by symmetry');
chk(norm(tau) < 1e-12*W*p.mass.halfL, 'no net torque with all quadrants equal');

%% ---- 2. the allocation matrix ----
fprintf('\n=== 2. allocation matrix B (Q = B*u.^2) ===\n');
fprintf('  %-6s %10s %10s %10s %10s %10s %10s\n','',p.idx.actPairName{:});
rn = {'Fx','Fy','Fz','taux','tauy','tauz'};
sc = [1e9 1e9 1e9 1e12 1e12 1e12];
un = {'nN/V^2','nN/V^2','nN/V^2','pNm/V^2','pNm/V^2','pNm/V^2'};
for i=1:6
    fprintf('  %-6s', rn{i});
    fprintf(' %10.3f', B(i,:)*sc(i));
    fprintf('   %s\n', un{i});
end
fprintf('\n  rank(B) at the centred pose = %d\n', rank(B));
chk(all(abs(B(1:2,5:6)) < 1e-15), 'Cp and TIA produce no lateral force, as geometry requires');
chk(all(abs(B(6,:)) < 1e-18), 'tauz is identically zero at the centred pose');
% is the deficiency structural, or only at the symmetric pose?
poseOff = [12e-6; -8e-6; 3e-6; 5e-4; -3e-4; 2e-3];
[~,Jo] = trapCapacitance(p, poseOff); [~,~,Bo] = trapWrench(p, Jo, zeros(6,1));
fprintf('  rank(B) at an off-centre pose = %d\n', rank(Bo));
fprintf('  |tauz| row off-centre: max %.3e N*m/V^2\n', max(abs(Bo(6,:))));
s = svd(B); so = svd(Bo);
fprintf('  singular values (centred) : %s\n', sprintf('%.2e ', s));
fprintf('  singular values (off-centre): %s\n', sprintf('%.2e ', so));

%% ---- 3. differential tilt actuation ----
fprintf('\n=== 3. differential drive produces tilt ===\n');
ud = uHover*[1.2 1.2 0.8 0.8 0 0]';          % quadrants 1,2 harder
[Fd,taud] = trapWrench(p, dCdq, ud);
fprintf('  u = [%.1f %.1f %.1f %.1f] V -> taux %.3f nN*m, tauy %.3f nN*m\n', ...
    ud(1:4), taud(1)*1e9, taud(2)*1e9);
chk(abs(taud(1)) > 1e-3*abs(Fd(3))*p.mass.halfW*1e-3 || abs(taud(1))>0, ...
    'unequal quadrant voltages generate roll torque');

%% ---- 4. open-loop instability ----
fprintf('\n=== 4. open-loop pole from the model itself ===\n');
d = 1e-9;
pz = pose0; pz(3) = d;
[~,Jp] = trapCapacitance(p,pz); [~,~,Bp] = trapWrench(p,Jp,zeros(6,1));
pm = pose0; pm(3) = -d;
[~,Jm] = trapCapacitance(p,pm); [~,~,Bm] = trapWrench(p,Jm,zeros(6,1));
dFzdrz = ((Bp(3,:)-Bm(3,:))*sel)*vHover/(2*d);
fprintf('  dFz/drz at hover = %+.3f N/m\n', dFzdrz);
pole = sqrt(abs(dFzdrz)/p.mass.m);
if dFzdrz > 0, kind = 'UNSTABLE'; else, kind = 'stable'; end
fprintf('  heave pole = %.1f rad/s (%.1f Hz), %s\n', pole, pole/2/pi, kind);
chk(dFzdrz > 0, 'heave stiffness is NEGATIVE (destabilising), as expected at fixed voltage');

%% ---- 5. constraint ----
fprintf('\n=== 5. checkPoseValid ===\n');
[ok0,m0] = checkPoseValid(p, pose0);
fprintf('  centred: ok=%d, margin %.2f um (should be the %.0f um gap)\n', ok0, m0*1e6, p.gap.z0*1e6);
chk(ok0 && abs(m0-p.gap.z0) < 1e-12, 'centred pose margin equals the nominal gap');

pTilt = pose0; pTilt(4) = 0.1289*pi/180;      % the established tilt limit
[okT,mT,dT] = checkPoseValid(p, pTilt);
fprintf('  at thx = 0.1289 deg: margin %.3f um, tightest = %s\n', mT*1e6, dT.worst);
chk(okT, 'the documented tilt limit is still inside the channel');

pBad = pose0; pBad(3) = 0.9*p.gap.z0; pBad(4) = 0.1289*pi/180;
[okB,mB] = checkPoseValid(p, pBad);
fprintf('  rz = 90%% of gap AND full tilt: ok=%d, margin %.3f um\n', okB, mB*1e6);
chk(~okB, 'combined translation + tilt correctly detected as contact');

pYaw = pose0; pYaw(1) = 0.9*p.gap.xy0; pYaw(6) = 0.8404*pi/180;
[okY,~,dY] = checkPoseValid(p, pYaw);
fprintf('  rx = 90%% lateral AND full yaw: ok=%d, tightest = %s\n', okY, dY.worst);
chk(~okY, 'yaw swing on top of lateral offset correctly detected');

%% ---- 6. dynamics wiring ----
fprintf('\n=== 6. trapDynamics ===\n');
x0 = [pose0; zeros(6,1)];
[xd, dg] = trapDynamics(p, x0, u);
fprintf('  at hover: acc = [%.3e %.3e %.3e] m/s^2\n', xd(7:9));
chk(norm(xd(7:12)) < 1e-6, 'hover voltage gives zero acceleration (gravity balanced)');
chk(numel(xd)==12 && all(isfinite(xd)), 'xdot is 12-by-1 and finite');
chk(dg.ok && abs(dg.margin-p.gap.z0)<1e-12, 'diagnostics report the clearance');

% free fall with the drive off
[xd0] = trapDynamics(p, x0, zeros(6,1));
fprintf('  drive off: az = %.4f m/s^2 (should be -g = %.4f)\n', xd0(9), -p.phys.g);
chk(abs(xd0(9)+p.phys.g) < 1e-9, 'with no drive the plate is in free fall');

%% ---- 7. sensing branch ----
fprintf('\n=== 7. sensing branch ===\n');
C0 = trapCapacitance(p, pose0);
Ce = trapSensing(p, C0);
fprintf('  Ceff = %.3f fF (measured 414.550 fF)\n', Ce(1)*1e15);
chk(abs(Ce(1)-414.5503e-15)/414.5503e-15 < 1e-4, 'Ceff matches the measured value at equilibrium');
chk(max(abs(Ce-mean(Ce)))/mean(Ce) < 1e-9, 'all four pairs agree at the symmetric pose');

fprintf('\n==================================================\n');
if nFail==0, fprintf('pipeline checks complete\n'); end
fprintf('==================================================\n');
