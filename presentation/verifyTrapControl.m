%% verifyTrapControl.m
% Checks the discrete-time regulator from trapControlDesign: structure,
% placement, and the two robustness questions that the design model does
% not answer by itself.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'));
p = trapDeviceGeometry();
K = trapControlDesign(p);
R = load(fullfile(root,'presentation','reportNumbers.mat')); Kstiff = R.R.K;
verdict = {'FAIL','PASS'};
chk = @(ok,msg) fprintf('  [%s] %s\n', verdict{logical(ok)+1}, msg);

%% ---- 1. controllability and observability ----
fprintf('=== 1. structure ===\n');
chk(K.ctrb.full,    sprintf('ctrb(A,B) rank %d of 6', K.ctrb.rank));
chk(K.obsv.full,    sprintf('obsv(A,Cs) rank %d of 6', K.obsv.rank));
chk(K.ctrb.fullAug, sprintf('augmented ctrb rank %d of 9 (integrators reachable)', K.ctrb.rankAug));
fprintf('  Cs rank %d of 6  -> the 4 sensors resolve 3 positions; velocities\n', rank(K.Cs));
fprintf('  come from the estimator, not from measurement.\n');

%% ---- 2. placement ----
fprintf('\n=== 2. placement ===\n');
sR = log(K.clR)/K.Tisr; sE = log(K.clE)/K.Tisr;
fprintf('  regulator: requested Re = %.1f, achieved %.1f .. %.1f rad/s\n', ...
    -K.lambdar, min(real(sR)), max(real(sR)));
fprintf('  estimator: requested Re = %.1f, achieved %.1f .. %.1f rad/s\n', ...
    -K.lambdae, min(real(sE)), max(real(sE)));
chk(max(abs(K.clR))<1 && max(abs(K.clE))<1, ...
    sprintf('both loops inside the unit circle (%.4f, %.4f)', max(abs(K.clR)), max(abs(K.clE))));

fNyq = 1/(2*K.Tisr);
fprintf('\n  Nyquist at Tisr = %.0f us is %.0f Hz\n', K.Tisr*1e6, fNyq);
fprintf('  estimator at %.0f Hz is %.0f%% of Nyquist; pole lands at z = %.3f\n', ...
    K.lambdae/2/pi, 100*(K.lambdae/2/pi)/fNyq, max(abs(K.clE)));
chk(K.lambdae/2/pi < 0.5*fNyq, 'estimator bandwidth below half Nyquist');

%% ---- 3. does it stabilise the REAL plant? ----
% The design model treats u as force, so it carries no electrostatic
% stiffness. The real device at fixed voltage has NEGATIVE stiffness. If
% the force-to-voltage inversion is exact and uses the current pose, that
% stiffness cancels and the design model is right. Test both extremes.
fprintf('\n=== 3. closed loop against the real plant ===\n');
Areal = K.A;
Areal(2,1) = -Kstiff(3,3)/p.mass.m;          % heave:  +37.76/m
Areal(4,3) = -Kstiff(4,4)/p.mass.I(1,1);     % tilt x
Areal(6,5) = -Kstiff(5,5)/p.mass.I(2,2);     % tilt y
ev0=sort(real(eig(Areal)),'descend'); fprintf('  open-loop poles with stiffness: %s rad/s\n', sprintf('%+.1f ',ev0));

for mode = 1:2
    if mode==1
        Au = K.A;  nm = 'stiffness CANCELLED by the voltage inversion';
    else
        Au = Areal; nm = 'stiffness NOT cancelled';
    end
    sysd = c2d(ss(Au,K.B,K.Cs,0), K.Tisr, 'zoh');
    Ad=sysd.A; Bd=sysd.B; Cd=sysd.C;
    % observer-based output feedback with integral action
    nx=6; ni=3;
    Acl = [Ad - Bd*K.Kx,            -Bd*K.Ki,        Bd*K.Kx;
           -K.Tisr*K.Cc,             eye(ni),        zeros(ni,nx);
           zeros(nx,nx),             zeros(nx,ni),   Ad - K.L*Cd];
    ev = eig(Acl);
    lab = {'UNSTABLE','STABLE'};
    fprintf('  %-45s max|z| = %.4f  %s\n', nm, max(abs(ev)), lab{(max(abs(ev))<1)+1});
end

%% ---- 4. delay, which the design ignores entirely ----
fprintf('\n=== 4. loop delay vs the unstable pole ===\n');
wu = sqrt(abs(Kstiff(3,3))/p.mass.m);
tSense = p.readout.tI2C + p.readout.dtChannel*3;
tTotal = tSense + K.Tisr;
fprintf('  sensing delay   %.0f us\n', tSense*1e6);
fprintf('  + one ISR       %.0f us\n', K.Tisr*1e6);
fprintf('  = loop delay    %.0f us\n', tTotal*1e6);
fprintf('  wu*td = %.2f  (needs to be well below 1)\n', wu*tTotal);
chk(wu*tTotal < 1, sprintf('loop delay within the unstable-pole budget (%.0f us)', 1e6/wu));

%% ---- 5. estimator noise gain ----
fprintf('\n=== 5. estimator noise amplification ===\n');
sigC = 1e-15;                                 % 1 fF of capacitance noise
dCdz = 19.7e-15/1e-6;                         % fF per m, at the operating point
sigZ = sigC/dCdz;
fprintf('  1 fF sensor noise -> %.1f nm height noise\n', sigZ*1e9);
nz = K.L*ones(4,1)*sigZ;
fprintf('  propagated to the state estimate per ISR:\n');
lbl={'rz','drz','thx','dthx','thy','dthy'};
un={'nm','nm/s','nrad','nrad/s','nrad','nrad/s'};
sc=[1e9 1e9 1e9 1e9 1e9 1e9];
for i=1:6, fprintf('    %-6s %10.2f %s\n', lbl{i}, abs(nz(i))*sc(i), un{i}); end
fprintf('\n');

