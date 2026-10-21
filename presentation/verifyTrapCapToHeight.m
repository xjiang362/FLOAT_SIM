%% verifyTrapCapToHeight.m
% Round-trips the sensing inversion and quantifies its two error sources.
%
%   pose -> trapCapacitance -> trapSensing -> trapCapToHeight -> pose
%
% An exact round trip would only prove the algebra, since forward and
% inverse share constants. The informative tests are the two that can fail:
% the lateral cross-sensitivity (the inverse assumes height alone, the
% forward model does not), and inverting the MEASURED COMSOL capacitance.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'estimator'), fullfile(root,'presentation'));
p = trapDeviceGeometry('z0', 15e-6);   % COMSOL data was taken at the 15 um gap
verdict = {'FAIL','PASS'};
chk = @(ok,msg) fprintf('  [%s] %s\n', verdict{logical(ok)+1}, msg);

%% ---- 1. round trip with NO lateral motion ----
fprintf('=== 1. round trip, vertical and tilt only ===\n');
poses = { zeros(6,1), [0;0;4e-6;0;0;0], [0;0;-6e-6;8e-4;0;0], ...
          [0;0;2e-6;5e-4;-7e-4;0], [0;0;-9e-6;-1e-3;6e-4;0] };
names = {'centred','rz=+4um','rz=-6um, thx','rz, thx, thy','rz=-9um, both tilts'};
worstZ = 0;
for k = 1:numel(poses)
    q = poses{k};
    [zc,pose3,info] = trapCapToHeight(p, trapSensing(p, trapCapacitance(p,q)));
    zTrue = q(3) + q(4)*info.point(:,2) - q(5)*info.point(:,1);
    eZ = max(abs(zc - zTrue)); worstZ = max(worstZ,eZ);
    fprintf('  %-20s height err %8.2e m   rz err %8.2e m\n', names{k}, eZ, abs(pose3(1)-q(3)));
end
chk(worstZ < 1e-8, sprintf('heights exact to %.1e m when the plate is laterally centred', worstZ));

%% ---- 2. lateral cross-sensitivity (the real limitation) ----
fprintf('\n=== 2. how much does lateral motion corrupt the height? ===\n');
base = [0;0;3e-6;2e-4;0;0];
z0 = trapCapToHeight(p, trapSensing(p, trapCapacitance(p,base)));
fprintf('  %10s %14s %14s\n','rx [um]','max dz [nm]','dz/drx [nm/um]');
for rxv = [5 10 20 30 45]*1e-6
    q = base; q(1) = rxv;
    z1 = trapCapToHeight(p, trapSensing(p, trapCapacitance(p,q)));
    fprintf('  %10.0f %14.1f %14.2f\n', rxv*1e6, max(abs(z1-z0))*1e9, max(abs(z1-z0))/rxv*1e3);
end
q = base; q(6) = 2e-3;
zy = trapCapToHeight(p, trapSensing(p, trapCapacitance(p,q)));
fprintf('  yaw 2 mrad alone: max dz = %.1f nm\n', max(abs(zy-z0))*1e9);
fprintf('  -> lateral motion is NOT observable, but it LEAKS into the height\n');
fprintf('     estimate through the overlap term. This is a real cross-axis\n');
fprintf('     error, not a numerical artifact.\n');

%% ---- 3. inverting the MEASURED capacitance ----
fprintf('\n=== 3. inverting COMSOL data, rz block ===\n');
D = load(fullfile(root,'presentation','capMatrixSweep63_parsed.mat'));
Cm = D.Cm; RES = 21;
bad = false(D.nPose,1);
for q=1:D.nPose, M=Cm(:,:,q); o=M-diag(diag(M)); bad(q)=any(o(:)<0); end
b = D.blk{2}; bg = b(~bad(b));
[~,ord] = sort(D.rz(bg)); bo = bg(ord);
rzT = zeros(numel(bo),1); rzE = zeros(numel(bo),1);
for ii = 1:numel(bo)
    q = bo(ii); Ce = zeros(4,1);
    for pr = 1:4
        Ce(pr) = effectiveMutualCapacitance(Cm(:,:,q), p.idx.sensePair(pr,1), ...
                                            p.idx.sensePair(pr,2), RES);
    end
    [~,pose3] = trapCapToHeight(p, Ce);
    rzT(ii) = D.rz(q); rzE(ii) = pose3(1);
end
err = rzE - rzT;
gain = rzT\rzE;
fprintf('  RMS error, full block      : %.3f um\n', rms(err)*1e6);
fprintf('  RMS error, |rz| <= 5 um    : %.3f um\n', rms(err(abs(rzT)<=5e-6))*1e6);
fprintf('  error at rz = 0            : %.4f um\n', interp1(rzT,err,0)*1e6);
fprintf('  best-fit scale factor      : %.4f  (1.0 would be unbiased)\n', gain);
fprintf('  -> the residual is a SCALE error, not an offset: exact at the\n');
fprintf('     calibration pose, growing to %.0f%% at the extremes.\n', ...
    100*max(abs(err./max(abs(rzT),1e-9))));
chk(rms(err(abs(rzT)<=5e-6)) < 0.5e-6, 'sub-half-micron recovery within +/-5 um');

%% ---- 4. conditioning of the plane fit ----
fprintf('\n=== 4. plane-fit conditioning ===\n');
[~,~,info] = trapCapToHeight(p, trapSensing(p, trapCapacitance(p,base)));
M = [ones(4,1), info.point(:,2), -info.point(:,1)];
Mn = M ./ vecnorm(M);                       % scale columns: cond must not
fprintf('  sensing points: %s\n', sprintf('(%+.0f,%+.0f) ', info.point'*1e6));
fprintf('  cond(M) raw            : %.1f   <- dominated by unit mismatch\n', cond(M));
fprintf('  cond(M) column-scaled  : %.4f  <- the real geometric conditioning\n', cond(Mn));
chk(cond(Mn) < 2, 'four corners of a square give an essentially perfect fit');
fprintf('\n');
