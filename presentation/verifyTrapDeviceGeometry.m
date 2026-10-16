%% verifyTrapDeviceGeometry.m
% Checks config/trapDeviceGeometry.m against everything established from
% the COMSOL capacitance matrix. Run this after ANY edit to the geometry
% struct -- a wrong terminal assignment or centroid is silent otherwise:
% the model still runs, it just fits badly.
%
% Checks performed
%   1. internal consistency : dX*dY == A for every footprint
%   2. terminal assignment  : the analytic C_iR must fit the measured
%                             Cm(i,21) at R^2 >= 0.99 for all 20 electrodes
%   3. Cp orientation       : +x electrodes must couple more to terminal 18
%   4. frame decomposition  : 4 sides + bottom must predict Cm(21,22)
%   5. mass properties      : finite, positive, consistent inertia
%
% Reference data: comsol_capMatrixSweep_rxRzThx63_results.txt via
% capMatrixSweep63_parsed.mat. Poses whose matrix has negative
% off-diagonals (mesh error in the thin gap) are excluded throughout.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir, fullfile(fileparts(scriptDir),'config'));
eps0 = 8.8541878128e-12;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

geo = trapDeviceGeometry();
D = load(fullfile(scriptDir,'capMatrixSweep63_parsed.mat'));
Cm = D.Cm; nPose = D.nPose; RES = geo.idx.resonator; FR = geo.idx.frame;
bad = false(nPose,1);
for p=1:nPose, M=Cm(:,:,p); o=M-diag(diag(M)); bad(p)=any(o(:)<0); end
good = ~bad;
rz = D.rz; thx = D.thx*pi/180; thy = D.thy*pi/180;

nFail = 0;
verdict = {'FAIL','PASS'};
chk = @(ok,msg) fprintf('  [%s] %s\n', verdict{logical(ok)+1}, msg);

%% ---- 1. dX*dY == A ----
fprintf('=== 1. internal consistency: dX*dY == A ===\n');
worst = 0;
for i = geo.idx.electrodes
    e = geo.term(i);
    rel = abs(e.dX*e.dY - e.A)/e.A;
    worst = max(worst, rel);
end
ok = worst < 1e-12; nFail = nFail + ~ok;
chk(ok, sprintf('worst relative mismatch %.2e over 20 electrodes', worst));

%% ---- 2. terminal assignment via the analytic C_iR ----
fprintf('\n=== 2. terminal assignment: analytic C_iR vs measured ===\n');
fprintf('  %4s %-9s %9s %9s %9s\n','term','name','A[mm^2]','k_i','R^2');
R2 = zeros(20,1); kAll = zeros(20,1);
for i = geo.idx.electrodes
    e = geo.term(i);
    g = geo.gap.z0 - rz - thx*e.Yc + thy*e.Xc;
    Cmod = eps0*e.A./g;
    meas = squeeze(Cm(i,RES,:));
    kAll(i) = Cmod(good)\meas(good);
    R2(i)   = r2fun(meas(good), Cmod(good)*kAll(i));
    fprintf('  %4d %-9s %9.4f %9.4f %9.5f\n', i, e.name, e.A*1e6, kAll(i), R2(i));
end
ok = min(R2) >= 0.99; nFail = nFail + ~ok;
chk(ok, sprintf('min R^2 = %.5f across all 20 electrodes (threshold 0.99)', min(R2)));
ok = all(kAll > 0.98 & kAll < 1.10); nFail = nFail + ~ok;
chk(ok, sprintf('gains in [%.4f, %.4f], all near 1 as fringing requires', min(kAll), max(kAll)));

% The two-parameter empirical fringing law k_i = a + b*(P/A) was removed
% from the geometry when the model was reduced to pure parallel plate, so
% there is nothing left to check it against: geo.fringe now carries only
% the conformal-mapping edge extension deltaOverGap = 2/pi. What matters
% is that the UNCORRECTED gains above are already near 1, which the
% previous check asserts. Report the residual spread for information.
PA = arrayfun(@(i) geo.term(i).P/geo.term(i).A, geo.idx.electrodes)';
fprintf('  no fringing correction applied (geo.fringe.deltaOverGap = %.4f = 2/pi)\n', ...
    geo.fringe.deltaOverGap);
fprintf('  residual gain spread %.2f%% peak-to-peak; correlation with P/A = %+.3f\n', ...
    100*(max(kAll)-min(kAll)), corr(PA, kAll));

%% ---- 3. Cp orientation ----
fprintf('\n=== 3. Cp orientation (17 = -x half, 18 = +x half) ===\n');
k0 = 11;                                  % equilibrium pose
plusX  = [1 2 3];                         % electrodes at +x
minusX = [5 6 7];                         % their quadrant-2 counterparts
c17p = mean(abs(Cm(plusX,17,k0)));  c18p = mean(abs(Cm(plusX,18,k0)));
c17m = mean(abs(Cm(minusX,17,k0))); c18m = mean(abs(Cm(minusX,18,k0)));
fprintf('  +x electrodes: to 17 = %.3f fF, to 18 = %.3f fF\n', c17p*1e15, c18p*1e15);
fprintf('  -x electrodes: to 17 = %.3f fF, to 18 = %.3f fF\n', c17m*1e15, c18m*1e15);
ok = (c18p > c17p) && (c17m > c18m); nFail = nFail + ~ok;
chk(ok, 'coupling pattern matches the stored orientation');
ok = geo.term(17).Xc < 0 && geo.term(18).Xc > 0; nFail = nFail + ~ok;
chk(ok, sprintf('stored centroids: Xc(17) = %+.1f um, Xc(18) = %+.1f um', ...
    geo.term(17).Xc*1e6, geo.term(18).Xc*1e6));

%% ---- 4. frame decomposition ----
fprintf('\n=== 4. frame decomposition (4 sides + bottom) ===\n');
Cpred = eps0*geo.frame.bottom.A/geo.frame.bottom.gap;
for kS=1:4, Cpred = Cpred + eps0*geo.frame.side(kS).A/geo.frame.side(kS).gap; end
Cmeas = Cm(RES,FR,k0);
fprintf('  predicted %.4f pF, measured %.4f pF, ratio %.4f\n', ...
    Cpred*1e12, Cmeas*1e12, Cpred/Cmeas);
ok = abs(Cpred/Cmeas - 1) < 0.05; nFail = nFail + ~ok;
chk(ok, 'parameter-free prediction within 5% of measurement');

%% ---- 5. mass properties ----
fprintf('\n=== 5. mass properties ===\n');
m = geo.mass;
fprintf('  m = %.4e kg (%.2f mg), Ixx = %.4e kg m^2\n', m.m, m.m*1e6, m.I(1,1));
ok = m.m > 0 && all(diag(m.I) > 0) && abs(m.I(3,3)-(m.I(1,1)+m.I(2,2)-m.m/6*m.t^2)) < 1e-18;
nFail = nFail + ~ok;
chk(ok, 'mass positive and inertia satisfies the plate perpendicular-axis relation');

%% ---- summary ----
fprintf('\n==================================================\n');
if nFail == 0
    fprintf('ALL CHECKS PASSED\n');
else
    fprintf('%d CHECK(S) FAILED -- do not use this geometry\n', nFail);
end
fprintf('==================================================\n');
