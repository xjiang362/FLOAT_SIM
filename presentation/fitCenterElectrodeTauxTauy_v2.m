%% fitCenterElectrodeTauxTauy_v2.m
% The R^2~0.6 ceiling was NOT a bad fitting equation and NOT noise -- the
% two source datasets are on different scales. diagnoseCenterElectrode-
% DataIntegrity.m showed the Fz model (fitted+validated on run 1 to
% R^2=0.998) predicts run 2 at R^2=-19.2 with a median ratio of 0.276,
% i.e. run 2's fields are uniformly weaker (a voltage/excitation
% difference: both Fz and tau scale as V^2, so a single factor covers
% both). Stitching rz=0 (run 1) to rz=+-7um (run 2) therefore fitted the
% run-to-run scale jump, not slope(rz).
%
% Fix: measure the scale factor directly where the two runs OVERLAP in
% pose -- run 1's rz-only block passes through rz=+-7um at zero tilt, and
% run 2's brackets sit at rz=+-7um -- then rescale run 2 onto run 1's
% (validated, V=30) footing and refit.
%
% Model: for a symmetric centered electrode, expanding
%   tau_x = INT y * eps*V^2/(2*(g - y*thx)^2) dA
% over a symmetric span gives tau_x = (4/3)*b^3*eps*V^2*thx / (2*g^3),
% i.e. tau_x = kappa * V^2 * thx / g^3 -- a pure 1/g^3 law with ONE
% coefficient, no a1/a2 pair and no non-monotonic shape.
% Equivalently in the project's established moment-arm language:
%   tau_x = Leq * Fz * thx  with Leq = (kappa/k) / g,
% so the effective lever arm grows as 1/gap. The exponent n is fitted
% rather than assumed, to verify the n=3 derivation against the data.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir,'centerElectrodeModel.mat'),'k_fit','a2','a4','a5','g0_z');

V = 30; L = 6000e-6;
r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
minGap = @(rz,thxD,thyD) (g0_z-rz) - (L/2)*(abs(sind(thxD))+abs(sind(thyD)));

raw1 = readmatrix(fullfile(scriptDir,'comsol_centerElectrodeSweep_rzThxThy150_results.txt'),'NumHeaderLines',5);
raw2 = readmatrix(fullfile(scriptDir,'comsol_centerElectrode_rzCrossTerm200_results.txt'),'NumHeaderLines',5);
rz1=raw1(:,3); thx1=raw1(:,4); thy1=raw1(:,5); Fz1=raw1(:,9); Taux1=raw1(:,10); Tauy1=raw1(:,11);
rz2=raw2(:,3); thx2=raw2(:,4); thy2=raw2(:,5); Fz2=raw2(:,9); Taux2=raw2(:,10); Tauy2=raw2(:,11);

blkRz = 1:50; blkThx0 = 51:100; blkThy0 = 101:150;

%% ===== STEP 1: measure the run1->run2 scale factor at overlapping poses =====
% Run 1's rz-only block is at thx=thy=0, so interpolate its Fz(rz) at the
% bracket values. Run 2's Fz varies only ~6% with tilt, so extrapolate
% each bracket's Fz to thx=0 with a quadratic (tilt enters as thx^2 for
% this symmetric electrode -- already established in the Fz fit).
[rzS,ord] = sort(rz1(blkRz)); FzS = Fz1(blkRz); FzS = FzS(ord);

fprintf('===== STEP 1: cross-run scale factor at matching poses =====\n');
fprintf('%-10s %14s %14s %10s\n','rz [um]','Fz run1 @th=0','Fz run2 @th=0','ratio');
scales = zeros(2,1); rzBr = [-7e-6, 7e-6];
for i = 1:2
    m = abs(rz2(1:100) - rzBr(i)) < 1e-12;      % thx-block only
    th_r = thx2(m)*pi/180;
    c = [ones(sum(m),1), th_r.^2] \ Fz2(m);     % Fz = A + B*thx^2 -> A is Fz at zero tilt
    Fz2_at0 = c(1);
    Fz1_at  = interp1(rzS, FzS, rzBr(i), 'pchip');
    scales(i) = Fz2_at0 / Fz1_at;
    fprintf('%-10.1f %14.5e %14.5e %10.4f\n', rzBr(i)*1e6, Fz1_at, Fz2_at0, scales(i));
end
scale = mean(scales);
fprintf('mean scale = %.4f   (spread %.2f%% -> %s)\n', scale, ...
    100*abs(diff(scales))/scale, ...
    string(abs(diff(scales))/scale < 0.05));
fprintf('implied voltage ratio sqrt(scale) = %.4f  => V_run2 ~ %.2f V if V_run1 = 30 V\n\n', ...
    sqrt(scale), 30*sqrt(scale));

%% ===== STEP 2: rescale run 2 onto run 1's footing, collect clean slopes =====
% tau scales with V^2 exactly as Fz does, so the same factor applies.
fprintf('===== STEP 2: per-slice slopes, run 2 rescaled onto run 1 =====\n');
rzPts = [-7e-6; 0; 7e-6];
slopeX = zeros(3,1); slopeY = zeros(3,1); nSafe = zeros(3,1);

% rz = 0 from run 1 (native scale, geometrically safe throughout)
slopeX(2) = [thx1(blkThx0)*pi/180] \ Taux1(blkThx0);   % through origin: symmetric electrode has no baseline
slopeY(2) = [thy1(blkThy0)*pi/180] \ Tauy1(blkThy0);
nSafe(2)  = numel(blkThx0);

for i = [1 3]
    rzv = rzPts(i);
    mX = abs(rz2(1:100)-rzv) < 1e-12;
    idxX = find(mX);
    safeX = minGap(rzv, thx2(idxX), 0) > 4e-6;          % drop near-collision pts
    idxX = idxX(safeX);
    slopeX(i) = (thx2(idxX)*pi/180) \ (Taux2(idxX)/scale);

    mY = false(200,1); mY(101:200) = abs(rz2(101:200)-rzv) < 1e-12;
    idxY = find(mY);
    safeY = minGap(rzv, 0, thy2(idxY)) > 4e-6;
    idxY = idxY(safeY);
    slopeY(i) = (thy2(idxY)*pi/180) \ (Tauy2(idxY)/scale);
    nSafe(i) = numel(idxX);
end

gapN = g0_z - rzPts;
fprintf('%-10s %9s %15s %15s %7s\n','rz [um]','gap[um]','slopeX [Nm/rad]','slopeY [Nm/rad]','n_safe');
for i=1:3
    fprintf('%-10.1f %9.2f %15.5e %15.5e %7d\n', rzPts(i)*1e6, gapN(i)*1e6, slopeX(i), slopeY(i), nSafe(i));
end
fprintf('\nMonotonic increasing with rz now? %s\n', string(all(diff(slopeX)>0)));

%% ===== STEP 3: fit the power law slope = kappa*V^2 / g^n =====
fprintf('\n===== STEP 3: power-law fit  slope = kappa*V^2 / g^n =====\n');
cx = [ones(3,1), log(gapN)] \ log(slopeX);
cy = [ones(3,1), log(gapN)] \ log(slopeY);
nX = -cx(2); nY = -cy(2);
kappaX = exp(cx(1))/V^2; kappaY = exp(cy(1))/V^2;
fprintf('Taux: n = %.4f   (analytic prediction: 3)   kappa = %.6e\n', nX, kappaX);
fprintf('Tauy: n = %.4f   (analytic prediction: 3)   kappa = %.6e\n', nY, kappaY);
fprintf('R^2 of the power law on the 3 slopes: X=%.4f, Y=%.4f\n', ...
    r2fun(log(slopeX), [ones(3,1),log(gapN)]*cx), r2fun(log(slopeY), [ones(3,1),log(gapN)]*cy));

% Locked-exponent version (n=3 from theory), single coefficient
kx3 = (V^2./gapN.^3) \ slopeX;
ky3 = (V^2./gapN.^3) \ slopeY;
fprintf('\nWith n LOCKED at the analytic 3:\n');
fprintf('  kappa_x = %.6e, kappa_y = %.6e  (x<->y symmetry: %.2f%% apart)\n', ...
    kx3, ky3, 100*abs(kx3-ky3)/mean([kx3 ky3]));

%% ===== STEP 4: evaluate the final model against EVERY raw point =====
% Final model (run-1 / V=30 footing):
%   Taux = kappa * V^2 * thx_rad / (g0 - rz)^n
Taux_fun = @(rz,thxr,kap,n) kap*V^2 .* thxr ./ (g0_z - rz).^n;
Tauy_fun = @(rz,thyr,kap,n) kap*V^2 .* thyr ./ (g0_z - rz).^n;

fprintf('\n===== STEP 4: model vs ALL raw points (run 2 rescaled) =====\n');
fprintf('%-26s %10s %10s\n','dataset','R2 (n fit)','R2 (n=3)');
evalSet = @(rzv,thv,tauv,fun,kap,nn) r2fun(tauv, fun(rzv, thv*pi/180, kap, nn));

% run 1, rz=0 tilt block
fprintf('%-26s %10.4f %10.4f\n','run1 thx-only (rz=0)', ...
    evalSet(rz1(blkThx0),thx1(blkThx0),Taux1(blkThx0),Taux_fun,kappaX,nX), ...
    evalSet(rz1(blkThx0),thx1(blkThx0),Taux1(blkThx0),Taux_fun,kx3,3));
fprintf('%-26s %10.4f %10.4f\n','run1 thy-only (rz=0)', ...
    evalSet(rz1(blkThy0),thy1(blkThy0),Tauy1(blkThy0),Tauy_fun,kappaY,nY), ...
    evalSet(rz1(blkThy0),thy1(blkThy0),Tauy1(blkThy0),Tauy_fun,ky3,3));

% run 2 brackets, rescaled, safe points only
for i=[1 3]
    rzv = rzPts(i);
    idxX = find(abs(rz2(1:100)-rzv)<1e-12);
    idxX = idxX(minGap(rzv,thx2(idxX),0) > 4e-6);
    fprintf('%-26s %10.4f %10.4f\n', sprintf('run2 thx rz=%+.0fum',rzv*1e6), ...
        evalSet(rz2(idxX),thx2(idxX),Taux2(idxX)/scale,Taux_fun,kappaX,nX), ...
        evalSet(rz2(idxX),thx2(idxX),Taux2(idxX)/scale,Taux_fun,kx3,3));
end

% Pooled
rzPool = [rz1(blkThx0); rz2(1:100)];
thPool = [thx1(blkThx0); thx2(1:100)];
tauPool= [Taux1(blkThx0); Taux2(1:100)/scale];
safePool = minGap(rzPool, thPool, 0) > 4e-6;
fprintf('%-26s %10.4f %10.4f\n','POOLED Taux (safe pts)', ...
    r2fun(tauPool(safePool), Taux_fun(rzPool(safePool),thPool(safePool)*pi/180,kappaX,nX)), ...
    r2fun(tauPool(safePool), Taux_fun(rzPool(safePool),thPool(safePool)*pi/180,kx3,3)));

save(fullfile(scriptDir,'centerElectrodeTauxTauy_v2.mat'), ...
    'kappaX','kappaY','nX','nY','kx3','ky3','scale','g0_z','V');

%% ===== Plot =====
figure('Color','w','Position',[60 60 1300 520]);
cols = [0.80 0.25 0.20; 0.20 0.30 0.75; 0.15 0.60 0.30];

subplot(1,2,1); hold on; grid on;
for i=1:3
    rzv = rzPts(i);
    if i==2
        th = thx1(blkThx0); tau = Taux1(blkThx0);
    else
        idx = find(abs(rz2(1:100)-rzv)<1e-12);
        idx = idx(minGap(rzv,thx2(idx),0) > 4e-6);
        th = thx2(idx); tau = Taux2(idx)/scale;
    end
    scatter(th, tau, 26, cols(i,:), 'filled', 'DisplayName', sprintf('data rz=%+.0f\\mum',rzv*1e6));
    tg = linspace(min(th),max(th),100)';
    plot(tg, Taux_fun(rzv, tg*pi/180, kx3, 3), '-', 'Color',cols(i,:), 'LineWidth',1.7, 'HandleVisibility','off');
end
xlabel('\theta_x [deg]'); ylabel('Taux [N\cdotm]  (run-1 / V=30 footing)');
title('Taux: data (dots) vs \kappa V^2\theta_x/g^3 (lines)');
legend('Location','best','FontSize',8);

subplot(1,2,2); hold on; grid on;
gg = linspace(7e-6, 24e-6, 200)';
loglog(gapN*1e6, slopeX, 'o', 'MarkerSize',10,'MarkerFaceColor',[0.2 0.3 0.75],'DisplayName','measured slope');
plot(gg*1e6, kx3*V^2./gg.^3, 'k-', 'LineWidth',1.7, 'DisplayName','\kappa V^2/g^3 (n locked = 3)');
plot(gg*1e6, kappaX*V^2./gg.^nX, 'r--', 'LineWidth',1.4, 'DisplayName',sprintf('fitted n = %.2f',nX));
set(gca,'XScale','log','YScale','log'); grid on;
xlabel('gap g = g_0 - rz [\mum]'); ylabel('d(Taux)/d\theta_x [N\cdotm/rad]');
title('Slope vs gap: monotonic 1/g^3 once runs are on one scale');
legend('Location','best','FontSize',8);

sgtitle('Center electrode torque model v2 (cross-run scale corrected)');
exportgraphics(gcf, fullfile(scriptDir,'CenterElectrode_TauxTauy_v2Fit.png'));
fprintf('\nSaved plot to CenterElectrode_TauxTauy_v2Fit.png\n');
