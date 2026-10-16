%% plotCeffVsRz_3models.m
% Ceff (between terminals 1 and 2) vs rz, three datasets:
%   1. Ground truth: effectiveMutualCapacitance.m, eliminating all
%      genuinely floating terminals (5,6,9,10,13,14,21 -- resonator is
%      21; bottom Si, terminal 22, is fixed-voltage/grounded and
%      dropped automatically since it's not in the floating list).
%   2. Resonator-only analytical approximation, now using the FITTED
%      physics-informed models for S_R, Cm1p, Cm2p (fitSR_physicsModel.m,
%      fitCm1p_fringingModel.m, fitCm2p_fringingModel.m) instead of raw
%      simulation lookups for those three quantities. Cm12, S1, S2 are
%      not yet modeled, so they still come from the raw simulation data.
%   3. Empirical fit k1 + k2/(k3+rz), fit to the ground-truth data.
% R^2 for models 2 and 3 are both computed against the ground truth.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'capSweep300_withGND_allMatrices.mat'));  % Cm_all, rx,ry,rz,thx,thy,thz, nPose
load(fullfile(scriptDir, 'SR_physicsModel.mat'));      % pfit, modelFun, d0,x0,y0 (S_R model)
SR_model = @(p,X) modelFun(p,X);
SR_pfit = pfit; SR_fun = modelFun;
load(fullfile(scriptDir, 'Cm1p_fringingFit.mat'));     % pfit, modelFun (Cm1p model)
Cm1p_pfit = pfit; Cm1p_fun = modelFun;
load(fullfile(scriptDir, 'Cm2p_fringingFit.mat'));     % pfit, modelFun (Cm2p model)
Cm2p_pfit = pfit; Cm2p_fun = modelFun;

thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;

%% Dataset 1: ground truth via effectiveMutualCapacitance.m
floatingTerms = [5,6,9,10,13,14,21];
Ceff_gt = zeros(nPose,1);
for p = 1:nPose
    Ceff_gt(p) = effectiveMutualCapacitance(Cm_all(:,:,p), 1, 2, floatingTerms);
end

%% Dataset 2: resonator-only analytical approximation, fitted S_R/Cm1p/Cm2p
Cm12 = squeeze(Cm_all(1,2,:));
S1 = squeeze(sum(Cm_all(1,:,:),2));
S2 = squeeze(sum(Cm_all(2,:,:),2));

Cm1R = Cm1p_fun(Cm1p_pfit, [rx, ry, rz, thx_rad, thy_rad, thz_rad]);
Cm2R = Cm2p_fun(Cm2p_pfit, [rx, ry, rz, thx_rad, thy_rad, thz_rad]);
SR = SR_fun(SR_pfit, [rx, ry, rz]);

C_M2 = Cm12 + Cm1R.*Cm2R./SR;
C1_2 = S1 - Cm12 - Cm1R.*(Cm1R+Cm2R)./SR;
C2_2 = S2 - Cm12 - Cm2R.*(Cm1R+Cm2R)./SR;
Ceff_approx = C_M2 + C1_2.*C2_2./(C1_2+C2_2);

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
r2_approx = r2fun(Ceff_gt, Ceff_approx);
fprintf('Dataset 2 (resonator-only analytical, fitted S_R/Cm1p/Cm2p) R^2 vs ground truth: %.4f\n', r2_approx);

%% Dataset 3: two-point calibration of Cp + k/(d0-rz), d0 fixed at 15um
% d0=15um matches the established nominal gap (g0_z) used throughout
% this project's other physics-informed models. With d0 fixed, the
% model is LINEAR in (Cp,k), so it can be solved exactly from the two
% endpoint values (rz=-14um and rz=+14um) -- no nonlinear optimizer
% needed at all.
d0 = 15e-6;
[~, iMin] = min(rz);
[~, iMax] = max(rz);
rz1 = rz(iMin); rz2 = rz(iMax);
Ceff1 = Ceff_gt(iMin); Ceff2 = Ceff_gt(iMax);

A1 = 1/(d0-rz1); A2 = 1/(d0-rz2);
k = (Ceff2-Ceff1)/(A2-A1);
Cp = Ceff1 - k*A1;

modelFun = @(z) Cp + k./(d0-z);
Ceff_fit = modelFun(rz);
r2_fit = r2fun(Ceff_gt, Ceff_fit);
fprintf('Dataset 3 (two-point calibrated Cp+k/(d0-rz), d0=%.3gum fixed) R^2 vs ground truth: %.4f\n', d0*1e6, r2_fit);
fprintf('Calibration points: rz=%.3gum -> Ceff=%.4e;  rz=%.3gum -> Ceff=%.4e\n', rz1*1e6, Ceff1, rz2*1e6, Ceff2);
fprintf('Cp=%.4e, k=%.4e\n', Cp, k);

%% Plot
[rzSorted, order] = sort(rz);
figure('Color','w','Name','Ceff vs rz: ground truth, analytical approx, empirical fit');
hold on; grid on;
scatter(rz*1e6, Ceff_gt, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Ground truth (effectiveMutualCapacitance)');
scatter(rz*1e6, Ceff_approx, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', sprintf('Resonator-only analytical, fitted models (R^2=%.4f)', r2_approx));
plot(rzSorted*1e6, Ceff_fit(order), 'b-', 'LineWidth', 2, 'DisplayName', sprintf('C_p+k/(d_0-z), d_0=15\\mum fixed (R^2=%.4f)', r2_fit));
xlabel('rz [\mum]'); ylabel('C_{eff} (terminals 1-2) [F]');
title('Effective capacitance between terminals 1 and 2 vs rz');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'Ceff_vs_rz_3models.png'));
fprintf('\nSaved plot to Ceff_vs_rz_3models.png\n');
