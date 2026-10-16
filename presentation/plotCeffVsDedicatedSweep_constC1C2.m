%% plotCeffVsDedicatedSweep_constC1C2.m
% Fully closed-form Ceff model, no per-pose simulation lookup except for
% Cm12 (not yet modeled) and the one-time equilibrium calibration:
%   - Cm1R(pose), Cm2R(pose), S_R(pose): FITTED models only (no raw
%     lookup at all -- Cm1p_fringingFit.mat, Cm2p_fringingFit.mat,
%     SR_physicsModel.mat).
%   - C1, C2: held CONSTANT at their equilibrium (near rx=ry=rz=thx=
%     thy=thz=0) values, computed once from raw simulation data and then
%     frozen -- not re-looked-up per pose. Physically reasonable since
%     terminals 1 and 2 are fixed (only the resonator moves), so their
%     capacitance to everything else besides the resonator/each other
%     shouldn't depend much on the resonator's pose.
%   - Cm12(pose): still raw simulation lookup (no model exists for it).
%
% Tested on the dedicated rz/thx/thy-only sweep (sensingSweep60), using
% its own near-equilibrium point (thx-only block, smallest |thx|) for
% the one-time C1/C2 calibration.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'sensingSweep60_allMatrices.mat'));   % Cm_all, rx,ry,rz,thx,thy,thz, nPose
load(fullfile(scriptDir, 'SR_physicsModel.mat'));      % pfit, modelFun
SR_pfit = pfit; SR_fun = modelFun;
load(fullfile(scriptDir, 'Cm1p_fringingFit.mat'));     % pfit, modelFun
Cm1p_pfit = pfit; Cm1p_fun = modelFun;
load(fullfile(scriptDir, 'Cm2p_fringingFit.mat'));     % pfit, modelFun
Cm2p_pfit = pfit; Cm2p_fun = modelFun;

thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;

%% Ground truth
floatingTerms = [5,6,9,10,13,14,21];
Ceff_gt = zeros(nPose,1);
for pp = 1:nPose
    Ceff_gt(pp) = effectiveMutualCapacitance(Cm_all(:,:,pp), 1, 2, floatingTerms);
end

%% One-time equilibrium calibration of C1, C2
[~, ieq] = min(rx.^2 + ry.^2 + rz.^2 + thx_rad.^2 + thy_rad.^2 + thz_rad.^2);
fprintf('Equilibrium calibration pose (idx %d): rx=%.3e ry=%.3e rz=%.3e thx=%.4gdeg thy=%.4gdeg thz=%.4gdeg\n', ...
    ieq, rx(ieq), ry(ieq), rz(ieq), thx(ieq), thy(ieq), thz(ieq));

Cm12_eq = Cm_all(1,2,ieq);
Cm1R_eq = Cm_all(1,21,ieq);
Cm2R_eq = Cm_all(2,21,ieq);
S1_eq = sum(Cm_all(1,:,ieq));
S2_eq = sum(Cm_all(2,:,ieq));
SR_eq = sum(Cm_all(21,:,ieq));

C1_const = S1_eq - Cm12_eq - Cm1R_eq*(Cm1R_eq+Cm2R_eq)/SR_eq;
C2_const = S2_eq - Cm12_eq - Cm2R_eq*(Cm1R_eq+Cm2R_eq)/SR_eq;
fprintf('Calibrated constants: C1 = %.6e F, C2 = %.6e F (held fixed for all poses)\n', C1_const, C2_const);

%% New model: Cm1R, Cm2R, S_R from FITTED models only; C1, C2 constant; Cm12 raw lookup
Cm12 = squeeze(Cm_all(1,2,:));   % still raw lookup -- no model for this yet
Cm1R_model = Cm1p_fun(Cm1p_pfit, [rx, ry, rz, thx_rad, thy_rad, thz_rad]);
Cm2R_model = Cm2p_fun(Cm2p_pfit, [rx, ry, rz, thx_rad, thy_rad, thz_rad]);
SR_model = SR_fun(SR_pfit, [rx, ry, rz]);

C_M = Cm12 + Cm1R_model.*Cm2R_model./SR_model;
Ceff_new = C_M + C1_const*C2_const/(C1_const+C2_const);

%% Plot per block
blocks = {1:20, 21:40, 41:60};
blockNames = {'rz-only', 'thx-only', 'thy-only'};
blockX = {rz*1e6, thx, thy};
blockXNames = {'rz [\mum]', '\theta_x [deg]', '\theta_y [deg]'};

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

figure('Color','w','Name','Ceff: fitted-model-only (const C1,C2) vs dedicated sweep', 'Position',[50 50 1600 500]);
for col = 1:3
    idx = blocks{col};
    r2 = r2fun(Ceff_gt(idx), Ceff_new(idx));
    fprintf('%8s block: R^2 = %.4f\n', blockNames{col}, r2);

    subplot(1,3,col); hold on; grid on;
    scatter(blockX{col}(idx), Ceff_gt(idx), 30, 'k', 'filled', 'DisplayName', 'Ground truth (COMSOL)');
    scatter(blockX{col}(idx), Ceff_new(idx), 30, 'r', 'filled', 'DisplayName', sprintf('New model (R^2=%.4f)', r2));
    xlabel(blockXNames{col}); ylabel('C_{eff} (terminals 1-2) [F]');
    title(['C_{eff} vs ' blockXNames{col} ' (' blockNames{col} ' sweep)']);
    legend('Location','best');
end
sgtitle('C_{eff}: fitted-model-only for Cm1R/Cm2R/S_R, C1/C2 held constant at equilibrium, Cm12 raw lookup');

exportgraphics(gcf, fullfile(scriptDir, 'Ceff_vs_dedicatedSweep_constC1C2.png'));
fprintf('Saved plot to Ceff_vs_dedicatedSweep_constC1C2.png\n');
