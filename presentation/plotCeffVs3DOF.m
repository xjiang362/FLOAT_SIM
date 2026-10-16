%% plotCeffVs3DOF.m
% Ceff (between terminals 1 and 2) vs rz, thx, thy in 3 subplots:
% ground truth (effectiveMutualCapacitance.m) vs the fully fitted
% physics model (S_R, Cm1p, Cm2p fitted models plugged into the
% resonator-only analytical formula; Cm12, S1, S2 still from raw
% simulation since not yet modeled). Uses the existing scattered 300-pt
% full-6D LHS capacitance dataset (preliminary, while the dedicated
% rz/thx/thy-only sweep runs).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'capSweep300_withGND_allMatrices.mat'));  % Cm_all, rx,ry,rz,thx,thy,thz, nPose
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

%% Fitted-model Ceff
Cm12 = squeeze(Cm_all(1,2,:));
S1 = squeeze(sum(Cm_all(1,:,:),2));
S2 = squeeze(sum(Cm_all(2,:,:),2));

Cm1R = Cm1p_fun(Cm1p_pfit, [rx, ry, rz, thx_rad, thy_rad, thz_rad]);
Cm2R = Cm2p_fun(Cm2p_pfit, [rx, ry, rz, thx_rad, thy_rad, thz_rad]);
SR = SR_fun(SR_pfit, [rx, ry, rz]);

C_M = Cm12 + Cm1R.*Cm2R./SR;
C1 = S1 - Cm12 - Cm1R.*(Cm1R+Cm2R)./SR;
C2 = S2 - Cm12 - Cm2R.*(Cm1R+Cm2R)./SR;
Ceff_model = C_M + C1.*C2./(C1+C2);

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
r2 = r2fun(Ceff_gt, Ceff_model);
fprintf('Ceff (fitted model) overall R^2 vs ground truth: %.6f\n', r2);

%% Plot
dofs = {rz*1e6, thx, thy};
dofNames = {'rz [\mum]', '\theta_x [deg]', '\theta_y [deg]'};

figure('Color','w','Name','Ceff vs rz, thx, thy: fitted model vs ground truth', 'Position',[50 50 1600 500]);
for i = 1:3
    subplot(1,3,i); hold on; grid on;
    scatter(dofs{i}, Ceff_gt, 15, 'k', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'Ground truth (COMSOL)');
    scatter(dofs{i}, Ceff_model, 15, 'r', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'Fitted model');
    xlabel(dofNames{i}); ylabel('C_{eff} (terminals 1-2) [F]');
    title(['C_{eff} vs ' dofNames{i}]);
    if i == 1
        legend('Location','best');
    end
end
sgtitle(sprintf('C_{eff}: fitted model vs COMSOL ground truth (overall R^2=%.4f)', r2));

exportgraphics(gcf, fullfile(scriptDir, 'Ceff_vs_3DOF.png'));
fprintf('Saved plot to Ceff_vs_3DOF.png\n');
