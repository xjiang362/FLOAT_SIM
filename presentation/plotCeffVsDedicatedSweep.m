%% plotCeffVsDedicatedSweep.m
% Ceff (terminals 1-2): fitted model vs COMSOL ground truth, vs rz, thx,
% thy, using the dedicated single-axis sensing sweep (rz-only, thx-only,
% thy-only blocks, 20 pts each, all other DOF held exactly at zero).

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

%% Plot per block
blocks = {1:20, 21:40, 41:60};
blockNames = {'rz-only', 'thx-only', 'thy-only'};
blockX = {rz*1e6, thx, thy};
blockXNames = {'rz [\mum]', '\theta_x [deg]', '\theta_y [deg]'};

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

figure('Color','w','Name','Ceff vs rz, thx, thy: fitted model vs dedicated sweep', 'Position',[50 50 1600 500]);
for col = 1:3
    idx = blocks{col};
    r2 = r2fun(Ceff_gt(idx), Ceff_model(idx));
    fprintf('%8s block: R^2 = %.6f\n', blockNames{col}, r2);

    subplot(1,3,col); hold on; grid on;
    scatter(blockX{col}(idx), Ceff_gt(idx), 30, 'k', 'filled', 'DisplayName', 'Ground truth (COMSOL)');
    scatter(blockX{col}(idx), Ceff_model(idx), 30, 'r', 'filled', 'DisplayName', sprintf('Fitted model (R^2=%.4f)', r2));
    xlabel(blockXNames{col}); ylabel('C_{eff} (terminals 1-2) [F]');
    title(['C_{eff} vs ' blockXNames{col} ' (' blockNames{col} ' sweep)']);
    legend('Location','best');
end

exportgraphics(gcf, fullfile(scriptDir, 'Ceff_vs_dedicatedSweep.png'));
fprintf('Saved plot to Ceff_vs_dedicatedSweep.png\n');
