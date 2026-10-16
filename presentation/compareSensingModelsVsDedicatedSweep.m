%% compareSensingModelsVsDedicatedSweep.m
% Real test (not the scattered LHS preliminary check): compares the
% fitted S_R, Cm1p, Cm2p models, and the fully-assembled Ceff model,
% against the dedicated single-axis sensing sweep (rz-only, thx-only,
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

%% Simulation (ground truth) values
S_R_sim  = squeeze(sum(Cm_all(21,:,:), 2));
Cm1p_sim = squeeze(Cm_all(1,21,:));
Cm2p_sim = squeeze(Cm_all(2,21,:));

floatingTerms = [5,6,9,10,13,14,21];
Ceff_gt = zeros(nPose,1);
for pp = 1:nPose
    Ceff_gt(pp) = effectiveMutualCapacitance(Cm_all(:,:,pp), 1, 2, floatingTerms);
end

%% Fitted-model values
S_R_model  = SR_fun(SR_pfit, [rx, ry, rz]);
Cm1p_model = Cm1p_fun(Cm1p_pfit, [rx, ry, rz, thx_rad, thy_rad, thz_rad]);
Cm2p_model = Cm2p_fun(Cm2p_pfit, [rx, ry, rz, thx_rad, thy_rad, thz_rad]);

Cm12 = squeeze(Cm_all(1,2,:));
S1 = squeeze(sum(Cm_all(1,:,:),2));
S2 = squeeze(sum(Cm_all(2,:,:),2));
C_M = Cm12 + Cm1p_model.*Cm2p_model./S_R_model;
C1 = S1 - Cm12 - Cm1p_model.*(Cm1p_model+Cm2p_model)./S_R_model;
C2 = S2 - Cm12 - Cm2p_model.*(Cm1p_model+Cm2p_model)./S_R_model;
Ceff_model = C_M + C1.*C2./(C1+C2);

%% R^2 per quantity per block
blocks = {1:20, 21:40, 41:60};
blockNames = {'rz-only', 'thx-only', 'thy-only'};
blockX = {rz*1e6, thx, thy};
blockXNames = {'rz [\mum]', '\theta_x [deg]', '\theta_y [deg]'};

quantities = {S_R_sim, S_R_model, 'S_R'; ...
              Cm1p_sim, Cm1p_model, 'Cm(1,p)'; ...
              Cm2p_sim, Cm2p_model, 'Cm(2,p)'; ...
              Ceff_gt, Ceff_model, 'C_{eff}'};

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

figure('Color','w','Name','Sensing models vs dedicated rz/thx/thy sweep', 'Position',[50 50 1600 1500]);
for row = 1:4
    ySim = quantities{row,1};
    yModel = quantities{row,2};
    qName = quantities{row,3};
    for col = 1:3
        idx = blocks{col};
        r2 = r2fun(ySim(idx), yModel(idx));
        fprintf('%8s, %8s block: R^2 = %.4f\n', qName, blockNames{col}, r2);

        subplot(4,3,(row-1)*3+col); hold on; grid on;
        scatter(blockX{col}(idx), ySim(idx), 30, 'k', 'filled', 'DisplayName', 'Simulation');
        scatter(blockX{col}(idx), yModel(idx), 30, 'r', 'filled', 'DisplayName', sprintf('Fitted model (R^2=%.4f)', r2));
        xlabel(blockXNames{col}); ylabel([qName ' [F]']);
        title([qName ' vs ' blockXNames{col} ' (' blockNames{col} ' sweep)']);
        legend('Location','best');
    end
end

exportgraphics(gcf, fullfile(scriptDir, 'SensingModels_vs_dedicatedSweep.png'));
fprintf('\nSaved plot to SensingModels_vs_dedicatedSweep.png\n');
