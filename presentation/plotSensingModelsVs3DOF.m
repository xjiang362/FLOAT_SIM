%% plotSensingModelsVs3DOF.m
% Preliminary check while the dedicated rz/thx/thy sensing sweep runs:
% S_R, Cm1p, Cm2p vs rz, thx, thy (3x3 grid), simulation vs fitted
% physics model, using the existing scattered 300-pt full-6D LHS
% capacitance dataset (capSweep300_withGND_allMatrices.mat) -- not a
% clean single-axis sweep, so scatter from the other DOF is expected in
% every panel, same caveat as the Fz-vs-3DOF LHS comparison.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'capSweep300_withGND_allMatrices.mat'));  % Cm_all, rx,ry,rz,thx,thy,thz, nPose
load(fullfile(scriptDir, 'SR_physicsModel.mat'));      % pfit, modelFun, d0,x0,y0
SR_pfit = pfit; SR_fun = modelFun;
load(fullfile(scriptDir, 'Cm1p_fringingFit.mat'));     % pfit, modelFun
Cm1p_pfit = pfit; Cm1p_fun = modelFun;
load(fullfile(scriptDir, 'Cm2p_fringingFit.mat'));     % pfit, modelFun
Cm2p_pfit = pfit; Cm2p_fun = modelFun;

thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;

S_R_sim  = squeeze(sum(Cm_all(21,:,:), 2));
Cm1p_sim = squeeze(Cm_all(1,21,:));
Cm2p_sim = squeeze(Cm_all(2,21,:));

S_R_model  = SR_fun(SR_pfit, [rx, ry, rz]);
Cm1p_model = Cm1p_fun(Cm1p_pfit, [rx, ry, rz, thx_rad, thy_rad, thz_rad]);
Cm2p_model = Cm2p_fun(Cm2p_pfit, [rx, ry, rz, thx_rad, thy_rad, thz_rad]);

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

quantities = {S_R_sim, S_R_model, 'S_R'; ...
              Cm1p_sim, Cm1p_model, 'Cm(1,p)'; ...
              Cm2p_sim, Cm2p_model, 'Cm(2,p)'};

dofs = {rz*1e6, thx, thy};
dofNames = {'rz [\mum]', '\theta_x [deg]', '\theta_y [deg]'};

figure('Color','w','Name','S_R, Cm1p, Cm2p vs rz, thx, thy (LHS dataset, preliminary)', 'Position',[50 50 1600 1200]);
for row = 1:3
    ySim = quantities{row,1};
    yModel = quantities{row,2};
    qName = quantities{row,3};
    r2 = r2fun(ySim, yModel);
    for col = 1:3
        subplot(3,3,(row-1)*3+col); hold on; grid on;
        scatter(dofs{col}, ySim, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Simulation');
        scatter(dofs{col}, yModel, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Fitted model');
        xlabel(dofNames{col}); ylabel([qName ' [F]']);
        if col == 1
            title(sprintf('%s vs %s (overall R^2=%.4f)', qName, dofNames{col}, r2));
        else
            title([qName ' vs ' dofNames{col}]);
        end
        if row == 1 && col == 1
            legend('Location','best');
        end
    end
    fprintf('%s: overall R^2 (full LHS dataset) = %.4f\n', qName, r2);
end

exportgraphics(gcf, fullfile(scriptDir, 'SensingModels_vs_3DOF_LHS.png'));
fprintf('\nSaved plot to SensingModels_vs_3DOF_LHS.png\n');
