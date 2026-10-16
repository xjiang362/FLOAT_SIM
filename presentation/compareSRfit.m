%% compareSRfit.m
% Compares the fitted S_R physics model against the raw simulation
% S_R data: predicted-vs-actual scatter, and both vs rz to see shape
% match directly.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'capSweep300_withGND_allMatrices.mat'));  % Cm_all, rx,ry,rz,...
load(fullfile(scriptDir, 'SR_physicsModel.mat'));  % pfit, modelFun, d0,x0,y0, r2, r2_cv

S_R_sim = squeeze(sum(Cm_all(21,:,:), 2));
S_R_model = modelFun(pfit, [rx, ry, rz]);

fprintf('In-sample R^2: %.4f, 10-fold CV R^2: %.4f\n', r2, r2_cv);

figure('Color','w','Name','S_R: fitted model vs simulation', 'Position',[50 50 1400 500]);

subplot(1,2,1); hold on; grid on;
scatter(S_R_sim, S_R_model, 20, 'filled', 'MarkerFaceAlpha', 0.4);
lims = [min(S_R_sim) max(S_R_sim)];
plot(lims, lims, 'r--', 'LineWidth', 1.5);
xlabel('S_R (simulation) [F]'); ylabel('S_R (fitted model) [F]');
title(sprintf('Predicted vs actual (R^2=%.4f)', r2));

subplot(1,2,2); hold on; grid on;
scatter(rz*1e6, S_R_sim, 15, 'b', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'Simulation');
scatter(rz*1e6, S_R_model, 15, 'r', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'Fitted model');
xlabel('rz [\mum]'); ylabel('S_R [F]');
title('S_R vs rz: model vs simulation');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'SR_fit_vs_sim.png'));
fprintf('Saved plot to SR_fit_vs_sim.png\n');

err = S_R_model - S_R_sim;
fprintf('\nError: mean=%.4e  std=%.4e  max|err|=%.4e\n', mean(err), std(err), max(abs(err)));
