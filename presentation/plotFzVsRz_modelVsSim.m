%% plotFzVsRz_modelVsSim.m
% Fz vs rz: MST simulation results vs the fitted fringing-capacitance
% nonlinear model, plus the error between them, using the full 2000-pt
% 6D LHS dataset (points scattered across all 6 DOF, not a clean 1D
% sweep, so some scatter from the other DOF is expected).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));  % rx,ry,rz,thx,thy,thz,Fz_MST,...

U = 30*ones(size(rx));
[~, ~, Fz_model] = trapElectrode1Model(rx, ry, rz, thx, thy, thz, U);

err = Fz_model - Fz_MST;
relErr = err ./ Fz_MST;

r2 = 1 - sum((Fz_MST-Fz_model).^2)/sum((Fz_MST-mean(Fz_MST)).^2);
fprintf('Fz: R^2 = %.4f\n', r2);
fprintf('Error: mean=%.4e  std=%.4e  max|err|=%.4e\n', mean(err), std(err), max(abs(err)));

figure('Color','w','Name','Fz vs rz: simulation vs model, with error', 'Position',[50 50 1400 450]);

subplot(1,2,1); hold on; grid on;
scatter(rz*1e6, Fz_MST, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Simulation');
scatter(rz*1e6, Fz_model, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Model');
xlabel('rz [\mum]'); ylabel('Fz [N]');
title(sprintf('Fz vs rz (R^2 = %.4f)', r2));
legend('Location','best');

subplot(1,2,2); hold on; grid on;
scatter(rz*1e6, err, 10, 'filled', 'MarkerFaceAlpha', 0.3);
yline(0, 'k--');
xlabel('rz [\mum]'); ylabel('Fz_{model} - Fz_{MST} [N]');
title('Error vs rz');

exportgraphics(gcf, fullfile(scriptDir, 'Fz_vs_rz_modelVsSim.png'));
fprintf('Saved plot to Fz_vs_rz_modelVsSim.png\n');
