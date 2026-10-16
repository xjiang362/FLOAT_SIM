%% plotFxVsRx_MSTvsEnergy.m
% Fx vs rx, MST and fringing-model energy-derived overlaid, across the
% full 2000-pt scattered LHS dataset (not binned -- just a direct
% ballpark/trend sanity check, not a precision validation, since Fx
% depends weakly on other DOF too in both methods).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));
load(fullfile(scriptDir, 'fringingModelForces.mat'));

figure('Color','w','Name','Fx vs rx: MST vs fringing-model energy-derived');
subplot(1,2,1);
scatter(rx*1e6, Fx_MST, 8, 'filled', 'MarkerFaceAlpha', 0.3);
grid on; xlabel('rx [\mum]'); ylabel('Fx [N]'); title('Fx (MST)');
ylim([-1e-5, 1e-5]);

subplot(1,2,2);
scatter(rx*1e6, Fx_energy, 8, 'filled', 'MarkerFaceAlpha', 0.3);
grid on; xlabel('rx [\mum]'); ylabel('Fx [N]'); title('Fx (fringing energy model)');

sgtitle('Fx vs rx: ballpark/trend comparison');
exportgraphics(gcf, fullfile(scriptDir, 'Fx_vs_rx_MSTvsEnergy.png'));
fprintf('Saved plot to Fx_vs_rx_MSTvsEnergy.png\n');

fprintf('\nFx MST:    range [%.3g, %.3g] N, median=%.3g\n', min(Fx_MST), max(Fx_MST), median(Fx_MST));
fprintf('Fx energy: range [%.3g, %.3g] N, median=%.3g\n', min(Fx_energy), max(Fx_energy), median(Fx_energy));
