%% plotMSTvsEnergy_rxOnly.m
% Plots the two independent simulation-derived Fx estimates (MST probe
% vs energy-method gradient of the real simulated W) on the dedicated
% uniform rx-only sweep, plus the raw W(rx) curve itself for context
% (to show just how flat/noisy-looking it is over this range).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_uniform_energy_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1);
Fx_MST = raw(:,7);
W = raw(:,13);
Fx_energy = gradient(W, rx);

figure('Color','w','Name','MST vs energy method: rx-only sweep', 'Position',[50 50 1400 500]);

subplot(1,2,1); hold on; grid on;
plot(rx*1e6, W, 'o-', 'Color', [0 0.6 0], 'LineWidth', 1.2, 'MarkerFaceColor', [0 0.6 0]);
xlabel('rx [\mum]'); ylabel('W [J]');
title('Raw simulated total electric energy vs rx');

subplot(1,2,2); hold on; grid on;
scatter(rx*1e6, Fx_MST, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(rx*1e6, Fx_energy, 30, [0 0.6 0], 'filled', 'DisplayName', 'Simulation (energy method)');
yline(0, 'k--', 'HandleVisibility', 'off');
xlabel('rx [\mum]'); ylabel('Fx [N]');
title('Fx vs rx: MST vs energy method');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'MSTvsEnergy_rxOnly.png'));
fprintf('Saved plot to MSTvsEnergy_rxOnly.png\n');
