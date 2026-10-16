%% plotEnergyVsThz.m
% Plots the raw energy W vs thz from the thz-primary sweep
% (comsol_sweep_thzOnly_rxPair_energy_results.txt): both rx-bracket
% tracks (rx=-1um, +1um) and their average (~rx=0), to visually inspect
% the shape discussed (is there really a dip at thz=0, or does it peak
% there due to noise?).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_thzOnly_rxPair_energy_results.txt'), 'NumHeaderLines', 5);
thz_all = raw(:,6);
W_all = raw(:,13);

thzBase = thz_all(1:2:end);
W_rxMinus = W_all(1:2:end);
W_rxPlus  = W_all(2:2:end);
W_avg = (W_rxMinus + W_rxPlus) / 2;

fprintf('W_avg range: [%.10e, %.10e] J\n', min(W_avg), max(W_avg));
fprintf('swing/mean = %.4f%%\n', 100*(max(W_avg)-min(W_avg))/mean(W_avg));

figure('Color','w','Name','Energy vs thz');
hold on; grid on;
plot(thzBase, W_rxMinus, 'Color', [0.6 0.6 1], 'LineWidth', 1, 'DisplayName', 'W (rx=-1\mum)');
plot(thzBase, W_rxPlus, 'Color', [1 0.6 0.6], 'LineWidth', 1, 'DisplayName', 'W (rx=+1\mum)');
plot(thzBase, W_avg, 'k-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'W (averaged, ~rx=0)');
xline(0, 'k:', 'HandleVisibility', 'off');
xlabel('\theta_z [deg]'); ylabel('W [J]');
title('Energy vs \theta_z (thz-primary sweep)');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'EnergyVsThz.png'));
fprintf('Saved plot to EnergyVsThz.png\n');
