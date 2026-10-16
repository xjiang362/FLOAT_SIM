%% plotEnergyVsRx_ryThzPairs.m
% Plots the raw simulated energy W vs rx from the voltage-corrected
% rx/ry/thz bracket sweep -- the baseline (ry-averaged, ~ry=0,thz=0)
% energy curve that Fx is differentiated from.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_ryThzPairs_results.txt'), 'NumHeaderLines', 5);
rx_all = raw(:,1);
W_all = raw(:,13);

rxBase = rx_all(1:4:end);
W_ryMinus = W_all(1:4:end);   % ry=-1um, thz=0
W_ryPlus  = W_all(2:4:end);   % ry=+1um, thz=0
W_avg = (W_ryMinus + W_ryPlus) / 2;

fprintf('W_avg range: [%.10e, %.10e] J\n', min(W_avg), max(W_avg));
fprintf('swing/mean = %.4f%%\n', 100*(max(W_avg)-min(W_avg))/mean(W_avg));

figure('Color','w','Name','Energy vs rx');
hold on; grid on;
plot(rxBase*1e6, W_ryMinus, 'Color', [0.6 0.6 1], 'LineWidth', 1, 'DisplayName', 'W (ry=-1\mum)');
plot(rxBase*1e6, W_ryPlus, 'Color', [1 0.6 0.6], 'LineWidth', 1, 'DisplayName', 'W (ry=+1\mum)');
plot(rxBase*1e6, W_avg, 'k-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'W (averaged, ~ry=0)');
xlabel('rx [\mum]'); ylabel('W [J]');
title('Energy vs rx (voltage-corrected sweep)');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'EnergyVsRx_ryThzPairs.png'));
fprintf('Saved plot to EnergyVsRx_ryThzPairs.png\n');
