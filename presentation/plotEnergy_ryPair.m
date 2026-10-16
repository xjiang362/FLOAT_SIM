%% plotEnergy_ryPair.m
% Plots the raw W(rx) at both ry brackets (-1um, +1um) to see directly
% whether there's a resolvable difference between them, or whether it's
% within noise -- the thing Fy_energy = [W(+1um)-W(-1um)]/2um is trying
% to extract.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_ryPair_results.txt'), 'NumHeaderLines', 5);
rx_all = raw(:,1);
W_all = raw(:,13);

rxBase = rx_all(1:2:end);
W_minus = W_all(1:2:end);   % ry = -1um
W_plus  = W_all(2:2:end);   % ry = +1um

fprintf('W_minus (ry=-1um) range: [%.10e, %.10e] J\n', min(W_minus), max(W_minus));
fprintf('W_plus  (ry=+1um) range: [%.10e, %.10e] J\n', min(W_plus), max(W_plus));
fprintf('W_plus - W_minus range: [%.4e, %.4e] J\n', min(W_plus-W_minus), max(W_plus-W_minus));
fprintf('For reference, W itself is ~%.4e J -- the difference is ~%.2e relative\n', mean(W_minus), mean(abs(W_plus-W_minus))/mean(W_minus));

figure('Color','w','Name','W(rx) at ry=-1um and +1um', 'Position',[50 50 1400 500]);

subplot(1,2,1); hold on; grid on;
plot(rxBase*1e6, W_minus, 'b-o', 'MarkerSize', 4, 'DisplayName', 'ry = -1\mum');
plot(rxBase*1e6, W_plus, 'r-o', 'MarkerSize', 4, 'DisplayName', 'ry = +1\mum');
xlabel('rx [\mum]'); ylabel('W [J]');
title('W(rx) at both ry brackets');
legend('Location','best');

subplot(1,2,2); hold on; grid on;
plot(rxBase*1e6, W_plus - W_minus, 'g-o', 'MarkerSize', 4);
yline(0, 'k--');
xlabel('rx [\mum]'); ylabel('W(+1\mum) - W(-1\mum) [J]');
title('Energy difference driving Fy_{energy}');

exportgraphics(gcf, fullfile(scriptDir, 'Energy_ryPair.png'));
fprintf('\nSaved plot to Energy_ryPair.png\n');
