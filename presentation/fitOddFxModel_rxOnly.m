%% fitOddFxModel_rxOnly.m
% Fits Fx(rx) = F0 + c1*rx + c3*rx^3, with the constant F0 FIXED at
% 1.4e-7 N (the offset the curve is centered around, not the origin) and
% only c1, c3 (the ODD/antisymmetric part) fit by linear least squares.
% Target: the 3-pt movmean-smoothed energy-method Fx (validated as the
% more trustworthy simulation-derived Fx estimate), on the dedicated
% uniform rx-only sweep.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_uniform_energy_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1);
Fx_MST = raw(:,7);
W = raw(:,13);
W_movmean = movmean(W, 3);
Fx_movmean = gradient(W_movmean, rx);

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

F0_fixed = -1.4e-7;
target = Fx_movmean - F0_fixed;
designMat = [rx, rx.^3];
c = designMat \ target;
c1 = c(1); c3 = c(2);
Fx_odd_fixedF0 = F0_fixed + c1*rx + c3*rx.^3;

% Also fit F0 freely (full 3-param linear least squares: F0, c1, c3 all fit)
designMatFree = [ones(size(rx)), rx, rx.^3];
cFree = designMatFree \ Fx_movmean;
F0 = cFree(1); c1free = cFree(2); c3free = cFree(3);
Fx_odd = designMatFree * cFree;

fprintf('Fixed-F0 fit:  F0=%.4e (fixed), c1=%.6e, c3=%.6e, R^2 vs movmean=%.4f\n', ...
    F0_fixed, c1, c3, r2fun(Fx_movmean, Fx_odd_fixedF0));
fprintf('Free-F0 fit:   F0=%.4e (fit),   c1=%.6e, c3=%.6e, R^2 vs movmean=%.4f\n', ...
    F0, c1free, c3free, r2fun(Fx_movmean, Fx_odd));
fprintf('\nR^2 (free-F0 model vs MST) = %.4f\n', r2fun(Fx_MST, Fx_odd));

%% Plot
figure('Color','w','Name','Odd Fx model centered at F0, vs rx', 'Position',[50 50 800 600]);
hold on; grid on;
scatter(rx*1e6, Fx_MST, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(rx*1e6, Fx_movmean, 25, 'b', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', '3-pt movmean energy method');
rxFine = linspace(min(rx), max(rx), 200)';
plot(rxFine*1e6, F0_fixed + c1*rxFine + c3*rxFine.^3, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'LineStyle', '--', ...
    'DisplayName', sprintf('Fixed F0=%.2e (R^2=%.4f)', F0_fixed, r2fun(Fx_movmean, Fx_odd_fixedF0)));
plot(rxFine*1e6, F0 + c1free*rxFine + c3free*rxFine.^3, 'Color', [0 0.6 0], 'LineWidth', 2, ...
    'DisplayName', sprintf('Free F0=%.2e (R^2=%.4f)', F0, r2fun(Fx_movmean, Fx_odd)));
yline(0, 'k--', 'HandleVisibility', 'off');
xlabel('rx [\mum]'); ylabel('Fx [N]');
title('Odd Fx model: fixed vs free F0, vs rx');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'OddFxModel_rxOnly.png'));
fprintf('\nSaved plot to OddFxModel_rxOnly.png\n');
