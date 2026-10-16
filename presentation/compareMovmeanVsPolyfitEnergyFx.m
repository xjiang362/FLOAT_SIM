%% compareMovmeanVsPolyfitEnergyFx.m
% Compares two smoothing-then-differentiate approaches for the
% energy-method Fx(rx) estimate on the dedicated uniform rx-only sweep:
%   1. Polyfit (order 4, chosen by LOO CV) then analytic derivative
%      (polyFitEnergyFx_rxOnly.m).
%   2. 3-point moving average low-pass filter (movmean) on W, then
%      numerical gradient of the smoothed curve.
% Both are compared against MST and against each other.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_uniform_energy_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1);
Fx_MST = raw(:,7);
W = raw(:,13);

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

%% Method 1: raw numerical gradient (baseline, no smoothing)
Fx_raw = gradient(W, rx);

%% Method 2: polyfit (order 4, from LOO CV) then analytic derivative
muRx = mean(rx); sigmaRx = std(rx);
xhat = (rx - muRx) / sigmaRx;
bestOrder = 4;
p = polyfit(xhat, W, bestOrder);
dp = polyder(p);
Fx_poly = polyval(dp, xhat) / sigmaRx;
W_poly = polyval(p, xhat);

%% Method 3: 3-point moving average low-pass filter, then gradient
W_movmean = movmean(W, 3);
Fx_movmean = gradient(W_movmean, rx);

%% Compare
fprintf('%-22s %12s\n', 'method', 'R^2 vs MST');
fprintf('%-22s %12.4f\n', 'raw gradient', r2fun(Fx_MST, Fx_raw));
fprintf('%-22s %12.4f\n', 'polyfit (order 4)', r2fun(Fx_MST, Fx_poly));
fprintf('%-22s %12.4f\n', '3-pt movmean', r2fun(Fx_MST, Fx_movmean));
fprintf('\npolyfit vs movmean, R^2 = %.4f\n', r2fun(Fx_poly, Fx_movmean));

%% Plot
figure('Color','w','Name','Movmean vs polyfit energy method', 'Position',[50 50 1400 500]);

subplot(1,2,1); hold on; grid on;
scatter(rx*1e6, W, 25, 'k', 'filled', 'DisplayName', 'Simulated W (raw)');
plot(rx*1e6, W_movmean, 'b-o', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', '3-pt movmean');
rxFine = linspace(min(rx), max(rx), 200)';
plot(rxFine*1e6, polyval(p, (rxFine-muRx)/sigmaRx), 'r-', 'LineWidth', 2, 'DisplayName', 'Polyfit (order 4)');
xlabel('rx [\mum]'); ylabel('W [J]');
title('W(rx): raw vs 3-pt movmean vs polyfit');
legend('Location','best');

subplot(1,2,2); hold on; grid on;
scatter(rx*1e6, Fx_MST, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(rx*1e6, Fx_movmean, 25, 'b', 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', sprintf('3-pt movmean (R^2=%.4f)', r2fun(Fx_MST, Fx_movmean)));
plot(rxFine*1e6, polyval(dp, (rxFine-muRx)/sigmaRx)/sigmaRx, 'Color', [0 0.6 0], 'LineWidth', 2, 'DisplayName', sprintf('Polyfit order 4 (R^2=%.4f)', r2fun(Fx_MST, Fx_poly)));
yline(0, 'k--', 'HandleVisibility', 'off');
xlabel('rx [\mum]'); ylabel('Fx [N]');
title('Fx vs rx: MST vs two smoothing methods');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'MovmeanVsPolyfit_energyFx.png'));
fprintf('\nSaved plot to MovmeanVsPolyfit_energyFx.png\n');
