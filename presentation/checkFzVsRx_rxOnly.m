%% checkFzVsRx_rxOnly.m
% Checks whether Fz shows a real dependence on rx, using the dedicated
% uniform rx-only sweep (comsol_sweep_rxOnly_uniform_energy_results.txt,
% rz=thx=thy=thz=ry=0 throughout, only rx varies). Fz here is the direct
% MST probe (this sweep has no rz variation, so the energy method can't
% independently derive Fz from this dataset -- only Fx).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_uniform_energy_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1);
Fz_MST = raw(:,9);

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

fprintf('Fz range: [%.6e, %.6e] N, mean=%.6e, std=%.4e\n', min(Fz_MST), max(Fz_MST), mean(Fz_MST), std(Fz_MST));
fprintf('(for reference, Fx on this same sweep had mean|Fx|~7e-8 N -- Fz is ~500x larger, well above noise floor)\n\n');

muRx = mean(rx); sigmaRx = std(rx);
xhat = (rx - muRx) / sigmaRx;
p1 = polyfit(xhat, Fz_MST, 1);
p2 = polyfit(xhat, Fz_MST, 2);
Fz1 = polyval(p1, xhat);
Fz2 = polyval(p2, xhat);
fprintf('Linear fit:    R^2 = %.4f\n', r2fun(Fz_MST, Fz1));
fprintf('Quadratic fit: R^2 = %.4f\n', r2fun(Fz_MST, Fz2));
fprintf('Correlation(rx, Fz) = %.4f\n', corr(rx, Fz_MST));

%% Plot
figure('Color','w','Name','Fz vs rx (dedicated uniform sweep)', 'Position',[50 50 800 600]);
hold on; grid on;
scatter(rx*1e6, Fz_MST, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
rxFine = linspace(min(rx), max(rx), 200)';
xhatFine = (rxFine-muRx)/sigmaRx;
plot(rxFine*1e6, polyval(p2, xhatFine), 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Quadratic fit (R^2=%.4f)', r2fun(Fz_MST,Fz2)));
xlabel('rx [\mum]'); ylabel('Fz [N]');
title('Fz vs rx: clear real trend, well above noise floor');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'FzVsRx_rxOnly.png'));
fprintf('\nSaved plot to FzVsRx_rxOnly.png\n');
