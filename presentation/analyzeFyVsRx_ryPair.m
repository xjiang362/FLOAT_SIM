%% analyzeFyVsRx_ryPair.m
% Computes Fy(rx) via the energy method using a two-point ry bracket
% (ry=-1um, +1um) at each rx grid point: Fy(rx) = [W(rx,+1um) -
% W(rx,-1um)] / 2um -- a genuine, cheap, clean central-difference
% estimate, independent of any fitted model. Compares against the
% direct MST Fy probe at the same points.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_ryPair_results.txt'), 'NumHeaderLines', 5);
rx_all = raw(:,1); ry_all = raw(:,2);
Fy_MST_all = raw(:,8);
W_all = raw(:,13);

% Rows alternate (rx_i,-1um),(rx_i,+1um)
rxBase = rx_all(1:2:end);
W_minus = W_all(1:2:end);   % ry = -1um
W_plus  = W_all(2:2:end);   % ry = +1um
Fy_MST_minus = Fy_MST_all(1:2:end);
Fy_MST_plus  = Fy_MST_all(2:2:end);

ry_bracket = mean(ry_all(2:2:end));   % should be 1e-6
Fy_energy = (W_plus - W_minus) / (2*ry_bracket);
Fy_MST_avg = (Fy_MST_minus + Fy_MST_plus) / 2;   % MST reference at same rx (averaged over the two ry brackets)

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
fprintf('Fy_energy range: [%.4e, %.4e] N\n', min(Fy_energy), max(Fy_energy));
fprintf('Fy_MST_avg range: [%.4e, %.4e] N\n', min(Fy_MST_avg), max(Fy_MST_avg));
fprintf('R^2 (Fy_energy vs Fy_MST_avg) = %.4f\n', r2fun(Fy_MST_avg, Fy_energy));
fprintf('corr(rx, Fy_energy) = %.4f\n', corr(rxBase, Fy_energy));

muRx = mean(rxBase); sigmaRx = std(rxBase);
xhat = (rxBase-muRx)/sigmaRx;
p1 = polyfit(xhat, Fy_energy, 1);
p3 = polyfit(xhat, Fy_energy, 3);
fprintf('Fy_energy vs rx: linear R^2=%.4f, cubic R^2=%.4f\n', ...
    r2fun(Fy_energy, polyval(p1,xhat)), r2fun(Fy_energy, polyval(p3,xhat)));

%% Plot
figure('Color','w','Name','Fy vs rx: energy method (ry bracket) vs MST', 'Position',[50 50 900 600]);
hold on; grid on;
scatter(rxBase*1e6, Fy_MST_avg, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST, avg of \pm1\mum bracket)');
scatter(rxBase*1e6, Fy_energy, 30, [0 0.6 0], 'filled', 'DisplayName', 'Simulation (energy method, ry bracket)');
yline(0, 'k--', 'HandleVisibility', 'off');
xlabel('rx [\mum]'); ylabel('Fy [N]');
title('Fy vs rx: energy method (2-pt ry bracket) vs MST');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'FyVsRx_energyMethod_ryPair.png'));
fprintf('\nSaved plot to FyVsRx_energyMethod_ryPair.png\n');
