%% analyzeRxSweep_ryThzPairs.m
% Analyzes the voltage-corrected rx sweep with ry and thz brackets
% (comsol_sweep_rxOnly_ryThzPairs_results.txt, 200 rows = 50 rx values x
% 4 rows each: (ry=-1um,thz=0), (ry=+1um,thz=0), (ry=0,thz=-thzB(rx)),
% (ry=0,thz=+thzB(rx))). Derives THREE energy-method relationships:
%   Fx(rx)   = d/drx of W_avg(rx), where W_avg = mean(W(ry=-1um), W(ry=+1um))
%              approximates W(rx,ry=0,thz=0) to O(ry^2)
%   Fy(rx)   = [W(ry=+1um) - W(ry=-1um)] / (2*1um)   (central diff in ry)
%   Tauz(rx) = [W(thz=+thzB) - W(thz=-thzB)] / (2*thzB_rad)  (central diff in thz)
% Tauz is undefined at the two rx endpoints (+-45um) where the safe thz
% bracket collapses to 0 -- excluded from that fit.
% Compares each against the corresponding direct MST probe (averaged
% over the relevant bracket pair) for validation.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_ryThzPairs_results.txt'), 'NumHeaderLines', 5);
rx_all = raw(:,1); ry_all = raw(:,2); thz_all = raw(:,6);
Fx_all = raw(:,7); Fy_all = raw(:,8); Tauz_all = raw(:,12);
W_all = raw(:,13);

nPtsRx = numel(rx_all)/4;
rxBase = rx_all(1:4:end);

W_ryMinus = W_all(1:4:end);   % ry=-1um, thz=0
W_ryPlus  = W_all(2:4:end);   % ry=+1um, thz=0
W_thzMinus = W_all(3:4:end);  % ry=0, thz=-thzB(rx)
W_thzPlus  = W_all(4:4:end);  % ry=0, thz=+thzB(rx)

Fx_MST_avg = (Fx_all(1:4:end) + Fx_all(2:4:end)) / 2;
Fy_MST_avg = (Fy_all(1:4:end) + Fy_all(2:4:end)) / 2;
Tauz_MST_avg = (Tauz_all(3:4:end) + Tauz_all(4:4:end)) / 2;

ry_bracket = 1e-6;
thz_bracket_deg = (thz_all(4:4:end) - thz_all(3:4:end)) / 2;   % half the total span, in deg
thz_bracket_rad = thz_bracket_deg * pi/180;

%% Fx(rx): differentiate the (ry-averaged) energy baseline
W_avg = (W_ryMinus + W_ryPlus) / 2;
Fx_energy = gradient(W_avg, rxBase);

%% Fy(rx): central difference across the ry bracket
Fy_energy = (W_ryPlus - W_ryMinus) / (2*ry_bracket);

%% Tauz(rx): central difference across the thz bracket (skip zero-bracket endpoints)
validThz = thz_bracket_rad > 0;
Tauz_energy = nan(nPtsRx,1);
Tauz_energy(validThz) = (W_thzPlus(validThz) - W_thzMinus(validThz)) ./ (2*thz_bracket_rad(validThz));

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
fprintf('W range: [%.4e, %.4e] J\n', min(W_all), max(W_all));
fprintf('R^2 (Fx_energy vs Fx_MST_avg) = %.4f\n', r2fun(Fx_MST_avg, Fx_energy));
fprintf('R^2 (Fy_energy vs Fy_MST_avg) = %.4f\n', r2fun(Fy_MST_avg, Fy_energy));
fprintf('R^2 (Tauz_energy vs Tauz_MST_avg, valid pts only) = %.4f\n', ...
    r2fun(Tauz_MST_avg(validThz), Tauz_energy(validThz)));

%% Plot 1: Fx vs rx
figure('Color','w','Name','Fx vs rx (energy method)');
hold on; grid on;
scatter(rxBase*1e6, Fx_MST_avg, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(rxBase*1e6, Fx_energy, 30, [0 0.6 0], 'filled', 'DisplayName', 'Simulation (energy method)');
xlabel('rx [\mum]'); ylabel('Fx [N]');
title(sprintf('Fx vs rx: energy method vs MST (R^2=%.4f)', r2fun(Fx_MST_avg, Fx_energy)));
legend('Location','best');
exportgraphics(gcf, fullfile(scriptDir, 'FxVsRx_energyMethod_voltageFixed.png'));
fprintf('Saved plot to FxVsRx_energyMethod_voltageFixed.png\n');

%% Plot 2: Fy vs rx
figure('Color','w','Name','Fy vs rx (energy method)');
hold on; grid on;
scatter(rxBase*1e6, Fy_MST_avg, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(rxBase*1e6, Fy_energy, 30, [0 0.6 0], 'filled', 'DisplayName', 'Simulation (energy method)');
xlabel('rx [\mum]'); ylabel('Fy [N]');
title(sprintf('Fy vs rx: energy method vs MST (R^2=%.4f)', r2fun(Fy_MST_avg, Fy_energy)));
legend('Location','best');
exportgraphics(gcf, fullfile(scriptDir, 'FyVsRx_energyMethod_voltageFixed.png'));
fprintf('Saved plot to FyVsRx_energyMethod_voltageFixed.png\n');

%% Plot 3: Tauz vs rx
figure('Color','w','Name','Tauz vs rx (energy method)');
hold on; grid on;
scatter(rxBase(validThz)*1e6, Tauz_MST_avg(validThz), 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(rxBase(validThz)*1e6, Tauz_energy(validThz), 30, [0 0.6 0], 'filled', 'DisplayName', 'Simulation (energy method)');
xlabel('rx [\mum]'); ylabel('Tauz [N*m]');
title(sprintf('Tauz vs rx: energy method vs MST (R^2=%.4f)', r2fun(Tauz_MST_avg(validThz), Tauz_energy(validThz))));
legend('Location','best');
exportgraphics(gcf, fullfile(scriptDir, 'TauzVsRx_energyMethod_voltageFixed.png'));
fprintf('Saved plot to TauzVsRx_energyMethod_voltageFixed.png\n');
