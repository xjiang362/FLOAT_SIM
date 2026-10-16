%% analyzeRxOnlySweep_uniform.m
% Analyzes the dedicated, uniformly-spaced rx-only sweep
% (comsol_sweep_rxOnly_uniform_energy_results.txt, 50 pts, rest of DOF
% held at exactly zero), which has BOTH MST force/torque probes and a
% genuine Total electric energy probe. Since this is a clean 1D uniform
% sweep (unlike the scattered 6D LHS data), a simple numerical gradient
% of W vs rx gives a trustworthy, model-free energy-method Fx estimate
% to compare against MST directly -- resolving whether the coupled-WL
% model's flat/wrong-sign Fx prediction is a real physics gap or a
% simulation-side sign/noise issue.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_uniform_energy_results.txt'), 'NumHeaderLines', 5);
if size(raw,2) ~= 13
    error('Expected 13 columns, got %d -- check the export format.', size(raw,2));
end

rx = raw(:,1); ry = raw(:,2); rz = raw(:,3);
thx = raw(:,4); thy = raw(:,5); thz = raw(:,6);
Fx_MST = raw(:,7); Fy_MST = raw(:,8); Fz_MST = raw(:,9);
Taux_MST = raw(:,10); Tauy_MST = raw(:,11); Tauz_MST = raw(:,12);
W = raw(:,13);
n = numel(rx);

fprintf('n = %d, rx range [%.4g, %.4g] m\n', n, min(rx), max(rx));
fprintf('W range: [%.6e, %.6e] J\n', min(W), max(W));

%% Energy method: clean central difference (uniform spacing, single axis)
Fx_energy = gradient(W, rx);

%% Coupled-WL model prediction
U = 30*ones(n,1);
[Fx_model, ~, ~, ~, ~, ~] = trapElectrode1Model_coupledWL(rx, ry, rz, thx, thy, thz, U);

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
fprintf('\nFx: energy-method (gradient of real W) vs MST, R^2 = %.4f\n', r2fun(Fx_MST, Fx_energy));
fprintf('Fx: coupled-WL model vs MST, R^2 = %.4f\n', r2fun(Fx_MST, Fx_model));
fprintf('Fx: coupled-WL model vs energy-method, R^2 = %.4f\n', r2fun(Fx_energy, Fx_model));

fprintf('\nMagnitudes: mean|Fx_MST|=%.4e, mean|Fx_energy|=%.4e, mean(Fx_model)=%.4e\n', ...
    mean(abs(Fx_MST)), mean(abs(Fx_energy)), mean(Fx_model));

%% Plot
figure('Color','w','Name','Fx vs rx: MST, energy method, coupled-WL model', 'Position',[50 50 900 600]);
hold on; grid on;
scatter(rx*1e6, Fx_MST, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(rx*1e6, Fx_energy, 30, [0 0.6 0], 'filled', 'DisplayName', 'Simulation (energy method, gradient of W)');
plot(rx*1e6, Fx_model, 'r-', 'LineWidth', 2, 'DisplayName', 'Coupled-WL model');
xlabel('rx [\mum]'); ylabel('Fx [N]');
title('Fx vs rx (uniform dedicated sweep)');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'Fx_vs_rx_uniform_MSTvsEnergyVsModel.png'));
fprintf('\nSaved plot to Fx_vs_rx_uniform_MSTvsEnergyVsModel.png\n');
