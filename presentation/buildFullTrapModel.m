%% buildFullTrapModel.m
% Assembles the full 6D model from the single-axis characterization:
%   Fx   ~ poly fit vs rx        (validated, LOO R^2=0.84)
%   Fy   ~ poly fit vs ry        (validated, LOO R^2=0.55)
%   Fz   ~ k_Fz   / (g0_z-rz)^2  (physics-informed, LOO R^2=0.99)
%   Taux ~ k_Taux / (g0_z-rz)^2  (physics-informed, LOO R^2=0.99)
%   Tauy ~ k_Tauy / (g0_z-rz)^2  (physics-informed, LOO R^2=0.99)
%   Tauz ~ best available single-axis fit (weakest link, R^2=0.41)
%
% Then validates this "dominant axis only" additive model against the
% REAL 6D LHS dataset (55 points, all 6 DOF genuinely varying
% together) -- this directly tests whether ignoring the small
% secondary cross-couplings we found is actually good enough, or
% whether real interaction terms need to be added on top.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

g0_z = 15e-6;

%% Load single-axis characterization data/models
loaded = load(fullfile(scriptDir, 'combined300_fittedModels.mat'));
sweepData = loaded.data;
axisBlocks = loaded.axisBlocks;

model = struct();
model.Fx = loaded.results.Fx;   % .model (fitlm), .xScale
model.Fy = loaded.results.Fy;

% Inverse-square gap coefficients, refit here for a self-contained script
rzRows = axisBlocks{strcmp(axisBlocks(:,2), 'rz'), 1};
rz = sweepData.rz(rzRows);
predictor = 1 ./ (g0_z - rz).^2;
model.k_Fz   = predictor \ sweepData.Fz(rzRows);
model.k_Taux = predictor \ sweepData.Taux(rzRows);
model.k_Tauy = predictor \ sweepData.Tauy(rzRows);

% Tauz: best available single-axis fit (rx gave the highest LOO R^2,
% 0.412, in the coupling matrix -- still the weakest part of this model)
rxRows = axisBlocks{strcmp(axisBlocks(:,2), 'rx'), 1};
rx_ = sweepData.rx(rxRows);
tauz_ = sweepData.Tauz(rxRows);
[rx_, ord] = sort(rx_); tauz_ = tauz_(ord);
model.Tauz.xScale = max(abs(rx_));
model.Tauz.model = fitlm(rx_ / model.Tauz.xScale, tauz_, 'poly3');

fprintf('Assembled model:\n');
fprintf('  Fx   ~ poly%d(rx)\n', model.Fx.bestOrder);
fprintf('  Fy   ~ poly%d(ry)\n', model.Fy.bestOrder);
fprintf('  Fz   ~ %.4g / (g0_z-rz)^2\n', model.k_Fz);
fprintf('  Taux ~ %.4g / (g0_z-rz)^2\n', model.k_Taux);
fprintf('  Tauy ~ %.4g / (g0_z-rz)^2\n', model.k_Tauy);
fprintf('  Tauz ~ poly3(rx)  [weakest link, LOO R^2=0.41 in isolation]\n\n');

%% Validate against the REAL 6D LHS dataset (all DOF genuinely coupled)
lhsData = loadAllComsolResults();
X = table2array(lhsData(:, {'rx','ry','rz','thx','thy','thz'}));
Ypred = predictTrapOutputs(X, model, g0_z);

outputNames = {'Fx', 'Fy', 'Fz', 'Taux', 'Tauy', 'Tauz'};
figure('Color', 'w', 'Name', 'Assembled model vs real 6D LHS data');
fprintf('%8s  %12s  %12s\n', 'output', 'RMSE', 'R^2');
for i = 1:numel(outputNames)
    name = outputNames{i};
    yActual = lhsData.(name);
    yModel = Ypred.(name);

    rmse = sqrt(mean((yModel - yActual).^2));
    ss_res = sum((yActual - yModel).^2);
    ss_tot = sum((yActual - mean(yActual)).^2);
    r2 = 1 - ss_res / ss_tot;
    fprintf('%8s  %12.4g  %12.3f\n', name, rmse, r2);

    subplot(2, 3, i); hold on; grid on;
    scatter(yActual, yModel, 30, 'filled');
    lims = [min([yActual; yModel]), max([yActual; yModel])];
    plot(lims, lims, 'r--');
    xlabel('actual (6D LHS)'); ylabel('predicted (assembled model)');
    title(sprintf('%s (R^2=%.2f)', name, r2));
end

exportgraphics(gcf, fullfile(scriptDir, 'assembledModel_vs_realLHS.png'));
fprintf('\nSaved plot to assembledModel_vs_realLHS.png\n');

save(fullfile(scriptDir, 'assembledTrapModel.mat'), 'model', 'g0_z');
fprintf('Saved model to assembledTrapModel.mat\n');


function Y = predictTrapOutputs(X, model, g0_z)
    % X: [N x 6] = [rx ry rz thx thy thz]
    rx = X(:,1); rz = X(:,3);
    Y = struct();
    Y.Fx   = predict(model.Fx.model, rx / model.Fx.xScale);
    Y.Fy   = predict(model.Fy.model, X(:,2) / model.Fy.xScale);
    Y.Fz   = model.k_Fz   ./ (g0_z - rz).^2;
    Y.Taux = model.k_Taux ./ (g0_z - rz).^2;
    Y.Tauy = model.k_Tauy ./ (g0_z - rz).^2;
    Y.Tauz = predict(model.Tauz.model, rx / model.Tauz.xScale);
end
