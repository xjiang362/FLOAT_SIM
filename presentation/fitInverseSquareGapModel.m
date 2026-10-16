%% fitInverseSquareGapModel.m
% Fits the physically-motivated inverse-square gap model
%   Fz   = k_Fz   / (g0_z - rz)^2
%   Taux = k_Taux / (g0_z - rz)^2
%   Tauy = k_Tauy / (g0_z - rz)^2
% to the rz-only single-axis sweep data, as an alternative to the
% polynomial fits that topped out around LOO R^2 ~ 0.40 for these
% three outputs. Since g0_z is already known (15um, the nominal
% vertical gap), k enters linearly once you transform the predictor to
% 1/(g0_z-rz)^2 -- so this is fit as simple linear regression through
% the origin (no intercept, matching the proposed form exactly), with
% the same leave-one-out validation used throughout for a fair
% comparison against the polynomial fits.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

loaded = load(fullfile(scriptDir, 'combined300_fittedModels.mat'));
data = loaded.data;
axisBlocks = loaded.axisBlocks;

g0_z = 15e-6;   % [m] nominal vertical gap, matches sampleTrapPoseLHS.m

rows = axisBlocks{strcmp(axisBlocks(:,2), 'rz'), 1};
rz = data.rz(rows);
[rz, order] = sort(rz);
n = numel(rz);

predictor = 1 ./ (g0_z - rz).^2;   % the nonlinear part, now a known regressor

outputNames = {'Fz', 'Taux', 'Tauy'};
figure('Color', 'w', 'Name', 'Inverse-square gap model vs rz');
fprintf('%8s  %12s  %12s  %12s\n', 'output', 'k (fit)', 'LOO R^2', 'poly5 R^2 (ref)');

% Reference poly5 LOO R^2 for comparison (same recipe as before)
refR2 = struct('Fz', 0.408, 'Taux', 0.407, 'Tauy', 0.404);

for i = 1:numel(outputNames)
    name = outputNames{i};
    y = data.(name)(rows);
    y = y(order);

    % Full-data fit (no intercept, matching the proposed proportional form)
    k = predictor \ y;   % least-squares solution of y = k*predictor

    % Leave-one-out cross-validation
    yLOO = zeros(n, 1);
    for j = 1:n
        keep = true(n, 1); keep(j) = false;
        k_j = predictor(keep) \ y(keep);
        yLOO(j) = k_j * predictor(j);
    end
    ss_res = sum((y - yLOO).^2);
    ss_tot = sum((y - mean(y)).^2);
    r2 = 1 - ss_res / ss_tot;

    fprintf('%8s  %12.4g  %12.3f  %12.3f\n', name, k, r2, refR2.(name));

    rzFine = linspace(min(rz), max(rz), 300)';
    yFine = k ./ (g0_z - rzFine).^2;

    subplot(1, 3, i); hold on; grid on;
    scatter(rz*1e6, y, 25, 'filled', 'DisplayName', 'Sim');
    plot(rzFine*1e6, yFine, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Inverse-square fit');
    xlabel('rz [\mum]'); ylabel(name);
    title(sprintf('%s (LOO R^2=%.3f)', name, r2));
    legend('Location', 'best');
end

exportgraphics(gcf, fullfile(scriptDir, 'inverseSquareGapModel_fit.png'));
fprintf('\nSaved plot to inverseSquareGapModel_fit.png\n');
