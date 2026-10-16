%% testFxGapModel_onGrid.m
% Tests the pure "shape(rx) * k/(g0_z-rz)^2" model (no thx/thy terms)
% against the rx-rz grid data, where thx=thy=0 by design -- this
% isolates whether the 1/(gap)^2 exponent is correct at all, cleanly
% separated from whether thx/thy coupling terms are also needed
% (which the full 55-point LHS data can't cleanly test, since it's
% missing those terms and has thx/thy varying too).
%
% NOTE: this export's column order is rx,ry,rz,thx,thy,thz,... (NOT
% rz,rx,ry,... like the earlier LHS exports), so it's parsed directly
% here rather than via parseComsolResults.m.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

g0_z = 15e-6;
xScale = 4.5e-5;
Fx_shape = @(rx) 2.20e-8 + 6.28e-7*(rx/xScale) - 9.80e-8*(rx/xScale).^2 ...
    - 1.080e-6*(rx/xScale).^3 + 8.20e-8*(rx/xScale).^4 + 5.30e-7*(rx/xScale).^5;

raw = readmatrix(fullfile(scriptDir, 'comsol_rxRzGrid_results.txt'), 'NumHeaderLines', 5);
if size(raw,2) ~= 14
    error('Expected 14 columns, got %d -- check export format.', size(raw,2));
end

rx = raw(:,1); rz = raw(:,3); thx = raw(:,4); thy = raw(:,5);
Fx = raw(:,7);
n = numel(rx);

fprintf('Loaded %d points. thx range: [%.3g %.3g], thy range: [%.3g %.3g] (should be ~0)\n', ...
    n, min(thx), max(thx), min(thy), max(thy));

%% Fit pure gap model: Fx = k * shape(rx) / (g0_z-rz)^2
predictor = Fx_shape(rx) ./ (g0_z - rz).^2;
k = predictor \ Fx;
Fx_pred = k * predictor;
r2 = 1 - sum((Fx-Fx_pred).^2) / sum((Fx-mean(Fx)).^2);

% LOO validation
yLOO = zeros(n,1);
for i = 1:n
    keep = true(n,1); keep(i) = false;
    k_i = predictor(keep) \ Fx(keep);
    yLOO(i) = k_i * predictor(i);
end
r2_loo = 1 - sum((Fx-yLOO).^2) / sum((Fx-mean(Fx)).^2);

fprintf('\nPure gap model: Fx = k * shape(rx) / (g0_z-rz)^2\n');
fprintf('k = %.4g\n', k);
fprintf('In-sample R^2: %.3f   LOO R^2: %.3f\n', r2, r2_loo);

%% Plot: Fx vs rz, one line per rx value, actual vs model
rxValues = unique(rx);
figure('Color','w','Name','Fx: pure gap model vs rx-rz grid data');
colors = lines(numel(rxValues));
hold on; grid on;
for i = 1:numel(rxValues)
    mask = rx == rxValues(i);
    [rzSorted, order] = sort(rz(mask));
    FxVals = Fx(mask); FxVals = FxVals(order);
    FxModelVals = Fx_pred(mask); FxModelVals = FxModelVals(order);

    scatter(rzSorted*1e6, FxVals, 25, colors(i,:), 'filled', 'DisplayName', sprintf('rx=%.0fum (sim)', rxValues(i)*1e6));
    plot(rzSorted*1e6, FxModelVals, '--', 'Color', colors(i,:), 'LineWidth', 1.2, 'HandleVisibility', 'off');
end
xlabel('rz [\mum]'); ylabel('Fx [N]');
title(sprintf('Fx vs rz at fixed rx (LOO R^2=%.3f)', r2_loo));
legend('Location', 'best');

exportgraphics(gcf, fullfile(scriptDir, 'FxGapModel_onGrid.png'));
fprintf('Saved plot to FxGapModel_onGrid.png\n');
