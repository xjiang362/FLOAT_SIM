%% plotFxVsRxWithSlope.m
% Plots Fx vs rx using the already-fitted, LOO-validated model from
% combined300_fittedModels.mat, and extracts the slope (dFx/drx) at
% rx=0 -- the local restoring-force stiffness at the equilibrium
% point. For a polynomial fit, the linear coefficient IS the exact
% derivative at x=0 (all higher-order terms vanish there), so this is
% read straight off the fitted coefficients, not a finite-difference
% approximation.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

loaded = load(fullfile(scriptDir, 'combined300_fittedModels.mat'));
results = loaded.results;
data = loaded.data;
axisBlocks = loaded.axisBlocks;

fxResult = results.Fx;
mdl = fxResult.model;
xScale = fxResult.xScale;

% Linear coefficient of the normalized-x fit is the exact slope at
% xn=0; converting back to physical rx units divides by xScale (chain
% rule: d/dx = (1/xScale) * d/dxn).
coeffs = mdl.Coefficients.Estimate;
slopeNorm = coeffs(2);          % coefficient of the x^1 term
slopeAtZero = slopeNorm / xScale;   % [N/m]

% Standard error on that coefficient, propagated the same way, for an
% honest uncertainty on the slope estimate
slopeSE = mdl.Coefficients.SE(2) / xScale;

fprintf('Fx vs rx (poly%d fit, LOO R^2=%.3f)\n', fxResult.bestOrder, fxResult.looR2);
fprintf('Slope at rx=0: dFx/drx = %.4g +/- %.4g N/m\n', slopeAtZero, slopeSE);

%% Plot
rowsRx = axisBlocks{strcmp(axisBlocks(:,2), 'rx'), 1};
x = data.rx(rowsRx);
y = data.Fx(rowsRx);
[x, order] = sort(x);
y = y(order);

xFine = linspace(min(x), max(x), 300)';
yFine = predict(mdl, xFine / xScale);

% Tangent line at rx=0, drawn over a modest span around the
% equilibrium point for visual reference (not the full data range,
% since the tangent is only a good local approximation there)
tangentSpan = 0.2 * max(abs(x));
xTangent = [-tangentSpan, tangentSpan];
yTangent = predict(mdl, 0) + slopeAtZero * xTangent;

figure('Color', 'w', 'Name', 'Fx vs rx');
hold on; grid on;
scatter(x, y, 25, 'filled', 'DisplayName', 'COMSOL data');
plot(xFine, yFine, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('poly%d fit', fxResult.bestOrder));
plot(xTangent, yTangent, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('tangent at rx=0 (%.3g N/m)', slopeAtZero));
plot(0, predict(mdl, 0), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8, 'DisplayName', 'rx=0');
xlabel('rx [m]'); ylabel('Fx [N]');
title(sprintf('Fx vs rx (LOO R^2=%.3f, slope at 0 = %.4g N/m)', fxResult.looR2, slopeAtZero));
legend('Location', 'best');

exportgraphics(gcf, fullfile(scriptDir, 'Fx_vs_rx_slope.png'));
fprintf('Saved plot to Fx_vs_rx_slope.png\n');
