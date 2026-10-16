%% analyzeSingleAxisSweeps.m
% Analyzes the three dedicated single-axis sweeps (rx-only, ry-only,
% thz-only, every other DOF held at zero) to check the reported
% "linear near zero, then turns over and decreases" shape for Fx/Fy,
% and the corresponding thz->Tauz relationship. Fits a cubic
% (poly3) model to each -- the right functional form for an odd,
% rise-then-fall restoring force/torque (Fx~k1*rx+k3*rx^3), which a
% plain quadratic model cannot represent at all since its extra terms
% are even-order and contribute nothing to an odd-symmetric curve.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

sweeps = {
    'comsol_sweep_rxOnly_results.txt',  'rx',  'Fx'
    'comsol_sweep_ryOnly_results.txt',  'ry',  'Fy'
    'comsol_sweep_thzOnly_results.txt', 'thz', 'Tauz'
};

figure('Color', 'w', 'Name', 'Single-axis sweeps: poly2 vs poly3 vs poly5');
fprintf('%10s  %10s  %10s  %10s\n', 'relation', 'R^2(quad)', 'R^2(cubic)', 'R^2(quintic)');
for i = 1:size(sweeps, 1)
    data = parseComsolResults(fullfile(scriptDir, sweeps{i,1}));
    x = data.(sweeps{i,2});
    y = data.(sweeps{i,3});

    [x, order] = sort(x);
    y = y(order);

    % Normalize x before fitting -- rx/ry in raw meters (~1e-5) makes
    % x^5 ~1e-25, which rank-deficiencies the design matrix and
    % silently zeros out the higher-order terms (as seen with thz,
    % which is O(1) in degrees and didn't have this problem).
    xScale = max(abs(x));
    xn = x / xScale;

    mdlQuad    = fitlm(xn, y, 'poly2');
    mdlCubic   = fitlm(xn, y, 'poly3');
    mdlQuintic = fitlm(xn, y, 'poly5');

    xFine = linspace(min(x), max(x), 200)';
    xFineN = xFine / xScale;
    yFineQuad    = predict(mdlQuad, xFineN);
    yFineCubic   = predict(mdlCubic, xFineN);
    yFineQuintic = predict(mdlQuintic, xFineN);

    fprintf('%4s->%-5s  %10.3f  %10.3f  %10.3f\n', sweeps{i,2}, sweeps{i,3}, ...
        mdlQuad.Rsquared.Ordinary, mdlCubic.Rsquared.Ordinary, mdlQuintic.Rsquared.Ordinary);

    subplot(1, 3, i); hold on; grid on;
    scatter(x, y, 30, 'filled', 'DisplayName', 'COMSOL data');
    plot(xFine, yFineQuad, 'b--', 'LineWidth', 1, 'DisplayName', 'quadratic');
    plot(xFine, yFineCubic, 'r-', 'LineWidth', 1.5, 'DisplayName', 'cubic');
    plot(xFine, yFineQuintic, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'quintic');
    xlabel(sweeps{i,2}); ylabel(sweeps{i,3});
    title(sprintf('%s vs %s', sweeps{i,3}, sweeps{i,2}));
    legend('Location', 'best');

    fprintf('  %s quintic coefficients: %s\n', sweeps{i,3}, mat2str(mdlQuintic.Coefficients.Estimate', 4));
end
exportgraphics(gcf, fullfile(scriptDir, 'singleAxisSweeps_cubicFit.png'));
