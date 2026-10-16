%% analyzeCombined300Sweep.m
% Processes the 300-point single-variable characterization sweep
% (comsol_sweep_combined300_results.txt): splits it back into its 6
% row blocks (rx, ry, rz, thx, thy, thz, 50 points each, all other DOF
% zero in each block), fits each axis's dominant output
% (Fx-rx, Fy-ry, Fz-rz, Taux-thx, Tauy-thy, Tauz-thz) with poly2/3/5,
% and validates with leave-one-out cross-validation -- not just
% in-sample R^2, since the earlier 20-point pilot's poly5 fit looked
% great in-sample but was overfitting (~3:1 data:parameter ratio).
% 50 points per axis gives a much safer ~8:1 ratio for poly5 (6
% coefficients), and LOO gives an honest read on whether that
% flexibility is actually earning its keep.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

resultsFile = fullfile(scriptDir, 'comsol_sweep_combined300_results.txt');
data = parseComsolResults(resultsFile);

if height(data) ~= 300
    error('Expected 300 rows, got %d -- check the export matches the 300-point combined sweep.', height(data));
end

% Row blocks, in the same order the sweep was generated
axisBlocks = {
    1:50,    'rx',  'Fx'
    51:100,  'ry',  'Fy'
    101:150, 'rz',  'Fz'
    151:200, 'thx', 'Taux'
    201:250, 'thy', 'Tauy'
    251:300, 'thz', 'Tauz'
};

polyOrders = [2, 3, 5];
results = struct();

figure('Color', 'w', 'Name', '300-point single-axis sweep: LOO-validated fits');
fprintf('%6s->%-5s', 'axis', 'out');
for p = polyOrders
    fprintf('  R^2(poly%d,LOO)', p);
end
fprintf('\n');

for b = 1:size(axisBlocks, 1)
    rows = axisBlocks{b, 1};
    inName = axisBlocks{b, 2};
    outName = axisBlocks{b, 3};

    x = data.(inName)(rows);
    y = data.(outName)(rows);
    [x, order] = sort(x);
    y = y(order);
    n = numel(x);

    xScale = max(abs(x));
    xn = x / xScale;

    fprintf('%6s->%-5s', inName, outName);
    bestR2 = -inf; bestOrder = polyOrders(1); bestYFit = [];
    for p = polyOrders
        modelSpec = sprintf('poly%d', p);

        % Leave-one-out cross-validation (cheap for OLS)
        yLOO = zeros(n, 1);
        for i = 1:n
            keep = true(n, 1); keep(i) = false;
            mdl_i = fitlm(xn(keep), y(keep), modelSpec);
            yLOO(i) = predict(mdl_i, xn(i));
        end
        ss_res = sum((y - yLOO).^2);
        ss_tot = sum((y - mean(y)).^2);
        r2 = 1 - ss_res / ss_tot;
        fprintf('  %14.3f', r2);

        if r2 > bestR2
            bestR2 = r2;
            bestOrder = p;
        end
    end
    fprintf('   (best: poly%d)\n', bestOrder);

    % Full-data fit at the best order, for plotting and later querying
    mdl = fitlm(xn, y, sprintf('poly%d', bestOrder));
    results.(outName) = struct('model', mdl, 'xScale', xScale, 'bestOrder', bestOrder, 'looR2', bestR2);

    xFine = linspace(min(x), max(x), 200)';
    yFine = predict(mdl, xFine / xScale);

    subplot(2, 3, b); hold on; grid on;
    scatter(x, y, 20, 'filled', 'DisplayName', 'COMSOL data');
    plot(xFine, yFine, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('poly%d fit', bestOrder));
    xlabel(inName); ylabel(outName);
    title(sprintf('%s vs %s (LOO R^2=%.3f)', outName, inName, bestR2));
    legend('Location', 'best');
end

exportgraphics(gcf, fullfile(scriptDir, 'combined300_dominantAxisFits.png'));

save(fullfile(scriptDir, 'combined300_fittedModels.mat'), 'results', 'data', 'axisBlocks');
fprintf('\nSaved fitted models and data to combined300_fittedModels.mat\n');
