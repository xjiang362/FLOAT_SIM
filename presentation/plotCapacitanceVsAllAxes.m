%% plotCapacitanceVsAllAxes.m
% Plots capacitance (capacitance_pF, "probe4" / es.Q0_5/es.V0_5) against
% each of the 6 DOF, using the 300-point single-variable sweep (each
% row block has only one DOF nonzero, so this shows capacitance's
% dependence on that axis in clean isolation, same row-block structure
% as analyzeCombined300Sweep.m).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

data = parseComsolResults(fullfile(scriptDir, 'comsol_sweep_combined300_results.txt'));

axisBlocks = {
    1:50,    'rx'
    51:100,  'ry'
    101:150, 'rz'
    151:200, 'thx'
    201:250, 'thy'
    251:300, 'thz'
};

% First pass: gather each axis's trimmed (endpoint-dropped) data and
% cubic fit, so the shared y-limits can be computed before plotting.
axisData = cell(size(axisBlocks, 1), 1);
yAllTrimmed = [];
for b = 1:size(axisBlocks, 1)
    rows = axisBlocks{b, 1};
    inName = axisBlocks{b, 2};

    x = data.(inName)(rows);
    y = data.capacitance_pF(rows);
    [x, order] = sort(x);
    y = y(order);

    % Discard the first and last point of each sweep (suspected
    % anomalies at the extreme ends of the range)
    x = x(2:end-1);
    y = y(2:end-1);

    mdl = fitlm(x / max(abs(x)), y, 'poly3');
    xFine = linspace(min(x), max(x), 200)';
    yFine = predict(mdl, xFine / max(abs(x)));

    axisData{b} = struct('inName', inName, 'x', x, 'y', y, 'xFine', xFine, 'yFine', yFine);
    yAllTrimmed = [yAllTrimmed; y]; %#ok<AGROW>
end

yLims = [min(yAllTrimmed), max(yAllTrimmed)];
pad = 0.05 * diff(yLims);
yLims = yLims + [-pad, pad];

figure('Color', 'w', 'Name', 'Capacitance vs each DOF');
for b = 1:size(axisBlocks, 1)
    s = axisData{b};

    subplot(2, 3, b); hold on; grid on;
    scatter(s.x, s.y, 20, 'filled', 'DisplayName', 'data (endpoints dropped)');
    plot(s.xFine, s.yFine, 'r-', 'LineWidth', 1.5, 'DisplayName', 'poly3 fit');
    xlabel(s.inName); ylabel('capacitance [pF]');
    ylim(yLims);
    title(sprintf('capacitance vs %s', s.inName));
    legend('Location', 'best');
end

exportgraphics(gcf, fullfile(scriptDir, 'capacitance_vs_allAxes.png'));
fprintf('Saved plot to capacitance_vs_allAxes.png\n');
