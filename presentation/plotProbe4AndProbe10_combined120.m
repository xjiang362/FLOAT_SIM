%% plotProbe4AndProbe10_combined120.m
% Plots Probe 4 and Probe 10 against each of the 6 DOF, using the
% 120-point single-variable sweep (comsol_sweep_combined120_results.txt,
% 20 points/axis, all other DOF zero in each block). Same recipe as
% the earlier capacitance-vs-all-axes plot: shared y-scale per figure,
% first/last point of each sweep dropped (suspected endpoint
% anomalies), poly3 fit overlaid instead of connecting raw points.
%
% NOTE: this export has a DIFFERENT column layout than the earlier
% 15-column format parseComsolResults.m expects -- 14 columns here,
% with the two Floating Potential probes replaced by a single Probe
% 10 at the end. Parsed directly here rather than forcing it through
% the standard parser.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

resultsFile = fullfile(scriptDir, 'comsol_sweep_combined120_results.txt');
raw = readmatrix(resultsFile, 'NumHeaderLines', 5);

if size(raw, 2) ~= 14
    error('Expected 14 columns in %s, got %d -- check the export format.', resultsFile, size(raw, 2));
end

data = table();
data.rz  = raw(:,1);
data.rx  = raw(:,2);
data.ry  = raw(:,3);
data.thx = raw(:,4);
data.thy = raw(:,5);
data.thz = raw(:,6);
data.Fx   = raw(:,7);
data.Fy   = raw(:,8);
data.Fz   = raw(:,9);
data.probe4  = raw(:,10);
data.Taux = raw(:,11);
data.Tauy = raw(:,12);
data.Tauz = raw(:,13);
data.probe10 = raw(:,14);

if height(data) ~= 120
    error('Expected 120 rows, got %d -- check the export matches the 120-point combined sweep.', height(data));
end

axisBlocks = {
    1:20,    'rx'
    21:40,   'ry'
    41:60,   'rz'
    61:80,   'thx'
    81:100,  'thy'
    101:120, 'thz'
};

% Per-axis display scaling and labels: rx/ry/rz shown in um, thx/thy/thz
% converted from degrees to radians and labeled with LaTeX theta symbols.
axisDisplay = struct( ...
    'rx',  struct('scale', 1e6,    'label', 'rx [\mum]',        'titleLabel', 'rx',           'latex', false), ...
    'ry',  struct('scale', 1e6,    'label', 'ry [\mum]',        'titleLabel', 'ry',           'latex', false), ...
    'rz',  struct('scale', 1e6,    'label', 'rz [\mum]',        'titleLabel', 'rz',           'latex', false), ...
    'thx', struct('scale', pi/180, 'label', '$\theta_x$ [rad]', 'titleLabel', '$\theta_x$',   'latex', true), ...
    'thy', struct('scale', pi/180, 'label', '$\theta_y$ [rad]', 'titleLabel', '$\theta_y$',   'latex', true), ...
    'thz', struct('scale', pi/180, 'label', '$\theta_z$ [rad]', 'titleLabel', '$\theta_z$',   'latex', true) ...
);

plotProbe(data, axisBlocks, axisDisplay, 'probe4', 'Cap [pF]', scriptDir);
plotProbe(data, axisBlocks, axisDisplay, 'probe10', 'Probe 10', scriptDir);


function plotProbe(data, axisBlocks, axisDisplay, colName, yLabelStr, scriptDir)

% First pass: trimmed data + quadratic fit per axis, to compute shared y-limits
axisData = cell(size(axisBlocks, 1), 1);
yAllTrimmed = [];
for b = 1:size(axisBlocks, 1)
    rows = axisBlocks{b, 1};
    inName = axisBlocks{b, 2};
    disp_ = axisDisplay.(inName);

    x = data.(inName)(rows) * disp_.scale;
    y = data.(colName)(rows);
    [x, order] = sort(x);
    y = y(order);

    x = x(2:end-1);   % drop first/last point (suspected anomalies)
    y = y(2:end-1);

    mdl = fitlm(x / max(abs(x)), y, 'poly3');
    xFine = linspace(min(x), max(x), 200)';
    yFine = predict(mdl, xFine / max(abs(x)));

    axisData{b} = struct('inName', inName, 'x', x, 'y', y, 'xFine', xFine, 'yFine', yFine, 'disp', disp_);
    yAllTrimmed = [yAllTrimmed; y]; %#ok<AGROW>
end

yLims = [min(yAllTrimmed), max(yAllTrimmed)];
pad = 0.05 * diff(yLims);
yLims = yLims + [-pad, pad];

figure('Color', 'w', 'Name', sprintf('%s vs each DOF', yLabelStr));
for b = 1:size(axisBlocks, 1)
    s = axisData{b};

    ax = subplot(2, 3, b); hold(ax, 'on'); grid(ax, 'on');
    scatter(ax, s.x, s.y, 20, 'filled', 'DisplayName', 'Sim');
    plot(ax, s.xFine, s.yFine, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Fit');
    if s.disp.latex
        xlabel(ax, s.disp.label, 'Interpreter', 'latex');
        title(ax, sprintf('%s vs %s', yLabelStr, s.disp.titleLabel), 'Interpreter', 'latex');
    else
        xlabel(ax, s.disp.label);
        title(ax, sprintf('%s vs %s', yLabelStr, s.disp.titleLabel));
    end
    ylabel(ax, yLabelStr);
    ylim(ax, yLims);
    legend(ax, 'Location', 'best');
end

outFile = fullfile(scriptDir, sprintf('%s_vs_allAxes.png', colName));
exportgraphics(gcf, outFile);
fprintf('Saved %s\n', outFile);

end
