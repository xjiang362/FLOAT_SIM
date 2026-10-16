%% plotAllOutputsPerAxis_combined300.m
% For the 300-point single-variable sweep, generates 6 figures (one
% per DOF: rx, ry, rz, thx, thy, thz), each with 6 subplots showing
% how all 6 outputs (Fx, Fy, Fz, Taux, Tauy, Tauz) vary with that one
% DOF (using that axis's own row block, where every other DOF is
% zero). This is the transpose of the earlier "1 output vs all 6
% inputs" plots -- here it's "1 input vs all 6 outputs" per figure.
%
% Translations (rx,ry,rz) shown in um; tilts (thx,thy,thz) shown in
% radians with LaTeX theta labels, matching the convention used for
% the capacitance/probe plots. Poly3 fit overlaid on each subplot.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

loaded = load(fullfile(scriptDir, 'combined300_fittedModels.mat'));
data = loaded.data;
axisBlocks = loaded.axisBlocks;

outputNames = {'Fx', 'Fy', 'Fz', 'Taux', 'Tauy', 'Tauz'};

axisDisplay = struct( ...
    'rx',  struct('scale', 1e6,    'label', 'rx [\mum]',        'titleLabel', 'rx',         'latex', false), ...
    'ry',  struct('scale', 1e6,    'label', 'ry [\mum]',        'titleLabel', 'ry',         'latex', false), ...
    'rz',  struct('scale', 1e6,    'label', 'rz [\mum]',        'titleLabel', 'rz',         'latex', false), ...
    'thx', struct('scale', pi/180, 'label', '$\theta_x$ [rad]', 'titleLabel', '$\theta_x$', 'latex', true), ...
    'thy', struct('scale', pi/180, 'label', '$\theta_y$ [rad]', 'titleLabel', '$\theta_y$', 'latex', true), ...
    'thz', struct('scale', pi/180, 'label', '$\theta_z$ [rad]', 'titleLabel', '$\theta_z$', 'latex', true) ...
);

for b = 1:size(axisBlocks, 1)
    rows = axisBlocks{b, 1};
    inName = axisBlocks{b, 2};
    disp_ = axisDisplay.(inName);

    x = data.(inName)(rows) * disp_.scale;
    [x, order] = sort(x);

    figure('Color', 'w', 'Name', sprintf('All outputs vs %s', inName));
    for k = 1:numel(outputNames)
        outName = outputNames{k};
        y = data.(outName)(rows);
        y = y(order);

        ax = subplot(2, 3, k); hold(ax, 'on'); grid(ax, 'on');
        scatter(ax, x, y, 20, 'filled');
        if disp_.latex
            xlabel(ax, disp_.label, 'Interpreter', 'latex');
            title(ax, sprintf('%s vs %s', outName, disp_.titleLabel), 'Interpreter', 'latex');
        else
            xlabel(ax, disp_.label);
            title(ax, sprintf('%s vs %s', outName, disp_.titleLabel));
        end
        ylabel(ax, outName);
    end

    outFile = fullfile(scriptDir, sprintf('combined300_allOutputs_vs_%s.png', inName));
    exportgraphics(gcf, outFile);
    fprintf('Saved %s\n', outFile);
end
