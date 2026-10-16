%% processInitialBatch40Results.m
% Loads the COMSOL-exported results for the 40-point initial LHS batch
% (comsol_lhs_initial40_results.txt), organizes them into a table with
% consistent naming, sanity-checks the ranges, saves a clean .mat for
% later use (e.g. GP fitting), and gives a quick first-look plot.
%
% The raw export's input columns are ordered [rz, rx, ry, thx, thy,
% thz] -- NOT rx,ry,rz -- reordered here to match the project's usual
% rx,ry,rz,thx,thy,thz convention. "probe4" is an unlabeled Global
% Variable Probe in the export; rename it below once you confirm what
% quantity it corresponds to in COMSOL.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
resultsFile = fullfile(scriptDir, 'comsol_lhs_initial40_results.txt');

raw = readmatrix(resultsFile, 'NumHeaderLines', 5);

if size(raw, 2) ~= 15
    error('Expected 15 columns in %s, got %d -- check the export format hasn''t changed.', ...
        resultsFile, size(raw, 2));
end

rz  = raw(:,1);
rx  = raw(:,2);
ry  = raw(:,3);
thx = raw(:,4);
thy = raw(:,5);
thz = raw(:,6);
Fx      = raw(:,7);
Fy      = raw(:,8);
Fz      = raw(:,9);
probe4  = raw(:,10);   % TODO: rename once you confirm what this Global Variable Probe evaluates
Taux    = raw(:,11);
Tauy    = raw(:,12);
Tauz    = raw(:,13);
Vfloat8 = raw(:,14);
Vfloat9 = raw(:,15);

data = table(rx, ry, rz, thx, thy, thz, Fx, Fy, Fz, probe4, Taux, Tauy, Tauz, Vfloat8, Vfloat9);

fprintf('Loaded %d combinations from %s\n\n', height(data), resultsFile);

%% Sanity checks
fprintf('%10s  %14s  %14s\n', 'var', 'min', 'max');
for v = data.Properties.VariableNames
    val = data.(v{1});
    fprintf('%10s  %14.4g  %14.4g\n', v{1}, min(val), max(val));
end

badRows = any(~isfinite(table2array(data)), 2);
if any(badRows)
    warning('%d row(s) contain NaN/Inf -- check combination(s): %s', ...
        sum(badRows), mat2str(find(badRows)'));
end

%% Save cleaned data for later use (e.g. GP fitting)
outFile = fullfile(scriptDir, 'comsol_lhs_initial40_results.mat');
save(outFile, 'data');
fprintf('\nSaved cleaned table to %s\n', outFile);

%% Quick first-look plot: each output against its most directly-coupled
% input (Fx-rx, Fy-ry, Fz-rz, Taux-thx, Tauy-thy, Tauz-thz). Only a
% single-axis slice through a 6-D input space, so scatter (not a
% clean curve) is expected -- useful as a sanity check, not a model.
inPairs  = {'rx','ry','rz','thx','thy','thz'};
outPairs = {'Fx','Fy','Fz','Taux','Tauy','Tauz'};

figure('Color', 'w', 'Name', 'Initial 40-point batch: dominant-axis coupling');
for i = 1:6
    subplot(2, 3, i);
    scatter(data.(inPairs{i}), data.(outPairs{i}), 30, 'filled');
    grid on;
    xlabel(inPairs{i});
    ylabel(outPairs{i});
    title(sprintf('%s vs %s', outPairs{i}, inPairs{i}));
end
exportgraphics(gcf, fullfile(scriptDir, 'initialBatch40_dominantAxisCoupling.png'));
