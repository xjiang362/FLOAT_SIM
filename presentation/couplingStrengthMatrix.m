%% couplingStrengthMatrix.m
% Systematically checks all 36 input-output pairs (6 axes x 6 outputs)
% from the 300-point single-variable sweep, to distinguish real
% cross-coupling from near-zero/noise-floor relationships. For each
% (axis, output) pair, computes the output's range (max-min) within
% that axis's single-variable block, then normalizes each output's row
% by its own largest range across all 6 axes (its dominant axis) --
% giving a 0-1 relative coupling strength: 1.0 = dominant axis, small
% values = likely negligible coupling (that axis barely moves this
% output at all, relative to how much its own dominant axis does).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

loaded = load(fullfile(scriptDir, 'combined300_fittedModels.mat'));
data = loaded.data;
axisBlocks = loaded.axisBlocks;

axisNames = axisBlocks(:,2)';
outputNames = {'Fx', 'Fy', 'Fz', 'Taux', 'Tauy', 'Tauz'};

nAxes = numel(axisNames);
nOut = numel(outputNames);
rangeMatrix = zeros(nOut, nAxes);   % rows = outputs, cols = axes

for a = 1:nAxes
    rows = axisBlocks{a, 1};
    for o = 1:nOut
        y = data.(outputNames{o})(rows);
        rangeMatrix(o, a) = max(y) - min(y);
    end
end

% Normalize each output's row by its own max (dominant-axis) range
dominantRange = max(rangeMatrix, [], 2);
normMatrix = rangeMatrix ./ dominantRange;

fprintf('Raw range (max-min) of each output within each single-axis sweep:\n');
fprintf('%8s', 'output');
for a = 1:nAxes, fprintf('  %10s', axisNames{a}); end
fprintf('\n');
for o = 1:nOut
    fprintf('%8s', outputNames{o});
    for a = 1:nAxes, fprintf('  %10.3g', rangeMatrix(o,a)); end
    fprintf('\n');
end

fprintf('\nNormalized coupling strength (relative to each output''s own dominant axis, 1.0 = dominant):\n');
fprintf('%8s', 'output');
for a = 1:nAxes, fprintf('  %10s', axisNames{a}); end
fprintf('\n');
for o = 1:nOut
    fprintf('%8s', outputNames{o});
    for a = 1:nAxes, fprintf('  %10.4f', normMatrix(o,a)); end
    fprintf('\n');
end

threshold = 0.01;   % below 1% of dominant-axis range -> likely negligible
fprintf('\nPairs below %.0f%% of dominant-axis range (likely negligible/noise-floor coupling):\n', threshold*100);
for o = 1:nOut
    for a = 1:nAxes
        if normMatrix(o,a) < threshold
            fprintf('  %s vs %s: %.4f%% of %s''s dominant range\n', ...
                outputNames{o}, axisNames{a}, normMatrix(o,a)*100, outputNames{o});
        end
    end
end

%% Heatmap visualization (log scale, since ratios span orders of magnitude)
figure('Color', 'w', 'Name', 'Coupling strength matrix');
imagesc(log10(max(normMatrix, 1e-6)));   % floor at 1e-6 to avoid log(0)
colormap(flipud(hot));
cb = colorbar;
cb.Label.String = 'log_{10}(relative coupling strength)';
set(gca, 'XTick', 1:nAxes, 'XTickLabel', axisNames, ...
         'YTick', 1:nOut, 'YTickLabel', outputNames);
xlabel('input axis'); ylabel('output');
title('Coupling strength: output range vs. axis, normalized to dominant axis');

for o = 1:nOut
    for a = 1:nAxes
        text(a, o, sprintf('%.3f', normMatrix(o,a)), ...
            'HorizontalAlignment', 'center', 'Color', [0.2 0.2 0.2], 'FontSize', 8);
    end
end

exportgraphics(gcf, fullfile(scriptDir, 'couplingStrengthMatrix.png'));
fprintf('\nSaved heatmap to couplingStrengthMatrix.png\n');
