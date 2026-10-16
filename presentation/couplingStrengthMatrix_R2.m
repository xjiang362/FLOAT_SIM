%% couplingStrengthMatrix_R2.m
% Rebuilds the coupling strength matrix using leave-one-out validated
% R^2 instead of raw range ratio -- range (max-min) can't distinguish
% a small-but-real signal from pure noise scatter of similar spread,
% since both produce a similar range. LOO R^2 directly tests whether
% there's a genuine, reproducible, fittable pattern (linear or
% nonlinear) regardless of its absolute magnitude: real signal scores
% meaningfully positive under honest cross-validation, pure noise
% scores ~0 or negative, same diagnostic already used to validate the
% 6 "dominant axis" pairs earlier.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

loaded = load(fullfile(scriptDir, 'combined300_fittedModels.mat'));
data = loaded.data;
axisBlocks = loaded.axisBlocks;

axisNames = axisBlocks(:,2)';
outputNames = {'Fx', 'Fy', 'Fz', 'Taux', 'Tauy', 'Tauz'};
polyOrders = [2, 3, 5];

nAxes = numel(axisNames);
nOut = numel(outputNames);
r2Matrix = zeros(nOut, nAxes);

for a = 1:nAxes
    rows = axisBlocks{a, 1};
    x = data.(axisNames{a})(rows);
    [x, order] = sort(x);
    n = numel(x);
    xn = x / max(abs(x));

    for o = 1:nOut
        y = data.(outputNames{o})(rows);
        y = y(order);

        bestR2 = -inf;
        for p = polyOrders
            yLOO = zeros(n, 1);
            for i = 1:n
                keep = true(n, 1); keep(i) = false;
                mdl_i = fitlm(xn(keep), y(keep), sprintf('poly%d', p));
                yLOO(i) = predict(mdl_i, xn(i));
            end
            ss_res = sum((y - yLOO).^2);
            ss_tot = sum((y - mean(y)).^2);
            r2 = 1 - ss_res / ss_tot;
            bestR2 = max(bestR2, r2);
        end
        r2Matrix(o, a) = bestR2;
    end
    fprintf('Done axis %s\n', axisNames{a});
end

fprintf('\nLeave-one-out R^2 (best of poly2/poly3/poly5), all 36 pairs:\n');
fprintf('%8s', 'output');
for a = 1:nAxes, fprintf('  %10s', axisNames{a}); end
fprintf('\n');
for o = 1:nOut
    fprintf('%8s', outputNames{o});
    for a = 1:nAxes, fprintf('  %10.3f', r2Matrix(o,a)); end
    fprintf('\n');
end

%% Heatmap
figure('Color', 'w', 'Name', 'Coupling strength matrix (LOO R^2)');
imagesc(r2Matrix, [min(r2Matrix(:), [], 'omitnan'), 1]);
colormap(flipud(hot));
cb = colorbar;
cb.Label.String = 'leave-one-out R^2';
set(gca, 'XTick', 1:nAxes, 'XTickLabel', axisNames, ...
         'YTick', 1:nOut, 'YTickLabel', outputNames);
xlabel('input axis'); ylabel('output');
title('Coupling strength: LOO R^2 per (axis, output) pair');

for o = 1:nOut
    for a = 1:nAxes
        text(a, o, sprintf('%.2f', r2Matrix(o,a)), ...
            'HorizontalAlignment', 'center', 'Color', [0.2 0.2 0.2], 'FontSize', 8);
    end
end

exportgraphics(gcf, fullfile(scriptDir, 'couplingStrengthMatrix_R2.png'));
fprintf('\nSaved heatmap to couplingStrengthMatrix_R2.png\n');
