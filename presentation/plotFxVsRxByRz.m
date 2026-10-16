%% plotFxVsRxByRz.m
% Raw Fx vs rx, one curve per fixed rz value (no de-scaling/model
% assumptions applied) -- direct look at how the rx-shape changes
% across rz.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_rxRzGrid_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1); rz = raw(:,3); Fx = raw(:,7);

rzValues = unique(rz);

figure('Color', 'w', 'Name', 'Fx vs rx, one curve per rz');
hold on; grid on;
colors = jet(numel(rzValues));
for j = 1:numel(rzValues)
    mask = abs(rz - rzValues(j)) < 1e-12;
    rxSlice = rx(mask);
    FxSlice = Fx(mask);
    [rxSlice, order] = sort(rxSlice);
    FxSlice = FxSlice(order);

    plot(rxSlice*1e6, FxSlice, '-o', 'Color', colors(j,:), ...
        'DisplayName', sprintf('rz=%.2fum', rzValues(j)*1e6));
end
xlabel('rx [\mum]'); ylabel('Fx [N]');
title('Fx vs rx at each fixed rz');
legend('Location', 'eastoutside', 'FontSize', 7);

exportgraphics(gcf, fullfile(scriptDir, 'Fx_vs_rx_by_rz.png'));
fprintf('Saved plot to Fx_vs_rx_by_rz.png\n');
