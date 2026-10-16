%% plotFxEnergyRawVsFit.m
% Raw energy-method Fx (point-to-point -gradient(W,rx)) vs the smoothed
% fit (LOO-order-selected polynomial fit to W, differentiated
% analytically), one subplot per rz slice (5 total).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_rxFineSweepEnergy_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1); rz = raw(:,3); W = raw(:,15);

rzValues = unique(rz);
nRz = numel(rzValues);
orders = 2:6;

figure('Color','w','Name','Fx (energy method): raw vs fitted, by rz slice', 'Position', [50 50 1800 380]);

for j = 1:nRz
    mask = abs(rz - rzValues(j)) < 1e-12;
    rxSlice = rx(mask);
    WSlice = W(mask);
    [rxSlice, order] = sort(rxSlice);
    WSlice = WSlice(order);
    n = numel(rxSlice);
    xScale = max(abs(rxSlice));
    xn = rxSlice / xScale;

    % Raw: point-to-point finite difference
    FxRaw = -gradient(WSlice, rxSlice);

    % Fitted: LOO-order-selected poly fit to W, differentiated analytically
    r2loo = zeros(numel(orders),1);
    for oi = 1:numel(orders)
        ord = orders(oi);
        yLOO = zeros(n,1);
        for i = 1:n
            keep = true(n,1); keep(i) = false;
            p = polyfit(xn(keep), WSlice(keep), ord);
            yLOO(i) = polyval(p, xn(i));
        end
        r2loo(oi) = 1 - sum((WSlice-yLOO).^2)/sum((WSlice-mean(WSlice)).^2);
    end
    [~, bi] = max(r2loo);
    bestOrd = orders(bi);
    p = polyfit(xn, WSlice, bestOrd);
    dpdxn = polyder(p);
    FxFit = -polyval(dpdxn, xn) / xScale;

    subplot(1,5,j); hold on; grid on;
    plot(rxSlice*1e6, FxRaw, 'o', 'Color', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerSize', 4, 'DisplayName', 'raw (point-to-point)');
    plot(rxSlice*1e6, FxFit, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('fit (poly%d)', bestOrd));
    xlabel('rx [\mum]'); ylabel('Fx [N]');
    title(sprintf('rz = %.1f um', rzValues(j)*1e6));
    legend('Location','best','FontSize',7);
end
sgtitle('Fx via energy method: raw point-to-point vs smoothed fit');
exportgraphics(gcf, fullfile(scriptDir, 'Fx_energyMethod_rawVsFit.png'));
fprintf('Saved plot to Fx_energyMethod_rawVsFit.png\n');
