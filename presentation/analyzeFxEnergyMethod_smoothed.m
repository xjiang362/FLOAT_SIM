%% analyzeFxEnergyMethod_smoothed.m
% Raw point-to-point differencing of W(rx) (analyzeFxEnergyMethod.m) was
% dominated by mesh rebuild-to-rebuild noise (implied Fx noise floor
% ~5e-7 N, larger than the true signal ~1e-8 to 1e-7 N) -- naive finite
% differencing amplifies uncorrelated per-point mesh noise. Fix: fit a
% smooth polynomial through ALL 41 points of each rz slice (order
% picked via LOO, same convention as analyzeCombined300Sweep.m), then
% differentiate the FIT analytically: Fx = -d/drx [poly(rx)].

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_rxFineSweepEnergy_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1); rz = raw(:,3);
Fx_MST = raw(:,7);
W = raw(:,15);

rzValues = unique(rz);
nRz = numel(rzValues);
orders = 2:6;

figure('Color','w','Name','Fx: smoothed energy-derivative vs MST, by rz slice');
colors = lines(nRz);
Fx_smooth_all = cell(nRz,1);
rx_all = cell(nRz,1);
bestOrders = zeros(nRz,1);

for j = 1:nRz
    mask = abs(rz - rzValues(j)) < 1e-12;
    rxSlice = rx(mask);
    WSlice = W(mask);
    FxMstSlice = Fx_MST(mask);
    [rxSlice, order] = sort(rxSlice);
    WSlice = WSlice(order);
    FxMstSlice = FxMstSlice(order);
    n = numel(rxSlice);
    xScale = max(abs(rxSlice));
    xn = rxSlice / xScale;

    % LOO order selection on W fit itself
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
    bestOrders(j) = bestOrd;

    p = polyfit(xn, WSlice, bestOrd);
    dpdxn = polyder(p);
    % dW/drx = dW/dxn * dxn/drx = polyval(dpdxn, xn) / xScale
    FxSmooth = -polyval(dpdxn, xn) / xScale;

    rx_all{j} = rxSlice;
    Fx_smooth_all{j} = FxSmooth;

    subplot(2,3,j); hold on; grid on;
    plot(rxSlice*1e6, FxSmooth, '-', 'Color', colors(j,:), 'LineWidth', 2, 'DisplayName', sprintf('Fx (energy, poly%d smooth)', bestOrd));
    plot(rxSlice*1e6, FxMstSlice, '--x', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Fx (MST)');
    xlabel('rx [\mum]'); ylabel('Fx [N]');
    title(sprintf('rz = %.1f um (LOO best order=%d)', rzValues(j)*1e6, bestOrd));
    legend('Location','best','FontSize',7);
end
sgtitle('Fx via SMOOTHED energy-derivative (poly fit + analytic diff) vs MST');
exportgraphics(gcf, fullfile(scriptDir, 'Fx_energyMethod_smoothed_vs_MST.png'));
fprintf('Saved plot to Fx_energyMethod_smoothed_vs_MST.png\n');
fprintf('LOO-selected poly orders per rz slice: %s\n', mat2str(bestOrders'));

%% Separability check on smoothed estimate
g0_z = 15e-6;
figure('Color','w','Name','Smoothed energy-method Fx: de-scaled f(rx), by rz');
hold on; grid on;
for j = 1:nRz
    fApparent = Fx_smooth_all{j} * (g0_z - rzValues(j))^2;
    plot(rx_all{j}*1e6, fApparent, '-', 'Color', colors(j,:), 'LineWidth', 1.5, 'DisplayName', sprintf('rz=%.1fum', rzValues(j)*1e6));
end
xlabel('rx [\mum]'); ylabel('f_{apparent}(rx) = Fx_{energy,smooth} \times (g_0-rz)^2');
title('Separability check (smoothed energy method): does de-scaled f(rx) collapse across rz?');
legend('Location','eastoutside','FontSize',7);
exportgraphics(gcf, fullfile(scriptDir, 'Fx_energyMethod_smoothed_separabilityCheck.png'));
fprintf('Saved plot to Fx_energyMethod_smoothed_separabilityCheck.png\n');

fprintf('\n--- CV of f_apparent across rz, at each rx (smoothed energy method) ---\n');
rxCommon = rx_all{1};
fApparentMat = zeros(numel(rxCommon), nRz);
for j = 1:nRz
    fApparentMat(:,j) = Fx_smooth_all{j} * (g0_z - rzValues(j))^2;
end
fprintf('%8s  %14s\n', 'rx(um)', 'CV across rz');
for i = 1:numel(rxCommon)
    cv = std(fApparentMat(i,:)) / abs(mean(fApparentMat(i,:)));
    fprintf('%8.1f  %14.3f\n', rxCommon(i)*1e6, cv);
end
