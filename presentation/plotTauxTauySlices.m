%% plotTauxTauySlices.m
% Taux/Tauy energy-derived vs MST didn't reconcile as a clean x=y
% relationship -- shows an organized but nonlinear curve instead,
% suggesting the two track real structure differently across thx/thy
% and rz. Since the full-6D data is scattered (LHS, not a grid), bin
% by rz into a handful of slices and plot Taux vs thx (MST vs
% energy-derived overlaid) per bin, same for Tauy vs thy, to see where
% the shapes agree/diverge.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));
load(fullfile(scriptDir, 'energyDerivedForces_v2.mat'));

nBins = 6;
edges = quantile(rz, linspace(0,1,nBins+1));
binIdx = discretize(rz, edges);

%% Taux vs thx, by rz bin
figure('Color','w','Name','Taux vs thx, by rz bin (MST vs energy-derived)', 'Position',[50 50 1600 700]);
for b = 1:nBins
    mask = binIdx == b;
    subplot(2,3,b); hold on; grid on;
    scatter(thx(mask), Taux_MST(mask), 10, 'b', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'MST');
    scatter(thx(mask), Taux_energy(mask), 10, 'r', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'energy-derived');
    xlabel('thx [deg]'); ylabel('Taux [N*m]');
    title(sprintf('rz in [%.1f, %.1f] um (n=%d)', edges(b)*1e6, edges(b+1)*1e6, sum(mask)));
    legend('Location','best','FontSize',7);
end
sgtitle('Taux vs thx, by rz bin');
exportgraphics(gcf, fullfile(scriptDir, 'Taux_vs_thx_byRzBin.png'));
fprintf('Saved plot to Taux_vs_thx_byRzBin.png\n');

%% Tauy vs thy, by rz bin
figure('Color','w','Name','Tauy vs thy, by rz bin (MST vs energy-derived)', 'Position',[50 50 1600 700]);
for b = 1:nBins
    mask = binIdx == b;
    subplot(2,3,b); hold on; grid on;
    scatter(thy(mask), Tauy_MST(mask), 10, 'b', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'MST');
    scatter(thy(mask), Tauy_energy(mask), 10, 'r', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'energy-derived');
    xlabel('thy [deg]'); ylabel('Tauy [N*m]');
    title(sprintf('rz in [%.1f, %.1f] um (n=%d)', edges(b)*1e6, edges(b+1)*1e6, sum(mask)));
    legend('Location','best','FontSize',7);
end
sgtitle('Tauy vs thy, by rz bin');
exportgraphics(gcf, fullfile(scriptDir, 'Tauy_vs_thy_byRzBin.png'));
fprintf('Saved plot to Tauy_vs_thy_byRzBin.png\n');
