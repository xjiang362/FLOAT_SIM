%% plotFxVsRxByRzBin.m
% Fx vs rx, MST and fringing-model energy-derived overlaid, split into
% rz bins -- same style as the earlier Taux/Tauy by-rz-bin slices.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));
load(fullfile(scriptDir, 'fringingModelForces.mat'));

nBins = 6;
edges = quantile(rz, linspace(0,1,nBins+1));
binIdx = discretize(rz, edges);

figure('Color','w','Name','Fx vs rx, by rz bin (MST vs energy-derived)', 'Position',[50 50 1600 700]);
for b = 1:nBins
    mask = binIdx == b;
    subplot(2,3,b); hold on; grid on;
    scatter(rx(mask)*1e6, Fx_MST(mask), 10, 'b', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'MST');
    scatter(rx(mask)*1e6, Fx_energy(mask), 10, 'r', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'energy-derived');
    xlabel('rx [\mum]'); ylabel('Fx [N]');
    title(sprintf('rz in [%.1f, %.1f] um (n=%d)', edges(b)*1e6, edges(b+1)*1e6, sum(mask)));
    legend('Location','best','FontSize',7);
end
sgtitle('Fx vs rx, by rz bin');
exportgraphics(gcf, fullfile(scriptDir, 'Fx_vs_rx_byRzBin.png'));
fprintf('Saved plot to Fx_vs_rx_byRzBin.png\n');
