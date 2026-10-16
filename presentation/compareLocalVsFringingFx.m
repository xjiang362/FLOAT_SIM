%% compareLocalVsFringingFx.m
% Cross-check between two INDEPENDENT energy-based Fx estimates at the
% same (rx,rz) points (ry=thx=thy=thz=0, from the dedicated fine rx
% sweep):
%   1. LOCAL: smoothed finite-difference on real per-point simulated W
%      (comsol_rxFineSweepEnergy_results.txt, analyzeFxEnergyMethod_smoothed.m)
%   2. GLOBAL: analytic derivative of the fitted fringing-capacitance
%      model (fit on the separate 2000-pt full-6D LHS sweep)
% If these agree, that's real independent validation of the fringing
% model -- not just a fit to the same data it was trained on. Also
% shows MST for reference (known noisy at this scale).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'fringingCapModel.mat'));  % modelFun, p, g0_z, W0, L0

raw = readmatrix(fullfile(scriptDir, 'comsol_rxFineSweepEnergy_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1); rz = raw(:,3);
Fx_MST = raw(:,7);
W = raw(:,15);

rzValues = unique(rz);
nRz = numel(rzValues);

figure('Color','w','Name','Local (data) vs Global (fringing model) vs MST, Fx vs rx', 'Position',[50 50 1800 380]);
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

    % Local: RAW individual data points -- point-to-point finite
    % difference on the actual simulated W (no smoothing/fitting), per
    % the user's request to see the raw energy-method points directly.
    FxLocal = gradient(WSlice, rxSlice);   % +dW/drx convention

    % Global: fringing model, evaluated at the same points (ry=thx=thy=thz=0)
    Xg = [rxSlice, zeros(n,1), rzValues(j)*ones(n,1), zeros(n,1), zeros(n,1), zeros(n,1)];
    h = 1e-6 * xScale;
    Xp = Xg; Xp(:,1) = Xp(:,1) + h;
    Xm = Xg; Xm(:,1) = Xm(:,1) - h;
    FxGlobal = (modelFun(p,Xp) - modelFun(p,Xm)) / (2*h);

    subplot(1,5,j); hold on; grid on;
    plot(rxSlice*1e6, FxLocal, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 4, 'DisplayName', 'local (raw data points)');
    plot(rxSlice*1e6, FxGlobal, 'm-', 'LineWidth', 2, 'DisplayName', 'global (fringing model)');
    plot(rxSlice*1e6, FxMstSlice, 'x', 'Color', [0.6 0.6 0.6], 'DisplayName', 'MST');
    xlabel('rx [\mum]'); ylabel('Fx [N]');
    title(sprintf('rz = %.1f um', rzValues(j)*1e6));
    legend('Location','best','FontSize',6);
end
sgtitle('Fx vs rx: local raw data points vs global fringing-model vs MST');
exportgraphics(gcf, fullfile(scriptDir, 'localVsFringingFx.png'));
fprintf('Saved plot to localVsFringingFx.png\n');
