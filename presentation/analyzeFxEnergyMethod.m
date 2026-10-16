%% analyzeFxEnergyMethod.m
% Computes Fx via the energy-derivative method (Fx = -dW/drx, central
% difference across the fine uniform 2um rx grid) and compares against
% the Maxwell Stress Tensor probe (Fx_MST), at each fixed rz slice.
% MST was confirmed non-convergent under mesh refinement (sign flips);
% this checks whether the energy-derivative estimate is better behaved
% (smooth, consistent sign/shape, roughly odd-symmetric in rx, and
% consistent in shape across rz slices).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_rxFineSweepEnergy_results.txt'), 'NumHeaderLines', 5);
% Columns: 1 rx, 2 ry, 3 rz, 4 thx, 5 thy, 6 thz, 7 Fx(MST), 8 Fy, 9 Fz,
% 10 Probe4(unused/Inf), 11 Taux, 12 Tauy, 13 Tauz, 14 Probe10(unused/Inf),
% 15 W (total electric energy, J)
rx = raw(:,1); rz = raw(:,3);
Fx_MST = raw(:,7);
W = raw(:,15);

rzValues = unique(rz);
nRz = numel(rzValues);

figure('Color','w','Name','Fx: energy-derivative vs MST, by rz slice');
colors = lines(nRz);
Fx_energy_all = cell(nRz,1);
rx_all = cell(nRz,1);

for j = 1:nRz
    mask = abs(rz - rzValues(j)) < 1e-12;
    rxSlice = rx(mask);
    WSlice = W(mask);
    FxMstSlice = Fx_MST(mask);

    [rxSlice, order] = sort(rxSlice);
    WSlice = WSlice(order);
    FxMstSlice = FxMstSlice(order);

    FxEnergySlice = -gradient(WSlice, rxSlice);

    rx_all{j} = rxSlice;
    Fx_energy_all{j} = FxEnergySlice;

    subplot(2,3,j); hold on; grid on;
    plot(rxSlice*1e6, FxEnergySlice, '-o', 'Color', colors(j,:), 'LineWidth', 1.5, 'DisplayName', 'Fx (energy method)');
    plot(rxSlice*1e6, FxMstSlice, '--x', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Fx (MST)');
    xlabel('rx [\mum]'); ylabel('Fx [N]');
    title(sprintf('rz = %.1f um', rzValues(j)*1e6));
    legend('Location','best','FontSize',7);
end
sgtitle('Fx via energy-derivative method vs MST, per rz slice');
exportgraphics(gcf, fullfile(scriptDir, 'Fx_energyMethod_vs_MST.png'));
fprintf('Saved plot to Fx_energyMethod_vs_MST.png\n');

%% Odd-symmetry check on the energy-derivative estimate (poly3 fit per slice)
fprintf('\n--- poly3 fit of Fx(energy) vs rx, per rz slice ---\n');
fprintf('%10s  %10s  %10s  %10s  %10s  %14s\n', 'rz(um)', 'c0', 'c1(x)', 'c2(x^2)', 'c3(x^3)', '(|c0|+|c2|)/(|c1|+|c3|)');
xScale = max(abs(rx_all{1}));
coeffs = zeros(nRz,4);
for j = 1:nRz
    mdl = fitlm(rx_all{j}/xScale, Fx_energy_all{j}, 'poly3');
    c = mdl.Coefficients.Estimate;
    coeffs(j,:) = c';
    evenMag = abs(c(1)) + abs(c(3));
    oddMag = abs(c(2)) + abs(c(4));
    fprintf('%10.1f  %10.3g  %10.3g  %10.3g  %10.3g  %14.3f\n', ...
        rzValues(j)*1e6, c(1), c(2), c(3), c(4), evenMag/oddMag);
end

%% Separability check (same de-scaling approach as Fx_separabilityCheck_v2.m)
g0_z = 15e-6;
figure('Color','w','Name','Energy-method Fx: de-scaled f(rx)=Fx*(g0_z-rz)^2, by rz');
hold on; grid on;
for j = 1:nRz
    fApparent = Fx_energy_all{j} * (g0_z - rzValues(j))^2;
    plot(rx_all{j}*1e6, fApparent, '-o', 'Color', colors(j,:), 'DisplayName', sprintf('rz=%.1fum', rzValues(j)*1e6));
end
xlabel('rx [\mum]'); ylabel('f_{apparent}(rx) = Fx_{energy} \times (g_0-rz)^2');
title('Separability check (energy method): does de-scaled f(rx) collapse across rz?');
legend('Location','eastoutside','FontSize',7);
exportgraphics(gcf, fullfile(scriptDir, 'Fx_energyMethod_separabilityCheck.png'));
fprintf('\nSaved plot to Fx_energyMethod_separabilityCheck.png\n');

fprintf('\n--- CV of f_apparent across rz, at each rx (energy method) ---\n');
rxCommon = rx_all{1};   % same grid at every rz slice
fApparentMat = zeros(numel(rxCommon), nRz);
for j = 1:nRz
    fApparentMat(:,j) = Fx_energy_all{j} * (g0_z - rzValues(j))^2;
end
fprintf('%8s  %14s\n', 'rx(um)', 'CV across rz');
for i = 1:numel(rxCommon)
    cv = std(fApparentMat(i,:)) / abs(mean(fApparentMat(i,:)));
    fprintf('%8.1f  %14.3f\n', rxCommon(i)*1e6, cv);
end
