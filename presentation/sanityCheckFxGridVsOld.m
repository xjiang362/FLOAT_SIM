%% sanityCheckFxGridVsOld.m
% Sanity check: does the new rx-rz grid's rz~=0 slice match the OLD
% single-axis poly5(rx) model it's supposed to be consistent with?
% If these disagree even at rz=0 (nominally the same physical
% condition as the original single-axis characterization), that
% points to a data/pipeline issue rather than genuine non-separable
% physics.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

xScale = 4.5e-5;
Fx_shape_OLD = @(rx) 2.20e-8 + 6.28e-7*(rx/xScale) - 9.80e-8*(rx/xScale).^2 ...
    - 1.080e-6*(rx/xScale).^3 + 8.20e-8*(rx/xScale).^4 + 5.30e-7*(rx/xScale).^5;

raw = readmatrix(fullfile(scriptDir, 'comsol_rxRzGrid_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1); rz = raw(:,3); Fx = raw(:,7);

rzValues = unique(rz);
[~, idxZero] = min(abs(rzValues));
rzNear0 = rzValues(idxZero);
fprintf('Closest rz to 0 in this grid: %.4g um\n', rzNear0*1e6);

mask = abs(rz - rzNear0) < 1e-12;
rxSlice = rx(mask);
FxSlice = Fx(mask);
[rxSlice, order] = sort(rxSlice);
FxSlice = FxSlice(order);

fprintf('\n%10s  %14s  %14s  %10s\n', 'rx(um)', 'Fx (grid, actual)', 'Fx (OLD model)', 'ratio');
for i = 1:numel(rxSlice)
    modelVal = Fx_shape_OLD(rxSlice(i));
    fprintf('%10.2f  %14.4g  %14.4g  %10.3f\n', rxSlice(i)*1e6, FxSlice(i), modelVal, FxSlice(i)/modelVal);
end

figure('Color','w','Name','Sanity check: grid rz~0 vs OLD poly5 model');
hold on; grid on;
xFine = linspace(min(rxSlice), max(rxSlice), 200)';
plot(xFine*1e6, Fx_shape_OLD(xFine), 'r-', 'LineWidth', 1.5, 'DisplayName', 'OLD poly5 model (single-axis fit)');
scatter(rxSlice*1e6, FxSlice, 40, 'filled', 'DisplayName', sprintf('NEW grid data at rz=%.2fum', rzNear0*1e6));
xlabel('rx [\mum]'); ylabel('Fx [N]');
title('Sanity check: does new grid data match old single-axis model at rz~0?');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'sanityCheck_FxGrid_vs_OLD.png'));
fprintf('\nSaved plot to sanityCheck_FxGrid_vs_OLD.png\n');
