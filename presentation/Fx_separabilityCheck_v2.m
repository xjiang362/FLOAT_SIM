%% Fx_separabilityCheck_v2.m
% Direct test of whether Fx = k/(g0_z-rz)^2 * f(rx) is separable:
% for each fixed rz, "undo" the assumed gap scaling by multiplying by
% (g0_z-rz)^2, giving an apparent f(rx) at that rz. If the model is
% truly separable, EVERY rz slice's de-scaled curve should collapse
% onto the same f(rx) shape (both in shape AND absolute magnitude,
% since k is a single shared constant) -- no further normalization
% applied, since deviation in scale is just as informative as
% deviation in shape here.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

g0_z = 15e-6;
raw = readmatrix(fullfile(scriptDir, 'comsol_rxRzGrid_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1); rz = raw(:,3); Fx = raw(:,7);

rzValues = unique(rz);   % same 15 values shared across all 6 rx groups
rxValues = unique(rx);

figure('Color', 'w', 'Name', 'De-scaled f(rx) = Fx*(g0_z-rz)^2, one curve per rz');
hold on; grid on;
colors = jet(numel(rzValues));
for j = 1:numel(rzValues)
    mask = abs(rz - rzValues(j)) < 1e-12;
    rxSlice = rx(mask);
    FxSlice = Fx(mask);
    [rxSlice, order] = sort(rxSlice);
    FxSlice = FxSlice(order);

    fApparent = FxSlice * (g0_z - rzValues(j))^2;

    plot(rxSlice*1e6, fApparent, '-o', 'Color', colors(j,:), ...
        'DisplayName', sprintf('rz=%.2fum', rzValues(j)*1e6));
end
xlabel('rx [\mum]'); ylabel('f_{apparent}(rx) = Fx \times (g_0-rz)^2');
title('Separability check: does de-scaled f(rx) collapse across all rz?');
legend('Location', 'eastoutside', 'FontSize', 7);

exportgraphics(gcf, fullfile(scriptDir, 'Fx_separabilityCheck_v2.png'));
fprintf('Saved plot to Fx_separabilityCheck_v2.png\n');

%% Quantify spread: at each rx value, compute the coefficient of
% variation (std/mean) of f_apparent across the 15 rz slices -- large
% values mean the curves are NOT collapsing at that rx.
fprintf('\n%8s  %14s\n', 'rx(um)', 'CV of f_apparent across rz');
for i = 1:numel(rxValues)
    mask = abs(rx - rxValues(i)) < 1e-12;
    fVals = Fx(mask) .* (g0_z - rz(mask)).^2;
    cv = std(fVals) / abs(mean(fVals));
    fprintf('%8.1f  %14.3f\n', rxValues(i)*1e6, cv);
end
