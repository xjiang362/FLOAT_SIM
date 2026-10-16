%% compareFxRxShapesByRz.m
% Fits poly3 to Fx vs rx independently at each of the 15 fixed-rz
% slices (6 points each), and reports the coefficients side by side --
% checks both whether each slice looks odd-symmetric (small
% intercept/x^2 terms relative to x/x^3) and whether the shape is
% consistent across rz or changes systematically.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_rxRzGrid_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1); rz = raw(:,3); Fx = raw(:,7);
xScale = max(abs(rx));

rzValues = unique(rz);
coeffs = zeros(numel(rzValues), 4);   % [intercept, x, x^2, x^3]

fprintf('%10s  %10s  %10s  %10s  %10s\n', 'rz(um)', 'c0', 'c1(x)', 'c2(x^2)', 'c3(x^3)');
for j = 1:numel(rzValues)
    mask = abs(rz - rzValues(j)) < 1e-12;
    rxSlice = rx(mask) / xScale;
    FxSlice = Fx(mask);
    mdl = fitlm(rxSlice, FxSlice, 'poly3');
    c = mdl.Coefficients.Estimate;
    coeffs(j,:) = c';
    fprintf('%10.2f  %10.3g  %10.3g  %10.3g  %10.3g\n', rzValues(j)*1e6, c(1), c(2), c(3), c(4));
end

fprintf('\n--- Odd-symmetry check: |even terms| / |odd terms| ---\n');
fprintf('%10s  %14s\n', 'rz(um)', '(|c0|+|c2|)/(|c1|+|c3|)');
for j = 1:numel(rzValues)
    evenMag = abs(coeffs(j,1)) + abs(coeffs(j,3));
    oddMag  = abs(coeffs(j,2)) + abs(coeffs(j,4));
    fprintf('%10.2f  %14.3f\n', rzValues(j)*1e6, evenMag/oddMag);
end

fprintf('\n--- Coefficient spread across rz (coefficient of variation) ---\n');
labels = {'c0','c1(x)','c2(x^2)','c3(x^3)'};
for k = 1:4
    cv = std(coeffs(:,k)) / abs(mean(coeffs(:,k)));
    fprintf('%8s: CV = %.3f\n', labels{k}, cv);
end

%% Plot: coefficients vs rz
figure('Color','w','Name','poly3(rx) coefficients vs rz');
for k = 1:4
    subplot(2,2,k);
    plot(rzValues*1e6, coeffs(:,k), '-o');
    grid on;
    xlabel('rz [\mum]'); ylabel(labels{k});
    title(sprintf('%s vs rz', labels{k}));
end
exportgraphics(gcf, fullfile(scriptDir, 'Fx_poly3coeffs_vs_rz.png'));
fprintf('\nSaved plot to Fx_poly3coeffs_vs_rz.png\n');
