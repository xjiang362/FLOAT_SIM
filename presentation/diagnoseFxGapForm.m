%% diagnoseFxGapForm.m
% Two diagnostics on the rx-rz grid data (thx=thy=0 throughout) to
% pin down why "shape(rx) * k/(g0_z-rz)^2" failed there:
%   1. Free-exponent fit: is some other power n != 2 the right form?
%   2. Separability check: does Fx(rz), normalized by its own value at
%      rz=0, have the SAME shape at every rx? If yes, Fx really is
%      separable (shape(rx)*g(rz)) and we just have the wrong g; if
%      the normalized curves differ across rx, Fx(rx,rz) doesn't
%      factor this way at all, regardless of what g is.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

g0_z = 15e-6;
xScale = 4.5e-5;
Fx_shape = @(rx) 2.20e-8 + 6.28e-7*(rx/xScale) - 9.80e-8*(rx/xScale).^2 ...
    - 1.080e-6*(rx/xScale).^3 + 8.20e-8*(rx/xScale).^4 + 5.30e-7*(rx/xScale).^5;

raw = readmatrix(fullfile(scriptDir, 'comsol_rxRzGrid_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1); rz = raw(:,3); Fx = raw(:,7);
n = numel(rx);

%% 1. Free-exponent fit: Fx = k * shape(rx) / (g0_z-rz)^p
modelFun = @(b,X) b(1) * Fx_shape(X(:,1)) ./ (g0_z - X(:,2)).^b(2);
X = [rx, rz];
b0 = [6.586e-11, 2];

opts = statset('MaxIter', 1000, 'Display', 'off');
mdl = fitnlm(X, Fx, modelFun, b0, 'Options', opts);
b = mdl.Coefficients.Estimate;

% LOO
yLOO = zeros(n,1);
for i = 1:n
    keep = true(n,1); keep(i) = false;
    try
        mdl_i = fitnlm(X(keep,:), Fx(keep), modelFun, b, 'Options', opts);
        yLOO(i) = predict(mdl_i, X(i,:));
    catch
        yLOO(i) = NaN;
    end
end
valid = ~isnan(yLOO);
r2_loo = 1 - sum((Fx(valid)-yLOO(valid)).^2) / sum((Fx(valid)-mean(Fx(valid))).^2);

fprintf('=== Free-exponent fit: Fx = k*shape(rx)/(g0_z-rz)^p ===\n');
fprintf('k=%.4g  p=%.4g\n', b(1), b(2));
fprintf('In-sample R^2: %.3f   LOO R^2: %.3f  (%d/%d converged)\n\n', ...
    mdl.Rsquared.Ordinary, r2_loo, sum(valid), n);

%% 2. Separability check: normalize each rx-slice by its value at rz=0
rxValues = unique(rx);
figure('Color','w','Name','Separability check: Fx(rz)/Fx(0) at each rx');
hold on; grid on;
colors = lines(numel(rxValues));
for i = 1:numel(rxValues)
    mask = rx == rxValues(i);
    rzSlice = rz(mask);
    FxSlice = Fx(mask);
    [rzSlice, order] = sort(rzSlice);
    FxSlice = FxSlice(order);

    [~, zeroIdx] = min(abs(rzSlice));   % index closest to rz=0
    FxNorm = FxSlice / FxSlice(zeroIdx);

    plot(rzSlice*1e6, FxNorm, '-o', 'Color', colors(i,:), ...
        'DisplayName', sprintf('rx=%.0fum', rxValues(i)*1e6));
end
xlabel('rz [\mum]'); ylabel('Fx(rz) / Fx(rz=0)  [normalized]');
title('Separability check: do normalized curves overlap across rx?');
legend('Location', 'best');
exportgraphics(gcf, fullfile(scriptDir, 'Fx_separabilityCheck.png'));
fprintf('Saved separability plot to Fx_separabilityCheck.png\n');
