%% fitFxFyGapScaled.m
% Fx = k_Fx / (g0_z-rz)^2 * poly5_shape(rx)
% Fy = k_Fy / (g0_z-rz)^2 * poly3_shape(ry)
% The poly5/poly3 SHAPES are the already-validated single-axis fits,
% held fixed (not refit); only k_Fx and k_Fy are new free parameters.
% Since the shape is fixed, k enters linearly (y = k * predictor), so
% this is simple 1-parameter linear regression -- much more robust
% than the earlier 6-7 parameter nonlinear coupled fits, which
% overfit catastrophically (LOO R^2 = -5907).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

g0_z = 15e-6;
xScale = 4.5e-5;
yScale = 4.5e-5;

Fx_shape = @(rx) 2.20e-8 + 6.28e-7*(rx/xScale) - 9.80e-8*(rx/xScale).^2 ...
    - 1.080e-6*(rx/xScale).^3 + 8.20e-8*(rx/xScale).^4 + 5.30e-7*(rx/xScale).^5;
Fy_shape = @(ry) 1.06e-8 + 3.061e-7*(ry/yScale) - 2.09e-8*(ry/yScale).^2 ...
    - 2.364e-7*(ry/yScale).^3;

lhsData = loadAllComsolResults();
rz = lhsData.rz; rx = lhsData.rx; ry = lhsData.ry;
n = height(lhsData);

fprintf('%8s  %12s  %12s  %12s\n', 'output', 'k (fit)', 'in-sample R^2', 'LOO R^2');

%% Fx
predFx = Fx_shape(rx) ./ (g0_z - rz).^2;
Fx_actual = lhsData.Fx;
k_Fx = predFx \ Fx_actual;
Fx_pred = k_Fx * predFx;
r2_fx = 1 - sum((Fx_actual-Fx_pred).^2)/sum((Fx_actual-mean(Fx_actual)).^2);

yLOO = zeros(n,1);
for i = 1:n
    keep = true(n,1); keep(i) = false;
    k_i = predFx(keep) \ Fx_actual(keep);
    yLOO(i) = k_i * predFx(i);
end
r2_fx_loo = 1 - sum((Fx_actual-yLOO).^2)/sum((Fx_actual-mean(Fx_actual)).^2);
fprintf('%8s  %12.4g  %12.3f  %12.3f\n', 'Fx', k_Fx, r2_fx, r2_fx_loo);

%% Fy
predFy = Fy_shape(ry) ./ (g0_z - rz).^2;
Fy_actual = lhsData.Fy;
k_Fy = predFy \ Fy_actual;
Fy_pred = k_Fy * predFy;
r2_fy = 1 - sum((Fy_actual-Fy_pred).^2)/sum((Fy_actual-mean(Fy_actual)).^2);

yLOO = zeros(n,1);
for i = 1:n
    keep = true(n,1); keep(i) = false;
    k_i = predFy(keep) \ Fy_actual(keep);
    yLOO(i) = k_i * predFy(i);
end
r2_fy_loo = 1 - sum((Fy_actual-yLOO).^2)/sum((Fy_actual-mean(Fy_actual)).^2);
fprintf('%8s  %12.4g  %12.3f  %12.3f\n', 'Fy', k_Fy, r2_fy, r2_fy_loo);

%% Plots
figure('Color','w','Name','Fx, Fy gap-scaled model vs real 6D LHS data');
subplot(1,2,1); hold on; grid on;
scatter(Fx_actual, Fx_pred, 30, 'filled');
lims = [min([Fx_actual;Fx_pred]), max([Fx_actual;Fx_pred])];
plot(lims, lims, 'r--');
xlabel('actual Fx'); ylabel('model Fx');
title(sprintf('Fx (R^2=%.2f, LOO R^2=%.2f)', r2_fx, r2_fx_loo));

subplot(1,2,2); hold on; grid on;
scatter(Fy_actual, Fy_pred, 30, 'filled');
lims = [min([Fy_actual;Fy_pred]), max([Fy_actual;Fy_pred])];
plot(lims, lims, 'r--');
xlabel('actual Fy'); ylabel('model Fy');
title(sprintf('Fy (R^2=%.2f, LOO R^2=%.2f)', r2_fy, r2_fy_loo));

exportgraphics(gcf, fullfile(scriptDir, 'FxFy_gapScaled_fit.png'));
fprintf('\nSaved plot to FxFy_gapScaled_fit.png\n');
