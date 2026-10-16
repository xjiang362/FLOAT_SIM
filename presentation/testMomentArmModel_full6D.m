%% testMomentArmModel_full6D.m
% Tests the moment-arm model for Taux/Tauy (fit on the 150-pt dedicated
% rz/thx/thy-only sweep) against the full 2000-pt scattered 6D LHS
% dataset (full6D_energyData.mat) -- a genuine out-of-sample test, since
% this dataset was never used to fit the moment-arm coefficients.
% Uses the REAL Fz_MST from this dataset (not a fitted Fz model), to
% isolate testing just the moment-arm relationship itself:
%   Leqx = a0 + a1*rz + a2*thx + a3*thy
%   Leqy = -(a0 + a1*rz + a2*thy + a3*thx)
%   Taux_pred = Leqx .* Fz_MST,  Tauy_pred = Leqy .* Fz_MST

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));   % rx,ry,rz,thx,thy,thz,Fx_MST...Tauz_MST,W
load(fullfile(scriptDir, 'tauxTauyMomentArm.mat'));   % coeffs, a0,a1,a2,a3

thx_rad = thx*pi/180; thy_rad = thy*pi/180;

Leqx = a0 + a1*rz + a2*thx_rad + a3*thy_rad;
Leqy = -(a0 + a1*rz + a2*thy_rad + a3*thx_rad);

Taux_pred = Leqx .* Fz_MST;
Tauy_pred = Leqy .* Fz_MST;

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
r2_taux = r2fun(Taux_MST, Taux_pred);
r2_tauy = r2fun(Tauy_MST, Tauy_pred);
fprintf('Full 2000-pt LHS dataset (out-of-sample test):\n');
fprintf('R^2 (Taux) = %.4f\n', r2_taux);
fprintf('R^2 (Tauy) = %.4f\n', r2_tauy);

%% Plot: predicted vs actual scatter
figure('Color','w','Name','Moment-arm model vs MST: full 2000-pt LHS dataset', 'Position',[50 50 1200 550]);

subplot(1,2,1); hold on; grid on;
scatter(Taux_MST, Taux_pred, 10, 'filled', 'MarkerFaceAlpha', 0.3);
lims = [min(Taux_MST) max(Taux_MST)];
plot(lims, lims, 'r--', 'LineWidth', 1.5);
xlabel('Taux (MST)'); ylabel('Taux (moment-arm model)');
title(sprintf('Taux: predicted vs actual (R^2=%.4f)', r2_taux));

subplot(1,2,2); hold on; grid on;
scatter(Tauy_MST, Tauy_pred, 10, 'filled', 'MarkerFaceAlpha', 0.3);
lims = [min(Tauy_MST) max(Tauy_MST)];
plot(lims, lims, 'r--', 'LineWidth', 1.5);
xlabel('Tauy (MST)'); ylabel('Tauy (moment-arm model)');
title(sprintf('Tauy: predicted vs actual (R^2=%.4f)', r2_tauy));

exportgraphics(gcf, fullfile(scriptDir, 'MomentArmModel_vs_full6D.png'));
fprintf('\nSaved plot to MomentArmModel_vs_full6D.png\n');
