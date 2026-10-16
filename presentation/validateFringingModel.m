%% validateFringingModel.m
% Differentiates the fitted fringing-capacitance W model via central
% finite-difference on the closed form itself (safe -- smooth analytic
% function, not raw simulation data, so no mesh-rebuild noise), then
% validates Fz, Taux, Tauy against trusted MST.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));
load(fullfile(scriptDir, 'fringingCapModel.mat'));

thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;
X = [rx, ry, rz, thx_rad, thy_rad, thz_rad];
n = numel(rz);

Wpred = modelFun(p, X);
r2 = 1 - sum((W-Wpred).^2)/sum((W-mean(W)).^2);
fprintf('Reconstructed W R^2: %.4f\n\n', r2);

% Central differences on the fitted closed form, per-column step sized
% to each variable's own range
steps = [1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9];  % base steps (m or rad)
scales = [max(abs(rx)), max(abs(ry)), max(abs(rz)), max(abs(thx_rad)), max(abs(thy_rad)), max(abs(thz_rad))];
h = 1e-6 * scales;  % relative step: 1e-6 of each variable's range

dWdq = zeros(n, 6);
for j = 1:6
    Xp = X; Xp(:,j) = Xp(:,j) + h(j);
    Xm = X; Xm(:,j) = Xm(:,j) - h(j);
    dWdq(:,j) = (modelFun(p, Xp) - modelFun(p, Xm)) / (2*h(j));
end

% +dW/dq convention (confirmed sign throughout)
Fx_energy = dWdq(:,1); Fy_energy = dWdq(:,2); Fz_energy = dWdq(:,3);
Taux_energy = dWdq(:,4); Tauy_energy = dWdq(:,5); Tauz_energy = dWdq(:,6);

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
fprintf('Fz:   R^2 vs MST = %.4f\n', r2fun(Fz_MST, Fz_energy));
fprintf('Taux: R^2 vs MST = %.4f\n', r2fun(Taux_MST, Taux_energy));
fprintf('Tauy: R^2 vs MST = %.4f\n', r2fun(Tauy_MST, Tauy_energy));

figure('Color','w','Name','Fringing model vs MST');
subplot(1,3,1); scatter(Fz_MST, Fz_energy, 6,'filled','MarkerFaceAlpha',0.3); hold on; grid on;
lims=[min(Fz_MST) max(Fz_MST)]; plot(lims,lims,'r--');
xlabel('Fz (MST)'); ylabel('Fz (fringing model)'); title(sprintf('Fz R^2=%.3f', r2fun(Fz_MST,Fz_energy)));
subplot(1,3,2); scatter(Taux_MST, Taux_energy, 6,'filled','MarkerFaceAlpha',0.3); hold on; grid on;
lims=[min(Taux_MST) max(Taux_MST)]; plot(lims,lims,'r--');
xlabel('Taux (MST)'); ylabel('Taux (fringing model)'); title(sprintf('Taux R^2=%.3f', r2fun(Taux_MST,Taux_energy)));
subplot(1,3,3); scatter(Tauy_MST, Tauy_energy, 6,'filled','MarkerFaceAlpha',0.3); hold on; grid on;
lims=[min(Tauy_MST) max(Tauy_MST)]; plot(lims,lims,'r--');
xlabel('Tauy (MST)'); ylabel('Tauy (fringing model)'); title(sprintf('Tauy R^2=%.3f', r2fun(Tauy_MST,Tauy_energy)));
exportgraphics(gcf, fullfile(scriptDir, 'fringingModel_vs_MST.png'));
fprintf('\nSaved plot to fringingModel_vs_MST.png\n');

save(fullfile(scriptDir, 'fringingModelForces.mat'), 'Fx_energy','Fy_energy','Fz_energy','Taux_energy','Tauy_energy','Tauz_energy');
