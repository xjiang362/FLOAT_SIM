%% validateTwoMechanismModel.m
% W = P1(thx,thy,thz)/gap + P2(rx,ry), analytic derivatives:
%   Fx = dP2/drx,  Fy = dP2/dry,  Fz = P1/gap^2
%   Taux = (1/gap)*dP1/dthx,  Tauy = (1/gap)*dP1/dthy,  Tauz = (1/gap)*dP1/dthz
% (all +dW/dq, confirmed sign convention)

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));
load(fullfile(scriptDir, 'twoMechanismModel.mat'));

g0_z = 15e-6;
thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;
gap = g0_z - rz;

% coefs order: invgap,thx_g,thy_g,thz_g,thx2_g,thy2_g,thz2_g,thxthy_g,thxthz_g,thythz_g,const,rx,ry,rx2,ry2,rxry
k0=coefs(1); a1=coefs(2); a2=coefs(3); a3=coefs(4); a4=coefs(5); a5=coefs(6); a6=coefs(7);
a7=coefs(8); a8=coefs(9); a9=coefs(10);
b0=coefs(11); b1=coefs(12); b2=coefs(13); b3=coefs(14); b4=coefs(15); b5=coefs(16);

P1 = k0 + a1*thx_rad + a2*thy_rad + a3*thz_rad + a4*thx_rad.^2 + a5*thy_rad.^2 + a6*thz_rad.^2 ...
   + a7*thx_rad.*thy_rad + a8*thx_rad.*thz_rad + a9*thy_rad.*thz_rad;
P2 = b0 + b1*rx + b2*ry + b3*rx.^2 + b4*ry.^2 + b5*rx.*ry;
Wpred = P1./gap + P2;
r2 = 1 - sum((W-Wpred).^2)/sum((W-mean(W)).^2);
fprintf('Reconstructed W R^2: %.4f\n\n', r2);

Fx_energy = b1 + 2*b3*rx + b5*ry;
Fy_energy = b2 + 2*b4*ry + b5*rx;
Fz_energy = P1 ./ gap.^2;
Taux_energy = (a1 + 2*a4*thx_rad + a7*thy_rad + a8*thz_rad) ./ gap;
Tauy_energy = (a2 + 2*a5*thy_rad + a7*thx_rad + a9*thz_rad) ./ gap;
Tauz_energy = (a3 + 2*a6*thz_rad + a8*thx_rad + a9*thy_rad) ./ gap;

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
fprintf('Fz:   R^2 vs MST = %.4f\n', r2fun(Fz_MST, Fz_energy));
fprintf('Taux: R^2 vs MST = %.4f\n', r2fun(Taux_MST, Taux_energy));
fprintf('Tauy: R^2 vs MST = %.4f\n', r2fun(Tauy_MST, Tauy_energy));

figure('Color','w','Name','Two-mechanism model vs MST');
subplot(1,3,1); scatter(Fz_MST, Fz_energy, 6,'filled','MarkerFaceAlpha',0.3); hold on; grid on;
lims=[min(Fz_MST) max(Fz_MST)]; plot(lims,lims,'r--');
xlabel('Fz (MST)'); ylabel('Fz (energy)'); title(sprintf('Fz R^2=%.3f', r2fun(Fz_MST,Fz_energy)));
subplot(1,3,2); scatter(Taux_MST, Taux_energy, 6,'filled','MarkerFaceAlpha',0.3); hold on; grid on;
lims=[min(Taux_MST) max(Taux_MST)]; plot(lims,lims,'r--');
xlabel('Taux (MST)'); ylabel('Taux (energy)'); title(sprintf('Taux R^2=%.3f', r2fun(Taux_MST,Taux_energy)));
subplot(1,3,3); scatter(Tauy_MST, Tauy_energy, 6,'filled','MarkerFaceAlpha',0.3); hold on; grid on;
lims=[min(Tauy_MST) max(Tauy_MST)]; plot(lims,lims,'r--');
xlabel('Tauy (MST)'); ylabel('Tauy (energy)'); title(sprintf('Tauy R^2=%.3f', r2fun(Tauy_MST,Tauy_energy)));
exportgraphics(gcf, fullfile(scriptDir, 'twoMechanismModel_vs_MST.png'));
fprintf('\nSaved plot to twoMechanismModel_vs_MST.png\n');
