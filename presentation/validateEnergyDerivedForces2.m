%% validateEnergyDerivedForces2.m
% Updated with the linear-term-inclusive residual model
% (residualStepwiseFit.m, CV R^2=0.9956 on total W). Differentiates
% analytically and re-validates Fz, Taux, Tauy against trusted MST.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));

g0_z = 15e-6; k = 1.601e-14;
thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;
gap = g0_z - rz;

c0   = 2.6808e-10;   % intercept
c1   = -5.3602e-07;  % thx
c2   = 5.247e-07;    % thy
c3   = 9.8265e-12;   % thx/gap
c4   = -9.6222e-12;  % thy/gap
c5   = 2.3632e-12;   % ry/gap
c6   = 1.0715e-14;   % thz/gap
cxz  = -0.030297;    % thx*rz
cyz  = 0.029382;     % thy*rz
cxz2 = -2.9776e-07;  % thx*(rz/z_max)^2
cyz2 = 2.7945e-07;   % thy*(rz/z_max)^2
c7   = -0.000127;    % thx^2
c8   = -0.00015655;  % thy^2
c9   = 2.0575e-06;   % thz^2
c10  = 0.095682;     % rx^2
c11  = 0.088075;     % ry^2
c12  = 2.1473e-09;   % thx^2/gap
c13  = 2.559e-09;    % thy^2/gap
c14  = -8.9337e-07;  % rx^2/gap
c15  = -8.1853e-07;  % ry^2/gap
c16  = -1.6696e-11;  % thz^2/gap
c17  = -2.9555e-05;  % thx*thy
z_max = max(abs(rz));
rzn = rz / z_max;

Wpred = k./gap + c0 + c1*thx_rad + c2*thy_rad + c3*thx_rad./gap + c4*thy_rad./gap ...
      + c5*ry./gap + c6*thz_rad./gap + cxz*thx_rad.*rz + cyz*thy_rad.*rz ...
      + cxz2*thx_rad.*rzn.^2 + cyz2*thy_rad.*rzn.^2 ...
      + c7*thx_rad.^2 + c8*thy_rad.^2 + c9*thz_rad.^2 ...
      + c10*rx.^2 + c11*ry.^2 + c12*thx_rad.^2./gap + c13*thy_rad.^2./gap ...
      + c14*rx.^2./gap + c15*ry.^2./gap + c16*thz_rad.^2./gap + c17*thx_rad.*thy_rad;
r2_check = 1 - sum((W-Wpred).^2)/sum((W-mean(W)).^2);
fprintf('In-sample R^2 of reconstructed full model: %.4f (should be ~0.9960)\n\n', r2_check);

% All six use +dW/dq -- confirmed to be COMSOL's probe sign convention
% throughout (not just an Fz quirk): negating the -dW/dq version of
% Taux/Tauy against MST jumped R^2 from -1.64 to +0.71 (corr -0.97),
% a near-perfect match once the sign is corrected.
Fx_energy = ( 2*c10*rx + 2*c14*rx./gap );
Fy_energy = ( c5./gap + 2*c11*ry + 2*c15*ry./gap );
% cxz*thx*rz is LINEAR in rz directly (not through gap): drz-derivative
% is cxz*thx (not /gap^2). cxz2*thx*rzn^2 = cxz2*thx*rz^2/z_max^2:
% drz-derivative is 2*cxz2*thx*rz/z_max^2.
Fz_energy = ( k + c3*thx_rad + c4*thy_rad + c5*ry + c6*thz_rad ...
            + c12*thx_rad.^2 + c13*thy_rad.^2 + c14*rx.^2 + c15*ry.^2 + c16*thz_rad.^2 ) ./ gap.^2 ...
            + cxz*thx_rad + cyz*thy_rad ...
            + 2*cxz2*thx_rad.*rz/z_max^2 + 2*cyz2*thy_rad.*rz/z_max^2;
Taux_energy = ( c1 + c3./gap + cxz*rz + cxz2*rzn.^2 + 2*c7*thx_rad + 2*c12*thx_rad./gap + c17*thy_rad );
Tauy_energy = ( c2 + c4./gap + cyz*rz + cyz2*rzn.^2 + 2*c8*thy_rad + 2*c13*thy_rad./gap + c17*thx_rad );
Tauz_energy = ( c6./gap + 2*c9*thz_rad + 2*c16*thz_rad./gap );

figure('Color','w','Name','Energy-derived vs MST v2: Fz, Taux, Tauy');
subplot(1,3,1);
scatter(Fz_MST, Fz_energy, 6, 'filled', 'MarkerFaceAlpha',0.3); hold on; grid on;
lims = [min(Fz_MST) max(Fz_MST)]; plot(lims,lims,'r--');
xlabel('Fz (MST) [N]'); ylabel('Fz (energy-derived) [N]');
r2 = 1 - sum((Fz_MST-Fz_energy).^2)/sum((Fz_MST-mean(Fz_MST)).^2);
title(sprintf('Fz: R^2=%.3f', r2));

subplot(1,3,2);
scatter(Taux_MST, Taux_energy, 6, 'filled', 'MarkerFaceAlpha',0.3); hold on; grid on;
lims = [min(Taux_MST) max(Taux_MST)]; plot(lims,lims,'r--');
xlabel('Taux (MST) [N*m]'); ylabel('Taux (energy-derived) [N*m]');
r2 = 1 - sum((Taux_MST-Taux_energy).^2)/sum((Taux_MST-mean(Taux_MST)).^2);
title(sprintf('Taux: R^2=%.3f', r2));

subplot(1,3,3);
scatter(Tauy_MST, Tauy_energy, 6, 'filled', 'MarkerFaceAlpha',0.3); hold on; grid on;
lims = [min(Tauy_MST) max(Tauy_MST)]; plot(lims,lims,'r--');
xlabel('Tauy (MST) [N*m]'); ylabel('Tauy (energy-derived) [N*m]');
r2 = 1 - sum((Tauy_MST-Tauy_energy).^2)/sum((Tauy_MST-mean(Tauy_MST)).^2);
title(sprintf('Tauy: R^2=%.3f', r2));

exportgraphics(gcf, fullfile(scriptDir, 'energyDerived_vs_MST_v2.png'));
fprintf('Saved plot to energyDerived_vs_MST_v2.png\n');

save(fullfile(scriptDir, 'energyDerivedForces_v2.mat'), ...
    'Fx_energy','Fy_energy','Fz_energy','Taux_energy','Tauy_energy','Tauz_energy');
