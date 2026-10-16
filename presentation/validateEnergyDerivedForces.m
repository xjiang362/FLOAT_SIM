%% validateEnergyDerivedForces.m
% Analytically differentiates the fitted global energy model
% W = k/(g0_z-rz) + c1*thx^2 + c2*thy^2 + c3*thz^2 + c4*rx^2 + c5*ry^2
%     + c6*thx^2/gap + c7*thy^2/gap + c8*rx^2/gap + c9*ry^2/gap + c10*thz^2/gap
%     + c11*thx*thy   (thx,thy,thz in RADIANS; gap = g0_z-rz)
% to get all six forces/torques, then cross-checks Fz, Taux, Tauy
% against the MST probes (already trusted/validated) as a consistency
% check on the whole energy-model approach before trusting the new
% Fx, Fy, Tauz predictions (MST was confirmed unreliable for those).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));

g0_z = 15e-6; k = 1.601e-14;
thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;
gap = g0_z - rz;

% Coefficients from the stepwise fit (residualStepwiseFit.m output)
c0  = 2.7093e-10;   % intercept (constant offset)
c1  = -0.00013943;  % thx^2
c2  = -0.0001005;   % thy^2
c3  = 1.8594e-06;   % thz^2
c4  = 0.10304;      % rx^2
c5  = 0.080201;     % ry^2
c6  = 2.3741e-09;   % thx^2/gap
c7  = 1.7672e-09;   % thy^2/gap
c8  = -9.1064e-07;  % rx^2/gap
c9  = -8.948e-07;   % ry^2/gap
c10 = -1.6413e-11;  % thz^2/gap
c11 = -2.8979e-05;  % thx*thy

% Predicted W (sanity check against combined CV R^2=0.9815)
Wpred = k./gap + c0 + c1*thx_rad.^2 + c2*thy_rad.^2 + c3*thz_rad.^2 + c4*rx.^2 + c5*ry.^2 ...
      + c6*thx_rad.^2./gap + c7*thy_rad.^2./gap + c8*rx.^2./gap + c9*ry.^2./gap + c10*thz_rad.^2./gap ...
      + c11*thx_rad.*thy_rad;
r2_check = 1 - sum((W-Wpred).^2)/sum((W-mean(W)).^2);
fprintf('In-sample R^2 of reconstructed full model: %.4f (should be close to 0.9815 CV)\n\n', r2_check);

% Analytic derivatives (forces/torques = -dW/dq)
Fx_energy = -( 2*c4*rx + 2*c8*rx./gap );
Fy_energy = -( 2*c5*ry + 2*c9*ry./gap );
Fz_energy = ( k + c6*thx_rad.^2 + c7*thy_rad.^2 + c8*rx.^2 + c9*ry.^2 + c10*thz_rad.^2 ) ./ gap.^2;
Taux_energy = -( 2*c1*thx_rad + 2*c6*thx_rad./gap + c11*thy_rad );
Tauy_energy = -( 2*c2*thy_rad + 2*c7*thy_rad./gap + c11*thx_rad );
Tauz_energy = -( 2*c3*thz_rad + 2*c10*thz_rad./gap );

%% Cross-check Fz, Taux, Tauy against trusted MST probes
figure('Color','w','Name','Energy-derived vs MST: Fz, Taux, Tauy (trusted axes)');
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

exportgraphics(gcf, fullfile(scriptDir, 'energyDerived_vs_MST_trustedAxes.png'));
fprintf('Saved plot to energyDerived_vs_MST_trustedAxes.png\n');

%% Show Fx, Fy, Tauz (previously untrustworthy via MST) predicted by energy model
figure('Color','w','Name','Energy-derived Fx, Fy, Tauz (MST was unreliable here)');
subplot(1,3,1);
scatter(Fx_MST, Fx_energy, 6, 'filled', 'MarkerFaceAlpha',0.3); hold on; grid on;
xlabel('Fx (MST) [N]'); ylabel('Fx (energy-derived) [N]');
title('Fx: MST vs energy-derived');

subplot(1,3,2);
scatter(Fy_MST, Fy_energy, 6, 'filled', 'MarkerFaceAlpha',0.3); hold on; grid on;
xlabel('Fy (MST) [N]'); ylabel('Fy (energy-derived) [N]');
title('Fy: MST vs energy-derived');

subplot(1,3,3);
scatter(Tauz_MST, Tauz_energy, 6, 'filled', 'MarkerFaceAlpha',0.3); hold on; grid on;
xlabel('Tauz (MST) [N*m]'); ylabel('Tauz (energy-derived) [N*m]');
title('Tauz: MST vs energy-derived');

exportgraphics(gcf, fullfile(scriptDir, 'energyDerived_FxFyTauz.png'));
fprintf('Saved plot to energyDerived_FxFyTauz.png\n');

fprintf('\n--- Magnitude sanity check ---\n');
fprintf('%8s  %14s  %14s\n', 'var', 'energy-derived range', 'MST range');
fprintf('%8s  [%9.3g,%9.3g]  [%9.3g,%9.3g]\n', 'Fx', min(Fx_energy), max(Fx_energy), min(Fx_MST), max(Fx_MST));
fprintf('%8s  [%9.3g,%9.3g]  [%9.3g,%9.3g]\n', 'Fy', min(Fy_energy), max(Fy_energy), min(Fy_MST), max(Fy_MST));
fprintf('%8s  [%9.3g,%9.3g]  [%9.3g,%9.3g]\n', 'Tauz', min(Tauz_energy), max(Tauz_energy), min(Tauz_MST), max(Tauz_MST));
