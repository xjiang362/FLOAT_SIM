%% plotFzVsRz_3models.m
% Fz vs rz: three datasets, using the full 2000-pt 6D LHS dataset
% (points scattered across all 6 DOF, not a clean 1D sweep, so some
% scatter from the other DOF is expected for all three series):
%   1. COMSOL/MST simulation (ground truth).
%   2. Fitted fringing-capacitance model (trapElectrode1Model.m).
%   3. Simple parallel-plate model Fz = eps_r*eps0*A*V^2/(2*(d0-z)^2),
%      with d0=15um and eps_r=1 fixed (matching the fringing model's
%      nominal gap) and A fit by linear least squares (the model is
%      linear in A, so this is a closed-form fit, not an optimizer).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));  % rx,ry,rz,thx,thy,thz,Fz_MST,...

V = 30;
U = V*ones(size(rx));

%% Dataset 2: fringing-capacitance model
[~, ~, Fz_fringe] = trapElectrode1Model(rx, ry, rz, thx, thy, thz, U);

%% Dataset 3: simple parallel-plate model, A fit by linear least squares
d0 = 15e-6;
eps0 = 8.854e-12;
eps_r = 1;

predictor = eps_r*eps0*V^2 ./ (2*(d0-rz).^2);
A_fit = (predictor' * Fz_MST) / (predictor' * predictor);
Fz_pp = A_fit * predictor;

fprintf('Parallel-plate model fit: A = %.6e m^2 (equivalent square side = %.1f um)\n', ...
    A_fit, sqrt(abs(A_fit))*1e6);

%% R^2 for both models
r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
r2_fringe = r2fun(Fz_MST, Fz_fringe);
r2_pp = r2fun(Fz_MST, Fz_pp);
fprintf('Fringing-capacitance model: R^2 = %.4f\n', r2_fringe);
fprintf('Parallel-plate model:       R^2 = %.4f\n', r2_pp);

err_fringe = Fz_fringe - Fz_MST;
err_pp = Fz_pp - Fz_MST;
fprintf('Fringing error: mean=%.4e std=%.4e max|err|=%.4e\n', mean(err_fringe), std(err_fringe), max(abs(err_fringe)));
fprintf('Parallel-plate error: mean=%.4e std=%.4e max|err|=%.4e\n', mean(err_pp), std(err_pp), max(abs(err_pp)));

%% Plot
figure('Color','w','Name','Fz vs rz: simulation, fringing model, parallel-plate model', 'Position',[50 50 1500 800]);

subplot(2,2,[1 2]); hold on; grid on;
scatter(rz*1e6, Fz_MST, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Simulation');
scatter(rz*1e6, Fz_fringe, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', sprintf('Fringing model (R^2=%.4f)', r2_fringe));
scatter(rz*1e6, Fz_pp, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', sprintf('Parallel-plate model (R^2=%.4f)', r2_pp));
xlabel('rz [\mum]'); ylabel('Fz [N]');
title('Fz vs rz');
legend('Location','best');

subplot(2,2,3); hold on; grid on;
scatter(rz*1e6, err_fringe, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.3);
yline(0, 'k--');
xlabel('rz [\mum]'); ylabel('Fz_{model} - Fz_{sim} [N]');
title('Fringing model error vs rz');

subplot(2,2,4); hold on; grid on;
scatter(rz*1e6, err_pp, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
yline(0, 'k--');
xlabel('rz [\mum]'); ylabel('Fz_{model} - Fz_{sim} [N]');
title('Parallel-plate model error vs rz');

exportgraphics(gcf, fullfile(scriptDir, 'Fz_vs_rz_3models.png'));
fprintf('\nSaved plot to Fz_vs_rz_3models.png\n');
