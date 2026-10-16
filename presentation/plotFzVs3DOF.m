%% plotFzVs3DOF.m
% Fz vs rz, thetax, thetay in 3 subplots, each showing all three
% datasets (simulation, fringing model, parallel-plate model), using the
% full 2000-pt 6D LHS dataset. Since this is a 6D LHS sample (not a
% clean 1D sweep along each axis), scatter from the other DOF is
% expected in every panel.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));  % rx,ry,rz,thx,thy,thz,Fz_MST,...

V = 30;
U = V*ones(size(rx));

%% Fringing-capacitance model
[~, ~, Fz_fringe] = trapElectrode1Model(rx, ry, rz, thx, thy, thz, U);

%% Parallel-plate model, A fit by linear least squares
d0 = 15e-6;
eps0 = 8.854e-12;
eps_r = 1;

predictor = eps_r*eps0*V^2 ./ (2*(d0-rz).^2);
A_fit = (predictor' * Fz_MST) / (predictor' * predictor);
Fz_pp = A_fit * predictor;

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
r2_fringe = r2fun(Fz_MST, Fz_fringe);
r2_pp = r2fun(Fz_MST, Fz_pp);

%% Plot
dofs = {rz*1e6, thx, thy};
dofNames = {'rz [\mum]', '\theta_x [deg]', '\theta_y [deg]'};

figure('Color','w','Name','Fz vs rz, thetax, thetay', 'Position',[50 50 1600 500]);
for i = 1:3
    subplot(1,3,i); hold on; grid on;
    scatter(dofs{i}, Fz_MST, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Simulation');
    scatter(dofs{i}, Fz_fringe, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', sprintf('Fringing model (R^2=%.4f)', r2_fringe));
    scatter(dofs{i}, Fz_pp, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', sprintf('Parallel-plate model (R^2=%.4f)', r2_pp));
    xlabel(dofNames{i}); ylabel('Fz [N]');
    title(['Fz vs ' dofNames{i}]);
    if i == 1
        legend('Location','best');
    end
end

exportgraphics(gcf, fullfile(scriptDir, 'Fz_vs_3DOF.png'));
fprintf('Saved plot to Fz_vs_3DOF.png\n');
