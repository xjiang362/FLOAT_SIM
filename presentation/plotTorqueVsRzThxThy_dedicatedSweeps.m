%% plotTorqueVsRzThxThy_dedicatedSweeps.m
% Same idea as plotFzVsRzThxThy_dedicatedSweeps.m, but for Taux and
% Tauy: 2 rows (Taux, Tauy) x 3 columns (rz, thx, thy), using the
% dedicated single-axis sweeps (comsol_sweep_combined300_results.txt):
%   rows 101-150: rz swept, rx=ry=thx=thy=thz=0
%   rows 151-200: thx swept, rx=ry=rz=thy=thz=0
%   rows 201-250: thy swept, rx=ry=rz=thx=thz=0
%
% Compares simulation vs the fitted fringing-capacitance model
% (trapElectrode1Model). No parallel-plate baseline here -- the
% Fz=eps_r*eps0*A*V^2/(2*(d0-z)^2) formula the user gave has no torque
% analog, so this is model-vs-simulation only.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

data = parseComsolResults(fullfile(scriptDir, 'comsol_sweep_combined300_results.txt'));
V = 30;

blocks = {101:150, 151:200, 201:250};
dofNames = {'rz', 'thx', 'thy'};
dofUnits = {'\mum', 'deg', 'deg'};
dofScale = {1e6, 1, 1};

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

tauOutputs = {'Taux', 'Tauy'};

figure('Color','w','Name','Taux, Tauy vs rz, thx, thy: dedicated single-axis sweeps', 'Position',[50 50 1600 900]);
for r = 1:2
    outName = tauOutputs{r};
    for i = 1:3
        b = data(blocks{i}, :);
        n = height(b);
        U = V*ones(n,1);

        [~, ~, ~, Taux_fringe, Tauy_fringe] = trapElectrode1Model(b.rx, b.ry, b.rz, b.thx, b.thy, b.thz, U);
        if r == 1
            tau_fringe = Taux_fringe;
            tau_sim = b.Taux;
        else
            tau_fringe = Tauy_fringe;
            tau_sim = b.Tauy;
        end

        r2_fringe = r2fun(tau_sim, tau_fringe);
        fprintf('%4s vs %3s-only sweep: Fringing R^2=%.4f\n', outName, dofNames{i}, r2_fringe);

        x = b.(dofNames{i}) * dofScale{i};

        subplot(2,3,(r-1)*3+i); hold on; grid on;
        scatter(x, tau_sim, 30, 'k', 'filled', 'DisplayName', 'Simulation');
        scatter(x, tau_fringe, 30, 'r', 'filled', 'DisplayName', sprintf('Fringing model (R^2=%.4f)', r2_fringe));
        xlabel([dofNames{i} ' [' dofUnits{i} ']']); ylabel([outName ' [N*m]']);
        title([outName ' vs ' dofNames{i} ' (dedicated sweep)']);
        legend('Location','best');
    end
end

exportgraphics(gcf, fullfile(scriptDir, 'Torque_vs_rzThxThy_dedicatedSweeps.png'));
fprintf('Saved plot to Torque_vs_rzThxThy_dedicatedSweeps.png\n');
