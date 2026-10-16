%% plotFzVsRzThxThy_dedicatedSweeps.m
% Fz vs rz, thx, thy in 3 subplots, using the DEDICATED single-axis
% characterization sweeps (comsol_sweep_combined300_results.txt, from
% dedicatedSingleAxisSweeps.m) instead of the scattered full6D LHS data:
%   rows 101-150: rz swept, rx=ry=thx=thy=thz=0
%   rows 151-200: thx swept, rx=ry=rz=thy=thz=0
%   rows 201-250: thy swept, rx=ry=rz=thx=thz=0
% Each block isolates that one DOF cleanly (no scatter from the other
% 5), unlike the full6D LHS comparison.
%
% Compares: simulation, fringing-capacitance model (trapElectrode1Model),
% and the parallel-plate model (Fz=eps_r*eps0*A*V^2/(2*(d0-z)^2), A fit
% by linear least squares on the full6D LHS rz data -- so it is only a
% function of rz and is flat/constant across the thx and thy blocks,
% which sit at rz=0).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

data = parseComsolResults(fullfile(scriptDir, 'comsol_sweep_combined300_results.txt'));
V = 30;

blocks = {101:150, 151:200, 201:250};
dofNames = {'rz', 'thx', 'thy'};
dofUnits = {'\mum', 'deg', 'deg'};
dofScale = {1e6, 1, 1};   % rz reported in um, thx/thy already in deg

%% Parallel-plate model: A fit by linear least squares on the full6D LHS rz data
load(fullfile(scriptDir, 'full6D_energyData.mat'));   % rz, Fz_MST (full6D LHS)
d0 = 15e-6; eps0 = 8.854e-12; eps_r = 1;
predictor_full = eps_r*eps0*V^2 ./ (2*(d0-rz).^2);
A_fit = (predictor_full' * Fz_MST) / (predictor_full' * predictor_full);
ppModel = @(z) A_fit * eps_r*eps0*V^2 ./ (2*(d0-z).^2);

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

figure('Color','w','Name','Fz vs rz, thx, thy: dedicated single-axis sweeps', 'Position',[50 50 1600 500]);
for i = 1:3
    b = data(blocks{i}, :);
    n = height(b);
    U = V*ones(n,1);

    [~, ~, Fz_fringe] = trapElectrode1Model(b.rx, b.ry, b.rz, b.thx, b.thy, b.thz, U);
    Fz_pp = ppModel(b.rz);   % rz=0 in the thx/thy blocks -> constant there

    r2_fringe = r2fun(b.Fz, Fz_fringe);
    r2_pp = r2fun(b.Fz, Fz_pp);
    fprintf('%3s-only sweep: Fringing R^2=%.4f, Parallel-plate R^2=%.4f\n', dofNames{i}, r2_fringe, r2_pp);

    x = b.(dofNames{i}) * dofScale{i};
    [xs, order] = sort(x);

    subplot(1,3,i); hold on; grid on;
    scatter(x, b.Fz, 30, 'k', 'filled', 'DisplayName', 'Simulation');
    scatter(x, Fz_fringe, 30, 'r', 'filled', 'DisplayName', sprintf('Fringing model (R^2=%.4f)', r2_fringe));
    plot(xs, Fz_pp(order), 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Parallel-plate model (R^2=%.4f)', r2_pp));
    xlabel([dofNames{i} ' [' dofUnits{i} ']']); ylabel('Fz [N]');
    title(['Fz vs ' dofNames{i} ' (dedicated sweep)']);
    legend('Location','best');
end

exportgraphics(gcf, fullfile(scriptDir, 'Fz_vs_rzThxThy_dedicatedSweeps.png'));
fprintf('Saved plot to Fz_vs_rzThxThy_dedicatedSweeps.png\n');
