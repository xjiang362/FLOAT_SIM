%% plotCm1p2p_rawShape.m
% Quick look at Cm(1,p) and Cm(2,p) (resonator coupling to terminals 1
% and 2) vs each DOF, to inform the fringing-model parametrization.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'capSweep300_withGND_allMatrices.mat'));

Cm1p = squeeze(Cm_all(1,21,:));
Cm2p = squeeze(Cm_all(2,21,:));

dofs = {rz, thx, thy, rx, ry, thz};
dofNames = {'rz','thx','thy','rx','ry','thz'};

figure('Color','w','Name','Cm(1,p) and Cm(2,p) vs each DOF', 'Position',[50 50 1600 700]);
for i = 1:6
    subplot(2,6,i); hold on; grid on;
    scatter(dofs{i}, Cm1p, 10, 'filled', 'MarkerFaceAlpha', 0.4);
    xlabel(dofNames{i}); ylabel('Cm(1,p) [F]');
    title(['Cm1p vs ' dofNames{i}]);

    subplot(2,6,i+6); hold on; grid on;
    scatter(dofs{i}, Cm2p, 10, 'filled', 'MarkerFaceAlpha', 0.4);
    xlabel(dofNames{i}); ylabel('Cm(2,p) [F]');
    title(['Cm2p vs ' dofNames{i}]);
end
exportgraphics(gcf, fullfile(scriptDir, 'Cm1p2p_rawShape.png'));
fprintf('Saved plot to Cm1p2p_rawShape.png\n');

fprintf('\nCm1p range: [%.4e, %.4e]\n', min(Cm1p), max(Cm1p));
fprintf('Cm2p range: [%.4e, %.4e]\n', min(Cm2p), max(Cm2p));
