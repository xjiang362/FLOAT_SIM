%% plotParallelPlateEnergyShape.m
% Illustrative shape only (not fit to any data): stored electrostatic
% energy between parallel plates vs position, W(x) = k/(d0-x), showing
% the characteristic divergence as the moving plate approaches the
% fixed plate (x -> d0).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

d0 = 1;     % normalized nominal gap
k = 1;      % normalized scale (eps0*A*V^2/2, arbitrary units here)

x = linspace(-0.9*d0, 0.9*d0, 500)';
W = k ./ (d0 - x);

figure('Color','w','Name','Parallel-plate energy shape');
plot(x, W, 'b-', 'LineWidth', 2.5);
xlabel('Position x');
ylabel('Energy W(x)');
title('Parallel-plate capacitor energy vs position');
grid on;
xlim([-0.9*d0, 0.9*d0]);
set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(scriptDir, 'ParallelPlateEnergyShape.png'));
fprintf('Saved plot to ParallelPlateEnergyShape.png\n');
