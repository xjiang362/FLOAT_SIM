%% plotOpticalPlusElectrostatic.m
% Illustrative shape only (not fit to any data): Gaussian optical dipole
% trap potential plus the parallel-plate electrostatic energy, shown
% both separately and summed, vs position.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

U0 = 1;      % normalized optical trap depth
w0 = 1;      % normalized beam waist
d0 = 3;      % normalized electrostatic gap (plate contact position)
k  = 1;      % normalized electrostatic energy scale

x = linspace(-0.9*d0, 0.9*d0, 500)';
U_optical = -U0 * exp(-2*x.^2/w0^2);
W_electrostatic = k ./ (d0 - x);
U_total = U_optical + W_electrostatic;

figure('Color','w','Name','Optical + electrostatic potential shape');
hold on; grid on;
plot(x, U_optical, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Optical (Gaussian well)');
plot(x, W_electrostatic, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Electrostatic (parallel-plate)');
plot(x, U_total, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Total');
xlabel('Position x');
ylabel('Potential energy U(x)');
title('Optical + electrostatic potential vs position');
legend('Location','best');
set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(scriptDir, 'OpticalPlusElectrostaticShape.png'));
fprintf('Saved plot to OpticalPlusElectrostaticShape.png\n');
