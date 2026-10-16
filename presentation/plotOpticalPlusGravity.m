%% plotOpticalPlusGravity.m
% Illustrative shape only (not fit to any data): Gaussian optical dipole
% trap potential plus a linear gravitational potential, shown both
% separately and summed, vs position.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

U0 = 1;      % normalized optical trap depth
w0 = 1;      % normalized beam waist
mg = 0.3;    % normalized gravitational "weight" (slope of m*g*x)

x = linspace(-3*w0, 3*w0, 500)';
U_optical = -U0 * exp(-2*x.^2/w0^2);
U_gravity = mg * x;
U_total = U_optical + U_gravity;

figure('Color','w','Name','Optical + gravity potential shape');
hold on; grid on;
plot(x, U_optical, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Optical (Gaussian well)');
plot(x, U_gravity, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Gravity (linear)');
plot(x, U_total, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Total');
xlabel('Position x');
ylabel('Potential energy U(x)');
title('Optical + gravitational potential vs position');
legend('Location','best');
set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(scriptDir, 'OpticalPlusGravityShape.png'));
fprintf('Saved plot to OpticalPlusGravityShape.png\n');
