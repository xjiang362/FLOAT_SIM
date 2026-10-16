%% plotGaussianOpticalPotentialShape.m
% Illustrative shape only (not fit to any data): a Gaussian-beam optical
% dipole trap potential well vs position, U(x) = -U0*exp(-2x^2/w0^2),
% for contrast against the parallel-plate electrostatic energy shape.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

U0 = 1;     % normalized trap depth
w0 = 1;     % normalized beam waist

x = linspace(-3*w0, 3*w0, 500)';
U = -U0 * exp(-2*x.^2/w0^2);

figure('Color','w','Name','Gaussian optical potential well shape');
plot(x, U, 'r-', 'LineWidth', 2.5);
xlabel('Position x');
ylabel('Potential U(x)');
title('Gaussian optical dipole trap potential vs position');
grid on;
xlim([-3*w0, 3*w0]);
set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(scriptDir, 'GaussianOpticalPotentialShape.png'));
fprintf('Saved plot to GaussianOpticalPotentialShape.png\n');
