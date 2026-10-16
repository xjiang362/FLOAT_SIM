%% plotElectrostaticPlusGravity.m
% Illustrative shape only (not fit to any data): parallel-plate
% electrostatic energy plus a linear gravitational potential, shown
% both separately and summed, vs position.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

d0 = 1;      % normalized nominal gap (plate contact position)
k  = 1;      % normalized electrostatic energy scale
mg = 2;      % normalized gravitational "weight" (slope of m*g*x) -- increased

x = linspace(-0.9*d0, 0.9*d0, 500)';
W_electrostatic = -k ./ (d0 - x);   % sign flipped: dives to -inf as x->d0 instead of +inf
U_gravity = mg * x;
U_total = W_electrostatic + U_gravity;

figure('Color','w','Name','Electrostatic + gravity potential shape');
hold on; grid on;
plot(x, W_electrostatic, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Electrostatic (parallel-plate)');
plot(x, U_gravity, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Gravity (linear)');
plot(x, U_total, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Total');
xlabel('Position x');
ylabel('Potential energy U(x)');
title('Electrostatic + gravitational potential vs position');
legend('Location','best');
set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(scriptDir, 'ElectrostaticPlusGravityShape.png'));
fprintf('Saved plot to ElectrostaticPlusGravityShape.png\n');
