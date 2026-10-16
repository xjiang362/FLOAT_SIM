%% compareSmallVsLargeP_impulse.m
% Compares the impulse response of G(s) = 1/(s+p)^2 (a repeated real
% pole at s=-p) for a small p (slow) vs a large p (fast), under an
% impulse disturbance. Closed-form impulse response:
%   g(t) = t*exp(-p*t)   (inverse Laplace of 1/(s+p)^2)

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

p_small = 1;    % slow system
p_large = 8;    % fast system

t = linspace(0, 6, 2000)';
g_small = t .* exp(-p_small*t);
g_large = t .* exp(-p_large*t);

figure('Color','w','Name','Impulse response: small p vs large p');
hold on; grid on;
plot(t, g_small, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('p = %g (slow)', p_small));
plot(t, g_large, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('p = %g (fast)', p_large));
yline(0, 'k--', 'HandleVisibility', 'off');
xlabel('Time [s]');
ylabel('y(t)');
title('Impulse response of G(s) = 1/(s+p)^2: small p vs large p');
legend('Location','best');
set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(scriptDir, 'SmallVsLargeP_impulse.png'));
fprintf('Saved plot to SmallVsLargeP_impulse.png\n');

fprintf('Peak time (t=1/p): small p -> %.4g s, large p -> %.4g s\n', 1/p_small, 1/p_large);
fprintf('Peak value (1/(p*e)): small p -> %.4g, large p -> %.4g\n', 1/(p_small*exp(1)), 1/(p_large*exp(1)));
