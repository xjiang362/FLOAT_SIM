%% compareRealVsComplexSmallLargeP.m
% Compares four first-order/complex-pole systems' impulse response (real
% part), pairing real vs complex pole at both a small and a large decay
% rate -- isolates the effect of ADDING oscillation (real -> complex,
% same decay rate) at each decay-rate regime:
%   1/(s+p)        -- small decay rate, real pole (no oscillation)
%   1/(s+p+j*W)     -- small decay rate, complex pole (adds oscillation)
%   1/(s+P)        -- large decay rate, real pole (no oscillation)
%   1/(s+P+j*W)     -- large decay rate, complex pole (adds oscillation)
% Re{y(t)} = exp(-a*(t-t0))*cos(b*(t-t0)) for t>=t0, 0 for t<t0 -- an
% impulse disturbance applied at t=t0 (=1) instead of at t=0, so the
% response starts flat at equilibrium and only kicks off once the
% impulse hits, rather than starting from a displaced initial condition.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

p = 1;    % small decay rate
P = 6;    % large decay rate
W = 18;   % angular frequency [rad/s], shared by both complex cases
t0 = 1;   % time of the impulse

t = linspace(0, 5, 5000)';
tau = t - t0;
u = tau >= 0;   % impulse response only kicks in for t>=t0

g_p      = u .* exp(-p*tau);
g_p_W    = u .* exp(-p*tau) .* cos(W*tau);
g_P      = u .* exp(-P*tau);
g_P_W    = u .* exp(-P*tau) .* cos(W*tau);

figure('Color','w','Name','Real vs complex pole at small vs large decay rate');
hold on; grid on;
plot(t, g_p,   'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('1/(s+p): p=%g', p));
plot(t, g_p_W, 'b-',  'LineWidth', 1.75, 'DisplayName', sprintf('1/(s+p+j\\Omega): p=%g, \\Omega=%g', p, W));
plot(t, g_P,   'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('1/(s+P): P=%g', P));
plot(t, g_P_W, 'r-',  'LineWidth', 1.75, 'DisplayName', sprintf('1/(s+P+j\\Omega): P=%g, \\Omega=%g', P, W));
xline(t0, 'k:', 'Impulse', 'HandleVisibility', 'off');
yline(0, 'k:', 'HandleVisibility', 'off');
xlabel('Time [s]');
ylabel('Re\{y(t)\}');
title(sprintf('Impulse response (impulse at t=%g s): real vs complex pole, small vs large decay rate', t0));
legend('Location','best');
set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(scriptDir, 'RealVsComplex_smallLargeP.png'));
fprintf('Saved plot to RealVsComplex_smallLargeP.png\n');
