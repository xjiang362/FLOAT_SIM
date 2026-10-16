%% compareSZeroSmallLargeP.m
% Compares the impulse response of two systems with a zero at the
% origin (relative degree 0 -- direct feedthrough):
%   s/(s+p+j*W)   -- small decay rate p
%   s/(s+P+j*W)   -- large decay rate P
% Using s/(s+a) = 1 - a/(s+a): the impulse response is
%   g(t) = delta(t-t0) - a*exp(-a*(t-t0))*u(t-t0)
% The delta term has weight exactly 1 in BOTH cases (an unavoidable
% instantaneous full-strength pass-through of the impulsive input,
% independent of p/P/W -- the "no high-frequency roll-off" signature of
% a zero at the origin) and is drawn as a vertical arrow, since a true
% Dirac delta can't be plotted as a finite curve. The REMAINING smooth
% tail, Re{-a*exp(-a*tau)} = -exp(-p*tau)*(p*cos(W*tau)+W*sin(W*tau)),
% is what actually depends on the pole location and is plotted as a
% continuous curve.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

p = 1;    % small decay rate
P = 6;    % large decay rate
W = 18;   % angular frequency [rad/s], shared by both cases
t0 = 1;   % time of the impulse

t = linspace(0, 5, 5000)';
tau = t - t0;
u = tau >= 0;

tail_p = u .* (-exp(-p*tau) .* (p*cos(W*tau) + W*sin(W*tau)));
tail_P = u .* (-exp(-P*tau) .* (P*cos(W*tau) + W*sin(W*tau)));

figure('Color','w','Name','s/(s+p+jW) vs s/(s+P+jW) impulse response');
hold on; grid on;
plot(t, tail_p, 'b-', 'LineWidth', 1.75, 'DisplayName', sprintf('s/(s+p+j\\Omega) tail: p=%g, \\Omega=%g', p, W));
plot(t, tail_P, 'r-', 'LineWidth', 1.75, 'DisplayName', sprintf('s/(s+P+j\\Omega) tail: P=%g, \\Omega=%g', P, W));
% Delta arrows at t0, weight 1 in both cases -- drawn at a fixed
% SYMBOLIC height (a true Dirac delta has infinite height/zero width,
% so this is not to scale against the tail, just a marker for "a unit
% impulse occurs here")
arrowHeight = 6;
quiver(t0, 0, 0, arrowHeight, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
quiver(t0+0.03, 0, 0, arrowHeight, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
text(t0+0.08, arrowHeight, '\delta(t-t_0), weight 1 (both -- symbolic height)', 'FontSize', 10);
yline(0, 'k:', 'HandleVisibility', 'off');
xlabel('Time [s]');
ylabel('y(t)');
title(sprintf('Impulse response (impulse at t=%g s): s/(s+p+j\\Omega) vs s/(s+P+j\\Omega)', t0));
legend('Location','best');
ylim([-20, 18]);
set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(scriptDir, 'SZero_smallLargeP.png'));
fprintf('Saved plot to SZero_smallLargeP.png\n');
