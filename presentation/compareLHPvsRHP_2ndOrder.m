%% compareLHPvsRHP_2ndOrder.m
% Compares the time-domain (initial-condition) response of two 2nd
% order systems with the same natural frequency and damping magnitude,
% but opposite pole locations:
%   LHP (stable):   G(s) = wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
%   RHP (unstable): G(s) = wn^2 / (s^2 - 2*zeta*wn*s + wn^2)
% Both start at y(0) = -1, y'(0) = 0. The LHP case settles to 0; the RHP
% case grows away from 0 instead.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

wn = 2*pi*1;     % natural frequency [rad/s] (1 Hz)
zeta = 0.15;     % damping ratio magnitude
x0 = [-1; 0];    % [y(0); y'(0)]

% State-space form: x1=y, x2=y'
odeLHP = @(t,x) [x(2); -wn^2*x(1) - 2*zeta*wn*x(2)];
odeRHP = @(t,x) [x(2); -wn^2*x(1) + 2*zeta*wn*x(2)];

tspan = linspace(0, 4, 2000);
[tL, xL] = ode45(odeLHP, tspan, x0);
[tR, xR] = ode45(odeRHP, tspan, x0);

figure('Color','w','Name','LHP vs RHP 2nd-order response');
hold on; grid on;
plot(tL, xL(:,1), 'b-', 'LineWidth', 2, 'DisplayName', 'LHP poles (stable) -- settles to 0');
plot(tR, xR(:,1), 'r-', 'LineWidth', 2, 'DisplayName', 'RHP poles (unstable) -- diverges');
yline(0, 'k--', 'HandleVisibility', 'off');
xlabel('Time [s]');
ylabel('y(t)');
title('2nd-order system response: LHP vs RHP poles, y(0)=-1');
legend('Location','best');
ylim([-2, 2]);
set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(scriptDir, 'LHPvsRHP_2ndOrder.png'));
fprintf('Saved plot to LHPvsRHP_2ndOrder.png\n');

poles_LHP = roots([1, 2*zeta*wn, wn^2]);
poles_RHP = roots([1, -2*zeta*wn, wn^2]);
fprintf('LHP poles: %s\n', mat2str(poles_LHP, 4));
fprintf('RHP poles: %s\n', mat2str(poles_RHP, 4));
