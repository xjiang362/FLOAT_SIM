%% compareSmallVsLargePW.m
% Compares the impulse response (real part) of three complex-pole
% systems, varying decay rate and oscillation frequency independently:
%   1/(s+p+j*w)  -- small decay rate p, small frequency w
%   1/(s+P+j*w)  -- LARGE decay rate P, same small frequency w
%   1/(s+P+j*W)  -- large decay rate P, LARGE frequency W
% so the middle case isolates the effect of decay rate alone (same
% oscillation as case 1, faster envelope), and the third isolates the
% added effect of frequency alone (same fast envelope as case 2, faster
% oscillation). Re{y(t)} = exp(-a*t)*cos(b*t) for pole -(a+jb).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

p = 1;    % small decay rate
P = 6;    % large decay rate
w = 2;    % small angular frequency [rad/s]
W = 18;   % large angular frequency [rad/s]

t = linspace(0, 4, 4000)';
g_p_w = exp(-p*t) .* cos(w*t);
g_P_w = exp(-P*t) .* cos(w*t);
g_P_W = exp(-P*t) .* cos(W*t);

figure('Color','w','Name','Small vs large decay rate and frequency');
hold on; grid on;
plot(t, g_p_w, 'b-', 'LineWidth', 1.75, 'DisplayName', sprintf('1/(s+p+j\\omega): p=%g, \\omega=%g', p, w));
plot(t, g_P_w, 'Color', [0.9 0.4 0], 'LineWidth', 1.5, 'DisplayName', sprintf('1/(s+P+j\\omega): P=%g, \\omega=%g', P, w));
plot(t, g_P_W, 'r-', 'LineWidth', 1, 'DisplayName', sprintf('1/(s+P+j\\Omega): P=%g, \\Omega=%g', P, W));
yline(0, 'k:', 'HandleVisibility', 'off');
xlabel('Time [s]');
ylabel('Re\{y(t)\}');
title('Impulse response: small vs large decay rate and frequency');
legend('Location','best');
set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(scriptDir, 'SmallVsLargePW_impulse.png'));
fprintf('Saved plot to SmallVsLargePW_impulse.png\n');
