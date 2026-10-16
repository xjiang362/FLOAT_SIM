%% compareRealVsComplexPole.m
% Compares the impulse response of three first-order systems sharing
% the same real part of their pole (same decay rate p), but different
% imaginary parts:
%   1/(s+p)        -- pure real pole, no oscillation
%   1/(s+p+j*w)     -- complex pole, small imaginary part (slow oscillation)
%   1/(s+p+j*W)     -- complex pole, large imaginary part (fast oscillation)
% For a real pole, g(t) = exp(-p*t). For a complex pole s=-(p+jw), the
% impulse response exp(-(p+jw)t) is complex-valued; the REAL PART,
% Re{g(t)} = exp(-p*t)*cos(w*t), is plotted (the physically meaningful
% part for a real system, as if paired with its conjugate pole) -- same
% decay envelope exp(-p*t) in all three cases, differing only in
% oscillation frequency.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

p = 1;       % decay rate (shared by all three)
w = 2;       % small angular frequency [rad/s]
W = 20;      % large angular frequency [rad/s]

t = linspace(0, 5, 4000)';
g_real     = exp(-p*t);
g_small_w  = exp(-p*t) .* cos(w*t);
g_large_W  = exp(-p*t) .* cos(W*t);

figure('Color','w','Name','Real pole vs complex pole (small/large imaginary part)');
hold on; grid on;
plot(t, g_real, 'k-', 'LineWidth', 2, 'DisplayName', '1/(s+p): real pole');
plot(t, g_small_w, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('1/(s+p+j\\omega): \\omega=%g (slow)', w));
plot(t, g_large_W, 'r-', 'LineWidth', 1, 'DisplayName', sprintf('1/(s+p+j\\Omega): \\Omega=%g (fast)', W));
plot(t, exp(-p*t), 'k--', 'LineWidth', 1, 'HandleVisibility','off');
plot(t, -exp(-p*t), 'k--', 'LineWidth', 1, 'HandleVisibility','off');
yline(0, 'k:', 'HandleVisibility', 'off');
xlabel('Time [s]');
ylabel('Re\{y(t)\}');
title('Impulse response: real pole vs complex pole (small vs large \omega)');
legend('Location','best');
set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(scriptDir, 'RealVsComplexPole_impulse.png'));
fprintf('Saved plot to RealVsComplexPole_impulse.png\n');
