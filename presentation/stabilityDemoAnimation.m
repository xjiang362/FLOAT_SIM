%% stabilityDemoAnimation.m
% Presentation demo: two 6x6 mm plates side by side, both released
% from z=-5mm, following explicit keyframe trajectories -- one
% oscillates with a growing envelope (unstable), the other oscillates
% with a decaying envelope that settles to z=0 (stable). Illustrates
% the practical consequence of the sign of the local restoring-force
% slope (dF/dx) extracted from the COMSOL sweeps: a stabilizing
% (restoring) force damps oscillations back to equilibrium, while an
% unstable/destabilizing one grows them without bound.
%
% Saves an animation (both plates together) and a static z(t) plot.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

%% Explicit keyframe trajectories [mm]
zUnstableMM = [-5, 1, 4, 0, -5, 5, -5, 5];
zStableMM   = [-5, -3, -2, -1, 0, 0, 0, 0];
nFrames = numel(zUnstableMM);

Ttotal = 10e-3;   % [s] total animation duration
t = linspace(0, Ttotal, nFrames)';

zUnstable = zUnstableMM(:) * 1e-3;   % [m]
zStable   = zStableMM(:) * 1e-3;     % [m]

%% Static plot: z(t) for both trajectories
figure('Color', 'w', 'Name', 'Stable vs unstable trajectories');
hold on; grid on;
plot(t*1e3, zStableMM, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'stable (dF/dx < 0)');
plot(t*1e3, zUnstableMM, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'DisplayName', 'unstable (dF/dx > 0)');
yline(0, 'k--', 'DisplayName', 'equilibrium');
xlabel('time [ms]'); ylabel('z [mm]');
title('Stable vs. unstable oscillation envelope');
legend('Location', 'best');
exportgraphics(gcf, fullfile(scriptDir, 'stabilityDemo_trajectories.png'));

%% Animation: two plates side by side, each moving per its own z(t)
L = 6e-3; W = 6e-3;      % [m] plate footprint
Lh = L/2; Wh = W/2;
faces = [1 2 3 4];
cx = [ Lh; -Lh; -Lh;  Lh];
cy = [ Wh;  Wh; -Wh; -Wh];

xOffset = 8e-3;   % [m] lateral separation between the two plates in the scene
zMax = 5e-3;      % [m] max height, for axis limits

fig = figure('Color', 'w', 'Name', 'Stability demo animation');
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); view(ax, 35, 20);
xlabel(ax, 'x [mm]'); ylabel(ax, 'y [mm]'); zlabel(ax, 'z [mm]');
xlim(ax, [-1 1] * (xOffset + Lh) * 1.3e3);
ylim(ax, [-1 1] * Wh * 1.6e3);
zlim(ax, [-zMax, zMax] * 1.2e3);

hStable = patch(ax, 'Faces', faces, 'Vertices', zeros(4,3), ...
    'FaceColor', [0.2 0.5 0.9], 'FaceAlpha', 0.85, 'EdgeColor', 'k', 'LineWidth', 1.2);
hUnstable = patch(ax, 'Faces', faces, 'Vertices', zeros(4,3), ...
    'FaceColor', [0.85 0.2 0.2], 'FaceAlpha', 0.85, 'EdgeColor', 'k', 'LineWidth', 1.2);
legend(ax, [hStable, hUnstable], {'stable', 'unstable'}, 'Location', 'northoutside');
titleH = title(ax, '');

videoFile = fullfile(scriptDir, 'stabilityDemo_animation.mp4');
vw = VideoWriter(videoFile, 'MPEG-4');
vw.FrameRate = 2;   % slow, presentation-paced -- only 8 keyframes
open(vw);

for k = 1:nFrames
    hStable.Vertices = [(cx - xOffset)*1e3, cy*1e3, repmat(zStableMM(k), 4, 1)];
    hUnstable.Vertices = [(cx + xOffset)*1e3, cy*1e3, repmat(zUnstableMM(k), 4, 1)];
    titleH.String = sprintf('t = %.2f ms   |   stable: %.1f mm   unstable: %.1f mm', ...
        t(k)*1e3, zStableMM(k), zUnstableMM(k));

    drawnow;
    writeVideo(vw, getframe(fig));
end

close(vw);
fprintf('Saved animation to %s\n', videoFile);
