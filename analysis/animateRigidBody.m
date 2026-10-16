function animateRigidBody(log, frameStep, boxSize, videoFile)
%ANIMATERIGIDBODY Simple 3D playback of the simulated trajectory: draws a
%box with body-frame axes at the rigid body's position/attitude.
%
%   frameStep  render every Nth logged sample (speed up playback)
%   boxSize    [lx ly lz] body dimensions [m] for the drawn box
%   videoFile  optional output path (e.g. 'out.mp4'); if given, every
%              rendered frame is also written to an MPEG-4 video at a
%              frame rate chosen so the video plays back at the same
%              speed as the simulated time (1 video-second = 1 sim-
%              second). Recording works from a headless `matlab -batch`
%              call too (no live display needed) -- useful for reviewing
%              a run without opening MATLAB interactively.

if nargin < 2 || isempty(frameStep), frameStep = 5; end
if nargin < 3 || isempty(boxSize),   boxSize = [0.3 0.2 0.15]; end
if nargin < 4, videoFile = ''; end

recording = ~isempty(videoFile);
if recording
    dt = log.t(2) - log.t(1);
    fps = min(60, max(5, 1 / (dt * frameStep)));   % keep playback speed ~real-time, clamped to a sane range
    vw = VideoWriter(videoFile, 'MPEG-4');
    vw.FrameRate = fps;
    open(vw);
end

t = log.t;
r = log.x(1:3, :);
q = log.x(7:10, :);

% Unit box vertices (centered), scaled by boxSize
hx = boxSize(1)/2; hy = boxSize(2)/2; hz = boxSize(3)/2;
V0 = [-hx -hy -hz; hx -hy -hz; hx hy -hz; -hx hy -hz; ...
      -hx -hy  hz; hx -hy  hz; hx hy  hz; -hx hy  hz]';
faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

fig = figure('Name', 'Rigid Body Animation', 'Color', 'w');
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
xlabel(ax, 'x [m]'); ylabel(ax, 'y [m]'); zlabel(ax, 'z [m]');
view(ax, 3);

pad = 1.2 * max([1e-3; abs(r(:))]);
xlim(ax, [-pad pad]); ylim(ax, [-pad pad]); zlim(ax, [-pad pad]);

plot3(ax, r(1,:), r(2,:), r(3,:), 'Color', [0.6 0.6 0.6]);

hBody = patch(ax, 'Faces', faces, 'Vertices', V0', ...
    'FaceColor', [0.2 0.5 0.9], 'FaceAlpha', 0.6, 'EdgeColor', 'k');
hX = plot3(ax, [0 0], [0 0], [0 0], 'r', 'LineWidth', 2);
hY = plot3(ax, [0 0], [0 0], [0 0], 'g', 'LineWidth', 2);
hZ = plot3(ax, [0 0], [0 0], [0 0], 'b', 'LineWidth', 2);
titleH = title(ax, '');

axLen = 1.5 * max(boxSize);

for k = 1:frameStep:numel(t)
    R = quat2rotmCustom(q(:, k));
    Vk = R * V0 + r(:, k);
    hBody.Vertices = Vk';

    origin = r(:, k);
    set(hX, 'XData', origin(1) + [0 R(1,1)]*axLen, 'YData', origin(2) + [0 R(2,1)]*axLen, 'ZData', origin(3) + [0 R(3,1)]*axLen);
    set(hY, 'XData', origin(1) + [0 R(1,2)]*axLen, 'YData', origin(2) + [0 R(2,2)]*axLen, 'ZData', origin(3) + [0 R(3,2)]*axLen);
    set(hZ, 'XData', origin(1) + [0 R(1,3)]*axLen, 'YData', origin(2) + [0 R(2,3)]*axLen, 'ZData', origin(3) + [0 R(3,3)]*axLen);

    titleH.String = sprintf('t = %.2f s', t(k));

    if recording
        drawnow;                      % render every frame (no throttling) before capture
        writeVideo(vw, getframe(fig));
    else
        drawnow limitrate;
    end
end

if recording
    close(vw);
    fprintf('Saved animation to %s\n', videoFile);
end

end
