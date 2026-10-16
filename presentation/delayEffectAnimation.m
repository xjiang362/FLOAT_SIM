%% delayEffectAnimation.m
% Presentation demo: shows what a sequential (multiplexed), per-channel
% sensor delay does to a rigid plate's reconstructed height/tilt, even
% when the plate itself only ever translates straight up (no real tilt).
%
% Setup: a 6x6 mm plate heaves at a constant rate (30 um / 10 ms) through
% a zigzag: up 4 ms, down 2 ms, up 2 ms, down 2 ms (0 -> 12 -> 6 -> 12 ->
% 6 um). Its true height is sensed at 4 corners, but the 4 channels are
% read one at a time off a multiplexed ADC rather than simultaneously:
% corner i's reading lags the nominal sample time by i*250us (corner 1
% read first/freshest, corner 4 last/stalest -- reorder chOrderDelay
% below if your hardware scans corners in a different sequence). The
% plate pose is then reconstructed every 1 ms from those 4 (differently
% stale) readings using the simple rigid-plate estimator
%   z    = (z1+z2+z3+z4)/4
%   roll = (z1+z2-z3-z4)/(2W)
%   pitch= (z2+z3-z1-z4)/(2L)
% Because the 4 corners are sampled at 4 different true times during
% translation, this estimator reports a spurious roll/pitch on top of a
% lagged z -- an artifact of the delay, not a real deformation.
%
% Saves two videos next to this script: the true (undelayed, flat)
% motion, and the delay-corrupted reconstruction (with the true height
% overlaid as a faint reference for comparison).

clear; clc; close all;

%% Parameters
L = 6e-3;  W = 6e-3;           % [m] plate footprint (x length, y width)
zMin = 0;  zMax = 30e-6;       % [m] heave travel range, bottom to top
chDelay = 250e-6;              % [s] per-channel sequential read delay

% Motion profile: a zigzag built from up/down phases, all moving at the
% same rate (zMax-zMin)/sum(segDur) -- so shorter phases cover less of
% the range rather than moving faster. segDir(i)=+1 is "up", -1 is "down".
segDur = [4, 2, 2, 2] * 1e-3;    % [s] phase durations
segDir = [+1, -1, +1, -1];       % up/down per phase
Ttotal = sum(segDur);            % [s] total motion duration

frameDt = 1e-3;                 % [s] time between plotted/animated frames
nFrames = 10;
frameTimes = (1:nFrames) * frameDt;   % 1..10 ms -> last frame is at the top

tiltExaggeration = 1;   % visual multiplier on roll/pitch in the rendered
                         % patch; 1 = physically accurate (the delay-
                         % induced tilt here is already a few um across
                         % the plate, so no exaggeration is needed to see
                         % it -- raise this if you shrink chDelay/Ttotal
                         % and the effect gets too small to read).
fps = 2;                 % slow, presentation-paced playback (not real-time)

scriptDir = fileparts(mfilename('fullpath'));
trueVideoFile  = fullfile(scriptDir, 'delayEffect_truePlate.mp4');
reconVideoFile = fullfile(scriptDir, 'delayEffect_reconstructedPlate.mp4');

%% True trajectory: rigid, flat plate; height is a piecewise-linear
% zigzag built from segDur/segDir at a constant rate
rate = (zMax - zMin) / Ttotal;   % [m/s] same speed for every phase
tBreaks = [0, cumsum(segDur)];
zBreaks = zMin + [0, cumsum(segDir .* segDur * rate)];
zTraj = @(t) interp1(tBreaks, zBreaks, min(max(t, 0), Ttotal), 'linear');

%% Corner layout (goes counter-clockwise: front-right, front-left,
% back-left, back-right). This order is what makes the roll/pitch
% formulas above come out simple -- see the algebra in the estimator
% loop below if you reorder the corners.
Lh = L/2; Wh = W/2;
cx = [ Lh; -Lh; -Lh;  Lh];
cy = [ Wh;  Wh; -Wh; -Wh];

chOrderDelay = (1:4)' * chDelay;   % corner i's read lags by i*chDelay

%% Sample the true trajectory once per corner per frame, then run the
% simple rigid-plate estimator on the 4 (staggered-in-time) readings
zTrueFrame = zeros(1, nFrames);
zMeasAll   = zeros(4, nFrames);
z0Est      = zeros(1, nFrames);
rollEst    = zeros(1, nFrames);
pitchEst   = zeros(1, nFrames);

for k = 1:nFrames
    tk = frameTimes(k);
    zTrueFrame(k) = zTraj(tk);

    zMeas = zTraj(tk - chOrderDelay);   % 4x1, one delayed sample per corner
    zMeasAll(:, k) = zMeas;

    z0Est(k)    = mean(zMeas);
    rollEst(k)  = (zMeas(1) + zMeas(2) - zMeas(3) - zMeas(4)) / (2*W);
    pitchEst(k) = (zMeas(2) + zMeas(3) - zMeas(1) - zMeas(4)) / (2*L);
end

fprintf('Frame   t[ms]   z_true[um]  z_est[um]  err[nm]   roll[urad]  pitch[urad]\n');
for k = 1:nFrames
    fprintf('%5d  %6.1f  %10.3f  %9.3f  %7.1f  %11.3f  %11.3f\n', ...
        k, frameTimes(k)*1e3, zTrueFrame(k)*1e6, z0Est(k)*1e6, ...
        (z0Est(k) - zTrueFrame(k))*1e9, rollEst(k)*1e6, pitchEst(k)*1e6);
end

%% Figure: the 4 delayed corner readings vs. the true (undelayed) height
figCorners = figure('Color', 'w', 'Name', 'Corner readings vs. true z');
cornerColors = lines(4);
hold on; grid on; box on;
plot(frameTimes*1e3, zTrueFrame*1e6, 'k-', 'LineWidth', 2.5, 'DisplayName', 'z real');
cornerLabels = {'z_1', 'z_2', 'z_3', 'z_4'};
for i = 1:4
    plot(frameTimes*1e3, zMeasAll(i,:)*1e6, '-o', 'Color', cornerColors(i,:), ...
        'LineWidth', 1.2, 'DisplayName', cornerLabels{i});
end
xlim([0 10]); xlabel('time [ms]'); ylabel('z [\mum]');
title('Delayed corner readings vs. true plate height');
legend('Location', 'best');
exportgraphics(figCorners, fullfile(scriptDir, 'delayEffect_cornerHeights.png'));

%% Figure: roll (theta_x) -- real (always 0, the plate never tilts) vs.
% the estimator's delay-corrupted measurement
figRoll = figure('Color', 'w', 'Name', 'theta_x real vs. measured');
hold on; grid on; box on;
plot(frameTimes*1e3, zeros(1, nFrames), 'k-', 'LineWidth', 2.5, 'DisplayName', '\theta_x real');
plot(frameTimes*1e3, rollEst*1e6, '-o', 'Color', [0.85 0.4 0.1], ...
    'LineWidth', 1.2, 'DisplayName', '\theta_x measured');
xlim([0 10]); xlabel('time [ms]'); ylabel('\theta_x [\murad]');
title('Roll (\theta_x): real vs. delay-corrupted estimate');
legend('Location', 'best');
exportgraphics(figRoll, fullfile(scriptDir, 'delayEffect_rollComparison.png'));

%% Shared z-axis range for both videos, with margin for the (possibly
% exaggerated) tilt so the patch never clips the axes
tiltMargin = tiltExaggeration * max([abs(rollEst)*Wh, abs(pitchEst)*Lh]);
zAxisRange = [zMin, zMax + tiltMargin] * 1e6;   % [um]
pad = 0.05 * (zAxisRange(2) - zAxisRange(1));
zAxisRange = zAxisRange + [-pad, pad];

%% Animation 1: true plate motion (flat, no tilt, no delay)
renderPlateVideo(trueVideoFile, fps, cx, cy, zAxisRange, frameTimes, ...
    zTrueFrame, zeros(1,nFrames), zeros(1,nFrames), 1, [], ...
    'True plate motion (no sensor delay)');

%% Animation 2: reconstructed plate from the sequentially-delayed reads
renderPlateVideo(reconVideoFile, fps, cx, cy, zAxisRange, frameTimes, ...
    z0Est, rollEst, pitchEst, tiltExaggeration, zTrueFrame, ...
    sprintf('Reconstructed plate (%.0f us/channel delay)', chDelay*1e6));


function renderPlateVideo(videoFile, fps, cx, cy, zAxisRangeUm, frameTimes, ...
    z0Arr, rollArr, pitchArr, tiltExaggeration, zRefArr, titleStr)
%RENDERPLATEVIDEO Draws + saves an N-frame video of a flat/tilted plate
% reconstructed from [z0, roll, pitch] at each frame. zRefArr (optional,
% pass [] to skip) overlays a faint flat reference patch at the true
% height, for a visual before/after comparison.

nFrames = numel(frameTimes);
faces = [1 2 3 4];

fig = figure('Color', 'w', 'Name', titleStr);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); view(ax, 35, 20);
xlabel(ax, 'x [mm]'); ylabel(ax, 'y [mm]'); zlabel(ax, 'z [\mum]');
xlim(ax, [-1 1] * max(abs(cx)) * 1.3e3);
ylim(ax, [-1 1] * max(abs(cy)) * 1.3e3);
zlim(ax, zAxisRangeUm);   % z not to physical scale vs. x/y (30 um travel on a 6 mm plate)

if ~isempty(zRefArr)
    hRef = patch(ax, 'Faces', faces, 'Vertices', zeros(4,3), ...
        'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
hPlate = patch(ax, 'Faces', faces, 'Vertices', zeros(4,3), ...
    'FaceColor', [0.85 0.4 0.1], 'FaceAlpha', 0.85, 'EdgeColor', 'k', 'LineWidth', 1.2);
if ~isempty(zRefArr)
    legend(ax, [hPlate hRef], {'reconstructed', 'true (no delay)'}, 'Location', 'northoutside');
end
titleH = title(ax, '');

vw = VideoWriter(videoFile, 'MPEG-4');
vw.FrameRate = fps;
open(vw);

for k = 1:nFrames
    zDraw = z0Arr(k) + tiltExaggeration * (rollArr(k)*cy - pitchArr(k)*cx);   % 4x1 [m]
    hPlate.Vertices = [cx*1e3, cy*1e3, zDraw*1e6];

    if ~isempty(zRefArr)
        hRef.Vertices = [cx*1e3, cy*1e3, repmat(zRefArr(k)*1e6, 4, 1)];
    end

    titleH.String = sprintf('%s\nt = %.1f ms', titleStr, frameTimes(k)*1e3);

    drawnow;
    writeVideo(vw, getframe(fig));
end

close(vw);
fprintf('Saved %s\n', videoFile);

end
