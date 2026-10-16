%% dedicatedSingleAxisSweeps.m
% Single-variable characterization sweeps for all 6 DOF: for each of
% rx, ry, rz, thx, thy, thz in turn, sweep that one variable while
% holding the other five at exactly zero. nPts points per axis,
% combined into one COMSOL sweep file (output filename reflects the
% total point count, e.g. comsol_sweep_combined120.txt for nPts=20).
%
% Non-uniform spacing per axis: generated from points uniform in u in
% [-1,1] via x = sign(u)*|u|^p*xmax (p=2), which concentrates
% resolution near the equilibrium point (0,0,0,0,0,0) where the
% controller actually operates, while still reaching the full range
% to capture turnover/softening behavior at larger displacement (as
% seen for Fx/Fy/Tauz). 50 points per axis gives a comfortable ~8:1
% data:parameter ratio for a poly5 fit, versus the original 20-point
% pilot's ~3:1 (overfitting territory).
%
% thx/thy/thz's max ranges depend on rz (or rx/ry) per
% sampleTrapPoseLHS.m's constraint formula; since every other DOF is
% zero in each of these single-axis sweeps, that reduces to the
% simplest case (full tilt_disp_factor budget against the nominal
% gap, no translation eating into it).
%
% NOTE: the geometry constants below must match sampleTrapPoseLHS.m --
% duplicated here rather than shared, since this is a one-off
% diagnostic sweep, not part of the regular LHS pipeline.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

nPts = 20;
p = 2;   % warp exponent: higher = more concentrated near zero
paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};

u = linspace(-1, 1, nPts)';
warpedSpacing = @(xmax) sign(u) .* abs(u).^p * xmax;
zero1 = zeros(nPts, 1);

% --- Nominal geometry (must match sampleTrapPoseLHS.m) ---
g0_z  = 15e-6;
g0_xy = 50e-6;
L = 6000e-6; W = L;
tilt_disp_factor = 0.9;

x_max = 45e-6;
z_max = 14e-6;
th_xy_max_at_zero = asind((tilt_disp_factor * g0_z) / (L/2 + W/2));       % thx, thy at rz=0
thz_max_at_zero   = asind((tilt_disp_factor * g0_xy) / max(L/2, W/2));    % thz at rx=ry=0

rxSweep  = warpedSpacing(x_max);
rySweep  = warpedSpacing(x_max);
rzSweep  = warpedSpacing(z_max);
thxSweep = warpedSpacing(th_xy_max_at_zero);
thySweep = warpedSpacing(th_xy_max_at_zero);
thzSweep = warpedSpacing(thz_max_at_zero);

axisSweeps = {rxSweep, rySweep, rzSweep, thxSweep, thySweep, thzSweep};
axisNames  = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};

fprintf('%6s  %14s  %14s\n', 'axis', 'min', 'max');
for i = 1:6
    fprintf('%6s  %14.4g  %14.4g\n', axisNames{i}, min(axisSweeps{i}), max(axisSweeps{i}));
end

%% Combined file: all 300 points in one sweep (single COMSOL run)
% Column i of a row is zero unless that row belongs to axis i's block.
nAxes = 6;
combined = zeros(nAxes*nPts, nAxes);
for i = 1:nAxes
    rows = (i-1)*nPts + (1:nPts);
    combined(rows, i) = axisSweeps{i};
end

paramData = num2cell(combined, 1);
outFile = fullfile(scriptDir, sprintf('comsol_sweep_combined%d.txt', nAxes*nPts));
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);

fprintf('\nSaved combined sweep (%d pts) to %s\n', nAxes*nPts, outFile);
fprintf('Row blocks: ');
for i = 1:nAxes
    fprintf('%s=%d-%d  ', axisNames{i}, (i-1)*nPts+1, i*nPts);
end
fprintf('\n');
