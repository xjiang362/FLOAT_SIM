%% dedicatedSingleAxisSweeps_uniform.m
% Same dedicated single-axis characterization sweep as
% dedicatedSingleAxisSweeps.m (for each of rx, ry, rz, thx, thy, thz in
% turn, sweep that one variable while holding the other five at exactly
% zero, 50 points per axis, 300 total), but with UNIFORM spacing
% (linspace) instead of the warped sign(u)*|u|^p*xmax spacing -- so
% points are evenly distributed across the full range rather than
% concentrated near zero.
%
% Same nominal geometry constants as dedicatedSingleAxisSweeps.m and
% sampleTrapPoseLHS.m (duplicated here, not shared, since this is a
% one-off diagnostic sweep).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

nPts = 50;
paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};

% --- Nominal geometry (must match sampleTrapPoseLHS.m) ---
g0_z  = 15e-6;
g0_xy = 50e-6;
L = 6000e-6; W = L;
tilt_disp_factor = 0.9;

x_max = 45e-6;
z_max = 14e-6;
th_xy_max_at_zero = asind((tilt_disp_factor * g0_z) / (L/2 + W/2));       % thx, thy at rz=0
thz_max_at_zero   = asind((tilt_disp_factor * g0_xy) / max(L/2, W/2));    % thz at rx=ry=0

uniformSpacing = @(xmax) linspace(-xmax, xmax, nPts)';

rxSweep  = uniformSpacing(x_max);
rySweep  = uniformSpacing(x_max);
rzSweep  = uniformSpacing(z_max);
thxSweep = uniformSpacing(th_xy_max_at_zero);
thySweep = uniformSpacing(th_xy_max_at_zero);
thzSweep = uniformSpacing(thz_max_at_zero);

axisSweeps = {rxSweep, rySweep, rzSweep, thxSweep, thySweep, thzSweep};
axisNames  = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};

fprintf('%6s  %14s  %14s\n', 'axis', 'min', 'max');
for i = 1:6
    fprintf('%6s  %14.4g  %14.4g\n', axisNames{i}, min(axisSweeps{i}), max(axisSweeps{i}));
end

%% Combined file: all 300 points in one sweep (single COMSOL run)
nAxes = 6;
combined = zeros(nAxes*nPts, nAxes);
for i = 1:nAxes
    rows = (i-1)*nPts + (1:nPts);
    combined(rows, i) = axisSweeps{i};
end

paramData = num2cell(combined, 1);
outFile = fullfile(scriptDir, sprintf('comsol_sweep_combined%d_uniform.txt', nAxes*nPts));
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);

fprintf('\nSaved combined uniform sweep (%d pts) to %s\n', nAxes*nPts, outFile);
fprintf('Row blocks: ');
for i = 1:nAxes
    fprintf('%s=%d-%d  ', axisNames{i}, (i-1)*nPts+1, i*nPts);
end
fprintf('\n');
fprintf('\nRun with the same probes as comsol_sweep_combined300_results.txt (MST force/torque\n');
fprintf('probes 1,2,3,5,6,7), and include Global Variable Probe 11 (Total electric energy)\n');
fprintf('if you want the energy-method comparison too.\n');
