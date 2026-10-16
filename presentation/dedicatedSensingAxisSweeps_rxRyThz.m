%% dedicatedSensingAxisSweeps_rxRyThz.m
% Dedicated single-variable characterization sweeps for the sensing
% (capacitance-matrix) model, covering the remaining 3 DOF not covered
% by dedicatedSensingAxisSweeps.m (rz, thx, thy): rx, ry, thz. Same
% approach -- sweep one variable while holding the other five at exactly
% zero -- to isolate each DOF's effect on the capacitance matrix
% cleanly.
%
% Non-uniform spacing per axis (same warp as the other sweeps): points
% uniform in u in [-1,1] via x = sign(u)*|u|^p*xmax (p=2), concentrating
% resolution near the equilibrium point while still reaching the full
% range.
%
% NOTE: geometry constants below must match sampleTrapPoseLHS.m --
% duplicated here rather than shared, since this is a one-off
% diagnostic sweep, not part of the regular LHS pipeline.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

nPts = 20;
p = 2;   % warp exponent: higher = more concentrated near zero

% --- Nominal geometry (must match sampleTrapPoseLHS.m) ---
g0_xy = 50e-6;
L = 6000e-6; W = L;
tilt_disp_factor = 0.9;

x_max = 45e-6;
thz_max_at_zero = asind((tilt_disp_factor * g0_xy) / max(L/2, W/2));   % thz at rx=ry=0

u = linspace(-1, 1, nPts)';
warpedSpacing = @(xmax) sign(u) .* abs(u).^p * xmax;

rxSweep  = warpedSpacing(x_max);
rySweep  = warpedSpacing(x_max);
thzSweep = warpedSpacing(thz_max_at_zero);

axisSweeps = {rxSweep, rySweep, thzSweep};
axisNames  = {'rx', 'ry', 'thz'};

fprintf('%6s  %14s  %14s\n', 'axis', 'min', 'max');
for i = 1:3
    fprintf('%6s  %14.4g  %14.4g\n', axisNames{i}, min(axisSweeps{i}), max(axisSweeps{i}));
end

%% Combined file: all 60 points in one sweep (single COMSOL run)
nAxes = 3;
combined = zeros(nAxes*nPts, 6);
colIdx = [1, 2, 6];   % rx,ry,thz positions in the 6-column (rx,ry,rz,thx,thy,thz) layout
for i = 1:nAxes
    rows = (i-1)*nPts + (1:nPts);
    combined(rows, colIdx(i)) = axisSweeps{i};
end

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData = num2cell(combined, 1);

outFile = fullfile(scriptDir, sprintf('comsol_sensingSweep_rxRyThz%d.txt', nAxes*nPts));
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);

fprintf('\nSaved combined sensing sweep (%d pts) to %s\n', nAxes*nPts, outFile);
fprintf('Row blocks: ');
for i = 1:nAxes
    fprintf('%s=%d-%d  ', axisNames{i}, (i-1)*nPts+1, i*nPts);
end
fprintf('\n');
fprintf('\nRun with the full 22-terminal mutual capacitance matrix export enabled at\n');
fprintf('each pose point (same Stationary Source Sweep setup used for the rz/thx/thy\n');
fprintf('sensing sweep), so it can be parsed the same way.\n');
