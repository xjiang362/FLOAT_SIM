%% dedicatedSensingAxisSweeps.m
% Dedicated single-variable characterization sweeps for the sensing
% (capacitance-matrix) model, mirroring dedicatedSingleAxisSweeps.m
% (used for the force/energy model): for rz, thx, thy in turn, sweep
% that one variable while holding the other five at exactly zero. This
% isolates each DOF's effect on the capacitance matrix cleanly, unlike
% the scattered 300-pt full-6D LHS sweep used so far for the S_R/Cm1p/
% Cm2p fits.
%
% Only rz, thx, thy are included (not rx, ry, thz) since those are the
% three DOF the fitted sensing models (S_R, Cm1p, Cm2p) actually depend
% on.
%
% Non-uniform spacing per axis (same warp as the force sweep): points
% uniform in u in [-1,1] via x = sign(u)*|u|^p*xmax (p=2), concentrating
% resolution near the equilibrium point (0,0,0) while still reaching the
% full range.
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
g0_z  = 15e-6;
L = 6000e-6; W = L;
tilt_disp_factor = 0.9;

z_max = 14e-6;
th_xy_max_at_zero = asind((tilt_disp_factor * g0_z) / (L/2 + W/2));   % thx, thy at rz=0

u = linspace(-1, 1, nPts)';
warpedSpacing = @(xmax) sign(u) .* abs(u).^p * xmax;
zero1 = zeros(nPts, 1);

rzSweep  = warpedSpacing(z_max);
thxSweep = warpedSpacing(th_xy_max_at_zero);
thySweep = warpedSpacing(th_xy_max_at_zero);

axisSweeps = {rzSweep, thxSweep, thySweep};
axisNames  = {'rz', 'thx', 'thy'};

fprintf('%6s  %14s  %14s\n', 'axis', 'min', 'max');
for i = 1:3
    fprintf('%6s  %14.4g  %14.4g\n', axisNames{i}, min(axisSweeps{i}), max(axisSweeps{i}));
end

%% Combined file: all 60 points in one sweep (single COMSOL run)
% Column layout matches the 6-terminal-pose convention used throughout
% (rx,ry,rz,thx,thy,thz); rx,ry,thz stay zero everywhere.
nAxes = 3;
combined = zeros(nAxes*nPts, 6);
colIdx = [3, 4, 5];   % rz,thx,thy positions in the 6-column layout
for i = 1:nAxes
    rows = (i-1)*nPts + (1:nPts);
    combined(rows, colIdx(i)) = axisSweeps{i};
end

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData = num2cell(combined, 1);

outFile = fullfile(scriptDir, sprintf('comsol_sensingSweep_rzThxThy%d.txt', nAxes*nPts));
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);

fprintf('\nSaved combined sensing sweep (%d pts) to %s\n', nAxes*nPts, outFile);
fprintf('Row blocks: ');
for i = 1:nAxes
    fprintf('%s=%d-%d  ', axisNames{i}, (i-1)*nPts+1, i*nPts);
end
fprintf('\n');
fprintf('\nRun with the full 22-terminal mutual/Maxwell/elastance capacitance matrix\n');
fprintf('export enabled at each pose point (same Stationary Source Sweep setup used\n');
fprintf('for the 300-pt capacitance sweep), so it can be parsed the same way.\n');
