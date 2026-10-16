%% centerElectrodeRzCrossTermSweep.m
% Tests whether Taux(thx), Tauy(thy) depend on rz (an rz*thx / rz*thy
% cross term, analogous to the quadrant electrodes' Leqx/Leqy moment-arm
% model) -- structurally invisible in the original rz-only/thx-only/
% thy-only sweep, since rz=0 throughout the thx-only and thy-only
% blocks there.
%
% Design: thx-primary sweep (50 pts, same range as before) at rz=+-7um,
% and thy-primary sweep (50 pts) at rz=+-7um -- an rz bracket at HALF
% the established +-14um range (not a tiny +-1um bracket like the
% energy-method sweeps needed): Taux/Tauy are large, clean MST signals
% here (R^2>0.999 in the rz=0 fit), so no noise-avoidance reason to use
% a small bracket -- a bigger one gives a clearer read on how the
% thx-slope changes with rz. Reads MST directly (no energy-method
% bracket-differencing needed for this).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

nPts = 50;
p = 2;
g0_z = 15e-6;
L = 6000e-6; W = L;
tilt_disp_factor = 0.9;
th_xy_max_at_zero = asind((tilt_disp_factor * g0_z) / (L/2 + W/2));

rz_bracket = 7e-6;   % [m] half of the established +-14um range

u = linspace(-1, 1, nPts)';
warpedSpacing = @(xmax) sign(u) .* abs(u).^p * xmax;
thxSweep = warpedSpacing(th_xy_max_at_zero);
thySweep = warpedSpacing(th_xy_max_at_zero);

% Block 1: thx-primary, rz=-7um and +7um (100 rows)
% Block 2: thy-primary, rz=-7um and +7um (100 rows)
thx1 = repelem(thxSweep, 2);
rz1  = repmat([-rz_bracket; rz_bracket], nPts, 1);
thy1 = zeros(2*nPts,1);

thy2 = repelem(thySweep, 2);
rz2  = repmat([-rz_bracket; rz_bracket], nPts, 1);
thx2 = zeros(2*nPts,1);

thx = [thx1; thx2];
thy = [thy1; thy2];
rz  = [rz1; rz2];
rx = zeros(size(rz)); ry = zeros(size(rz)); thz = zeros(size(rz));

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData = {rx, ry, rz, thx, thy, thz};

outFile = fullfile(scriptDir, 'comsol_centerElectrode_rzCrossTerm200.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);

fprintf('Saved %d points to %s\n', numel(rz), outFile);
fprintf('Row blocks: thx-primary (rz=+-%.4gm)=1-100, thy-primary (rz=+-%.4gm)=101-200\n', rz_bracket, rz_bracket);
fprintf('thx/thy range: [%.4f, %.4f] deg\n', min(thxSweep), max(thxSweep));
fprintf('\nVoltage should already be applied from the last run -- just re-run this\n');
fprintf('parameter file with the same probe setup (MST is sufficient, energy optional).\n');
