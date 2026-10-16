%% centerElectrodeTauRzFineSweep.m
% The rz=-7/0/+7um 3-point check showed Taux(thx)'s SLOPE is non-monotonic
% in rz (peaks near rz=0, falls off much faster toward rz<0 than rz>0) --
% a shape no linear-in-rz model (a2+a1*rz) or 1/Geff^n power law can fit.
% Need the actual slope(rz) curve at finer resolution across the full
% established +-14um domain to identify the real functional form.
%
% Design: fine rz sweep (25 pts, -14 to +14um) at a FIXED nonzero thx
% (thy=0), and a mirrored fine rz sweep (25 pts) at a fixed nonzero thy
% (thx=0). At each rz, Taux/thx_fixed (Tauy/thy_fixed) is then a direct,
% densely-sampled reading of the slope(rz) shape (Taux already confirmed
% ~linear in thx at each of the 3 known rz's, so dividing by a single
% fixed thx is valid here). MST is sufficient (large, clean signal,
% R^2>0.999 at rz=0).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

rz_max = 14e-6;   % [m] established full domain
nPts = 25;
rzSweep = linspace(-rz_max, rz_max, nPts)';

% Fixed tilt: reuse a mid-range value from the original thx/thy sweep
% domain (not the extreme, to stay well clear of any collision-clearance
% edge effects at the same time as sweeping rz to its extreme).
g0_z = 15e-6; L = 6000e-6; W = L; tilt_disp_factor = 0.9;
th_xy_max_at_zero = asind((tilt_disp_factor * g0_z) / (L/2 + W/2));
thx_fixed = 0.6 * th_xy_max_at_zero;
thy_fixed = 0.6 * th_xy_max_at_zero;

% Block 1: rz sweep at fixed thx, thy=0 (25 rows)
rzB1 = rzSweep; thxB1 = thx_fixed*ones(nPts,1); thyB1 = zeros(nPts,1);
% Block 2: rz sweep at fixed thy, thx=0 (25 rows)
rzB2 = rzSweep; thyB2 = thy_fixed*ones(nPts,1); thxB2 = zeros(nPts,1);

rz  = [rzB1; rzB2];
thx = [thxB1; thxB2];
thy = [thyB1; thyB2];
rx = zeros(size(rz)); ry = zeros(size(rz)); thz = zeros(size(rz));

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData = {rx, ry, rz, thx, thy, thz};

outFile = fullfile(scriptDir, 'comsol_centerElectrodeTauRzFine50.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);

fprintf('Saved %d points to %s\n', numel(rz), outFile);
fprintf('Row blocks: rz-sweep @ thx=%.4fdeg (thy=0) = 1-%d; rz-sweep @ thy=%.4fdeg (thx=0) = %d-%d\n', ...
    thx_fixed, nPts, thy_fixed, nPts+1, 2*nPts);
fprintf('rz range: [%.1f, %.1f] um, %d pts, linear spacing\n', -rz_max*1e6, rz_max*1e6, nPts);
fprintf('\nVoltage should already be applied from the last run -- just re-run this\n');
fprintf('parameter file with the same probe setup (MST is sufficient).\n');
