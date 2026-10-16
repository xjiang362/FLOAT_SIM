%% thzSweepWithRxPair.m
% Dedicated thz sweep with an rx bracket (rx=-1um, +1um) at every thz
% grid point, so Fx(thz) can be derived via the energy method: for each
% thz, a central difference in rx,
%   Fx(thz) ~= [W(thz,+1um) - W(thz,-1um)] / 2um
% ry=rz=thx=thy=0 throughout -- mirrors rxSweepWithRyThzPairs.m's
% ry-bracket idea, but with the roles of rx and thz swapped (thz is now
% the primary swept variable, rx is the small bracket).
%
% Safe thz range: uses sampleTrapPoseLHS.m's yaw constraint,
% thz_max = asind(min(x_extra,y_extra)/max(L/2,W/2)) with
% x_extra = max(tilt_disp_factor*g0_xy - |rx|, 0). Since |rx|=1um at
% every point here (both bracket values), the SAME thz_max applies to
% every thz row (no per-row shrinking needed, unlike the earlier
% rx-primary sweep where |rx| ranged up to 45um) -- evaluated at the
% worst case rx=1um: thz_max = 0.8404 deg.
%
% IMPORTANT: make sure the electrode voltage is applied (matching the
% rest of this project, e.g. 30V) before running this sweep.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

nPtsThz = 50;
rx_bracket = 1e-6;   % [m] +-1um rx bracket

g0_xy = 50e-6; L = 6000e-6; W = 6000e-6; tilt_disp_factor = 0.9;
xy_budget_total = tilt_disp_factor * g0_xy;
x_extra = max(xy_budget_total - rx_bracket, 0);
thz_max_safe = asind(x_extra / max(L/2, W/2));   % deg
fprintf('thz_max_safe (accounting for rx=+-%.4g m bracket) = %.6f deg\n', rx_bracket, thz_max_safe);

thzBase = linspace(-thz_max_safe, thz_max_safe, nPtsThz)';

% 2 rows per thz value: (rx=-1um), (rx=+1um)
thz = repelem(thzBase, 2);
rx = repmat([-rx_bracket; rx_bracket], nPtsThz, 1);
zeroCol = zeros(2*nPtsThz, 1);

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData  = {rx, zeroCol, zeroCol, zeroCol, zeroCol, thz};

outFile = fullfile(scriptDir, 'comsol_sweep_thzOnly_rxPair.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);
fprintf('\nSaved %d parameters x %d combinations to %s\n', numel(paramNames), numel(thz), outFile);
fprintf('thz: %d uniformly-spaced values in [%.4g, %.4g] deg, each paired with rx = -%.4g m and +%.4g m\n', ...
    nPtsThz, min(thzBase), max(thzBase), rx_bracket, rx_bracket);
fprintf('\nMAKE SURE THE ELECTRODE VOLTAGE IS APPLIED before running this sweep --\n');
fprintf('export with Total electric energy (Global Variable Probe, matching the energy-method model),\n');
fprintf('MST Fx probe optional for comparison.\n');
