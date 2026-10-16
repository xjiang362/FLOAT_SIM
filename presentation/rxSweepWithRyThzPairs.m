%% rxSweepWithRyThzPairs.m
% Extends the dedicated uniform rx sweep with TWO independent two-point
% brackets at every rx grid point:
%   - ry = -1um, +1um (thz=0)   -> Fy(rx) via energy method:
%       Fy(rx) ~= [W(rx,+1um,0) - W(rx,-1um,0)] / 2um
%   - thz = -thz_bracket(rx), +thz_bracket(rx) (ry=0) -> Tauz(rx) via
%     energy method:
%       Tauz(rx) ~= [W(rx,0,+thzB) - W(rx,0,-thzB)] / (2*thzB_rad)
% rz=thx=thy=0 throughout (only rx, and one bracket variable, ever
% nonzero in a given row).
%
% thz_bracket is NOT a single fixed value -- it's rx-DEPENDENT, using
% the full combined-usage constraint from sampleTrapPoseLHS.m:
%   xy_budget_total = tilt_disp_factor*g0_xy = 45um
%   x_extra(rx) = max(xy_budget_total - |rx|, 0)
%   thz_max(rx) = asind(x_extra(rx) / max(L/2,W/2))
% (y_extra is just xy_budget_total here since ry=0 in the thz-bracket
% rows). Using the rx=0 value of thz_max (0.8595deg) at EVERY rx --
% what the first version of this script did -- caused a real COMSOL
% meshing failure at rx=-41.3um, thz=-0.172deg: the safe thz_max at
% that rx is only ~0.070deg, so that bracket pushed a plate corner past
% the wall. Each row now uses thz_bracket(rx) = thz_max(rx)/5, so it
% shrinks safely as |rx| grows toward the 45um edge.
%
% IMPORTANT: this sweep needs the electrode VOLTAGE actually applied in
% the COMSOL model (the previous rx/ry-bracket run had no voltage set,
% which is why W was essentially flat/noise-floor -- there was no field
% energy to speak of). Set V (e.g. 30V, matching the rest of this
% project) before running this sweep.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

nPtsRx = 50;
x_max = 45e-6;      % [m] match the established rx range
ry_bracket = 1e-6;  % [m] +-1um ry bracket

g0_xy = 50e-6; L = 6000e-6; W = 6000e-6; tilt_disp_factor = 0.9;
xy_budget_total = tilt_disp_factor * g0_xy;

rxBase = linspace(-x_max, x_max, nPtsRx)';
x_extra = max(xy_budget_total - abs(rxBase), 0);
thz_max_rx = asind(x_extra / max(L/2, W/2));   % deg, one value per rx point
thz_bracket_rx = thz_max_rx / 5;               % deg, one value per rx point
fprintf('thz_max(rx) range: [%.6f, %.6f] deg -> thz_bracket(rx) range: [%.6f, %.6f] deg\n', ...
    min(thz_max_rx), max(thz_max_rx), min(thz_bracket_rx), max(thz_bracket_rx));

% 4 rows per rx value: (ry=-1um,thz=0), (ry=+1um,thz=0), (ry=0,thz=-B(rx)), (ry=0,thz=+B(rx))
rx = repelem(rxBase, 4);
ry = repmat([-ry_bracket; ry_bracket; 0; 0], nPtsRx, 1);
thz = reshape([zeros(1,nPtsRx); zeros(1,nPtsRx); -thz_bracket_rx'; thz_bracket_rx'], [], 1);
zeroCol = zeros(4*nPtsRx, 1);

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData  = {rx, ry, zeroCol, zeroCol, zeroCol, thz};

outFile = fullfile(scriptDir, 'comsol_sweep_rxOnly_ryThzPairs.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);
fprintf('\nSaved %d parameters x %d combinations to %s\n', numel(paramNames), numel(rx), outFile);
fprintf('rx: %d uniform values in [%.4g, %.4g] m, each with 4 rows:\n', nPtsRx, min(rxBase), max(rxBase));
fprintf('  (ry=-%.4g m, thz=0), (ry=+%.4g m, thz=0), (ry=0, thz=-thzB(rx)), (ry=0, thz=+thzB(rx))\n', ...
    ry_bracket, ry_bracket);
fprintf('thz bracket shrinks safely near the rx edges (see range printed above).\n');
fprintf('\nMAKE SURE THE ELECTRODE VOLTAGE IS APPLIED before running this sweep --\n');
fprintf('export with Total electric energy (Global Variable Probe 11 or equivalent),\n');
fprintf('MST Fy/Tauz probes optional for comparison.\n');
