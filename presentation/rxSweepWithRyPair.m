%% rxSweepWithRyPair.m
% Extends the dedicated uniform rx-only sweep with a two-point ry
% "bracket" (ry=-1um and ry=+1um) at EVERY rx grid point, so Fy(rx) can
% be derived via the energy method: for each rx, a central difference
% in ry, Fy(rx) ~= [W(rx,+1um) - W(rx,-1um)] / (2um), using two real
% simulated energy points rather than a numerical gradient across a
% full ry sweep. rz=thx=thy=thz=0 throughout, same as the original
% rx-only sweep.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

nPtsRx = 50;
x_max = 45e-6;   % [m] match the established rx range
ry_bracket = 1e-6;   % [m] +-1um ry bracket for the central difference

rxBase = linspace(-x_max, x_max, nPtsRx)';

% Two rows per rx value: (rx, ry=-1um) and (rx, ry=+1um)
rx = repelem(rxBase, 2);
ry = repmat([-ry_bracket; ry_bracket], nPtsRx, 1);
zeroCol = zeros(2*nPtsRx, 1);

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData  = {rx, ry, zeroCol, zeroCol, zeroCol, zeroCol};

outFile = fullfile(scriptDir, 'comsol_sweep_rxOnly_ryPair.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);
fprintf('Saved %d parameters x %d combinations to %s\n', numel(paramNames), numel(rx), outFile);
fprintf('rx: %d uniformly-spaced values in [%.4g, %.4g] m, each paired with ry = -%.4g m and +%.4g m\n', ...
    nPtsRx, min(rxBase), max(rxBase), ry_bracket, ry_bracket);
fprintf('\nExport with Total electric energy (Global Variable Probe 11 or equivalent) --\n');
fprintf('MST Fy probe optional/for comparison, but the energy bracket is the point here.\n');
