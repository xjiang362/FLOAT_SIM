%% rzOnlySweep10_gndCheck.m
% Quick 10-pt rz-only sweep (all other DOF = 0) to check the resonator's
% capacitance-to-ground behavior after updating the ground plane
% location -- uniform spacing across the full rz range for a clean
% look at the trend.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

n = 10;
z_max = 14e-6;   % [m] match the established rz range
rz = linspace(-z_max, z_max, n)';
zeroCol = zeros(n,1);

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData  = {zeroCol, zeroCol, rz, zeroCol, zeroCol, zeroCol};

outFile = fullfile(scriptDir, 'comsol_rzOnlySweep10_gndCheck.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);
fprintf('Saved %d parameters x %d combinations to %s\n', numel(paramNames), n, outFile);
fprintf('rz values [um]: %s\n', mat2str(rz'*1e6, 4));
