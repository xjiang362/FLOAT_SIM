%% rxOnlySweep_uniform.m
% Dedicated rx-only sweep (rest of DOF held at exactly zero), UNIFORMLY
% spaced (linspace), to diagnose the Fx model mismatch found via
% trapElectrode1Model_coupledWL (flat, wrong-sign prediction vs noisy
% MST data) without paying for the full 300-pt 6-axis sweep.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

nPts = 50;
x_max = 45e-6;   % [m] match the established rx range

rx = linspace(-x_max, x_max, nPts)';
zeroCol = zeros(nPts,1);

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData  = {rx, zeroCol, zeroCol, zeroCol, zeroCol, zeroCol};

outFile = fullfile(scriptDir, 'comsol_sweep_rxOnly_uniform.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);
fprintf('Saved %d parameters x %d combinations to %s\n', numel(paramNames), nPts, outFile);
fprintf('rx range: [%.4g, %.4g] m, uniformly spaced\n', min(rx), max(rx));
fprintf('\nExport with the MST force/torque probes (1,2,3,5,6,7) and, if possible,\n');
fprintf('Global Variable Probe 11 (Total electric energy) for the energy-method check.\n');
