%% uniformRxSweep20.m
% Uniformly-spaced rx-only sweep: 20 points, linearly spaced from
% -45um to +45um, every other DOF held at exactly zero.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

nPts = 20;
rx = linspace(-45e-6, 45e-6, nPts)';
zero1 = zeros(nPts, 1);

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};

outFile = fullfile(scriptDir, 'comsol_sweep_rxUniform20.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, ...
    {rx, zero1, zero1, zero1, zero1, zero1});

fprintf('Saved uniform rx sweep (%d pts, rx in [%.3g, %.3g] m) to %s\n', ...
    nPts, min(rx), max(rx), outFile);
