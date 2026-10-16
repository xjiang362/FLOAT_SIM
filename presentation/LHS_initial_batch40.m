%% LHS_initial_batch40.m
% Bootstrap batch for adaptive GP sampling: a fresh, independent
% 40-point Latin Hypercube sample (NOT a subset of the 200-point
% comsol_lhs_sweep.txt design -- subsetting an LHS breaks its per-axis
% stratification, so this calls sampleTrapPoseLHS fresh at n=40).
%
% Run this batch in COMSOL first. Once you have the corresponding
% output data (e.g. force/capacitance per combination), a GP surrogate
% can be fit on it and used to propose the next, uncertainty-guided
% batch rather than continuing to draw more LHS points blindly.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

[rx, ry, rz, thx, thy, thz] = sampleTrapPoseLHS(40);

fprintf('%6s  %12s  %12s\n', 'var', 'min', 'max');
for v = {'rx','ry','rz','thx','thy','thz'}
    name = v{1};
    val = eval(name);
    fprintf('%6s  %12.4g  %12.4g\n', name, min(val), max(val));
end

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData  = {rx, ry, rz, thx, thy, thz};

outFile = fullfile(scriptDir, 'comsol_lhs_initial40.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);
fprintf('Saved %d parameters x %d combinations to %s\n', numel(paramNames), numel(rx), outFile);
