%% capacitanceSweep300.m
% 300-pt full-6D LHS pose sweep for mapping the capacitance
% matrix/elastance matrix to plate pose. Reuses the same collision-safe
% sampler as the force/energy sweeps. Smaller than the 2000-pt force
% sweep since the target output (effective capacitance between two
% terminals) is a single scalar with a much simpler expected
% pose-dependence, and because each point here is more expensive
% (nested terminal/source sweep instead of a single stationary solve).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

n = 300;
[rx, ry, rz, thx, thy, thz] = sampleTrapPoseLHS(n);

fprintf('%6s  %12s  %12s\n', 'var', 'min', 'max');
for v = {'rx','ry','rz','thx','thy','thz'}
    name = v{1};
    val = eval(name);
    fprintf('%6s  %12.4g  %12.4g\n', name, min(val), max(val));
end

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData  = {rx, ry, rz, thx, thy, thz};

outFile = fullfile(scriptDir, 'comsol_capacitanceSweep300.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);
fprintf('\nSaved %d parameters x %d combinations to %s\n', numel(paramNames), numel(rx), outFile);
fprintf('Run with the full elastance/Maxwell/mutual matrix export enabled at each pose point\n');
fprintf('(saving the full matrix as a safety check, even though only S_11,S_22,S_12 are needed).\n');
