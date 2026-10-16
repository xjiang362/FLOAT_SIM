%% LHS_COMSOL_Sim_Points_Generation.m
% Generates a 200-point Latin Hypercube sample over the levitated
% plate's 6 rigid-body DOF (rx, ry, rz, thx, thy, thz) -- see
% sampleTrapPoseLHS.m for the sampling/constraint logic -- and exports
% it as a COMSOL parametric-sweep input file, plus a 3-point smoke-test
% subset for a quick pipeline check before committing to the full run.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

[rx, ry, rz, thx, thy, thz] = sampleTrapPoseLHS(200);

% Sanity-check ranges before exporting
fprintf('%6s  %12s  %12s\n', 'var', 'min', 'max');
for v = {'rx','ry','rz','thx','thy','thz'}
    name = v{1};
    val = eval(name);
    fprintf('%6s  %12.4g  %12.4g\n', name, min(val), max(val));
end

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData  = {rx, ry, rz, thx, thy, thz};

outFile = fullfile(scriptDir, 'comsol_lhs_sweep.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);
fprintf('Saved %d parameters x %d combinations to %s\n', numel(paramNames), numel(rx), outFile);

% Smoke-test subset: 3 combinations that exercise all 6 DOF at once
% (translation AND tilt together, not either alone), for a quick COMSOL
% run to confirm the parametric sweep mechanics/mesh hold up before
% committing to the full 200-point sweep. Selected by a combined
% "stress score" (each DOF normalized by its own max magnitude across
% all 200 rows, then summed) so the 3 picks span the space: the single
% most extreme combination, a middling/typical one, and the mildest one.
score = abs(rx)/max(abs(rx)) + abs(ry)/max(abs(ry)) + abs(rz)/max(abs(rz)) + ...
        abs(thx)/max(abs(thx)) + abs(thy)/max(abs(thy)) + abs(thz)/max(abs(thz));
[~, idxMax] = max(score);
[~, idxMin] = min(score);
[~, idxMed] = min(abs(score - median(score)));
testIdx = [idxMax, idxMed, idxMin];

fprintf('\nSmoke-test subset (row index, then rx ry rz [m] / thx thy thz [deg]):\n');
for k = 1:numel(testIdx)
    i = testIdx(k);
    fprintf('  #%d (row %3d): % .4g % .4g % .4g   % .4g % .4g % .4g\n', ...
        k, i, rx(i), ry(i), rz(i), thx(i), thy(i), thz(i));
end

testData = {rx(testIdx), ry(testIdx), rz(testIdx), thx(testIdx), thy(testIdx), thz(testIdx)};
testFile = fullfile(scriptDir, 'comsol_test_3pts.txt');
writeComsolParamFile(testFile, paramNames, paramUnits, testData);
fprintf('Saved 3-point smoke test to %s\n', testFile);
