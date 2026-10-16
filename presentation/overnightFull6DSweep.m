%% overnightFull6DSweep.m
% ~2000pt overnight sweep across the FULL 6-DOF space (rx,ry,rz,
% thx,thy,thz), replacing the earlier rx-rz-only plan. Instead of
% per-axis finite-differencing (only practical on dense isolated
% line-scans), the intent here is to fit ONE smooth global energy
% surrogate W(rx,ry,rz,thx,thy,thz) via regression over all sampled
% points, then differentiate that surrogate analytically to get
% Fx,Fy,Fz,Taux,Tauy,Tauz all at once, self-consistently -- this never
% touches the Maxwell Stress Tensor probes at all (confirmed
% non-convergent for Fx under mesh refinement), and avoids the
% mesh-rebuild point-to-point noise that dominated naive finite
% differencing on the fine rx-rz grid, since it's a global least-squares
% fit over thousands of points rather than a difference between two
% adjacent runs.
%
% Uses the existing collision-safe LHS sampler (sampleTrapPoseLHS.m).
%
% NOTE: in COMSOL, keep exporting the total electric field energy
% Global Evaluation (W) as an extra column alongside the usual
% Fx/Fy/Fz/capacitance/Taux/Tauy/Tauz probes, at the mesh resolution
% validated as converged.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

n = 2000;
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

outFile = fullfile(scriptDir, 'comsol_overnightFull6DSweep.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);
fprintf('\nSaved %d parameters x %d combinations to %s\n', numel(paramNames), numel(rx), outFile);
fprintf('At ~50min/205pts from the rx-rz fine sweep, expect roughly %.1f hours.\n', n/205*50/60);
fprintf('(Note: full 6D points solve the same geometry-rebuild pipeline as before, so the\n');
fprintf('per-point cost should be comparable -- this replaces, not adds to, the earlier\n');
fprintf('rx-rz-only overnight file given the ~8hr budget.)\n');
