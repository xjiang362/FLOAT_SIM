%% rxFineSweepEnergy.m
% Fine, uniformly-spaced rx sweep (2um step) at a handful of fixed rz
% slices, for computing Fx via the energy-derivative method
% (Fx = -dW/drx, central difference across the grid) instead of the
% Maxwell Stress Tensor surface integral -- MST was confirmed
% non-convergent under mesh refinement for Fx (sign flips), while the
% total electric field energy W converges cleanly. 2um step chosen for
% ~10x margin over the energy noise floor seen between normal/fine mesh
% at a single fixed geometry.
%
% ry, thx, thy, thz held at zero throughout.
%
% NOTE: in COMSOL, add a Global Evaluation (Derived Values) for the
% total electric field energy in the air domain (e.g. intop(es.We) or
% the equivalent energy variable exposed by your Electrostatics
% interface) so it's exported as an extra column alongside the
% existing Fx/Fy/Fz/capacitance/Taux/Tauy/Tauz probes. Use the SAME
% mesh resolution validated as converged in the mesh convergence check.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

x_max = 40e-6;      % [m] match rxRzGridSweep.m range
dx = 2e-6;           % [m] fixed 2um step, uniform (not warped -- need
                      % even spacing for clean finite differencing)
rxSweep = (-x_max:dx:x_max)';   % 41 points

z_max = 12e-6;       % [m] match rxRzGridSweep.m range
rzValues = [-1, -0.5, 0, 0.5, 1]' * z_max;   % 5 fixed slices

nRx = numel(rxSweep);
nRz = numel(rzValues);
rx = repmat(rxSweep, nRz, 1);            % [nRx*nRz x 1]
rz = repelem(rzValues, nRx);             % [nRx*nRz x 1]
zeroCol = zeros(nRx*nRz, 1);

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};

outFile = fullfile(scriptDir, 'comsol_rxFineSweepEnergy.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, ...
    {rx, zeroCol, rz, zeroCol, zeroCol, zeroCol});

fprintf('Saved fine rx sweep (%d pts: %d rx values x %d rz slices) to %s\n', ...
    nRx*nRz, nRx, nRz, outFile);
fprintf('rx range: [%.0f, %.0f] um, uniform step %.1fum\n', -x_max*1e6, x_max*1e6, dx*1e6);
fprintf('rz slices: %s um\n', mat2str(rzValues*1e6, 4));
fprintf('Row blocks: ');
for i = 1:nRz
    fprintf('rz=%.1fum: rows %d-%d  ', rzValues(i)*1e6, (i-1)*nRx+1, i*nRx);
end
fprintf('\n');
