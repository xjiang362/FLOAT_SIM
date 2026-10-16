%% overnightRxRzFineSweep.m
% ~2000pt overnight sweep for the Fx energy-derivative method. Keeps
% the already-validated fine rx grid (2um step, +-40um, 41 pts -- gave
% clean low-order LOO fits for the well-behaved rz slices in
% comsol_rxFineSweepEnergy_results.txt) but massively increases rz
% resolution (5 -> 49 slices, warped/denser near rz=0) so the
% rz-dependence of the fitted Fx(rx) shape can actually be
% characterized instead of eyeballed across 5 points. Total: 49*41 =
% 2009 points.
%
% ry, thx, thy, thz held at zero throughout.
%
% Same requirement as before: export W (total electric field energy,
% Global Evaluation) as an extra column alongside the usual probes, at
% the mesh resolution validated as converged.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

x_max = 40e-6;      % [m] match rxFineSweepEnergy.m range
dx = 2e-6;           % [m] uniform 2um step (validated: clean smoothed fits)
rxSweep = (-x_max:dx:x_max)';   % 41 points

z_max = 12e-6;       % [m] match rxRzGridSweep.m range (meshing-safe)
nRz = 49;
p = 2;               % warp exponent, denser near rz=0
u = linspace(-1, 1, nRz)';
rzValues = sign(u) .* abs(u).^p * z_max;

nRx = numel(rxSweep);
rx = repmat(rxSweep, nRz, 1);            % [nRx*nRz x 1]
rz = repelem(rzValues, nRx);             % [nRx*nRz x 1]
zeroCol = zeros(nRx*nRz, 1);

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};

outFile = fullfile(scriptDir, 'comsol_overnightRxRzFineSweep.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, ...
    {rx, zeroCol, rz, zeroCol, zeroCol, zeroCol});

fprintf('Saved overnight rx-rz fine sweep (%d pts: %d rx values x %d rz slices) to %s\n', ...
    nRx*nRz, nRx, nRz, outFile);
fprintf('rx range: [%.0f, %.0f] um, uniform step %.1fum\n', -x_max*1e6, x_max*1e6, dx*1e6);
fprintf('rz range: [%.3g, %.3g] um, warped (dense near 0), %d slices\n', min(rzValues)*1e6, max(rzValues)*1e6, nRz);
fprintf('At ~50min/205pts from tonight''s run, expect roughly %.1f hours.\n', (nRx*nRz)/205*50/60);
