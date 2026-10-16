%% rxRzGridSweep.m
% A 2D grid crossing rx and rz, to resolve Fx's true coupling with rz
% -- the single-axis rz-only sweep couldn't see this since it has
% rx=0 throughout (Fx near its small poly5 intercept there, a
% near-zero signal), and the real 6D LHS data is too sparse to isolate
% a clean rz trend at any fixed rx. This gives dense, clean slices:
% for each of 6 fixed rx values (nonzero, since rx=0 is already
% covered by existing single-axis data), a full rz sweep.
%
% ry, thx, thy, thz held at zero throughout.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

x_max = 40e-6;                                          % [m] dialed back from 45um
rxValues = [-1, -2/3, -1/3, 1/3, 2/3, 1] * x_max;       % [m], 6 fixed values
nRz = 15;
z_max = 12e-6;   % [m] dialed back from 14um -- meshing issues at the tighter gap
p = 2;   % warp exponent, same convention as before: concentrate near rz=0

u = linspace(-1, 1, nRz)';
rzSweep = sign(u) .* abs(u).^p * z_max;

nRx = numel(rxValues);
rx = repelem(rxValues(:), nRz);          % [nRx*nRz x 1]
rz = repmat(rzSweep, nRx, 1);            % [nRx*nRz x 1]
zeroCol = zeros(nRx*nRz, 1);

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};

outFile = fullfile(scriptDir, 'comsol_rxRzGrid.txt');
writeComsolParamFile(outFile, paramNames, paramUnits, ...
    {rx, zeroCol, rz, zeroCol, zeroCol, zeroCol});

fprintf('Saved rx-rz grid (%d pts: %d rx values x %d rz points) to %s\n', ...
    nRx*nRz, nRx, nRz, outFile);
fprintf('rx values: %s\n', mat2str(rxValues*1e6, 4));
fprintf('rz range: [%.3g, %.3g] um, warped (dense near 0)\n', min(rzSweep)*1e6, max(rzSweep)*1e6);
fprintf('Row blocks: ');
for i = 1:nRx
    fprintf('rx=%.0fum: rows %d-%d  ', rxValues(i)*1e6, (i-1)*nRz+1, i*nRz);
end
fprintf('\n');
