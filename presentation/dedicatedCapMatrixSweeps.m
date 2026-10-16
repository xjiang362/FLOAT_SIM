%% dedicatedCapMatrixSweeps.m
% Dedicated SINGLE-AXIS sweeps for the 22-terminal capacitance matrix,
% aimed at the resonator <-> frame mutual capacitance Cm(21,22).
%
% WHY: in the 300-pose full-6D LHS (comsol_capacitanceSweep300_withGND_*)
% Cm(21,22) is by far the largest coupling to the frame -- mean 39.8 fF,
% ranging 12.5 to 281 fF -- but it correlates -0.72 with rz, so the gap
% dependence swamps everything else. rx shows only corr -0.098 (p=0.09)
% there, which is NOT evidence of absence: an LHS cannot separate a weak
% rx effect from a dominant rz one. Sweeping one axis at a time with the
% other five pinned at exactly zero resolves each dependence cleanly.
%
% Three axes, one block each, all other DOF exactly zero:
%   rx   +/- 45 um      lateral  -> tests overlap/edge coupling to the frame
%   rz   +/- 14 um      vertical -> the known dominant dependence
%   thx  +/- 0.1289 deg tilt     -> tests tilt sensitivity of the frame term
%
% SPACING: uniform, not the p=2 warped spacing used by the force sweeps.
% The warped spacing concentrates points near zero, which suits a
% known-smooth 1/g law; here the SHAPE is the unknown (the rx dependence
% may well have a one-sided corner at rx=0, as the lateral force does), so
% even coverage of the full range is more informative and also permits
% clean finite differencing. 21 points per axis rather than 20 so that
% exactly zero is included as a shared reference pose in every block --
% three more solves than a round 60, in exchange for a common reference
% and a free cross-block repeatability check.
%
% This 63-point single-axis set is intended to REPLACE a re-run of the
% 300-pose LHS for this question: one-at-a-time sweeps buy far more
% information per solve when the goal is to separate individual
% dependences, and each point here is expensive (nested terminal/source
% sweep rather than a single stationary solve).
%
% COMSOL SETUP -- this must be run on the 22-TERMINAL model with the
% bottom Si included (the same setup that produced
% comsol_capacitanceSweep300_withGND_term_*), using the nested
% terminal/source sweep, and exporting the MUTUAL capacitance matrix
% (Maxwell and inverse/elastance too, for the usual cross-validation).
% A 21-terminal export has no frame terminal and cannot answer this.
%
% SYMBOLS: rx, ry, rz [m] plate translations; thx, thy, thz [deg] plate
%   rotations; g0_z [m] nominal gap; L, W [m] plate length and width;
%   tilt_disp_factor [dimensionless] fraction of the gap budget allowed
%   to be consumed by tilt before collision.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

nPts = 21;                  % odd, so rx=rz=thx=0 is sampled exactly

% --- Nominal geometry (must match sampleTrapPoseLHS.m) ---
g0_z = 15e-6;
L = 6000e-6; W = L;
tilt_disp_factor = 0.9;

x_max  = 45e-6;
z_max  = 14e-6;
thx_max = asind((tilt_disp_factor * g0_z) / (L/2 + W/2));

u = linspace(-1, 1, nPts)';          % uniform
rxSweep  = u * x_max;
rzSweep  = u * z_max;
thxSweep = u * thx_max;

axisSweeps = {rxSweep, rzSweep, thxSweep};
axisNames  = {'rx', 'rz', 'thx'};
axisCols   = [1, 3, 4];              % position in the (rx,ry,rz,thx,thy,thz) layout
axisUnits  = {'um', 'um', 'deg'};
axisScale  = [1e6, 1e6, 1];

fprintf('%6s  %14s  %14s  %10s  %6s\n','axis','min','max','step','unit');
for i = 1:3
    s = axisSweeps{i}*axisScale(i);
    fprintf('%6s  %14.4g  %14.4g  %10.4g  %6s\n', axisNames{i}, min(s), max(s), s(2)-s(1), axisUnits{i});
end

%% Combined file: all three blocks in one COMSOL run
nAxes = 3;
combined = zeros(nAxes*nPts, 6);
for i = 1:nAxes
    rows = (i-1)*nPts + (1:nPts);
    combined(rows, axisCols(i)) = axisSweeps{i};
end

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData  = num2cell(combined, 1);

outFile = fullfile(scriptDir, sprintf('comsol_capMatrixSweep_rxRzThx%d.txt', nAxes*nPts));
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);

fprintf('\nSaved %d points (%d axes x %d pts) to %s\n', nAxes*nPts, nAxes, nPts, outFile);
fprintf('Row blocks: ');
for i = 1:nAxes
    fprintf('%s = %d-%d   ', axisNames{i}, (i-1)*nPts+1, i*nPts);
end
fprintf('\n');
fprintf('\nNote: rows %d, %d and %d are all the SAME pose (all DOF zero) --\n', ...
    (nPts+1)/2, nPts+(nPts+1)/2, 2*nPts+(nPts+1)/2);
fprintf('a free repeatability check on the solver across the three blocks.\n');

fprintf('\n--- BEFORE RUNNING ---\n');
fprintf('1. Use the 22-TERMINAL model (bottom Si included), the same setup as\n');
fprintf('   comsol_capacitanceSweep300_withGND_term_*. A 21-terminal export has\n');
fprintf('   no frame terminal and cannot answer the resonator<->frame question.\n');
fprintf('2. Export the MUTUAL capacitance matrix (plus Maxwell and inverse for\n');
fprintf('   cross-validation), 22 rows per pose, same column layout as the 300-pose\n');
fprintf('   file so parseCapSweep300_withGND.m can be reused with nPose=%d.\n', nAxes*nPts);
fprintf('3. Keep the mesh at the resolution already validated as converged.\n');
