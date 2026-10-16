%% dedicatedCenterElectrodeSweeps.m
% Dedicated single-variable characterization sweep for the CENTRAL
% electrode (unlike quadrants 1-4, this one is symmetric about both
% x=0 and y=0, so rx, ry, thz are expected to have no meaningful
% effect -- confirmed by the same reasoning already validated for the
% sensing model, not modeled here). Only rz, thx, thy are swept, one at
% a time, 50 pts each (150 total) -- reuses the exact pattern already
% proven for the quadrant-1 Fz/Taux/Tauy fit and the sensing model.
%
% Fx, Fy may still show up nonzero vs thx/thy in the exported data --
% expected, since COMSOL's force probes report GLOBAL-frame components
% and the sweep tilts the plate geometry itself (not a co-rotating
% measurement frame), so a symmetric electrode's dominant Fz gets
% partly resolved onto the global x/y axes as the plate tilts. Not
% necessarily a separate physical coupling to fit independently.
%
% Lessons applied from the quadrant-1 fitting work (avoid repeating
% mistakes that cost extra sweep rounds there):
%   1. APPLY THE ELECTRODE VOLTAGE before running this sweep (30V,
%      matching the rest of this project) -- an earlier sweep run
%      without voltage had to be discarded entirely.
%   2. Export BOTH the MST force/torque probes AND the Total Electric
%      Energy probe from this SAME run -- needed a second sweep round
%      last time because energy wasn't included initially.
%   3. Translation range capped at what's safe -- thx,thy range here
%      uses the same collision-clearance formula as before; if the
%      center electrode's footprint (L,W) differs from the quadrant
%      electrodes' 6mm x 6mm, update L/W below before running.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

nPts = 50;
p = 2;   % warp exponent: higher = more concentrated near zero (matches established convention)

% --- Nominal geometry (must match sampleTrapPoseLHS.m unless the
% center electrode's own footprint differs -- UPDATE L,W if so) ---
g0_z  = 15e-6;
L = 6000e-6; W = L;
tilt_disp_factor = 0.9;

z_max = 14e-6;
th_xy_max_at_zero = asind((tilt_disp_factor * g0_z) / (L/2 + W/2));   % thx, thy at rz=0

u = linspace(-1, 1, nPts)';
warpedSpacing = @(xmax) sign(u) .* abs(u).^p * xmax;

rzSweep  = warpedSpacing(z_max);
thxSweep = warpedSpacing(th_xy_max_at_zero);
thySweep = warpedSpacing(th_xy_max_at_zero);

axisSweeps = {rzSweep, thxSweep, thySweep};
axisNames  = {'rz', 'thx', 'thy'};

fprintf('%6s  %14s  %14s\n', 'axis', 'min', 'max');
for i = 1:3
    fprintf('%6s  %14.4g  %14.4g\n', axisNames{i}, min(axisSweeps{i}), max(axisSweeps{i}));
end

%% Combined file: all 150 points in one sweep (single COMSOL run)
nAxes = 3;
combined = zeros(nAxes*nPts, 6);
colIdx = [3, 4, 5];   % rz,thx,thy positions in the 6-column (rx,ry,rz,thx,thy,thz) layout
for i = 1:nAxes
    rows = (i-1)*nPts + (1:nPts);
    combined(rows, colIdx(i)) = axisSweeps{i};
end

paramNames = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
paramUnits = {'m', 'm', 'm', 'deg', 'deg', 'deg'};
paramData = num2cell(combined, 1);

outFile = fullfile(scriptDir, sprintf('comsol_centerElectrodeSweep_rzThxThy%d.txt', nAxes*nPts));
writeComsolParamFile(outFile, paramNames, paramUnits, paramData);

fprintf('\nSaved combined center-electrode sweep (%d pts) to %s\n', nAxes*nPts, outFile);
fprintf('Row blocks: ');
for i = 1:nAxes
    fprintf('%s=%d-%d  ', axisNames{i}, (i-1)*nPts+1, i*nPts);
end
fprintf('\n\nMAKE SURE THE CENTER ELECTRODE VOLTAGE IS APPLIED before running --\n');
fprintf('export with BOTH the MST force/torque probes AND Total Electric Energy\n');
fprintf('(Global Variable Probe, same as the energy-method model) in this SAME run.\n');
