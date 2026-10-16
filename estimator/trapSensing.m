function [Ceff, info] = trapSensing(plant, CiR)
%TRAPSENSING Sensed capacitance of each quadrant readout pair.
%
%   [Ceff, info] = trapSensing(plant, CiR)
%
% Sensing stage of the levitation pipeline: takes the per-terminal
% capacitances produced by the capacitance stage and returns what the
% readout actually measures across each differential sensing pair.
%
% INPUTS
%   plant   geometry/calibration struct from trapDeviceGeometry()
%   CiR     [F] 20-by-1 terminal-to-resonator capacitances, FRINGING
%           CORRECTED (i.e. k_i*eps0*A_i/g_i, not the raw eps0*A_i/g_i).
%           Only the eight sensing terminals are read; the rest are
%           ignored and may be anything.
%
% OUTPUTS
%   Ceff    [F] 4-by-1, one per quadrant, Ceff(q) across pair q
%   info    struct with the intermediate branches, for diagnostics:
%           .Ca,.Cb [F] 4-by-2 each electrode's to-ground capacitance
%           .series [F] 4-by-1 the through-ground path
%           .pair   [-] 4-by-2 the terminal indices used
%
% MODEL
% Ports i and j are the two sensing electrodes of one quadrant; the
% resonator floats; every other terminal is grounded. Eliminating the
% floating resonator (the Y-delta / Schur step) and combining with the
% direct branches gives
%
%   Ceff = C12 + (Ca*Cb)/(Ca+Cb),   Ca = Cg + Cf + C_iR,  Cb = Cg + Cf + C_jR
%
%   C12  direct coplanar coupling between the two electrodes
%   Cg   each electrode's own direct capacitance to ground
%   Cf   per-electrode correction for the fringing the capacitance stage
%        does not model
%
% Cf sits INSIDE the series, alongside Cg, because it belongs to one
% electrode's branch. For matched electrodes the series halves it, which
% is why the equivalent closed form reads
%     Ceff = C12 + Cg/2 + Cf/2 + C_iR/2
%
% IMPORTANT -- Cf is calibrated against the UN-FRINGED eps0*A/g that
% trapCapacitance returns. If a fringing correction is ever reinstated
% upstream, Cf collapses from 31.8 fF to about 0.3 fF and must be
% re-extracted; using both at once double-counts fringing and costs ~4%.
%
% The through-ground series term dominates (~90% of Ceff) because the
% resonator is so strongly tied to ground -- C_rg/C_iR ~ 66 -- that it
% behaves as an AC ground rather than as a floating transfer node. That
% is also why the resonator-to-frame capacitance, despite being 41 pF and
% swinging 324% over rz, perturbs Ceff by only ~2%.
%
% ACCURACY, against the exact 22-terminal Schur reduction on the 63-pose
% dedicated sweep (poses with |rz| <= ~10 um, where the extracted matrix
% is clean): median 0.04%, RMS 1.76%, max 5.8%, R^2 = 0.990. The worst
% case sits at the largest gaps. Validated by verifyTrapSensing.m.
%
% LIMITATIONS
%   * No rx/ry dependence beyond whatever CiR carries. The sensed signal
%     does vary ~2% laterally, mostly through plate-electrode overlap
%     truncation, which the capacitance stage must supply.
%   * Constants are shared across all four quadrants by symmetry. The
%     eight sensing electrodes' fitted gains agree to 0.4%, so this is
%     well supported, but it has not been checked pair-by-pair against
%     a measurement of all four Ceff values.
%   * Calibrated at 30 V from one COMSOL mesh; re-extract if either
%     changes.
%
% See also TRAPDEVICEGEOMETRY, VERIFYTRAPSENSING.

pair = plant.idx.sensePair;             % 4-by-2 terminal indices
nPair = size(pair,1);

CiR = CiR(:);
if numel(CiR) < max(pair(:))
    error('trapSensing:shortInput', ...
        'CiR has %d entries; sensing pairs reference terminal %d.', ...
        numel(CiR), max(pair(:)));
end

Cg  = plant.ceff.Cg;
C12 = plant.ceff.C12;
Cf  = plant.ceff.Cf;

Ca = Cg + Cf + CiR(pair(:,1));
Cb = Cg + Cf + CiR(pair(:,2));
series = (Ca.*Cb)./(Ca + Cb);
Ceff = C12 + series;

if nargout > 1
    info.Ca     = Ca;
    info.Cb     = Cb;
    info.series = series;
    info.pair   = pair;
    info.names  = {plant.term(pair(:,1)).name; plant.term(pair(:,2)).name};
end

Ceff = reshape(Ceff, nPair, 1);
end
