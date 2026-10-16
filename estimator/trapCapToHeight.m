function [zc, pose3, info] = trapCapToHeight(plant, Ceff)
%TRAPCAPTOHEIGHT Local plate height at each sensing pair, from capacitance.
%
%   [zc, pose3, info] = trapCapToHeight(plant, Ceff)
%
% Inverts the sensing model. Each quadrant readout pair measures a single
% effective capacitance, and that capacitance is a one-to-one function of
% the plate height directly above that pair, so it can be solved in closed
% form with no iteration.
%
% INPUTS
%   plant  struct from trapDeviceGeometry()
%   Ceff   [F] 4-by-1 sensed capacitances, one per quadrant pair, in the
%          order given by plant.idx.sensePair
%
% OUTPUTS
%   zc     [m] 4-by-1 plate height above nominal at each pair's sensing
%          point. Positive zc CLOSES the gap (the plate has risen toward
%          the electrodes), matching the sign convention used everywhere
%          else in this model.
%   pose3  [m; rad; rad] 3-by-1 least-squares [rz; thx; thy] fitted to the
%          four heights. Optional; see below.
%   info   .point   [m] 4-by-2 sensing-point coordinates (x,y)
%          .Cgap    [F] 4-by-1 the gap-dependent part of each measurement
%          .gap     [m] 4-by-1 local gap, g0 - zc
%          .resid   [m] 4-by-1 fit residual, zc - model(pose3)
%          .cond    [-] condition number of the 4-by-3 fit
%
% FORWARD MODEL BEING INVERTED
%
%   Ceff = C12 + Cg/2 + Cf/2 + (eps0/2) * dX*dY / (g0 - zc)
%
% The first three terms are pose-independent, so collecting them as
% C0 = C12 + Cg/2 + Cf/2 leaves a single hyperbola in zc, inverted by
%
%   zc = g0 - (eps0 * dX * dY) / (2 * (Ceff - C0))
%
% SENSING POINT. zc is the height at the geometric centre of the pair, the
% midpoint of the two triangular electrodes' area centroids, which lands
% at (+/-2250, +/-2250) um. That is not the centroid of either electrode
% (those sit at 2503.3 and 1996.7 um); it is where their contributions
% balance, which is why a single scalar suffices per pair.
%
% RECOVERING THE POSE. The plate surface is planar to the small-angle
% approximation this device operates in, so
%       zc(x,y) = rz + thx*y - thy*x
% Four measurements at four known points overdetermine the three unknowns,
% and pose3 is the least-squares solution. The remaining three degrees of
% freedom (rx, ry, thz) are NOT observable from these four capacitances:
% they are the passively trapped axes, and the sensing model has no
% lateral dependence to invert.
%
% CONDITIONING. The four points sit at the corners of a square, so the fit
% is well conditioned and the residual is a direct check on the planar
% assumption: a residual much above the capacitance noise floor means
% either a bad channel or a pose outside the model's validity.
%
% VALIDITY. Ceff must exceed C0, since the gap term is strictly positive.
% Approaching C0 from above drives zc to minus infinity, so a measurement
% at or below C0 is rejected rather than returned as a wild height.
%
% See also TRAPSENSING, TRAPCAPACITANCE, TRAPSAMPLERSTEP.

Ceff = Ceff(:);
pair = plant.idx.sensePair;
nPair = size(pair,1);
if numel(Ceff) ~= nPair
    error('trapCapToHeight:size','Ceff must be %d-by-1; got %d.', nPair, numel(Ceff));
end

% ---- constants of the forward model ----
e1   = plant.term(pair(1,1));
A    = e1.dX * e1.dY;                       % = footprint area, by the dX convention
C0   = plant.ceff.C12 + plant.ceff.Cg/2 + plant.ceff.Cf/2;
half = plant.phys.eps0 * A / 2;
g0   = plant.gap.z0;

% ---- invert, one pair at a time ----
Cgap = Ceff - C0;
if any(Cgap <= 0)
    bad = find(Cgap <= 0);
    error('trapCapToHeight:belowFloor', ...
        ['Pair %s measured at or below the pose-independent floor C0 = %.3f fF. ', ...
         'The gap term cannot be negative, so no height exists. ', ...
         'Check the channel or the calibration constants.'], ...
        mat2str(bad'), C0*1e15);
end
zc = g0 - half./Cgap;

% ---- sensing points, from the geometry rather than hard-coded ----
pt = zeros(nPair,2);
for q = 1:nPair
    a = plant.term(pair(q,1));  b = plant.term(pair(q,2));
    pt(q,:) = [(a.Xc+b.Xc)/2, (a.Yc+b.Yc)/2];
end

% ---- least-squares plane: zc = rz + thx*y - thy*x ----
if nargout > 1
    M = [ones(nPair,1), pt(:,2), -pt(:,1)];
    pose3 = M \ zc;
end

if nargout > 2
    info.point = pt;
    info.Cgap  = Cgap;
    info.gap   = g0 - zc;
    M = [ones(nPair,1), pt(:,2), -pt(:,1)];
    info.resid = zc - M*(M\zc);
    info.cond  = cond(M);
end
end
