function [rx, ry, rz, thx, thy, thz] = sampleTrapPoseLHS(n)
%SAMPLETRAPPOSELHS Latin Hypercube sample of n rigid-body poses
%(rx,ry,rz,thx,thy,thz) for the levitated plate, geometrically
%constrained so no sampled pose collides with the trap walls/floor.
%
%   Units: rx, ry, rz in meters; thx, thy, thz in degrees.
%
%   All 6 DOF are sampled independently over their own full symmetric
%   range (positive or negative). Both translation and rotation are now
%   applied as geometry-level transforms (rebuild + remesh per
%   combination) rather than via a deforming domain, so bounds here
%   reflect actual physical/collision limits rather than moving-mesh
%   robustness -- hence the much larger tilt_disp_factor (0.9) than the
%   earlier ALE-constrained version (0.2) used.
%
%   Calling this with a different n produces an independent LHS design
%   (not a subset of any other call) -- a subset of an LHS loses its
%   per-axis stratification/space-filling guarantee, so always call
%   this fresh for each batch size you need rather than truncating a
%   larger design.

% --- Nominal geometry ---
g0_z  = 15e-6;      % [m] nominal vertical gap (top or bottom) when rz = 0,
                     % i.e. the plate sits centered in a 2*g0_z-tall channel
g0_xy = 50e-6;       % [m] nominal lateral gap
L = 6000e-6;         % [m] platform length, x (update to your actual dimension)
W = L;               % [m] platform width, y (square plate)
tilt_disp_factor = 0.9;   % max corner displacement induced by tilt alone,
                           % as a fraction of whatever gap remains after
                           % translation -- raised from 0.2 now that
                           % rotation is also geometry-based (no ALE
                           % robustness margin needed, just real collision
                           % clearance)

% 1. LHS Probabilities
rng('default');
P = lhsdesign(n, 6, 'criterion', 'maximin');

% 2. Calculate Translations First (constrained by lateral gap)
% 45 um lateral limit -- translation is now a geometry-level transform
% (rebuild + remesh per combination), so this reflects real collision
% clearance against the 50 um lateral gap, not moving-mesh robustness.
x_max = 45e-6;                                 % [m]
rx = (P(:,1) * 2 - 1) * x_max;                 % rx spans [-x_max, +x_max]
ry = (P(:,2) * 2 - 1) * x_max;                 % ry spans [-x_max, +x_max]

z_max = 14e-6;                                 % [m] 1 um margin below g0_z=15um -- similarly now
                                                % a real collision clearance, not an ALE limit
rz = (P(:,3) * 2 - 1) * z_max;                 % rz spans [-z_max, +z_max]

% 3. Vertical Constraint (Roll/Pitch)
% thx (roll, about x) moves a corner vertically by +/-(W/2)*sin(thx) and
% thy (pitch, about y) by +/-(L/2)*sin(thy); since thx and thy are
% sampled independently, the worst case has BOTH pushing the same
% corner the same way at once, so their combined vertical excursion,
% ON TOP OF rz, is capped so the TOTAL vertical usage never exceeds
% tilt_disp_factor*g0_z. Capping the *remainder* after rz (i.e.
% tilt_disp_factor*(g0_z-|rz|)) instead would let total usage drift
% toward 100% of g0_z as |rz| grows -- e.g. at 90%, |rz|+0.9*(g0_z-|rz|)
% = 0.9*g0_z + 0.1*|rz|, which creeps toward g0_z for large |rz| instead
% of holding at a flat 90% -- so the total budget is capped directly,
% then whatever rz already used is subtracted off before splitting the
% rest between thx/thy's lever arms.
z_budget_total = tilt_disp_factor * g0_z;
z_extra = max(z_budget_total - abs(rz), 0);
th_xy_max = asind( z_extra / (L/2 + W/2) );

thx = (P(:,4) * 2 - 1) .* th_xy_max;   % thx spans [-th_xy_max, +th_xy_max]
thy = (P(:,5) * 2 - 1) .* th_xy_max;   % thy spans [-th_xy_max, +th_xy_max]

% 4. Lateral Constraint (Yaw)
% thz (yaw, about z) swings a corner laterally by ~sin(thz) times the
% plate's in-plane lever arm, in BOTH x and y at once. Same fix as
% above: cap TOTAL lateral usage (translation + tilt swing) at
% tilt_disp_factor*g0_xy directly, rather than capping the remainder
% after rx/ry -- the remainder-based version let combined usage creep
% toward 100% of g0_xy for large rx/ry, which is exactly what produced
% a real COMSOL meshing failure at rx=23.8um, ry=29.6um, thz=-0.34deg
% (that combination landed a corner within ~3um of the outer wall).
xy_budget_total = tilt_disp_factor * g0_xy;
x_extra = max(xy_budget_total - abs(rx), 0);
y_extra = max(xy_budget_total - abs(ry), 0);
thz_max = asind( min(x_extra, y_extra) / max(L/2, W/2) );

thz = (P(:,6) * 2 - 1) .* thz_max;   % thz spans [-thz_max, +thz_max]

end
