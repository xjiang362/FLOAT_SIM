function [ok, margin, detail] = checkPoseValid(plant, pose)
%CHECKPOSEVALID Is the resonator still clear of the channel walls?
%
%   [ok, margin, detail] = checkPoseValid(plant, pose)
%
% The levitated mass floats inside a closed channel, and the constraint is
% NOT a box on the pose variables: translation and tilt consume the same
% clearance, so they must be checked together at the plate CORNERS.
%
% INPUTS
%   plant   struct from trapDeviceGeometry()
%   pose    6-by-1 [rx; ry; rz; thx; thy; thz], metres and RADIANS
%
% OUTPUTS
%   ok      logical, true if every corner is clear of every wall
%   margin  [m] smallest clearance over all corners and all walls.
%           Negative means penetration. Returning the margin rather than
%           just a flag lets an integrator back off BEFORE contact instead
%           of detecting it afterwards.
%   detail  struct .gapTop .gapBot [m] 4-by-1 vertical clearances per
%                  corner (top = toward the electrodes, bottom = floor)
%           .gapX .gapY [m] 4-by-1 lateral clearances
%           .corner [m] 4-by-2 corner positions after the pose transform
%           .worst  char, which constraint is tightest
%
% VERTICAL. The local gap under the plate is the same expression the
% capacitance model uses,
%       g(x,y) = g0 - rz - thx*y + thy*x
% and the plate sits centred in a channel of total height 2*g0, so the
% clearance above is g and below is 2*g0 - g. Evaluating at the four
% corners is sufficient because g is linear in x and y, so its extrema
% over a rectangle always occur at corners.
%
% LATERAL. Yaw swings the corners outward, so the corners are rotated by
% thz and translated by (rx,ry), then checked against walls at
% halfL + xy0 and halfW + xy0. A pure box test on rx, ry would miss this:
% at the yaw limit a corner reaches further than |rx| alone suggests.
%
% See also TRAPDEVICEGEOMETRY, TRAPDYNAMICS.

pose = pose(:);
if numel(pose) ~= 6
    error('checkPoseValid:pose','pose must be 6-by-1 [rx;ry;rz;thx;thy;thz].');
end
rx = pose(1); ry = pose(2); rz = pose(3);
thx = pose(4); thy = pose(5); thz = pose(6);

a  = plant.mass.halfL;
b  = plant.mass.halfW;
g0 = plant.gap.z0;
xy0 = plant.gap.xy0;

% ---- plate corners in the body frame ----
cx = [ a;  a; -a; -a];
cy = [ b; -b; -b;  b];

% ---- vertical clearance at each corner ----
% g is linear in (x,y), so corner evaluation bounds it everywhere.
gTop = g0 - rz - thx*cy + thy*cx;        % clearance to the electrode plane
gBot = 2*g0 - gTop;                      % clearance to the floor

% ---- lateral clearance, corners rotated by yaw then translated ----
c = cos(thz); s = sin(thz);
X = rx + cx*c - cy*s;
Y = ry + cx*s + cy*c;
wallX = a + xy0;
wallY = b + xy0;
gX = wallX - abs(X);
gY = wallY - abs(Y);

all4 = [gTop; gBot; gX; gY];
margin = min(all4);
ok = margin > 0;

if nargout > 2
    detail.gapTop = gTop;
    detail.gapBot = gBot;
    detail.gapX   = gX;
    detail.gapY   = gY;
    detail.corner = [X, Y];
    names = [repmat({'top'},4,1); repmat({'bottom'},4,1); ...
             repmat({'x-wall'},4,1); repmat({'y-wall'},4,1)];
    [~,iw] = min(all4);
    detail.worst = names{iw};
end
end
