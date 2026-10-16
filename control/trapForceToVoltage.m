function [u, info] = trapForceToVoltage(plant, f, poseEst, varargin)
%TRAPFORCETOVOLTAGE Invert the V^2 actuator, cancelling its stiffness.
%
%   [u, info] = trapForceToVoltage(plant, f, poseEst)
%
% This is the step that lets the controller be designed on a plant with no
% electrostatic stiffness. It is feedback linearisation, and it is worth
% being explicit about why it works and when it stops working.
%
% THE PROBLEM. At fixed voltage the force on the plate is, per pair,
%       F_p(z,V) = b_p(z) * V^2,        b_p(z) = (1/2) dC_p/dz
% and b_p grows as the gap closes. So moving toward the electrode
% increases the attraction: dF/dz > 0, a NEGATIVE stiffness. Linearised
% at hover that is -37.76 N/m in heave, giving the open-loop unstable pole
% at 1144 rad/s. A controller commanding voltage must fight it.
%
% THE CANCELLATION. Command FORCE instead, and solve for the voltage that
% produces it at the CURRENT pose:
%       V_p = sqrt( F_p,cmd / b_p(zhat) )
% Then F_p,actual = b_p(z) * F_p,cmd / b_p(zhat) = F_p,cmd whenever
% zhat = z. The position dependence divides out, and the plant seen by the
% outer loop is a pure double integrator. That is the A matrix the
% regulator is designed on.
%
% WHEN IT STOPS WORKING. The cancellation is only as good as zhat. With a
% stale estimate, zhat ~ z - zdot*td, and
%       F_actual = F_cmd * b(z)/b(z-zdot*td)
%                ~ F_cmd * [1 + (dln b/dz) * zdot * td]
% The residual is proportional to VELOCITY, not position, so a delayed
% inversion does not leave a stiffness behind -- it leaves ANTI-DAMPING.
% With b ~ 1/g^2, dln b/dz = 2/g, so the effective damping picks up
%       dc = -F_cmd * (2/g) * td
% which at hover (F = mg = 283 uN, g = 15 um, td = 1130 us) is -0.043
% N*s/m, five orders of magnitude larger than the squeeze-film damping and
% with the wrong sign. That is the mechanism by which sensing delay
% destabilises this loop, and why the delay budget is not negotiable.
%
% INPUTS
%   plant    struct from trapDeviceGeometry()
%   f        [N] 5-by-1 commanded force per pair, in the controller's
%            order [fp; f1; f2; f3; f4] = [Cp; quad1..quad4]
%   poseEst  6-by-1 pose estimate used for the inversion. Pass the LIVE
%            estimate; passing a constant reinstates the negative stiffness.
%
% NAME/VALUE
%   'Vtia'   [V] fixed bias on the TIA pair (default 40)
%   'Vmax'   [V] clamp (default 2*hover, generous)
%
% OUTPUTS
%   u        [V] 6-by-1 voltages in plant.idx.actPair order
%            [quad1..quad4, Cp, TIA]
%   info     .Fbias   [N] 5-by-1 gravity offset distributed by area
%            .Ftia    [N] force the TIA bias already supplies
%            .Ftot    [N] 5-by-1 total commanded force
%            .b       [N/V^2] 5-by-1 force coefficient at poseEst
%            .clipped logical, true if any channel hit a limit
%
% GRAVITY OFFSET. The controller is linearised about hover, so its output
% is the force INCREMENT. The weight still has to be carried: the TIA pair
% contributes whatever its fixed bias gives, and the remainder is split
% across the other five pairs in proportion to electrode area. Because
% F = b*V^2 and b is itself proportional to area, an area-weighted split
% puts the SAME voltage on every pair, which is the natural bias point.
%
% See also TRAPCONTROLDESIGN, TRAPWRENCH.

ip = inputParser;
ip.addParameter('Vtia', 40);
ip.addParameter('Vmax', 80);
ip.addParameter('clampFrac', 0.9);   % shrink factor per step, used ONLY on
                                     % a pose that is already invalid
ip.parse(varargin{:});
Vtia = ip.Results.Vtia;
Vmax = ip.Results.Vmax;
clampFrac = ip.Results.clampFrac;

f = f(:);
if numel(f) ~= 5
    error('trapForceToVoltage:size','f must be 5-by-1 [fp; f1; f2; f3; f4].');
end

% ---- keep the ESTIMATE inside the model domain ----
% poseEst comes from the estimator, not from the plant, so it can leave the
% physical envelope during a transient even while the plate itself is well
% clear of the walls. Evaluating the capacitance there would ask for a
% non-positive gap and throw, and a real controller cannot crash on its own
% state estimate: it saturates.
%
% This must leave any VALID pose bit-for-bit untouched. An earlier version
% clamped each axis to a fixed fraction of its travel, which silently moved
% legitimate poses -- at rz = -14 um it inverted at -13.5 um instead, and
% since b ~ 1/g^2 that alone put a 3.4% error on every delivered force.
% So: shrink toward the centred pose ONLY when the pose is actually outside
% the envelope, and otherwise do nothing at all.
poseEst = poseEst(:);
if ~checkPoseValid(plant, poseEst)
    for iShrink = 1:60
        poseEst = clampFrac * poseEst;
        if checkPoseValid(plant, poseEst), break; end
    end
end

% ---- force coefficients at the ESTIMATED pose ----
[~, dCdq] = trapCapacitance(plant, poseEst);
[~,~,B]   = trapWrench(plant, dCdq, zeros(6,1));
bAll = B(3,:)';                       % [N/V^2] per pair, in actPair order
idxQ  = 1:4;  idxCp = 5;  idxTia = 6;

% controller order is [Cp, quad1..4]; map to actPair order
map = [idxCp, idxQ];
b   = bAll(map);                      % 5-by-1 matching f

% ---- gravity offset ----
W    = plant.mass.m * plant.phys.g;
Ftia = bAll(idxTia) * Vtia^2;
Wrem = W - Ftia;
A5   = arrayfun(@(i) plant.term(plant.idx.actPair(i,1)).A + ...
                     plant.term(plant.idx.actPair(i,2)).A, map)';
Fbias = Wrem * A5/sum(A5);

% ---- total, then invert ----
Ftot = Fbias + f;
clipped = false;
Ftot(Ftot < 0) = 0;                   % V^2 cannot be negative
if any(Fbias + f < 0), clipped = true; end

V5 = sqrt(Ftot ./ b);
if any(V5 > Vmax), clipped = true; V5 = min(V5, Vmax); end

u = zeros(6,1);
u(map)    = V5;
u(idxTia) = Vtia;

if nargout > 1
    info.Fbias = Fbias; info.Ftia = Ftia; info.Ftot = Ftot;
    info.b = b; info.clipped = clipped; info.map = map;
    info.poseEstUsed = poseEst;
end
end
