function [u, info] = trapWrenchToVoltage(plant, uw, poseEst, varargin)
%TRAPWRENCHTOVOLTAGE Allocate a commanded WRENCH across the 5 drive pairs.
%
%   [u, info] = trapWrenchToVoltage(plant, uw, poseEst)
%
% The control input is the wrench increment
%       uw = [dFz; dtaux; dtauy]            [N; N*m; N*m]
% and the weight is added here, so the total wrench to deliver is
%       Q = [W + dFz; dtaux; dtauy].
%
% WHY THIS REPLACES THE PER-PAIR INVERSION. trapForceToVoltage commanded a
% force per PAIR and inverted each pair's own z-force coefficient,
%       v_p = Ftot_p / b_p,     v = V^2.
% That pins each pair's total z-force and nothing else. It is exact for
% heave, because both members of a pair share a gap. It is NOT exact for
% tilt: a pair's two members sit at different coordinates (quadrant 1
% drives (2250,750) AND (750,2250) off one voltage), so a tilt makes them
% see different gaps, force redistributes between them, and the effective
% arm drifts,
%       arm_p(th) = mean_p(Y) + (2*th/g0)*Var_p(Y).
% The surviving stiffness is -(2/g0)*sum_p F_p*Var_p(Y), which left only
% beta = 1 - K_res/K_tot = 0.800 (thx) and 0.724 (thy) of it cancelled.
%
% Allocating the WRENCH removes that entirely. Solving
%       B3(pose) * v = Q,       B3 = rows [Fz; taux; tauy] of the
%                               allocation at the CURRENT pose
% delivers all three components exactly, so beta = 1 in every degree of
% freedom and there is no residual stiffness to model. Measured:
% K_res = 2.7e-11 N/m, -0.0 and -1.2e-14 N*m/rad, i.e. machine zero.
%
% REDUNDANCY. B3 is 3-by-5: five pairs for three controlled directions, so
% the solution is a 2-parameter family. The minimum-norm solution
% v = pinv(B3)*Q is used, which spreads the demand over all five pairs and
% keeps the voltages as low as the geometry allows. v is then floored at
% zero, since V^2 cannot be negative -- the actuator is UNILATERAL, it can
% only pull toward the electrodes. info.clipped reports when that floor or
% the Vmax ceiling is active, which is the signal that the commanded
% wrench was not physically deliverable.
%
% INPUTS
%   plant    struct from trapDeviceGeometry()
%   uw       3-by-1 [dFz; dtaux; dtauy], the regulator output
%   poseEst  6-by-1 pose to evaluate the allocation at. Pass the LIVE
%            estimate; a frozen pose reinstates the negative stiffness.
%
% NAME/VALUE
%   'Vtia'   [V] fixed bias on the TIA pair (default 0)
%   'Vmax'   [V] clamp (default 80)
%   'clampFrac'  shrink factor applied ONLY to an invalid poseEst
%
% OUTPUTS
%   u        [V] 6-by-1 in plant.idx.actPair order [quad1..4, Cp, TIA]
%   info     .Q .v .B3 .clipped .poseEstUsed
%
% See also TRAPWRENCH, TRAPCONTROLDESIGN, TRAPFORCETOVOLTAGE.

ip = inputParser;
ip.addParameter('Vtia', 0);
ip.addParameter('Vmax', 80);
ip.addParameter('clampFrac', 0.9);
ip.parse(varargin{:});
Vtia = ip.Results.Vtia;
Vmax = ip.Results.Vmax;
clampFrac = ip.Results.clampFrac;

uw = uw(:);
if numel(uw) ~= 3
    error('trapWrenchToVoltage:size','uw must be 3-by-1 [dFz; dtaux; dtauy].');
end

% ---- keep the estimate inside the model domain (valid poses untouched) ----
poseEst = poseEst(:);
if ~checkPoseValid(plant, poseEst)
    for iShrink = 1:60
        poseEst = clampFrac * poseEst;
        if checkPoseValid(plant, poseEst), break; end
    end
end

% ---- allocation at the estimated pose ----
[~, dCdq] = trapCapacitance(plant, poseEst);
[~,~,B]   = trapWrench(plant, dCdq, zeros(6,1));
map = [5 1 2 3 4];                       % actPair rows: Cp, quad1..4
B3  = B([3 4 5], map);                   % 3-by-5, v -> [Fz; taux; tauy]

W = plant.mass.m * plant.phys.g;
Ftia = B(3, 6) * Vtia^2;                 % the TIA bias already carries this
Q = [W - Ftia + uw(1); uw(2); uw(3)];

v = pinv(B3) * Q;

clipped = false;
if any(v < 0), clipped = true; v = max(v, 0); end

V5 = sqrt(v);
if any(V5 > Vmax), clipped = true; V5 = min(V5, Vmax); end

u = zeros(6,1);
u(map)  = V5;
u(6)    = Vtia;

if nargout > 1
    info.Q = Q; info.v = v; info.B3 = B3;
    info.clipped = clipped; info.poseEstUsed = poseEst;
end
end
