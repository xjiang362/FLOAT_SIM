function [Kst, info] = trapOpenLoopStiffness(plant, varargin)
%TRAPOPENLOOPSTIFFNESS Electrostatic stiffness at hover, at FIXED voltage.
%
%   [Kst, info] = trapOpenLoopStiffness(plant)
%   [Kst, info] = trapOpenLoopStiffness(plant, 'Vtia', 0)
%
% Returns the 6-by-6 matrix Kst = -d(wrench)/d(pose) evaluated at the hover
% pose with the hover voltages HELD CONSTANT. Sign convention: Kst is a
% stiffness, so a stable (restoring) degree of freedom has Kst > 0 and an
% unstable (anti-restoring) one has Kst < 0.
%
% WHY FIXED VOLTAGE. This is the plant an observer actually sees. The
% force-to-voltage inversion in trapForceToVoltage cancels this stiffness
% only where it is evaluated -- at the ESTIMATED pose. The estimation error
% therefore sees the stiffness uncancelled, so an observer designed on the
% stiffness-free double integrator is modelling the wrong plant. Use this
% matrix to build the model the observer is designed on.
%
% HEAVE HAS A CLOSED FORM. Each pair gives F = b*V^2 with b = (1/2)dC/dz
% and C = eps0*A/g, so b ~ 1/g^2 and
%       dF/dz = F * (2/g).
% Summed over the pairs at hover, where the total force carries the weight,
%       dFz/drz = W*(2/g0)      ->   Kst(3,3) = -2*W/g0
% which for W = m*g = 283.22 uN and g0 = 15 um is -37.76 N/m, the pole at
% w_u = sqrt(2*W/(g0*m)) = 1144 rad/s = 182 Hz. info.KzzClosedForm carries
% this value so the finite difference can be checked against it.
%
% INPUTS
%   plant   struct from trapDeviceGeometry()
%
% NAME/VALUE
%   'Vtia'  [V] fixed bias on the TIA pair (default 0, matching the
%           simulations run with the TIA grounded)
%   'pose'  6-by-1 pose to linearise about (default zeros, i.e. hover)
%   'h'     finite-difference step, [m m m rad rad rad] (default 1e-9/1e-9)
%
% OUTPUTS
%   Kst     [N/m, N/rad; N*m/m, N*m/rad] 6-by-6 stiffness
%   info    .u      [V] hover voltages held during the differentiation
%           .wrench [N;N*m] wrench at the linearisation pose
%           .KzzClosedForm  [N/m] -2*W/g0, for cross-checking Kst(3,3)
%           .poles  [rad/s] sqrt(-Kst_ii/m_ii) for the diagonal, signed
%
% See also TRAPFORCETOVOLTAGE, TRAPCONTROLDESIGN, TRAPWRENCH.

ip = inputParser;
ip.addParameter('Vtia', 0);
ip.addParameter('pose', zeros(6,1));
ip.addParameter('h', [1e-9 1e-9 1e-9 1e-9 1e-9 1e-9]);
ip.addParameter('afterInversion', false);
ip.parse(varargin{:});
Vtia = ip.Results.Vtia;
pose = ip.Results.pose(:);
h    = ip.Results.h(:);
afterInv = ip.Results.afterInversion;

% ---- hover voltages: zero commanded increment, so the bias carries W ----
u = trapForceToVoltage(plant, zeros(5,1), pose, 'Vtia', Vtia);

if afterInv
    % RESIDUAL stiffness: recompute the voltage at every probe pose, i.e.
    % let the inversion track the pose perfectly. Whatever stiffness
    % survives is what the inversion CANNOT cancel. For heave this is zero
    % (both members of a pair share a gap); for tilt it is not, because a
    % pair's two members sit at different coordinates and a tilt
    % redistributes force between them while their total z-force is held.
    wr = @(q) wrenchAt(plant, q, ...
        trapForceToVoltage(plant, zeros(5,1), q, 'Vtia', Vtia, 'Vmax', 1e9));
else
    wr = @(q) wrenchAt(plant, q, u);
end

Kst = zeros(6,6);
for j = 1:6
    qp = pose; qp(j) = qp(j) + h(j);
    qm = pose; qm(j) = qm(j) - h(j);
    Kst(:,j) = -(wr(qp) - wr(qm)) / (2*h(j));
end

if nargout > 1
    W = plant.mass.m * plant.phys.g;
    info.u = u;
    info.wrench = wr(pose);
    info.KzzClosedForm = -2*W/plant.gap.z0;
    mm = [plant.mass.m; plant.mass.m; plant.mass.m; ...
          plant.mass.I(1,1); plant.mass.I(2,2); plant.mass.I(3,3)];
    d  = diag(Kst);
    info.poles = sign(-d) .* sqrt(abs(d)./mm);
end
end

% =====================================================================
function w = wrenchAt(plant, q, u)
[~, dCdq] = trapCapacitance(plant, q);
[F, tau]  = trapWrench(plant, dCdq, u);
w = [F; tau];
end
