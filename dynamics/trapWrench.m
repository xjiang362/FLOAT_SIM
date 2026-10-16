function [F, tau, B] = trapWrench(plant, dCdq, u)
%TRAPWRENCH Electrostatic wrench on the resonator from the drive voltages.
%
%   [F, tau, B] = trapWrench(plant, dCdq, u)
%
% Force stage of the levitation pipeline. Every conductor is held at a
% fixed potential, so the generalised force is
%
%       Q_j = +dW/dq_j = 0.5 * sum_i V_i^2 * dC_i/dq_j
%
% with the PLUS sign (see the project's energy-method convention: at fixed
% potential the source does work, and F = +dW/dq, not -dW/dq).
%
% INPUTS
%   plant  struct from trapDeviceGeometry()
%   dCdq   [F/m, F/rad] 20-by-6 Jacobian from trapCapacitance()
%   u      [V] 6-by-1 drive voltages, one per actuation PAIR:
%            u(1) -> terminals  3, 4   quadrant 1
%            u(2) -> terminals  7, 8   quadrant 2
%            u(3) -> terminals 11,12   quadrant 3
%            u(4) -> terminals 15,16   quadrant 4
%            u(5) -> terminals 17,18   Cp
%            u(6) -> terminals 19,20   TIA
%
% OUTPUTS
%   F      [N]    3-by-1 force  [Fx; Fy; Fz]
%   tau    [N*m]  3-by-1 torque [taux; tauy; tauz]
%   B      [N/V^2, N*m/V^2] 6-by-6 allocation matrix, Q = B*(u.^2)
%
% DIFFERENTIAL DRIVE. Each pair is driven at +u and -u, so BOTH members
% see the same u^2 and contribute with the same sign. The proof mass
% therefore floats near zero potential while still being pulled by both
% sub-electrodes -- which is why a pair produces force at all rather than
% cancelling.
%
% THE ALLOCATION IS LINEAR IN u^2, NOT u. That is the reason the legacy
% thruster mixer cannot be reused: its allocation matrix assumes force
% linear in each actuator command. Working in v = u.^2 restores linearity,
% so B is a genuine allocation matrix and v = pinv(B)*Q_desired is a valid
% inverse -- subject to v >= 0, since u^2 cannot be negative. That
% unilateral constraint is the real difference from a bidirectional
% actuator and has to be handled by the allocator, not here.
%
% The eight SENSING electrodes are treated as grounded and contribute no
% force. If a sense excitation is added later it enters exactly the same
% way, as extra terms in the sum.
%
% Gravity is NOT included; it belongs to the dynamics stage.
%
% See also TRAPCAPACITANCE, TRAPDYNAMICS.

pair = plant.idx.actPair;              % 6-by-2 terminal indices
nPair = size(pair,1);

u = u(:);
if numel(u) ~= nPair
    error('trapWrench:badU','u must be %d-by-1 (one per actuation pair); got %d.', ...
        nPair, numel(u));
end
if size(dCdq,1) < max(pair(:)) || size(dCdq,2) ~= 6
    error('trapWrench:badJacobian', ...
        'dCdq must be at least %d-by-6; got %d-by-%d.', ...
        max(pair(:)), size(dCdq,1), size(dCdq,2));
end

% Allocation matrix: column p is the wrench produced by unit u(p)^2.
% Both members of the pair add, since each sees the same squared voltage.
B = zeros(6, nPair);
for p = 1:nPair
    B(:,p) = 0.5*( dCdq(pair(p,1),:) + dCdq(pair(p,2),:) )';
end

Q = B*(u.^2);
F   = Q(1:3);
tau = Q(4:6);
end
