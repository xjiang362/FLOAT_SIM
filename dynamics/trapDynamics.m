function [xdot, diag] = trapDynamics(plant, x, u)
%TRAPDYNAMICS Rigid-body dynamics of the levitated resonator.
%
%   [xdot, diag] = trapDynamics(plant, x, u)
%
% Closes the levitation pipeline: given the state and the drive voltages
% it evaluates the capacitance, the electrostatic wrench and gravity, and
% returns the state derivative for an ODE integrator.
%
% STATE (12-by-1)
%   x(1:6)   pose  [rx; ry; rz; thx; thy; thz]   metres and RADIANS
%   x(7:12)  rates [vx; vy; vz; wx; wy; wz]      m/s and rad/s
%
% INPUTS
%   plant  struct from trapDeviceGeometry()
%   u      [V] 6-by-1 drive voltages, one per actuation pair
%
% OUTPUTS
%   xdot   12-by-1 state derivative
%   diag   struct .F .tau .C .margin .minGap, for logging and for a
%          contact-aware integrator
%
% SMALL-ANGLE FORMULATION. Euler angles are used directly rather than
% quaternions because the device cannot rotate far: the collision
% constraint caps tilt at about 0.13 deg and yaw at 0.84 deg. Over that
% range Euler rates and body rates agree to parts in 1e4, so the usual
% kinematic transform is the identity to well within the model's own ~4%
% accuracy. This keeps the state in the SAME coordinates the capacitance
% model consumes, which removes a whole class of conversion bugs -- the
% legacy trapDeviceEOM.m carried a 13-state quaternion form and had to
% convert back to Euler angles on every call just to evaluate the
% electrode models.
%
% The gyroscopic term is retained even though it is small (order 0.7% of a
% typical torque at the fastest rates the unstable modes reach), because
% it costs one cross product.
%
% OPEN LOOP AND UNSTABLE. This plant is open-loop unstable in heave and
% tilt -- the electrostatic stiffness is negative, with poles near
% 1200-1300 rad/s at the nominal gap, rising above 1e4 rad/s as the gap
% closes. Integrating it without a stabilising controller will diverge to
% the wall; that is physics, not a bug. Use diag.margin to stop.
%
% See also TRAPCAPACITANCE, TRAPWRENCH, CHECKPOSEVALID.

x = x(:);
if numel(x) ~= 12
    error('trapDynamics:state','x must be 12-by-1 [pose(6); rates(6)].');
end
pose = x(1:6);
rate = x(7:12);

% ---- electrostatics ----
[C, dCdq, capInfo] = trapCapacitance(plant, pose);
[F, tau] = trapWrench(plant, dCdq, u);

% ---- gravity (world -z) ----
m = plant.mass.m;
F = F + [0; 0; -m*plant.phys.g];

% ---- gas damping ----
% Squeeze-film damping opposes the rate in every axis. It is small here
% (it shifts the unstable pole by one part in 4e5) but it is what finally
% settles the passively trapped lateral modes, which are otherwise
% undamped and ring indefinitely.
Qd  = plant.damp.C * rate;
F   = F   - Qd(1:3);
tau = tau - Qd(4:6);

% ---- rigid body ----
I = plant.mass.I;
w = rate(4:6);
acc    = F/m;
angAcc = I \ (tau - cross(w, I*w));

xdot = [rate; acc; angAcc];

% ---- clearance ----
[ok, margin, detail] = checkPoseValid(plant, pose);

if nargout > 1
    diag.F      = F;
    diag.tau    = tau;
    diag.C      = C;
    diag.ok     = ok;
    diag.margin = margin;
    diag.worst  = detail.worst;
    diag.minGap = capInfo.minGap;
end
end


