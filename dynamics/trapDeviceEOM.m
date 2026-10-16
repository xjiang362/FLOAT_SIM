function xdot = trapDeviceEOM(t, x, params, T, uTrap)
%TRAPDEVICEEOM 6-DOF equations of motion for the trap device, using the
%REAL nonlinear 4-quadrant electrostatic actuation model (fit from
%dedicated COMSOL energy-method sweeps -- see FLOATS_SIM/presentation)
%in place of the earlier placeholder (4 generic body-z thrusters + an
%idealized linear x/y/yaw "trap" spring).
%
%   State x (13x1): [r(3); v(3); q(4); w(3)]
%   T (4x1): held over from the placeholder actuator interface for
%     signature compatibility with simulateTrapDevice.m/
%     trapDeviceController.m (NOT yet reworked -- see note below).
%     Linearly mapped to quadrant electrode VOLTAGES here, in [0,Vmax].
%   uTrap (3x1): accepted but UNUSED. The old linear "trap" spring
%     (F=-k*x) is superseded by the real electrostatic model, which
%     provides its own (highly nonlinear) restoring force/torque -- there
%     is no separate trap-stiffness knob anymore.
%
%   INTERIM WIRING NOTE: T's linear thrust->voltage mapping and its
%   corner->quadrant assignment (see below) are a placeholder so this
%   file can be swapped without also rewriting the controller/mixer/
%   sensor stack in the same change. thrusterMixer.m's allocation matrix
%   assumes forces linear in each actuator, which is NOT true for
%   electrostatic force (~V^2, strongly pose-dependent) -- a real
%   controller for this actuator needs its own voltage-allocation logic,
%   not the linear thruster mixer. Revisit once that rework happens.
%
%   Quadrant electrode <-> corner correspondence (from
%   params.mixer.rx/ry vs. the mirror/rotation transforms in
%   trapElectrode{2,3,4}Model_v2.m): thruster 1 (+L,+W) -> electrode 1,
%   thruster 2 (+L,-W) -> electrode 4, thruster 3 (-L,-W) -> electrode 3,
%   thruster 4 (-L,+W) -> electrode 2.

r = x(1:3);
q = x(7:10);
w = x(11:13);

R = quat2rotmCustom(q);

% --- Pose in the electrostatic model's convention ---
% rx,ry,rz [m]: plate translation from the nominal trap center. Treated
% here as equal to the world-frame position r -- valid to the same
% small-angle/small-displacement approximation the whole fringing-model
% dataset was built under (translations of tens of um, tilts of a
% fraction of a degree), so body-frame and world-frame lateral
% displacement are not meaningfully different.
rx = r(1); ry = r(2); rz = r(3);

% thx,thy,thz [deg]: roll/pitch/yaw. quat2eulCustom returns
% [yaw;pitch;roll] in RADIANS (ZYX convention).
eul = quat2eulCustom(q);
thx_deg = eul(3) * 180/pi;   % roll
thy_deg = eul(2) * 180/pi;   % pitch
thz_deg = eul(1) * 180/pi;   % yaw

% --- Interim actuator command: T (4x1, [0,Tmax]) -> quadrant voltages ---
Vmax = params.electrostatic.Vmax;
Tmax = params.thrust.Tmax;
Vquad = zeros(4,1);
Vquad(1) = Vmax * T(1) / Tmax;   % electrode 1 <- thruster 1 (+L,+W)
Vquad(4) = Vmax * T(2) / Tmax;   % electrode 4 <- thruster 2 (+L,-W)
Vquad(3) = Vmax * T(3) / Tmax;   % electrode 3 <- thruster 3 (-L,-W)
Vquad(2) = Vmax * T(4) / Tmax;   % electrode 2 <- thruster 4 (-L,+W)

% --- Sum the 4 quadrant electrodes' electrostatic wrench (body frame) ---
[Fx1,Fy1,Fz1,Taux1,Tauy1,Tauz1] = trapElectrode1Model_v2(rx,ry,rz,thx_deg,thy_deg,thz_deg,Vquad(1));
[Fx2,Fy2,Fz2,Taux2,Tauy2,Tauz2] = trapElectrode2Model_v2(rx,ry,rz,thx_deg,thy_deg,thz_deg,Vquad(2));
[Fx3,Fy3,Fz3,Taux3,Tauy3,Tauz3] = trapElectrode3Model_v2(rx,ry,rz,thx_deg,thy_deg,thz_deg,Vquad(3));
[Fx4,Fy4,Fz4,Taux4,Tauy4,Tauz4] = trapElectrode4Model_v2(rx,ry,rz,thx_deg,thy_deg,thz_deg,Vquad(4));

F_es_body   = [Fx1+Fx2+Fx3+Fx4;  Fy1+Fy2+Fy3+Fy4;  Fz1+Fz2+Fz3+Fz4];
tau_es_body = [Taux1+Taux2+Taux3+Taux4;  Tauy1+Tauy2+Tauy3+Tauy4;  Tauz1+Tauz2+Tauz3+Tauz4];

% --- Assemble total wrench (electrostatic force rotated body->world,
% same convention the old thruster force used; gravity unchanged) ---
F_total_world  = R * F_es_body + [0; 0; -params.m * params.g];
tau_total_body = tau_es_body;

% --- Rigid body dynamics ---
rdot = x(4:6);
vdot = F_total_world / params.m;
qdot = 0.5 * quatMultiply(q, [0; w]);
wdot = params.Iinv * (tau_total_body - cross(w, params.I * w));

xdot = [rdot; vdot; qdot; wdot];

end
