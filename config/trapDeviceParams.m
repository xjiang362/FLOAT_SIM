function params = trapDeviceParams()
%TRAPDEVICEPARAMS Physical properties, actuator geometry, and control
%gains for the free-flying 6-DOF "trap device":
%   - 4 corner thrusters, each firing along the BODY z-axis only
%     -> actively controls heave (z) and roll/pitch torque; cannot
%        produce yaw torque (r x [0;0;T] always has zero z-component).
%   - x, y, and yaw are passively "trapped" by a linear restoring
%     force/torque F = -k*x with k = c*u, where u is a controller
%     output that sets the trap stiffness (not a direct force/torque).

% --- Geometry: uniform rectangular plate, body frame centered on its
% centroid (x = length, y = width, z = thickness) ---
params.geom.L   = 5950e-6;  % [m] plate length (x) -- matches trapDeviceGeometry.m
params.geom.W   = 5950e-6;  % [m] plate width  (y) -- was 6e-3, a rounded value;
params.geom.t   = 350e-6;   % [m] plate thickness (z)
params.geom.rho = 2330;     % [kg/m^3] density (silicon)

% --- Mass properties, derived from geometry (uniform rectangular plate,
% inertia about the centroid): I_a = m/12*(sum of the squares of the two
% dimensions perpendicular to axis a) ---
params.m = params.geom.rho * params.geom.L * params.geom.W * params.geom.t;   % [kg]

Ixx = params.m / 12 * (params.geom.W^2 + params.geom.t^2);
Iyy = params.m / 12 * (params.geom.L^2 + params.geom.t^2);
Izz = params.m / 12 * (params.geom.L^2 + params.geom.W^2);
params.I = diag([Ixx, Iyy, Izz]);
params.Iinv = inv(params.I);

params.g = 9.81;                                % [m/s^2] world +z is up

% --- Thruster geometry: 4 corners of the plate itself (body frame, z=0) ---
L = params.geom.L / 2;   % [m] half-length (x)
W = params.geom.W / 2;   % [m] half-width  (y)
rx = [ L;  L; -L; -L];   % corner x-offsets, thrusters 1..4
ry = [ W; -W; -W;  W];   % corner y-offsets, thrusters 1..4
params.mixer.rx = rx;
params.mixer.ry = ry;

% Allocation matrix: [Fz; taux; tauy] = B * [T1;T2;T3;T4]
%   Fz   row: sum of thrusts
%   taux row: sum(ry_i * T_i)      (r x [0;0;T] -> x-component = ry*T)
%   tauy row: -sum(rx_i * T_i)     (r x [0;0;T] -> y-component = -rx*T)
B = [ones(1,4); ry'; -rx'];
params.mixer.B    = B;
params.mixer.Binv = pinv(B);   % minimum-norm allocation (redundant actuation)

% --- Electrostatic quadrant-electrode actuation (real physics, from
% FLOATS_SIM/presentation) -- Vmax is the nominal drive voltage used
% throughout that model-fitting work, not yet a real actuator spec ---
params.electrostatic.Vmax = 30;   % [V]

params.thrust.Tmin = 0;     % [N] thrusters are unilateral (push only)

% PLACEHOLDER pending a real actuator spec: sized off the mass this
% device actually has (hover needs m*g/4 per thruster), with a healthy
% margin for maneuvering rather than an arbitrary round number. Replace
% thrustMarginFactor / Tmax directly once real actuator numbers exist.
thrustMarginFactor = 5;
params.thrust.Tmax = thrustMarginFactor * (params.m * params.g / 4);   % [N] per-thruster limit

% --- Corner height sensing: z_i = z + ry_i*roll - rx_i*pitch (small angle) ---
M = [ones(4,1), ry, -rx];      % 4x3, maps [z;roll;pitch] -> 4 corner heights
params.cornerSensor.M    = M;
params.cornerSensor.Minv = pinv(M);   % least-squares inverse (yaw is unobservable, not included)

% --- ISR sample period for the discrete controller/estimator (must
% divide the fine integration dt evenly -- see simulateTrapDevice.m) ---
params.control.Tisr = 0.01;   % [s] 100 Hz -- typical embedded control loop rate

% --- Discrete state-space regulator + observer design for heave, roll,
% pitch (each a SISO double integrator -- see doubleIntegratorSSDesign.m).
% lambdar/lambdae are desired closed-loop BANDWIDTHS [rad/s] for a
% repeated real continuous pole (critically damped): the regulator
% places its poles at exp(-lambdar*Tisr), the observer at
% exp(-lambdae*Tisr). The observer is designed faster than the regulator
% (lambdae > lambdar, ~5x here) -- standard practice so the separation
% principle holds up well even though this decoupled-axis approximation
% isn't an exact linear model of the real coupled rigid-body dynamics.
params.lambdar = 3;    % [rad/s] regulator bandwidth
params.lambdae = 15;   % [rad/s] estimator bandwidth

params.ssHeave = doubleIntegratorSSDesign(1 / params.m,    params.control.Tisr, params.lambdar, params.lambdae);
params.ssRoll  = doubleIntegratorSSDesign(1 / params.I(1,1), params.control.Tisr, params.lambdar, params.lambdae);
params.ssPitch = doubleIntegratorSSDesign(1 / params.I(2,2), params.control.Tisr, params.lambdar, params.lambdae);

% --- Passive trap: k = c*u, bang-bang stiffness scheduling for damping ---
% PLACEHOLDER pending a real trap stiffness spec: c_pos/c_yaw are sized
% so the trap's natural frequency at u=1 is wn_trap, i.e. k = m*wn^2 (or
% Izz*wn^2 for yaw) -- the same "pick a bandwidth, derive the constant"
% approach used for lambdar/lambdae above, rather than a value carried
% over from a completely different mass scale. Replace wn_trap directly
% once a real trap stiffness range exists.
wn_trap = 10;   % [rad/s] nominal trap natural frequency at u=1
params.trap.c_pos = params.m * wn_trap^2;      % [N/m per unit u]     x, y stiffness scale
params.trap.c_yaw = Izz * wn_trap^2;           % [N m/rad per unit u] yaw stiffness scale
params.trap.u_lo  = 0.3;    % soft phase (swinging away from center)
params.trap.u_hi  = 2.0;    % stiff phase (swinging toward center)

% Yaw has no sensor (4 corner z-readings can't observe it -- see
% cornerHeightEstimator.m), so its trap can't run the bang-bang damping
% schedule; it's just held at a fixed constant stiffness (open loop).
params.trap.u_yaw_const = 0.5 * (params.trap.u_lo + params.trap.u_hi);

end

