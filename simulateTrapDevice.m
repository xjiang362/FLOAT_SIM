function log = simulateTrapDevice(x0, setpoint, params, tSpan, dt)
%SIMULATETRAPDEVICE Two-rate (sampled-data) simulation of the trap
%device:
%   - The TRUE plant (trapDeviceEOM) is integrated continuously via
%     fixed-step RK4 at the fine step dt.
%   - The controller + discrete state-space estimator only run once per
%     ISR sample period params.control.Tisr (a multiple of dt), and their
%     output [T; uTrap] is held constant (zero-order hold) across every
%     fine step until the next tick.
%
%   log fields: t, x (13 x N, true state), T (4 x N thruster commands),
%   u (3 x N trap gains), xhat (6 x N logged [zHat;vzHat;rollHat;wxHat;
%   pitchHat;wyHat] estimate, held between ticks)

t0 = tSpan(1);
tf = tSpan(2);
tVec = t0:dt:tf;
N = numel(tVec);

Tisr = params.control.Tisr;
stepsPerTick = max(1, round(Tisr / dt));
if abs(Tisr - stepsPerTick * dt) > 1e-9
    warning('simulateTrapDevice:TisrNotMultipleOfDt', ...
        'params.control.Tisr is not an integer multiple of dt; using %d fine steps/tick (Tisr_eff = %.6g s).', ...
        stepsPerTick, stepsPerTick * dt);
end

x = x0(:);
x(7:10) = quatNormalize(x(7:10));

X    = zeros(13, N);
Th   = zeros(4, N);
Uh   = zeros(3, N);
Xhat = zeros(6, N);

% Bootstrap: initialize the observer's position states from the first
% instantaneous measurement, rate states at zero (unknown) -- the
% observer's own bandwidth is what converges these to the truth.
zCorners0 = cornerHeightSensor(x, params);
meas0 = cornerHeightEstimator(zCorners0, params);
if isa(setpoint, 'function_handle')
    sp0 = setpoint(t0);
else
    sp0 = setpoint;
end
xhat.heave = [meas0.z - sp0.z; 0];
xhat.roll  = [meas0.roll; 0];
xhat.pitch = [meas0.pitch; 0];

T = zeros(4, 1);
uTrap = zeros(3, 1);

for k = 1:N
    t = tVec(k);

    if k == 1 || mod(k - 1, stepsPerTick) == 0
        [T, uTrap, xhat] = trapDeviceController(t, x, xhat, params, setpoint);
    end

    X(:, k)  = x;
    Th(:, k) = T;
    Uh(:, k) = uTrap;
    Xhat(:, k) = [xhat.heave; xhat.roll; xhat.pitch];

    if k < N
        odeFun = @(tt, xx) trapDeviceEOM(tt, xx, params, T, uTrap);
        x = rk4Step(odeFun, t, x, dt);
        x(7:10) = quatNormalize(x(7:10));
    end
end

log.t    = tVec;
log.x    = X;
log.T    = Th;
log.u    = Uh;
log.xhat = Xhat;

end
