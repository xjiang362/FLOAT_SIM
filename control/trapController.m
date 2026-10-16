function u = trapController(x, params)
%TRAPCONTROLLER Stiffness-only "control" of the passive x/y/yaw trap.
%
%   The trap law is F = -k*x with k = c*u: u cannot apply a direct
%   force/torque, only reshape the spring constant. A constant u gives
%   an undamped spring (x, y, yaw would oscillate forever). Real damping
%   is obtained by modulating k over the cycle instead: at a switching
%   instant (x continuous, k jumps) the stored energy 1/2*k*x^2 jumps by
%   dE = 1/2*(k_new - k_old)*x^2. So:
%     - LOWER k right at the outer turning points (|x| large): x^2 is
%       large there, so dropping k discards a lot of energy.
%     - RAISE k back near the center (|x| ~ 0): x^2 is small there, so
%       stiffening costs almost nothing, ready to catch the next swing.
%   Net effect over a cycle: amplitude decays. In terms of x and xdot,
%   swinging away from center (x*xdot > 0) should stiffen, swinging back
%   toward center (x*xdot < 0) should soften -- the opposite assignment
%   was verified (analytically and by a 1-DOF numerical test) to pump
%   energy IN instead and blow up.
%
%   This is one reasonable damping strategy for a stiffness-only
%   actuator, not a unique/optimal one -- replace pickGain below with a
%   different schedule if your hardware calls for it.
%
%   Yaw has no sensor (the 4 corner z-readings can't observe it -- see
%   cornerHeightEstimator.m), so it can't run this bang-bang schedule at
%   all: there's no feedback to switch on. It's held at a fixed constant
%   stiffness instead, which means an excited yaw offset will oscillate
%   undamped rather than settle out. That's an expected consequence of
%   having no yaw feedback, not a bug.

r = x(1:3);
v = x(4:6);

u = zeros(3, 1);
u(1) = pickGain(r(1), v(1), params.trap.u_lo, params.trap.u_hi);
u(2) = pickGain(r(2), v(2), params.trap.u_lo, params.trap.u_hi);
u(3) = params.trap.u_yaw_const;

end

function ui = pickGain(pos, rate, u_lo, u_hi)
if pos * rate > 0
    ui = u_hi;   % swinging away from center: stiffen (catch it early)
else
    ui = u_lo;   % swinging toward center (or at rest): soften (coast in)
end
end
