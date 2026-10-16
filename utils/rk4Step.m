function xNext = rk4Step(odeFun, t, x, dt)
%RK4STEP One fixed-step classical Runge-Kutta (RK4) integration step.
%   xNext = rk4Step(odeFun, t, x, dt) advances the state x by dt using
%   odeFun(t, x) -> xdot.

k1 = odeFun(t,        x);
k2 = odeFun(t + dt/2, x + dt/2 * k1);
k3 = odeFun(t + dt/2, x + dt/2 * k2);
k4 = odeFun(t + dt,   x + dt   * k3);

xNext = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

end
