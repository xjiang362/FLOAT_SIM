function est = cornerHeightEstimator(zCorners, params)
%CORNERHEIGHTESTIMATOR Recovers [z; roll; pitch] from 4 corner height
%measurements via the small-angle geometric model
%   z_i = z + ry_i*roll - rx_i*pitch
%(same rx_i, ry_i corner geometry as the thruster mixer), solved by
%least squares since 4 measurements over-determine 3 unknowns.
%
%   Yaw is NOT recoverable from z-only corner sensing: a pure rotation
%   about the vertical axis doesn't change any corner's height, so it
%   never appears in the measurement model above. Callers must get yaw
%   from elsewhere (or do without it).
%
%   The small-angle model is only accurate near level flight; the
%   tilt controller is what's expected to keep roll/pitch small enough
%   for this to hold.

sol = params.cornerSensor.Minv * zCorners;   % [z; roll; pitch]

est.z     = sol(1);
est.roll  = sol(2);
est.pitch = sol(3);

end
