function T = thrusterMixer(wrench, params)
%THRUSTERMIXER Allocate a desired [Fz; taux; tauy] wrench across the 4
%corner thrusters via the precomputed minimum-norm pseudo-inverse, then
%saturate to each thruster's unilateral [Tmin, Tmax] range. Saturation
%can distort the realized wrench if the command is infeasible.

T = params.mixer.Binv * wrench;
T = min(max(T, params.thrust.Tmin), params.thrust.Tmax);

end
