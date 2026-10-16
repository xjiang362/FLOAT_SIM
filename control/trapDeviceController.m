function [T, uTrap, xhatNext] = trapDeviceController(t, xTrue, xhatPrev, params, setpoint)
%TRAPDEVICECONTROLLER One ISR-period update, combining:
%   - cornerHeightSensor + cornerHeightEstimator: instantaneous [z; roll;
%     pitch] measurement from the 4 raw corner readings.
%   - heaveController + tiltController: discrete state-space regulator +
%     observer (see estimator/doubleIntegratorSSDesign.m) for heave,
%     roll, pitch. xhatPrev/xhatNext are structs with fields heave, roll,
%     pitch, each a 2x1 [position; rate] estimate -- NEITHER position NOR
%     rate is read from the true state for these three channels anymore;
%     both come entirely from the observer. Their combined output is the
%     desired [Fz; taux; tauy] wrench, allocated to 4 thrusters by
%     thrusterMixer.
%   - trapController: passive trap stiffness command uTrap. x/y still
%     read the true state directly (not corner-sensed in this model);
%     yaw runs open loop (no sensor for it at all).
%
%   setpoint may be a struct (static hover target, fields z/vz) or a
%   function handle setpoint(t) -> struct for a time-varying reference.

if isa(setpoint, 'function_handle')
    sp = setpoint(t);
else
    sp = setpoint;
end

zCorners = cornerHeightSensor(xTrue, params);
meas = cornerHeightEstimator(zCorners, params);

[Fz_cmd, xhatNext.heave] = heaveController(xhatPrev.heave, meas.z, sp, params);
[tau_xy, xhatNext.roll, xhatNext.pitch] = tiltController( ...
    xhatPrev.roll, xhatPrev.pitch, meas.roll, meas.pitch, params);

wrench = [Fz_cmd; tau_xy];
T = thrusterMixer(wrench, params);
uTrap = trapController(xTrue, params);

end
