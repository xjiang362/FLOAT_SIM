function [Fz_cmd, xhatNext] = heaveController(xhatPrev, zMeas, setpoint, params)
%HEAVECONTROLLER One ISR-period state-space regulator+observer update for
%heave. xhatPrev/xhatNext are [z - setpoint.z; vz - setpoint.vz]
%perturbation states (both position AND rate are now estimated by
%params.ssHeave -- see doubleIntegratorSSDesign.m -- not read from the
%true state at all).
%
%   Assumes a static setpoint (setpoint.vz is only used to define the
%   zero-perturbation reference, not tracked as a moving target); a
%   moving reference would need an explicit feedforward term, which is
%   out of scope here.

y = zMeas - setpoint.z;

[uPert, xhatNext] = ssRegulatorEstimatorStep(xhatPrev, y, params.ssHeave);

Fz_cmd = params.m * params.g + uPert;   % gravity feedforward + regulator perturbation

end
