function [tau_xy, xhatRollNext, xhatPitchNext] = tiltController(xhatRollPrev, xhatPitchPrev, rollMeas, pitchMeas, params)
%TILTCONTROLLER One ISR-period state-space regulator+observer update for
%roll and pitch (two independent SISO double-integrator channels, both
%regulated to zero -- level flight). xhatRoll/xhatPitch are [angle; rate]
%estimates, both now entirely sourced from params.ssRoll/params.ssPitch
%(see doubleIntegratorSSDesign.m) driven off the corner-sensor angle
%measurement -- no true-state rate is read here anymore.

[tauxCmd, xhatRollNext]  = ssRegulatorEstimatorStep(xhatRollPrev,  rollMeas,  params.ssRoll);
[tauyCmd, xhatPitchNext] = ssRegulatorEstimatorStep(xhatPitchPrev, pitchMeas, params.ssPitch);

tau_xy = [tauxCmd; tauyCmd];

end
