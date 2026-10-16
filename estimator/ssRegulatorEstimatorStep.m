function [u, xhatNext] = ssRegulatorEstimatorStep(xhat, y, ss)
%SSREGULATORESTIMATORSTEP One ISR-period update of a discrete predictor-
%form observer-based compensator for a single channel (see
%doubleIntegratorSSDesign.m for how ss.K, ss.L were designed):
%
%   u          = -K*xhat[k]                                (regulator, from the CURRENT estimate)
%   xhat[k+1]  = Ad*xhat[k] + Bd*u + L*(y[k] - Cd*xhat[k])  (observer update, using u just computed)
%
%   This is the standard separation-principle structure: the regulator
%   never sees the true state or the raw measurement, only xhat; the
%   observer never sees the true state either, only the measurement y
%   and the control it itself just issued.

u = -ss.K * xhat;
xhatNext = ss.Ad * xhat + ss.Bd * u + ss.L * (y - ss.Cd * xhat);

end
