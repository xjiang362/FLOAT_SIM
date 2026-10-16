%% fitSR_physicsModel.m
% Fits S_R (resonator's full mutual row sum) to a 6-pole physics-
% informed model across all 6 DOF, replacing the raw-simulation-data
% lookup used in the "resonator-only analytical" Ceff model:
%
%   S_R = A1/(d0-rz+(d0+rz)*(k1*thx+k2*thy))
%       + A2/(d0+rz+(d0+rz)*(k3*thx+k4*thy))
%       + A3/(x0-rx) + A4/(x0+rx) + A5/(y0-ry) + A6/(y0+ry)
%
% x0=y0=50um, d0=15um fixed (matching the established nominal gaps).
% thx,thy in RADIANS. Denominators are floored to a small positive
% value during fitting to avoid complex contamination (same issue
% found with the fringing capacitance model fit).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'capSweep300_withGND_allMatrices.mat'));  % Cm_all, rx,ry,rz,thx,thy,thz, nPose

S_R = squeeze(sum(Cm_all(21,:,:), 2));

d0 = 15e-6; x0 = 50e-6; y0 = 50e-6;
floorVal = 1e-8;

X = [rx, ry, rz];

% k1-k4 (tilt correction) dropped -- not identifiable via nonlinear
% fit, negligible contribution beyond the A1-A6 pole terms. Added a
% constant C0 instead. params: [C0,A1,A2,A3,A4,A5,A6]
modelFun = @(p, X) p(1) + ...
    p(2) ./ max(d0 - X(:,3), floorVal) + ...
    p(3) ./ max(d0 + X(:,3), floorVal) + ...
    p(4) ./ max(x0 - X(:,1), floorVal) + ...
    p(5) ./ max(x0 + X(:,1), floorVal) + ...
    p(6) ./ max(y0 - X(:,2), floorVal) + ...
    p(7) ./ max(y0 + X(:,2), floorVal);

p0 = [0, 2e-16, 2e-16, 5e-16, 5e-16, 5e-16, 5e-16];

opts = statset('MaxIter', 2000, 'Display', 'off', 'TolFun', 1e-30, 'TolX', 1e-20);
[pfit, resid_fit] = nlinfit(X, S_R, modelFun, p0, opts);

fprintf('Fitted parameters:\n');
names = {'C0','A1','A2','A3','A4','A5','A6'};
for i = 1:numel(names)
    fprintf('  %s = %.6e\n', names{i}, pfit(i));
end

r2 = 1 - sum(resid_fit.^2)/sum((S_R-mean(S_R)).^2);
fprintf('\nIn-sample R^2: %.4f\n', r2);
fprintf('max|imag(pfit)| = %.4e (should be 0)\n', max(abs(imag(pfit))));

%% 10-fold CV
rng(0);
n = numel(S_R);
kfold = 10;
cvIdx = mod((0:n-1)', kfold) + 1;
cvIdx = cvIdx(randperm(n));
yCV = zeros(n,1);
for f = 1:kfold
    testMask = (cvIdx == f);
    trainMask = ~testMask;
    try
        pf = nlinfit(X(trainMask,:), S_R(trainMask), modelFun, pfit, opts);
        yCV(testMask) = modelFun(pf, X(testMask,:));
    catch
        yCV(testMask) = NaN;
    end
end
valid = ~isnan(yCV);
r2_cv = 1 - sum((S_R(valid)-yCV(valid)).^2)/sum((S_R(valid)-mean(S_R(valid))).^2);
fprintf('10-fold CV R^2: %.4f (%d/%d valid)\n', r2_cv, sum(valid), n);

save(fullfile(scriptDir, 'SR_physicsModel.mat'), 'pfit', 'modelFun', 'd0', 'x0', 'y0', 'r2', 'r2_cv');
