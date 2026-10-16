%% fitCm12S1S2_models.m
% Fits closed-form models for Cm12 (direct terminal 1<->2 coupling), S1,
% and S2 (Maxwell self-capacitance of terminals 1 and 2) as functions of
% rz, thx, thy ONLY -- the last three quantities needed to make the
% Ceff(terminals 1,2) model fully closed-form (no simulation lookups at
% all). Uses the dedicated single-axis sweep (sensingSweep60), where
% rz, thx, thy are each varied independently with everything else held
% at exactly zero -- so unlike the scattered LHS fits, thx/thy
% dependence is cleanly identifiable here, and a LINEAR term is included
% (not just quadratic), since the force/torque and S_R models on this
% same kind of dedicated sweep turned out to have real linear tilt
% dependence that a quadratic-only model missed.
%
% Model (same form for all three quantities):
%   Q(rz,thx,thy) = C0 + A1/(d0-rz) + b1*thx + b2*thx^2 + b3*thy + b4*thy^2
% (thx,thy in RADIANS). d0=15um fixed (established nominal gap).
%
% Only a SINGLE pole (1/(d0-rz)) is used, not the two-pole form used for
% S_R: the raw data shows S1,S2,Cm12 vary monotonically across the full
% rz range with no sign of the second divergence near rz=-14um that
% S_R showed (that U-shape was specific to the resonator's coupling to
% the bottom Si ground plane, which terminals 1/2 don't exhibit here).
% A first attempt WITH both poles hit severe collinearity (RCOND~1e-23):
% since rz=0 throughout the thx-only and thy-only blocks, 1/(d0-rz) and
% 1/(d0+rz) collapse to the SAME constant (1/d0) on 40 of the 60 points,
% leaving only the 20-point rz-only block to separate C0, A1, A2 -- too
% little independent information for 3 parameters. Dropping A2 removes
% the collinearity entirely and matches the observed monotonic shape.
%
% This model is LINEAR in all 6 parameters (C0,A1,b1,b2,b3,b4), so it is
% solved by ordinary linear least squares (backslash) -- no iterative
% optimizer, no risk of unidentifiable/rank-deficient parameters.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'sensingSweep60_allMatrices.mat'));   % Cm_all, rx,ry,rz,thx,thy,thz, nPose

thx_rad = thx*pi/180; thy_rad = thy*pi/180;
d0 = 15e-6;

Cm12 = squeeze(Cm_all(1,2,:));
S1 = squeeze(sum(Cm_all(1,:,:),2));
S2 = squeeze(sum(Cm_all(2,:,:),2));

designMat = [ones(nPose,1), 1./(d0-rz), thx_rad, thx_rad.^2, thy_rad, thy_rad.^2];
names = {'C0','A1','b1(thx)','b2(thx2)','b3(thy)','b4(thy2)'};

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

quantities = {Cm12, 'Cm12'; S1, 'S1'; S2, 'S2'};
results = struct();
for q = 1:3
    y = quantities{q,1};
    qName = quantities{q,2};

    % Column-scale the design matrix before solving: raw columns span
    % ~10 orders of magnitude (pole term ~1e4-1e5 vs thx^2 ~1e-6), which
    % triggers spurious "near singular" warnings from backslash even
    % though the underlying regression is well-posed. Normalize each
    % column to unit norm, solve, then rescale coefficients back.
    colScale = sqrt(sum(designMat.^2, 1));
    designMatScaled = designMat ./ colScale;
    coeffsScaled = designMatScaled \ y;
    coeffs = (coeffsScaled ./ colScale')';
    coeffs = coeffs(:);
    yhat = designMat * coeffs;
    r2 = r2fun(y, yhat);

    % Leave-one-out CV (closed-form via hat matrix, exact for linear regression)
    H = designMatScaled * ((designMatScaled'*designMatScaled) \ designMatScaled');
    h = diag(H);
    looResid = (y - yhat) ./ (1 - h);
    r2_loo = 1 - sum(looResid.^2)/sum((y-mean(y)).^2);

    fprintf('--- %s ---\n', qName);
    for i = 1:numel(names)
        fprintf('  %s = %.6e\n', names{i}, coeffs(i));
    end
    fprintf('  In-sample R^2 = %.6f, LOO R^2 = %.6f\n\n', r2, r2_loo);

    results.(qName).coeffs = coeffs;
    results.(qName).r2 = r2;
    results.(qName).r2_loo = r2_loo;
end

modelFun = @(coeffs, X) coeffs(1) + coeffs(2)./(d0-X(:,1)) + ...
    coeffs(3)*X(:,2) + coeffs(4)*X(:,2).^2 + coeffs(5)*X(:,3) + coeffs(6)*X(:,3).^2;
% X columns: [rz, thx_rad, thy_rad]

save(fullfile(scriptDir, 'Cm12_S1_S2_models.mat'), 'results', 'modelFun', 'd0');
fprintf('Saved to Cm12_S1_S2_models.mat\n');
