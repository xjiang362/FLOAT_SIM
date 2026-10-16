%% fitRz2CrossTerms.m
% Reverts to the correct single top-electrode pole (no bottom
% electrode -- confirmed by the user; the earlier two-pole R^2 bump was
% just curve-fitting flexibility from a free extra term, not real
% physics). Directly tests the hypothesis that thx, thy have
% individual/bare terms in W not properly scaled by rz (or rz^2),
% causing them to dominate incorrectly at extreme (very negative) rz --
% by FORCING thx*rz^2, thy*rz^2 into the residual regression alongside
% the previously stepwise-selected terms, rather than leaving stepwise
% free to drop them.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));

g0_z = 15e-6;
n = numel(rz);
thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;
gap = g0_z - rz;

k = 1.601e-14; p = 1;   % single top-pole, as originally fit
Wfit1 = k ./ gap.^p;
resid1 = W - Wfit1;

% rz^2 in SI units is ~1e-10 (rz~1e-5), so thx*rz^2 columns end up
% ~1e-13 -- numerically indistinguishable from zero next to the other
% predictors, causing exact rank deficiency. Normalize rz by its range
% (z_max=14um) before squaring so these terms are O(1)-scaled.
z_max = max(abs(rz));
rzn = rz / z_max;

% Previously stepwise-selected terms (without rz^2 cross terms), PLUS
% forced-in thx*rzn^2, thy*rzn^2
predictors = [thx_rad, thy_rad, thx_rad./gap, thy_rad./gap, ry./gap, thz_rad./gap, ...
              thx_rad.*rz, thy_rad.*rz, thx_rad.*rzn.^2, thy_rad.*rzn.^2, ...
              thx_rad.^2, thy_rad.^2, thz_rad.^2, rx.^2, ry.^2, ...
              thx_rad.^2./gap, thy_rad.^2./gap, rx.^2./gap, ry.^2./gap, thz_rad.^2./gap, ...
              thx_rad.*thy_rad];
predNames = {'thx','thy','thx_g','thy_g','ry_g','thz_g','thxrz','thyrz','thxrz2','thyrz2', ...
             'thx2','thy2','thz2','rx2','ry2','thx2_g','thy2_g','rx2_g','ry2_g','thz2_g','thxthy'};

mdl = fitlm(predictors, resid1, 'VarNames', [predNames, {'resid'}]);
disp(mdl);
fprintf('In-sample R^2: %.4f\n', mdl.Rsquared.Ordinary);

rng(0);
kfold = 10;
cvIdx = mod((0:n-1)', kfold) + 1;
cvIdx = cvIdx(randperm(n));
yCV = zeros(n,1);
for f = 1:kfold
    testMask = (cvIdx == f);
    trainMask = ~testMask;
    mf = fitlm(predictors(trainMask,:), resid1(trainMask));
    yCV(testMask) = predict(mf, predictors(testMask,:));
end
r2_cv_resid = 1 - sum((resid1-yCV).^2)/sum((resid1-mean(resid1)).^2);
yCV_W = Wfit1 + yCV;
r2_cv_W = 1 - sum((W-yCV_W).^2)/sum((W-mean(W)).^2);
fprintf('Residual 10-fold CV R^2: %.4f\n', r2_cv_resid);
fprintf('Combined 10-fold CV R^2 on W: %.4f\n', r2_cv_W);

c = mdl.Coefficients.Estimate;
fprintf('\nthxrz2 coeff: %.4g (p=%.4g)\n', c(9), mdl.Coefficients.pValue(9));
fprintf('thyrz2 coeff: %.4g (p=%.4g)\n', c(10), mdl.Coefficients.pValue(10));

save(fullfile(scriptDir, 'rz2CrossTermModel.mat'), 'k','p','g0_z','mdl','predNames');
