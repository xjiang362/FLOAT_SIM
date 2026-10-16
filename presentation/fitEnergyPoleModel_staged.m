%% fitEnergyPoleModel_staged.m
% Staged version: the joint nonlinear fit (k,p,a,b together) kept
% getting stuck with a,b showing zero sensitivity regardless of
% starting value -- a parameter-identifiability problem (rz's dominant
% pole term absorbs the fitting capacity, leaving a,b collinear with
% k,p). Fix: fit k,p on rz alone first (2-param, well-conditioned),
% then fit the tilt correction as a separate LINEAR regression on the
% residual (already shown to be ~356x above the mesh-noise floor, so
% real signal, correlating with thx^2/thy^2 at r~0.18).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));

g0_z = 15e-6;
n = numel(rz);
thx_rad = thx * pi/180;
thy_rad = thy * pi/180;
thz_rad = thz * pi/180;

%% Stage 1: k, p on rz alone
modelFun1 = @(b,rz) b(1) ./ (g0_z - rz).^b(2);
mdl1 = fitnlm(rz, W, modelFun1, [1.5e-14, 1]);
k = mdl1.Coefficients.Estimate(1);
p = mdl1.Coefficients.Estimate(2);
fprintf('Stage 1 (rz only): k=%.4g  p=%.4g  R^2=%.4f\n', k, p, mdl1.Rsquared.Ordinary);

Wfit1 = k ./ (g0_z - rz).^p;
resid1 = W - Wfit1;

%% Stage 2: linear regression of residual on tilt/lateral terms
% Candidate predictors: thx^2, thy^2 (even, expected from symmetry),
% plus interaction with rz (tilt effect might scale with how close to
% the wall you already are), plus rx^2,ry^2,thz^2 as a check (expected
% negligible based on earlier near-zero correlation).
predictors = [thx_rad.^2, thy_rad.^2, thx_rad.^2./(g0_z-rz), thy_rad.^2./(g0_z-rz), ...
              rx.^2, ry.^2, thz_rad.^2];
predNames = {'thx^2','thy^2','thx^2/(g0-rz)','thy^2/(g0-rz)','rx^2','ry^2','thz^2'};

mdl2 = fitlm(predictors, resid1, 'VarNames', [predNames, {'resid'}]);
disp(mdl2);

%% 10-fold CV for the combined staged model
rng(0);
kfold = 10;
cvIdx = mod((0:n-1)', kfold) + 1;
cvIdx = cvIdx(randperm(n));
yCV = zeros(n,1);
for f = 1:kfold
    testMask = (cvIdx == f);
    trainMask = ~testMask;
    m1 = fitnlm(rz(trainMask), W(trainMask), modelFun1, [k, p]);
    kf = m1.Coefficients.Estimate(1); pf = m1.Coefficients.Estimate(2);
    Wfit1_train = kf ./ (g0_z - rz(trainMask)).^pf;
    resid1_train = W(trainMask) - Wfit1_train;
    m2 = fitlm(predictors(trainMask,:), resid1_train);

    Wfit1_test = kf ./ (g0_z - rz(testMask)).^pf;
    resid_pred_test = predict(m2, predictors(testMask,:));
    yCV(testMask) = Wfit1_test + resid_pred_test;
end
r2_cv = 1 - sum((W-yCV).^2) / sum((W-mean(W)).^2);
fprintf('\nStaged model 10-fold CV R^2: %.4f\n', r2_cv);
fprintf('(Stage 1 alone CV R^2 was 0.9495 -- compare to see how much the tilt correction adds)\n');

save(fullfile(scriptDir, 'stagedEnergyModel.mat'), 'k', 'p', 'g0_z', 'mdl2', 'r2_cv');
