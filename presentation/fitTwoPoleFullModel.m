%% fitTwoPoleFullModel.m
% Rebuilds the full W model with a TWO-POLE baseline (top AND bottom
% electrode, since the plate sits centered in a 2*g0_z channel) instead
% of the single top-pole used before -- confirmed via
% k1=1.57e-14, k2=1.15e-15 (R^2 0.9496->0.9602 on rz alone). Then
% refits the residual with the same broad candidate set as before, PLUS
% explicit thx*rz^2, thy*rz^2 terms per the hypothesis that individual
% (non-rz-scaled) thx/thy terms were dominating incorrectly at extreme
% negative rz where the bottom-pole term matters.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));

g0_z = 15e-6;
n = numel(rz);
thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;
gapTop = g0_z - rz;
gapBot = g0_z + rz;

k1 = 1.574e-14; k2 = 1.153e-15;  % from the confirmed two-pole fit
Wfit1 = k1./gapTop + k2./gapBot;
resid1 = W - Wfit1;
fprintf('Two-pole base R^2: %.4f\n', 1 - sum(resid1.^2)/sum((W-mean(W)).^2));

predictors = [thx_rad, thy_rad, thz_rad, rx, ry, ...
              thx_rad./gapTop, thy_rad./gapTop, rx./gapTop, ry./gapTop, thz_rad./gapTop, ...
              thx_rad.*rz, thy_rad.*rz, thx_rad.*rz.^2, thy_rad.*rz.^2, ...
              thx_rad.^2, thy_rad.^2, thz_rad.^2, rx.^2, ry.^2, ...
              thx_rad.^2./gapTop, thy_rad.^2./gapTop, rx.^2./gapTop, ry.^2./gapTop, thz_rad.^2./gapTop, ...
              rx.*ry, rx.*thz_rad, ry.*thz_rad, thx_rad.*thy_rad];
predNames = {'thx','thy','thz','rx','ry','thx_gT','thy_gT','rx_gT','ry_gT','thz_gT', ...
             'thxrz','thyrz','thxrz2','thyrz2', ...
             'thx2','thy2','thz2','rx2','ry2','thx2_gT','thy2_gT','rx2_gT','ry2_gT','thz2_gT', ...
             'rxry','rxthz','rythz','thxthy'};

tbl = array2table(predictors, 'VariableNames', predNames);
tbl.resid = resid1;
mdl = stepwiselm(tbl, 'ResponseVar', 'resid', 'Upper', 'linear', 'Verbose', 0);
disp(mdl);
fprintf('Stepwise in-sample R^2: %.4f\n', mdl.Rsquared.Ordinary);
selectedTerms = mdl.PredictorNames;
fprintf('Selected: %s\n', strjoin(selectedTerms, ', '));

rng(0);
kfold = 10;
cvIdx = mod((0:n-1)', kfold) + 1;
cvIdx = cvIdx(randperm(n));
yCV = zeros(n,1);
Xsel = predictors(:, ismember(predNames, selectedTerms));
for f = 1:kfold
    testMask = (cvIdx == f);
    trainMask = ~testMask;
    mf = fitlm(Xsel(trainMask,:), resid1(trainMask));
    yCV(testMask) = predict(mf, Xsel(testMask,:));
end
r2_cv_resid = 1 - sum((resid1-yCV).^2)/sum((resid1-mean(resid1)).^2);
yCV_W = Wfit1 + yCV;
r2_cv_W = 1 - sum((W-yCV_W).^2)/sum((W-mean(W)).^2);
fprintf('Residual 10-fold CV R^2: %.4f\n', r2_cv_resid);
fprintf('Combined 10-fold CV R^2 on W: %.4f\n', r2_cv_W);

save(fullfile(scriptDir, 'twoPoleModel.mat'), 'k1','k2','g0_z','mdl','selectedTerms','predNames');
