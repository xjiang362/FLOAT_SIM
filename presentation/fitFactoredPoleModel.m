%% fitFactoredPoleModel.m
% Physically cleaner structure: W = P(rx,ry,thx,thy,thz) / gap, where
% gap=g0_z-rz and P is a polynomial in the OTHER 5 DOF only (no
% separate rz-dependence inside P -- the single pole from gap already
% carries all the rz-dependence, matching a parallel-plate capacitor
% picture W~Area_eff/gap with Area_eff varying with lateral position
% and tilt). This turns into an exact LINEAR regression: multiply both
% sides by gap, giving W*gap = P(rx,ry,thx,thy,thz), fit via ordinary
% least squares -- no nonlinear optimizer, no staged pole+residual
% fitting, no numerical identifiability issues.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));

g0_z = 15e-6;
n = numel(rz);
thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;
gap = g0_z - rz;

y = W .* gap;   % = P(rx,ry,thx,thy,thz), if the hypothesis holds

predictors = [thx_rad, thy_rad, thz_rad, rx, ry, ...
              thx_rad.^2, thy_rad.^2, thz_rad.^2, rx.^2, ry.^2, ...
              thx_rad.*thy_rad, thx_rad.*thz_rad, thy_rad.*thz_rad, ...
              rx.*ry, rx.*thx_rad, rx.*thy_rad, rx.*thz_rad, ...
              ry.*thx_rad, ry.*thy_rad, ry.*thz_rad];
predNames = {'thx','thy','thz','rx','ry','thx2','thy2','thz2','rx2','ry2', ...
             'thxthy','thxthz','thythz','rxry','rxthx','rxthy','rxthz','rythx','rythy','rythz'};

tbl = array2table(predictors, 'VariableNames', predNames);
tbl.y = y;
mdl = stepwiselm(tbl, 'ResponseVar', 'y', 'Upper', 'linear', 'Verbose', 0);
disp(mdl);
fprintf('Stepwise in-sample R^2 (on W*gap): %.4f\n', mdl.Rsquared.Ordinary);
selectedTerms = mdl.PredictorNames;
fprintf('Selected: %s\n', strjoin(selectedTerms, ', '));

%% 10-fold CV, checked on W directly (not W*gap, to match earlier comparisons)
rng(0);
kfold = 10;
cvIdx = mod((0:n-1)', kfold) + 1;
cvIdx = cvIdx(randperm(n));
yCV_Wgap = zeros(n,1);
Xsel = predictors(:, ismember(predNames, selectedTerms));
for f = 1:kfold
    testMask = (cvIdx == f);
    trainMask = ~testMask;
    mf = fitlm(Xsel(trainMask,:), y(trainMask));
    yCV_Wgap(testMask) = predict(mf, Xsel(testMask,:));
end
Wpred_CV = yCV_Wgap ./ gap;
r2_cv_W = 1 - sum((W-Wpred_CV).^2)/sum((W-mean(W)).^2);
fprintf('\n10-fold CV R^2 on W (factored model): %.4f\n', r2_cv_W);
fprintf('(compare: additive-terms model was 0.9964 in-sample recon / ~0.996 CV)\n');

save(fullfile(scriptDir, 'factoredPoleModel.mat'), 'mdl', 'selectedTerms', 'predNames', 'g0_z');
