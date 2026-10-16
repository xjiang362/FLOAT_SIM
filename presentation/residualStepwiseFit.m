%% residualStepwiseFit.m
% The hand-picked even-quadratic candidates (thx^2, thy^2, rx^2, ry^2,
% thz^2, and /(g0-rz) variants) only explained ~3-5% of the residual
% variance after removing the dominant rz pole term W=k/(g0_z-rz), on
% both the full set and the near-collision-outlier-excluded subset.
% Rather than keep guessing forms by hand, let stepwiselm pick from a
% broader candidate set (adds cross terms between rx,ry,thz and
% thx,thy), then validate the FINAL selected model via 10-fold CV
% (stepwise's own in-sample selection criterion isn't a substitute for
% out-of-sample validation).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));

g0_z = 15e-6; k = 1.601e-14; p = 1;
n = numel(rz);
thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;

Wfit1 = k ./ (g0_z - rz).^p;
resid1 = W - Wfit1;

gap = g0_z - rz;
% Include LINEAR terms too (thx, thy, rx, ry, thz and their /gap
% versions) -- the mirror-symmetry assumption (even terms only) was
% contradicted by the earlier finding of a persistent nonzero Fx at
% rx=0 on the dense rx-rz grid, so let the data decide rather than
% assuming symmetry forbids them. Also added genuine thx*rz/thy*rz
% (linear, not /gap) cross terms: the Taux/Tauy vs MST comparison
% showed matching slope but a mismatched, mirror-antisymmetric
% intercept across rz bins -- a sign that thx/gap was standing in for
% a genuinely different functional form (plain rz cross term) rather
% than truly explaining the rz-dependence.
predictors = [thx_rad, thy_rad, thz_rad, rx, ry, ...
              thx_rad./gap, thy_rad./gap, rx./gap, ry./gap, thz_rad./gap, ...
              thx_rad.*rz, thy_rad.*rz, thx_rad.*rz.^2, thy_rad.*rz.^2, ...
              thx_rad.^2, thy_rad.^2, thz_rad.^2, rx.^2, ry.^2, ...
              thx_rad.^2./gap, thy_rad.^2./gap, rx.^2./gap, ry.^2./gap, thz_rad.^2./gap, ...
              rx.*ry, rx.*thz_rad, ry.*thz_rad, thx_rad.*thy_rad];
predNames = {'thx','thy','thz','rx','ry','thx_g','thy_g','rx_g','ry_g','thz_g', ...
             'thxrz','thyrz','thxrz2','thyrz2', ...
             'thx2','thy2','thz2','rx2','ry2','thx2_g','thy2_g','rx2_g','ry2_g','thz2_g', ...
             'rxry','rxthz','rythz','thxthy'};

tbl = array2table(predictors, 'VariableNames', predNames);
tbl.resid = resid1;

mdl = stepwiselm(tbl, 'ResponseVar', 'resid', 'Upper', 'linear', 'Verbose', 0);
disp(mdl);
fprintf('Stepwise in-sample R^2: %.4f\n', mdl.Rsquared.Ordinary);

%% 10-fold CV of the stepwise-selected model structure
selectedTerms = mdl.PredictorNames;
fprintf('\nSelected terms: %s\n', strjoin(selectedTerms, ', '));

rng(0);
kfold = 10;
cvIdx = mod((0:n-1)', kfold) + 1;
cvIdx = cvIdx(randperm(n));
yCV = zeros(n,1);
if isempty(selectedTerms)
    fprintf('Stepwise selected NO terms -- residual is indistinguishable from noise given this candidate set.\n');
else
    Xsel = predictors(:, ismember(predNames, selectedTerms));
    for f = 1:kfold
        testMask = (cvIdx == f);
        trainMask = ~testMask;
        mf = fitlm(Xsel(trainMask,:), resid1(trainMask));
        yCV(testMask) = predict(mf, Xsel(testMask,:));
    end
    r2_cv_resid = 1 - sum((resid1-yCV).^2) / sum((resid1-mean(resid1)).^2);
    fprintf('Residual model 10-fold CV R^2: %.4f\n', r2_cv_resid);

    % Combined W-space CV R^2
    yCV_W = Wfit1 + yCV;
    r2_cv_W = 1 - sum((W-yCV_W).^2) / sum((W-mean(W)).^2);
    fprintf('Combined (pole+residual) 10-fold CV R^2 on W: %.4f\n', r2_cv_W);
end
