%% fitTwoMechanismModel.m
% Two-mechanism physical model:
%   W = P1(thx,thy,thz) / gap        -- parallel-plate/near-field term:
%                                        tilt changes the local gap, so
%                                        it scales with proximity (1/gap)
%     + P2(rx,ry)                    -- fringing-field term: lateral
%                                        overlap loss, roughly
%                                        gap-independent (confirmed by
%                                        the pure-1/gap factoring test,
%                                        which dropped rx,ry entirely)
% Both P1/gap and P2 are linear in their coefficients once the /gap
% division is precomputed as a feature transform, so the whole thing
% fits as ONE ordinary least-squares regression -- no staging, no
% nonlinear optimizer.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));

g0_z = 15e-6;
n = numel(rz);
thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;
gap = g0_z - rz;

% Parallel-plate group (divided by gap): constant, tilt linear, tilt
% quadratic, tilt cross terms
P1terms = [1./gap, thx_rad./gap, thy_rad./gap, thz_rad./gap, ...
           thx_rad.^2./gap, thy_rad.^2./gap, thz_rad.^2./gap, ...
           thx_rad.*thy_rad./gap, thx_rad.*thz_rad./gap, thy_rad.*thz_rad./gap];
P1names = {'invgap','thx_g','thy_g','thz_g','thx2_g','thy2_g','thz2_g','thxthy_g','thxthz_g','thythz_g'};

% Fringing group (NOT divided by gap): constant, rx/ry linear, quadratic,
% cross term
P2terms = [ones(n,1), rx, ry, rx.^2, ry.^2, rx.*ry];
P2names = {'const','rx','ry','rx2','ry2','rxry'};

predictors = [P1terms, P2terms];
predNames = [P1names, P2names];

% invgap ~6.7e4 vs rx^2 ~2e-9 -- a 13-order-of-magnitude spread that
% crashes fitlm's rank check (silently zeroed rx2,ry2,rxry last run).
% Normalize each column by its own max abs value, fit, then rescale
% coefficients back to original units afterward.
scale = max(abs(predictors), [], 1);
predictors_n = predictors ./ scale;

mdl = fitlm(predictors_n, W, 'Intercept', false, 'VarNames', [predNames, {'W'}]);
disp(mdl);
fprintf('In-sample R^2: %.4f\n', mdl.Rsquared.Ordinary);

coefs_n = mdl.Coefficients.Estimate;
coefs = coefs_n(:) ./ scale(:);   % rescale back to original units
fprintf('\n--- Coefficients (original units) ---\n');
for i = 1:numel(predNames)
    fprintf('%10s: %.4g\n', predNames{i}, coefs(i));
end

%% 10-fold CV (on normalized columns, same scale factors reused)
rng(0);
kfold = 10;
cvIdx = mod((0:n-1)', kfold) + 1;
cvIdx = cvIdx(randperm(n));
yCV = zeros(n,1);
for f = 1:kfold
    testMask = (cvIdx == f);
    trainMask = ~testMask;
    mf = fitlm(predictors_n(trainMask,:), W(trainMask), 'Intercept', false);
    yCV(testMask) = predict(mf, predictors_n(testMask,:));
end
r2_cv = 1 - sum((W-yCV).^2)/sum((W-mean(W)).^2);
fprintf('\n10-fold CV R^2: %.4f\n', r2_cv);
fprintf('(compare: additive-terms model was ~0.996 CV; pure 1/gap-only model was 0.78 CV)\n');

save(fullfile(scriptDir, 'twoMechanismModel.mat'), 'mdl', 'predNames', 'g0_z', 'coefs');
