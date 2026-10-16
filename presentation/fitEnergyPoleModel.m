%% fitEnergyPoleModel.m
% Physics-informed fit for total electric energy W over the full 6D
% LHS sweep: W ~ k / (g0_z - rz - a*thx^2 - b*thy^2)^p -- reuses the
% same pole structure already validated for Fz=k/(g0_z-rz)^2 (since
% F=-dW/dz implies W~k/(g0_z-rz) to leading order), with an EVEN
% (mirror-symmetric) quadratic tilt correction to the effective gap
% rather than a linear one, since W(thx)=W(-thx) and W(thy)=W(-thy) by
% the trap's lateral symmetry. rx, ry, thz are checked afterward as
% secondary corrections on the residual, since they aren't expected to
% much affect the dominant vertical-gap energy term.
%
% Validated via 10-fold cross-validation (not plain in-sample R^2) --
% full LOO on a 5-parameter nonlinear fit over 2000 pts is unnecessary
% and slow; 10-fold gives a robust out-of-sample check.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));

g0_z = 15e-6;
n = numel(rz);

% thx, thy in degrees as loaded -- convert to radians for a physically
% meaningful gap correction (a,b will then be in meters/rad^2, order of
% plate dimension L/2 ~ 3e-3 if the small-angle corner-dip picture holds)
thx_rad = thx * pi/180;
thy_rad = thy * pi/180;

X = [rz, thx_rad, thy_rad];
% Clamp denominator to a small positive floor -- LM steps during the
% search can otherwise push (g0_z-rz-a*thx^2-b*thy^2) negative for
% large-tilt points even from a valid starting point, producing
% complex/NaN residuals that break the fit.
gapFloor = 1e-7;
modelFun = @(b,X) b(1) ./ max(g0_z - X(:,1) - b(2)*X(:,2).^2 - b(3)*X(:,3).^2, gapFloor).^b(4);
% a,b sized so a*thx_rad^2 ~ 10% of g0_z at the extreme tilt (thx_rad~0.0023,
% thx_rad^2~5e-6): a ~ 0.1*1.5e-5/5e-6 ~ 0.3. The earlier guess of 2 was
% ~7x too large, pushing the denominator negative for large-tilt points.
b0 = [1.5e-14, 0.3, 0.3, 1];

opts = statset('MaxIter', 500, 'Display', 'off');
mdl = fitnlm(X, W, modelFun, b0, 'Options', opts);
b = mdl.Coefficients.Estimate;
fprintf('=== Pole model: W = k / (g0_z - rz - a*thx_rad^2 - b*thy_rad^2)^p ===\n');
fprintf('k=%.4g  a=%.4g  b=%.4g  p=%.4g\n', b(1), b(2), b(3), b(4));
fprintf('In-sample R^2: %.4f\n', mdl.Rsquared.Ordinary);

%% 10-fold CV
rng(0);
kfold = 10;
cvIdx = crossvalind_manual(n, kfold);
yCV = zeros(n,1);
for f = 1:kfold
    testMask = (cvIdx == f);
    trainMask = ~testMask;
    try
        mdl_f = fitnlm(X(trainMask,:), W(trainMask), modelFun, b, 'Options', opts);
        yCV(testMask) = predict(mdl_f, X(testMask,:));
    catch ME
        fprintf('Fold %d failed: %s\n', f, ME.message);
        yCV(testMask) = NaN;
    end
end
valid = ~isnan(yCV);
r2_cv = 1 - sum((W(valid)-yCV(valid)).^2) / sum((W(valid)-mean(W(valid))).^2);
fprintf('10-fold CV R^2: %.4f  (%d/%d valid)\n', r2_cv, sum(valid), n);

%% Residual check: does anything remain vs rx, ry, thz?
Wfit = predict(mdl, X);
resid = W - Wfit;
fprintf('\n--- Residual correlation with rx, ry, thz ---\n');
fprintf('corr(resid, rx)   = %.4f\n', corr(resid, rx));
fprintf('corr(resid, ry)   = %.4f\n', corr(resid, ry));
fprintf('corr(resid, thz)  = %.4f\n', corr(resid, thz));
fprintf('corr(resid, rx^2) = %.4f\n', corr(resid, rx.^2));
fprintf('corr(resid, ry^2) = %.4f\n', corr(resid, ry.^2));
fprintf('corr(resid, thz^2)= %.4f\n', corr(resid, thz.^2));

figure('Color','w','Name','Pole model fit: predicted vs actual W');
scatter(W, Wfit, 8, 'filled', 'MarkerFaceAlpha', 0.3);
hold on; grid on;
plot([min(W) max(W)], [min(W) max(W)], 'r--', 'LineWidth', 1.5);
xlabel('W actual [J]'); ylabel('W predicted [J]');
title(sprintf('Pole model: in-sample R^2=%.4f, 10-fold CV R^2=%.4f', mdl.Rsquared.Ordinary, r2_cv));
exportgraphics(gcf, fullfile(scriptDir, 'poleModel_predVsActual.png'));
fprintf('\nSaved plot to poleModel_predVsActual.png\n');

save(fullfile(scriptDir, 'poleModelFit.mat'), 'mdl', 'b', 'g0_z', 'r2_cv');
fprintf('Saved fit to poleModelFit.mat\n');

%% local function
function idx = crossvalind_manual(n, k)
    idx = mod((0:n-1)', k) + 1;
    idx = idx(randperm(n));
end
