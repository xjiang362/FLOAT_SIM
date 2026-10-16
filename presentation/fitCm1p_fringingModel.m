%% fitCm1p_fringingModel.m
% Fits Cm(1,p) (resonator <-> terminal 1 coupling) using the same 3D
% fringing-capacitance formula as the electrode force model, treating
% terminal 1 and the resonator as the two "plates". Nominal gap (G0) and
% a shared nominal width/length (WL0) are free fit parameters here
% (unlike the force model, where they were known physical electrode
% dimensions) since this is a different physical relationship. rx,ry
% dependence was dropped (unidentifiable, see below); the gap G varies
% with rz, thx, thy only.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'capSweep300_withGND_allMatrices.mat'));

Cm1p = squeeze(Cm_all(1,21,:));
thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;

eps0 = 8.854e-12;
floorVal = 1e-9;

Cf = @(G,Wd,Ld) (eps0*Wd.*Ld./G) .* ...
     (1 + G./(pi*Wd) + G./(pi*Wd).*log(2*pi*Wd./G)) .* ...
     (1 + G./(pi*Ld) + G./(pi*Ld).*log(2*pi*Ld./G));

% params: [k, G0, WL0, a1(thx), a2(thx2), a3(thy), a4(thy2), a5(rzscale)]
% b2,b3 (rx->W) and c2,c3 (ry->L) DROPPED: raw-shape plot showed no real
% rx,ry dependence for Cm1p, and the first fit attempt confirmed they are
% unidentifiable (rank deficient, nonsensical fitted magnitudes like
% b3=18.2, c3=-1.3e5) -- same failure mode as the S_R tilt terms. W and L
% were then fit separately as fixed constants and converged to nearly
% identical values (~12.5mm) -- the formula can't tell them apart once
% both are constants, so they are merged into one shared WL0 parameter.
modelFun = @(p, X) p(1) * Cf( ...
    max(p(2) - p(8)*X(:,3) - p(4)*X(:,4) - p(5)*X(:,4).^2 - p(6)*X(:,5) - p(7)*X(:,5).^2, floorVal), ...
    p(3) * ones(size(X,1),1), ...
    p(3) * ones(size(X,1),1) );

X = [rx, ry, rz, thx_rad, thy_rad, thz_rad];

% Initial guesses: G0,WL0 order ~10-100um / ~mm (unknown geometry), k
% calibrated to roughly match typical Cm1p scale
G0_0 = 30e-6; WL0_0 = 500e-6;
Cf0 = Cf(G0_0, WL0_0, WL0_0);
k0 = median(Cm1p) / Cf0;
p0 = [k0, G0_0, WL0_0, 1e-3, 0.3, 1e-3, 0.3, 1];

opts = statset('MaxIter', 2000, 'Display', 'off', 'TolFun', 1e-30, 'TolX', 1e-20);
[pfit, resid_fit] = nlinfit(X, Cm1p, modelFun, p0, opts);

names = {'k','G0','WL0','a1(thx)','a2(thx2)','a3(thy)','a4(thy2)','a5(rzscale)'};
fprintf('Fitted parameters:\n');
for i = 1:numel(names)
    fprintf('  %s = %.6e\n', names{i}, pfit(i));
end

r2 = 1 - sum(resid_fit.^2)/sum((Cm1p-mean(Cm1p)).^2);
fprintf('\nIn-sample R^2: %.4f\n', r2);
fprintf('max|imag(pfit)| = %.4e\n', max(abs(imag(pfit))));

%% 10-fold CV
rng(0);
n = numel(Cm1p);
kfold = 10;
cvIdx = mod((0:n-1)', kfold) + 1;
cvIdx = cvIdx(randperm(n));
yCV = zeros(n,1);
for f = 1:kfold
    testMask = (cvIdx == f);
    trainMask = ~testMask;
    try
        pf = nlinfit(X(trainMask,:), Cm1p(trainMask), modelFun, pfit, opts);
        yCV(testMask) = modelFun(pf, X(testMask,:));
    catch
        yCV(testMask) = NaN;
    end
end
valid = ~isnan(yCV);
r2_cv = 1 - sum((Cm1p(valid)-yCV(valid)).^2)/sum((Cm1p(valid)-mean(Cm1p(valid))).^2);
fprintf('10-fold CV R^2: %.4f (%d/%d valid)\n', r2_cv, sum(valid), n);

save(fullfile(scriptDir, 'Cm1p_fringingFit.mat'), 'pfit', 'modelFun', 'r2', 'r2_cv');
