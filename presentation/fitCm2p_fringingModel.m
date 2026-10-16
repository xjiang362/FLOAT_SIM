%% fitCm2p_fringingModel.m
% Fits Cm(2,p) (resonator <-> terminal 2 coupling) using the same reduced
% fringing-capacitance model that worked for Cm(1,p): shared nominal
% width/length (WL0), gap G varying with rz, thx, thy only (rx,ry terms
% dropped as unidentifiable/no real dependence).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'capSweep300_withGND_allMatrices.mat'));

Cm2p = squeeze(Cm_all(2,21,:));
thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;

eps0 = 8.854e-12;
floorVal = 1e-9;

Cf = @(G,Wd,Ld) (eps0*Wd.*Ld./G) .* ...
     (1 + G./(pi*Wd) + G./(pi*Wd).*log(2*pi*Wd./G)) .* ...
     (1 + G./(pi*Ld) + G./(pi*Ld).*log(2*pi*Ld./G));

% params: [k, G0, WL0, a1(thx), a2(thx2), a3(thy), a4(thy2), a5(rzscale)]
modelFun = @(p, X) p(1) * Cf( ...
    max(p(2) - p(8)*X(:,3) - p(4)*X(:,4) - p(5)*X(:,4).^2 - p(6)*X(:,5) - p(7)*X(:,5).^2, floorVal), ...
    p(3) * ones(size(X,1),1), ...
    p(3) * ones(size(X,1),1) );

X = [rx, ry, rz, thx_rad, thy_rad, thz_rad];

G0_0 = 30e-6; WL0_0 = 500e-6;
Cf0 = Cf(G0_0, WL0_0, WL0_0);
k0 = median(Cm2p) / Cf0;
p0 = [k0, G0_0, WL0_0, 1e-3, 0.3, 1e-3, 0.3, 1];

opts = statset('MaxIter', 2000, 'Display', 'off', 'TolFun', 1e-30, 'TolX', 1e-20);
[pfit, resid_fit] = nlinfit(X, Cm2p, modelFun, p0, opts);

names = {'k','G0','WL0','a1(thx)','a2(thx2)','a3(thy)','a4(thy2)','a5(rzscale)'};
fprintf('Fitted parameters:\n');
for i = 1:numel(names)
    fprintf('  %s = %.6e\n', names{i}, pfit(i));
end

r2 = 1 - sum(resid_fit.^2)/sum((Cm2p-mean(Cm2p)).^2);
fprintf('\nIn-sample R^2: %.4f\n', r2);
fprintf('max|imag(pfit)| = %.4e\n', max(abs(imag(pfit))));

%% 10-fold CV
rng(0);
n = numel(Cm2p);
kfold = 10;
cvIdx = mod((0:n-1)', kfold) + 1;
cvIdx = cvIdx(randperm(n));
yCV = zeros(n,1);
for f = 1:kfold
    testMask = (cvIdx == f);
    trainMask = ~testMask;
    try
        pf = nlinfit(X(trainMask,:), Cm2p(trainMask), modelFun, pfit, opts);
        yCV(testMask) = modelFun(pf, X(testMask,:));
    catch
        yCV(testMask) = NaN;
    end
end
valid = ~isnan(yCV);
r2_cv = 1 - sum((Cm2p(valid)-yCV(valid)).^2)/sum((Cm2p(valid)-mean(Cm2p(valid))).^2);
fprintf('10-fold CV R^2: %.4f (%d/%d valid)\n', r2_cv, sum(valid), n);

save(fullfile(scriptDir, 'Cm2p_fringingFit.mat'), 'pfit', 'modelFun', 'r2', 'r2_cv');
