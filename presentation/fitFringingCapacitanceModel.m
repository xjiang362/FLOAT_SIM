%% fitFringingCapacitanceModel.m
% Physics-grounded fringing capacitance model (3D fringing formula from
% the literature, e.g. "A Fringing Capacitance Model for Electrostatic
% Microstructure"):
%   C_f(G,W,L) = (eps0*W*L/G) * (1+G/(pi*W)+G/(pi*W)*ln(2*pi*W/G))
%                              * (1+G/(pi*L)+G/(pi*L)*ln(2*pi*L/G))
% with effective geometry:
%   G_eff  = g0_z - a5*rz - a1*thx - a2*(rz*thx) - a3*thy - a4*(rz*thy)
%   W_eff  = W0 + a*rx + c*thz + d*(rx*thz) + e*(ry*thz)
%   L_eff  = L0 + b*ry - c*thz - d*(rx*thz) - e*(ry*thz)
% G_eff uses rz*thx/rz*thy cross terms (tilting shifts the gap by an
% amount scaling with how far the plate already is in rz -- validated
% against the dedicated single-axis sweeps: fixed Taux/Tauy from
% R^2=-0.27/-0.55 to +0.87/+0.79 versus the original same-axis
% thx^2/thy^2 version).
%
% W_eff/L_eff are now COUPLED: the thz and rx*thz/ry*thz terms enter
% with the SAME coefficients (c,d,e) but opposite sign in W vs L, since
% yawing the plate (thz) trades width for length in a complementary way
% rather than each varying independently -- W0=L0 (nominal square
% electrode), so this is a single shared yaw-coupling mechanism, not two
% separate ones.
%
% W_total = k * C_f(G_eff,W_eff,L_eff), k absorbing eps0*V^2/2 (unknown
% voltage/scale) as one free parameter. thx,thy,thz in RADIANS. Fit
% directly on the energy method's W data (full6D_energyData.mat) -- not
% MST -- consistent with how this model has always been fit.
%
% Forces/torques obtained by central-difference on the FITTED closed
% form itself (safe here -- smooth analytic function, not raw
% simulation data, so no mesh-rebuild noise).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));

g0_z = 15e-6; W0 = 6000e-6; L0 = 6000e-6;
eps0 = 8.854e-12;
n = numel(rz);
thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;

Cf = @(G,Wd,Ld) (eps0*Wd.*Ld./G) .* ...
     (1 + G./(pi*Wd) + G./(pi*Wd).*log(2*pi*Wd./G)) .* ...
     (1 + G./(pi*Ld) + G./(pi*Ld).*log(2*pi*Ld./G));

% params: [k, a1,a2,a3,a4, a5(rz scale), a(rx->W), b(ry->L), c(thz), d(rx*thz), e(ry*thz)]
% Floor G_eff,W_eff,L_eff to a small positive value -- during
% optimization the LM search can try parameter values that push these
% negative (log() then goes complex), and that complex contamination
% sticks in the fitted coefficients even after the search returns to
% valid territory. Clamping keeps modelFun strictly real everywhere.
%
% a5 (rz scale, free instead of fixed at 1) added to correct the
% systematic high-Fz overestimate seen with a plain "g0_z - rz" gap --
% letting the direct rz->gap mapping be a free scale rather than
% assumed exactly 1:1.
floorVal = 1e-8;
modelFun = @(p, X) p(1) * Cf( ...
    max(g0_z - p(6)*X(:,3) - p(2)*X(:,4) - p(3)*(X(:,3).*X(:,4)) - p(4)*X(:,5) - p(5)*(X(:,3).*X(:,5)), floorVal), ...
    max(W0 + p(7)*X(:,1) + p(9)*X(:,6) + p(10)*(X(:,1).*X(:,6)) + p(11)*(X(:,2).*X(:,6)), floorVal), ...
    max(L0 + p(8)*X(:,2) - p(9)*X(:,6) - p(10)*(X(:,1).*X(:,6)) - p(11)*(X(:,2).*X(:,6)), floorVal) );

X = [rx, ry, rz, thx_rad, thy_rad, thz_rad];

% Initial guess: k from matching nominal (all DOF=0) point; geometry
% correction coefficients start near zero (no correction) except a
% small nonzero seed so the optimizer has a gradient to follow.
Cf0 = Cf(g0_z, W0, L0);
k0 = mean(W) / Cf0;
p0 = [k0, 1e-6, 1e-3, 1e-6, 1e-3, 1, 1e-3, 1e-3, 1e-3, 1e3, 1e3];

% Use nlinfit directly (not fitnlm) -- fitnlm's post-fit statistics
% machinery crashes when some parameters are weakly identified, even
% though the fit itself converges fine (SSE dropped from 8.16e-10 to
% 1.76e-16 in testing). nlinfit skips that and just returns coefficients.
opts = statset('MaxIter', 1000, 'Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);
[p, resid_fit] = nlinfit(X, W, modelFun, p0, opts);
fprintf('\n=== Fringing capacitance model ===\n');
names = {'k','a1(thx)','a2(rz*thx)','a3(thy)','a4(rz*thy)','a5(rz scale)','a(rx->W)','b(ry->L)','c(thz)','d(rx*thz)','e(ry*thz)'};
for i = 1:numel(p)
    fprintf('%14s: %.6g\n', names{i}, p(i));
end
r2_insample = 1 - sum(resid_fit.^2)/sum((W-mean(W)).^2);
fprintf('In-sample R^2: %.4f\n', r2_insample);

%% 10-fold CV
rng(0);
kfold = 10;
cvIdx = mod((0:n-1)', kfold) + 1;
cvIdx = cvIdx(randperm(n));
yCV = zeros(n,1);
for f = 1:kfold
    testMask = (cvIdx == f);
    trainMask = ~testMask;
    try
        pf = nlinfit(X(trainMask,:), W(trainMask), modelFun, p, opts);
        yCV(testMask) = modelFun(pf, X(testMask,:));
    catch ME
        fprintf('Fold %d failed: %s\n', f, ME.message);
        yCV(testMask) = NaN;
    end
end
valid = ~isnan(yCV);
r2_cv = 1 - sum((W(valid)-yCV(valid)).^2)/sum((W(valid)-mean(W(valid))).^2);
fprintf('10-fold CV R^2: %.4f (%d/%d valid)\n', r2_cv, sum(valid), n);

save(fullfile(scriptDir, 'fringingCapModel.mat'), 'p', 'modelFun', 'g0_z', 'W0', 'L0', 'r2_cv');
