%% fitSimpleFzModel.m
% Simplified Fz model: drop the W_eff/L_eff fringing-formula machinery
% entirely (fix them at the nominal W0, L0 -- their rx/ry/thz
% corrections were already agreed negligible for Fz/Taux/Tauy) and use
% the direct parallel-plate-style force law against the effective gap:
%
%   G_eff = g0_z - a5*rz - a1*thx - a2*(rz*thx) - a3*thy - a4*(rz*thy)
%   Fz = k*V^2 / G_eff^2
%
% This is a direct algebraic force formula (no energy differentiation
% needed), so it's fit directly against Fz_MST -- using the dedicated
% rz-only, thx-only, thy-only blocks of comsol_sweep_combined300_results.txt
% (150 points total), which isolate each DOF cleanly.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

data = parseComsolResults(fullfile(scriptDir, 'comsol_sweep_combined300_results.txt'));
idx = [101:150, 151:200, 201:250];   % rz-only, thx-only, thy-only blocks
b = data(idx,:);
n = height(b);
V = 30;

thx_rad = b.thx*pi/180; thy_rad = b.thy*pi/180;
X = [b.rz, thx_rad, thy_rad];
Fz_MST = b.Fz;

% params: [k, a1(thx), a3(thy), a5(rz scale)]
% NOTE: the rz*thx / rz*thy cross terms (a2, a4) are dropped here -- in
% this single-axis-only dataset, every row has either rz=0 or thx=thy=0,
% so rz*thx and rz*thy are EXACTLY zero at every point (not just weakly
% identified, structurally invisible). They would need data where rz
% and thx/thy vary simultaneously (e.g. the scattered full6D LHS set) to
% be fit at all.
%
% Fz=k*V^2/G^2 with G=g0_z-a5*rz-a1*thx-a3*thy is badly scaled for a
% direct nlinfit (k~1e-17 vs a5~1 vs g0_z~1e-5 span ~20 orders of
% magnitude, causing the optimizer to move k and freeze everything else
% -- confirmed by the rank-deficient warnings and by a1/a3/a5 coming
% back EXACTLY at their initial guesses). Fixed with an exact
% linearizing transform instead of iterative NLS:
%   1/sqrt(Fz) = G/(V*sqrt(k)) = A - B*rz - C*thx - D*thy
% which is LINEAR in (A,B,C,D) -- solved by ordinary least squares
% (backslash), then k,a1,a3,a5 recovered algebraically from A,B,C,D.
g0_z = 15e-6;
y = 1 ./ sqrt(Fz_MST);
designMat = [ones(n,1), -X(:,1), -X(:,2), -X(:,3)];   % [1, -rz, -thx, -thy]
coeffs = designMat \ y;
A = coeffs(1); B = coeffs(2); C = coeffs(3); D = coeffs(4);

k_fit = (g0_z / (V*A))^2;
a5 = B * g0_z / A;
a1 = C * g0_z / A;
a3 = D * g0_z / A;
pfit = [k_fit, a1, a3, a5];

floorVal = 1e-8;
modelFun = @(p, X) (p(1)*V^2) ./ max(g0_z - p(4)*X(:,1) - p(2)*X(:,2) - p(3)*X(:,3), floorVal).^2;
Fz_pred = modelFun(pfit, X);
resid_fit = Fz_MST - Fz_pred;

names = {'k','a1(thx)','a3(thy)','a5(rz scale)'};
fprintf('Fitted parameters (via linearizing transform, closed-form):\n');
for i = 1:numel(names)
    fprintf('  %s = %.6e\n', names{i}, pfit(i));
end

r2 = 1 - sum(resid_fit.^2)/sum((Fz_MST-mean(Fz_MST)).^2);
fprintf('\nIn-sample R^2 (rz+thx+thy dedicated blocks combined): %.4f\n', r2);
fprintf('R^2 in the LINEARIZED (1/sqrt(Fz)) space: %.4f\n', 1 - sum((y-designMat*coeffs).^2)/sum((y-mean(y)).^2));

%% Per-block R^2
blockNames = {'rz-only','thx-only','thy-only'};
blockIdx = {1:50, 51:100, 101:150};
Fz_model = modelFun(pfit, X);
r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
for i = 1:3
    bi = blockIdx{i};
    fprintf('%9s block R^2 = %.4f\n', blockNames{i}, r2fun(Fz_MST(bi), Fz_model(bi)));
end

save(fullfile(scriptDir, 'simpleFzModel.mat'), 'pfit', 'modelFun', 'g0_z', 'r2');

%% Plot
figure('Color','w','Name','Simple Fz=kV^2/Geff^2 model vs MST', 'Position',[50 50 1500 500]);
xVals = {b.rz(blockIdx{1})*1e6, b.thx(blockIdx{2}), b.thy(blockIdx{3})};
xNames = {'rz [\mum]', '\theta_x [deg]', '\theta_y [deg]'};
for i = 1:3
    bi = blockIdx{i};
    subplot(1,3,i); hold on; grid on;
    scatter(xVals{i}, Fz_MST(bi), 25, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
    scatter(xVals{i}, Fz_model(bi), 25, 'r', 'filled', 'MarkerFaceAlpha', 0.6, ...
        'DisplayName', sprintf('Model (R^2=%.4f)', r2fun(Fz_MST(bi), Fz_model(bi))));
    xlabel(xNames{i}); ylabel('Fz [N]');
    title(['Fz vs ' xNames{i} ' (' blockNames{i} ')']);
    legend('Location','best');
end

exportgraphics(gcf, fullfile(scriptDir, 'SimpleFzModel_vs_dedicatedSweeps.png'));
fprintf('\nSaved plot to SimpleFzModel_vs_dedicatedSweeps.png\n');
