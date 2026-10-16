%% fitCoupledPhysicsModel.m
% Fits the proposed coupled physics models against the REAL 6D LHS
% data (55 points) -- required here, not the single-axis sweeps,
% since identifying the cross-coupling terms (k2,k3 for Fx; the
% rx/ry/thz dependence for Tauz) needs data where these variables
% genuinely vary together.
%
%   Fx(rx,rz,thx,thy)   = k   / (g0_z - rz - k2*thx - k3*thy)^2 * poly3(rx)
%   Tauz(rz,rx,ry,thz)  = k4  / (g0_z - rz)^2 * [poly3(rx) + poly3(ry) + poly3(thz)]
%                          (additive combination, one shared intercept --
%                           not a full multivariate poly3 with cross terms,
%                           which would need ~20 parameters for only 55 points)
%
% thx/thy converted to radians inside the gap term so k2/k3 come out
% as physical lever-arm lengths (m/rad), consistent with the rotation
% geometry used throughout this project.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

g0_z = 15e-6;
lhsData = loadAllComsolResults();

rx = lhsData.rx; ry = lhsData.ry; rz = lhsData.rz;
thx_rad = lhsData.thx * pi/180;
thy_rad = lhsData.thy * pi/180;
thz = lhsData.thz;
n = height(lhsData);

xScale = max(abs(rx));
yScale = max(abs(ry));
zScale = max(abs(thz));

opts = statset('MaxIter', 2000, 'Display', 'off');

%% Fx model
% k fixed to k_Taux (2.08502e-17), reused from the earlier validated
% Taux fit, rather than fit fresh -- cuts Fx's free parameters from 7
% to 6 to address the non-identifiability/overfitting seen when k was
% left free (LOO R^2 = -5907).
k_Taux_fixed = 2.08502e-17;

Fx = lhsData.Fx;
Xfx = [rz, thx_rad, thy_rad, rx/xScale];
fxFun = @(b,X) k_Taux_fixed ./ (g0_z - X(:,1) - b(1)*X(:,2) - b(2)*X(:,3)).^2 .* ...
    (b(3) + b(4)*X(:,4) + b(5)*X(:,4).^2 + b(6)*X(:,4).^3);
b0_fx = [0, 0, 2e-8, 6e-7, -1e-7, -1e-6];

mdlFx = fitnlm(Xfx, Fx, fxFun, b0_fx, 'Options', opts);
bFx = mdlFx.Coefficients.Estimate;

% Leave-one-out validation (refit with each point held out)
yLOO_fx = zeros(n,1);
for i = 1:n
    keep = true(n,1); keep(i) = false;
    try
        mdl_i = fitnlm(Xfx(keep,:), Fx(keep), fxFun, bFx, 'Options', opts);
        yLOO_fx(i) = predict(mdl_i, Xfx(i,:));
    catch
        yLOO_fx(i) = NaN;
    end
end
valid = ~isnan(yLOO_fx);
r2_fx_loo = 1 - sum((Fx(valid)-yLOO_fx(valid)).^2) / sum((Fx(valid)-mean(Fx(valid))).^2);

fprintf('=== Fx model ===\n');
fprintf('k=k_Taux=%.4g (fixed)  k2=%.4g  k3=%.4g\n', k_Taux_fixed, bFx(1), bFx(2));
fprintf('poly3(rx) coeffs: %.4g %.4g %.4g %.4g\n', bFx(3), bFx(4), bFx(5), bFx(6));
fprintf('In-sample R^2: %.3f   LOO R^2: %.3f  (%d/%d LOO refits converged)\n', ...
    mdlFx.Rsquared.Ordinary, r2_fx_loo, sum(valid), n);

%% Tauz model
Tauz = lhsData.Tauz;
Xtz = [rz, rx/xScale, ry/yScale, thz/zScale];
tzFun = @(b,X) b(1) ./ (g0_z - X(:,1)).^2 .* ...
    (b(2) + b(3)*X(:,2) + b(4)*X(:,2).^2 + b(5)*X(:,2).^3 + ...
            b(6)*X(:,3) + b(7)*X(:,3).^2 + b(8)*X(:,3).^3 + ...
            b(9)*X(:,4) + b(10)*X(:,4).^2 + b(11)*X(:,4).^3);
b0_tz = [1e-17, 1e-8, 1e-8, -1e-9, -1e-8, 1e-8, -1e-9, -1e-8, 1e-8, -1e-9, -1e-8];

mdlTz = fitnlm(Xtz, Tauz, tzFun, b0_tz, 'Options', opts);
bTz = mdlTz.Coefficients.Estimate;

yLOO_tz = zeros(n,1);
for i = 1:n
    keep = true(n,1); keep(i) = false;
    try
        mdl_i = fitnlm(Xtz(keep,:), Tauz(keep), tzFun, bTz, 'Options', opts);
        yLOO_tz(i) = predict(mdl_i, Xtz(i,:));
    catch
        yLOO_tz(i) = NaN;
    end
end
valid = ~isnan(yLOO_tz);
r2_tz_loo = 1 - sum((Tauz(valid)-yLOO_tz(valid)).^2) / sum((Tauz(valid)-mean(Tauz(valid))).^2);

fprintf('\n=== Tauz model ===\n');
fprintf('k4=%.4g\n', bTz(1));
fprintf('poly3(rx) coeffs: %.4g %.4g %.4g %.4g\n', bTz(2), bTz(3), bTz(4), bTz(5));
fprintf('poly3(ry) coeffs: %.4g %.4g %.4g\n', bTz(6), bTz(7), bTz(8));
fprintf('poly3(thz) coeffs: %.4g %.4g %.4g\n', bTz(9), bTz(10), bTz(11));
fprintf('In-sample R^2: %.3f   LOO R^2: %.3f  (%d/%d LOO refits converged)\n', ...
    mdlTz.Rsquared.Ordinary, r2_tz_loo, sum(valid), n);

%% Plots
figure('Color','w','Name','Coupled physics models vs real 6D LHS data');
subplot(1,2,1); hold on; grid on;
scatter(Fx, predict(mdlFx, Xfx), 30, 'filled');
lims = [min([Fx; predict(mdlFx,Xfx)]), max([Fx; predict(mdlFx,Xfx)])];
plot(lims, lims, 'r--');
xlabel('actual Fx'); ylabel('model Fx');
title(sprintf('Fx (in-sample R^2=%.2f, LOO R^2=%.2f)', mdlFx.Rsquared.Ordinary, r2_fx_loo));

subplot(1,2,2); hold on; grid on;
scatter(Tauz, predict(mdlTz, Xtz), 30, 'filled');
lims = [min([Tauz; predict(mdlTz,Xtz)]), max([Tauz; predict(mdlTz,Xtz)])];
plot(lims, lims, 'r--');
xlabel('actual Tauz'); ylabel('model Tauz');
title(sprintf('Tauz (in-sample R^2=%.2f, LOO R^2=%.2f)', mdlTz.Rsquared.Ordinary, r2_tz_loo));

exportgraphics(gcf, fullfile(scriptDir, 'coupledPhysicsModel_fit.png'));
fprintf('\nSaved plot to coupledPhysicsModel_fit.png\n');
