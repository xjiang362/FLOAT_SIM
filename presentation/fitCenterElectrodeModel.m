%% fitCenterElectrodeModel.m
% Fits the center electrode's Fz(rz,thx,thy), Taux(thx), Tauy(thy) from
% the dedicated 150-pt sweep (comsol_centerElectrodeSweep_rzThxThy150_
% results.txt: rz-only rows 1-50, thx-only 51-100, thy-only 101-150).
% Also checks whether Fx, Fy are explained by the global-frame
% projection of the dominant Fz under tilt (per the "no fringing,
% centered electrode" physical reasoning), rather than fit as
% independent physics.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_centerElectrodeSweep_rzThxThy150_results.txt'), 'NumHeaderLines', 5);
rx=raw(:,1); ry=raw(:,2); rz=raw(:,3); thx=raw(:,4); thy=raw(:,5); thz=raw(:,6);
Fx=raw(:,7); Fy=raw(:,8); Fz=raw(:,9); Taux=raw(:,10); Tauy=raw(:,11); Tauz=raw(:,12); W=raw(:,13);

blk.rz  = 1:50;
blk.thx = 51:100;
blk.thy = 101:150;

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
V = 30;
g0_z = 15e-6;

%% --- Fz(rz,thx,thy): parallel-plate form, but QUADRATIC tilt terms ---
% Raw data check showed Fz(thx) is genuinely quadratic (R^2=0.99 for a
% pure thx^2 fit), not linear -- expected for a symmetric electrode
% (tilt reduces the effective gap symmetrically, an EVEN function of
% thx/thy), unlike the quadrant electrodes' linear-dominated G_eff.
X = [rz, (thx*pi/180).^2, (thy*pi/180).^2];   % [rz, thx_rad^2, thy_rad^2]
y = 1 ./ sqrt(Fz);
designMat = [ones(150,1), -X(:,1), -X(:,2), -X(:,3)];
coeffs = designMat \ y;
A = coeffs(1); B = coeffs(2); C = coeffs(3); D = coeffs(4);
k_fit = (g0_z / (V*A))^2;
a5 = B * g0_z / A;
a2 = C * g0_z / A;
a4 = D * g0_z / A;

Geff_fun = @(rz,thxr2,thyr2) g0_z - a5*rz - a2*thxr2 - a4*thyr2;
Fz_model = (k_fit*V^2) ./ max(Geff_fun(X(:,1),X(:,2),X(:,3)), 1e-9).^2;

fprintf('=== Fz(rz,thx,thy) [quadratic tilt terms] ===\n');
fprintf('k=%.6e, a2(thx^2)=%.6e, a4(thy^2)=%.6e, a5(rz)=%.6e\n', k_fit, a2, a4, a5);
for f = fieldnames(blk)'
    idx = blk.(f{1});
    fprintf('  %s-only R^2 vs MST = %.4f\n', f{1}, r2fun(Fz(idx), Fz_model(idx)));
end

%% --- Taux(thx), Tauy(thy): check raw shape, fit simple polynomial ---
fprintf('\n=== Taux(thx) raw shape (thx-only block) ===\n');
thxB = thx(blk.thx); TauxB = Taux(blk.thx);
p1 = polyfit(thxB, TauxB, 1);
p3 = polyfit(thxB, TauxB, 3);
fprintf('linear R^2=%.4f, cubic R^2=%.4f\n', r2fun(TauxB,polyval(p1,thxB)), r2fun(TauxB,polyval(p3,thxB)));
fprintf('linear coeff (slope, intercept): %.6e, %.6e\n', p1(1), p1(2));

fprintf('\n=== Tauy(thy) raw shape (thy-only block) ===\n');
thyB = thy(blk.thy); TauyB = Tauy(blk.thy);
q1 = polyfit(thyB, TauyB, 1);
q3 = polyfit(thyB, TauyB, 3);
fprintf('linear R^2=%.4f, cubic R^2=%.4f\n', r2fun(TauyB,polyval(q1,thyB)), r2fun(TauyB,polyval(q3,thyB)));
fprintf('linear coeff (slope, intercept): %.6e, %.6e\n', q1(1), q1(2));

%% --- Fx, Fy: check the global-frame tilt-projection hypothesis ---
% For a small rotation by (thx,thy) about the centroid, a body-frame
% force purely along z, [0;0;Fz_body], resolves onto the GLOBAL x,y axes
% approximately as Fx_global ~ Fz_body*thy_rad, Fy_global ~ -Fz_body*thx_rad
% (small-angle rotation-matrix projection; sign depends on rotation
% convention, checked empirically below rather than assumed).
fprintf('\n=== Fx, Fy vs tilt-projection hypothesis (thx-only, thy-only blocks) ===\n');
idx = blk.thx;
Fx_pred_a = Fz(idx) .* (thy(idx)*pi/180);       % should be ~0 here (thy=0 in this block)
Fy_pred_a = -Fz(idx) .* (thx(idx)*pi/180);
fprintf('thx-only: corr(Fy, -Fz*thx_rad) = %.4f, corr(Fx, Fz*thy_rad)=n/a (thy=0 here)\n', corr(Fy(idx), Fy_pred_a));
idx = blk.thy;
Fx_pred_b = Fz(idx) .* (thy(idx)*pi/180);
fprintf('thy-only: corr(Fx, Fz*thy_rad) = %.4f\n', corr(Fx(idx), Fx_pred_b));

fprintf('\nRaw Fx,Fy magnitude check:\n');
fprintf('thx-only block: mean|Fx|=%.4e, mean|Fy|=%.4e, mean Fz=%.4e\n', mean(abs(Fx(blk.thx))), mean(abs(Fy(blk.thx))), mean(Fz(blk.thx)));
fprintf('thy-only block: mean|Fx|=%.4e, mean|Fy|=%.4e, mean Fz=%.4e\n', mean(abs(Fx(blk.thy))), mean(abs(Fy(blk.thy))), mean(Fz(blk.thy)));

save(fullfile(scriptDir, 'centerElectrodeModel.mat'), 'k_fit','a2','a4','a5','g0_z','p1','q1');
