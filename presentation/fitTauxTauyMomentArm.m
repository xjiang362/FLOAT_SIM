%% fitTauxTauyMomentArm.m
% Equivalent-moment-arm model: Taux = Leqx*Fz, Tauy = Leqy*Fz, with
%   Leqx = a0 + a1*rz + a2*thx + a3*thy + a4*(rz*thx) + a5*(rz*thy)
%   Leqy = a0 + a1*rz + a2*thy + a3*thx + a4*(rz*thy) + a5*(rz*thx)
% Leqx and Leqy SHARE the same 6 coefficients (a0-a5), just with thx/thy
% (and their rz-cross terms) swapped between "own axis" and "other axis"
% roles -- reflecting the physical x/y symmetry of the trap. Fit jointly
% by stacking the Leqx and Leqy equations (with the appropriate variable
% swap) into one combined linear least-squares problem.
%
% Dataset: same as the Fz fit -- the dedicated rz-only/thx-only/thy-only
% blocks of comsol_sweep_combined300_results.txt (150 points), using the
% real Taux_MST/Tauy_MST/Fz_MST directly (Leqx_data = Taux_MST/Fz_MST,
% Leqy_data = Tauy_MST/Fz_MST), not the fitted Fz model.
%
% NOTE: a4 (rz*thx) and a5 (rz*thy) are exactly unidentifiable from this
% dataset for the same structural reason found in the Fz fit -- every
% point has either rz=0 or thx=thy=0, so both cross-term columns are
% identically zero throughout. Reported here anyway to show this
% explicitly (expect them to come back exactly at whatever initial
% value/collinear combination, not to be trusted) -- the a0-a3-only
% fit is the one to actually use.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

data = parseComsolResults(fullfile(scriptDir, 'comsol_sweep_combined300_results.txt'));
idx = [101:150, 151:200, 201:250];   % rz-only, thx-only, thy-only blocks
b = data(idx,:);
n = height(b);

rz = b.rz; thx = b.thx*pi/180; thy = b.thy*pi/180;
Fz_MST = b.Fz; Taux_MST = b.Taux; Tauy_MST = b.Tauy;

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

%% Fit a0-a3 only (a4,a5 dropped -- unidentifiable from this dataset)
% Regress Taux/Tauy DIRECTLY against Fz*[1,rz,thx,thy] rather than
% fitting the ratio Taux/Fz: Fz swings ~700x across the rz-only block
% (near a pull-in-like divergence at large rz), so dividing by it
% amplifies a few extreme points into huge outliers that wreck an
% unweighted least-squares fit of the ratio. Multiplying instead of
% dividing keeps everything on Taux's own natural scale.
% Leqx row: [1, rz, thx, thy]*Fz, Leqy row: -[1, rz, thy, thx]*Fz (swapped
% AND sign-flipped -- raw data shows Leqy's range is the negative mirror
% of Leqx's, i.e. Leqy=-(a0+a1*rz+a2*thy+a3*thx), not the plain +swap
% originally assumed; a plain-sign shared fit forced one coefficient to
% be both positive and negative at once, which is why it failed badly).
designX = [ones(n,1), rz, thx, thy] .* Fz_MST;
designY = -[ones(n,1), rz, thy, thx] .* Fz_MST;
designStacked = [designX; designY];
targetStacked = [Taux_MST; Tauy_MST];

coeffs = designStacked \ targetStacked;
a0 = coeffs(1); a1 = coeffs(2); a2 = coeffs(3); a3 = coeffs(4);
fprintf('Shared moment-arm coefficients (a0-a3 only, a4/a5 dropped as unidentifiable):\n');
fprintf('  a0 = %.6e\n  a1(rz) = %.6e\n  a2(own-axis thx/thy) = %.6e\n  a3(cross-axis thy/thx) = %.6e\n', a0, a1, a2, a3);

Taux_fit4 = designX * coeffs;
Tauy_fit4 = designY * coeffs;

fprintf('\nR^2 (Taux, a0-a3 model) = %.4f\n', r2fun(Taux_MST, Taux_fit4));
fprintf('R^2 (Tauy, a0-a3 model) = %.4f\n', r2fun(Tauy_MST, Tauy_fit4));

%% Also try full a0-a5 (for reference/diagnosis -- expect degeneracy)
designX6 = [ones(n,1), rz, thx, thy, rz.*thx, rz.*thy] .* Fz_MST;
designY6 = -[ones(n,1), rz, thy, thx, rz.*thy, rz.*thx] .* Fz_MST;
designStacked6 = [designX6; designY6];
coeffs6 = designStacked6 \ targetStacked;
fprintf('\nFull a0-a5 fit (diagnostic -- a4,a5 expected unreliable):\n');
names6 = {'a0','a1(rz)','a2(own)','a3(cross)','a4(rz*own)','a5(rz*cross)'};
for i = 1:6
    fprintf('  %s = %.6e\n', names6{i}, coeffs6(i));
end
Taux_fit6 = designX6*coeffs6;
Tauy_fit6 = designY6*coeffs6;
fprintf('R^2 (Taux, a0-a5 model) = %.4f\n', r2fun(Taux_MST, Taux_fit6));
fprintf('R^2 (Tauy, a0-a5 model) = %.4f\n', r2fun(Tauy_MST, Tauy_fit6));

save(fullfile(scriptDir, 'tauxTauyMomentArm.mat'), 'coeffs', 'a0','a1','a2','a3');

%% Plot (a0-a3 model)
blockNames = {'rz-only','thx-only','thy-only'};
blockIdx = {1:50, 51:100, 101:150};
xVals = {rz(blockIdx{1})*1e6, thx(blockIdx{2})*180/pi, thy(blockIdx{3})*180/pi};
xNames = {'rz [\mum]', '\theta_x [deg]', '\theta_y [deg]'};

figure('Color','w','Name','Taux, Tauy: moment-arm model vs MST', 'Position',[50 50 1700 900]);
for i = 1:3
    bi = blockIdx{i};
    subplot(2,3,i); hold on; grid on;
    scatter(xVals{i}, Taux_MST(bi), 25, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
    scatter(xVals{i}, Taux_fit4(bi), 25, 'r', 'filled', 'MarkerFaceAlpha', 0.6, ...
        'DisplayName', sprintf('Moment-arm model (R^2=%.4f)', r2fun(Taux_MST(bi), Taux_fit4(bi))));
    xlabel(xNames{i}); ylabel('Taux [N*m]');
    title(['Taux vs ' xNames{i} ' (' blockNames{i} ')']);
    legend('Location','best');

    subplot(2,3,i+3); hold on; grid on;
    scatter(xVals{i}, Tauy_MST(bi), 25, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
    scatter(xVals{i}, Tauy_fit4(bi), 25, 'r', 'filled', 'MarkerFaceAlpha', 0.6, ...
        'DisplayName', sprintf('Moment-arm model (R^2=%.4f)', r2fun(Tauy_MST(bi), Tauy_fit4(bi))));
    xlabel(xNames{i}); ylabel('Tauy [N*m]');
    title(['Tauy vs ' xNames{i} ' (' blockNames{i} ')']);
    legend('Location','best');
end

exportgraphics(gcf, fullfile(scriptDir, 'TauxTauy_momentArmModel.png'));
fprintf('\nSaved plot to TauxTauy_momentArmModel.png\n');
