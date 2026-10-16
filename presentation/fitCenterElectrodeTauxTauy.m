%% fitCenterElectrodeTauxTauy.m
% Fits Taux = Fz_model(rz,thx,thy)*thx*(a2+a1*rz), Tauy = Fz_model*thy*
% (a2+a1*rz) [same a1,a2 by the established x<->y symmetry] -- NOT a
% standalone a1*rz term (that would wrongly add a nonzero torque at
% thx=0, which the center electrode's data rules out), but rz MODULATING
% the thx-proportionality itself. Combines the original thx-only/
% thy-only sweep (rz=0) with the new rz=+-7um bracket sweep for a joint,
% well-identified fit -- linear in (a1,a2) given the already-fitted
% Fz_model, so solved in closed form.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'centerElectrodeModel.mat'), 'k_fit','a2','a4','a5','g0_z');
k_c = k_fit; a2_fz = a2; a4_fz = a4; a5_fz = a5;   % Fz model coefficients (renamed to avoid clashing with the a2 we're fitting here)

V = 30;
Geff_fun = @(rz,thxr,thyr) g0_z - a5_fz*rz - a2_fz*thxr.^2 - a4_fz*thyr.^2;
Fz_model = @(rz,thxr,thyr) (k_c*V^2) ./ max(Geff_fun(rz,thxr,thyr), 1e-9).^2;

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

%% --- Load original rz/thx/thy sweep (rz=0 for the thx-only/thy-only blocks) ---
raw1 = readmatrix(fullfile(scriptDir, 'comsol_centerElectrodeSweep_rzThxThy150_results.txt'), 'NumHeaderLines', 5);
rz1=raw1(:,3); thx1=raw1(:,4); thy1=raw1(:,5); Taux1=raw1(:,10); Tauy1=raw1(:,11);
blkThx0 = 51:100; blkThy0 = 101:150;

%% --- Load rz=+-7um bracket sweep ---
raw2 = readmatrix(fullfile(scriptDir, 'comsol_centerElectrode_rzCrossTerm200_results.txt'), 'NumHeaderLines', 5);
rz2=raw2(:,3); thx2=raw2(:,4); thy2=raw2(:,5); Taux2=raw2(:,10); Tauy2=raw2(:,11);
blkThxBr = 1:100; blkThyBr = 101:200;

%% --- Assemble combined Taux dataset: [rz, thx (rad)] -> Taux ---
rz_all_x  = [rz1(blkThx0);  rz2(blkThxBr)];
thx_all   = [thx1(blkThx0); thx2(blkThxBr)] * pi/180;
Taux_all  = [Taux1(blkThx0); Taux2(blkThxBr)];
thy_all_forx = [thy1(blkThx0); zeros(numel(blkThxBr),1)] * pi/180;   % thy=0 throughout both

% Regress Taux DIRECTLY against Fz*thx and Fz*rz*thx (NOT the ratio
% Taux/(Fz*thx) -- that blows up numerically wherever thx is near zero,
% which the warped sweep spacing concentrates points at; same pitfall
% found and fixed earlier in this project for Taux/Tauy vs Fz).
Fz_at_x = Fz_model(rz_all_x, thx_all, thy_all_forx);
designX = [Fz_at_x.*thx_all, Fz_at_x.*rz_all_x.*thx_all];
cX = designX \ Taux_all;
a2_leq = cX(1); a1_leq = cX(2);

Taux_model = Fz_at_x .* thx_all .* (a2_leq + a1_leq*rz_all_x);
fprintf('=== Taux model: Fz*thx*(a2+a1*rz) ===\n');
fprintf('a2=%.6e, a1=%.6e, R^2=%.4f\n', a2_leq, a1_leq, r2fun(Taux_all, Taux_model));

%% --- Same for Tauy (shared a1,a2 by x<->y symmetry -- verify, don't assume) ---
rz_all_y  = [rz1(blkThy0);  rz2(blkThyBr)];
thy_all   = [thy1(blkThy0); thy2(blkThyBr)] * pi/180;
Tauy_all  = [Tauy1(blkThy0); Tauy2(blkThyBr)];
thx_all_fory = [thx1(blkThy0); zeros(numel(blkThyBr),1)] * pi/180;

Fz_at_y = Fz_model(rz_all_y, thx_all_fory, thy_all);
designY = [Fz_at_y.*thy_all, Fz_at_y.*rz_all_y.*thy_all];
cY = designY \ Tauy_all;
a2_leqy = cY(1); a1_leqy = cY(2);

Tauy_model_own = Fz_at_y .* thy_all .* (a2_leqy + a1_leqy*rz_all_y);
fprintf('\n=== Tauy model (own fit): Fz*thy*(a2+a1*rz) ===\n');
fprintf('a2=%.6e, a1=%.6e, R^2=%.4f\n', a2_leqy, a1_leqy, r2fun(Tauy_all, Tauy_model_own));

Tauy_model_shared = Fz_at_y .* thy_all .* (a2_leq + a1_leq*rz_all_y);
fprintf('\nUsing Taux''s shared (a1,a2) on Tauy data: R^2=%.4f (check for x<->y symmetry)\n', r2fun(Tauy_all, Tauy_model_shared));

save(fullfile(scriptDir, 'centerElectrodeTauxTauy.mat'), 'a1_leq','a2_leq','Fz_model');
