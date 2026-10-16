%% fitCenterElectrodeTauxTauy_smoothed.m
% Same target model as fitCenterElectrodeTauxTauy.m (Taux = Fz_model*thx*
% (a2+a1*rz), Tauy analogous, shared a1/a2 by x<->y symmetry), but fixes
% the R^2~0.6 weakness by smoothing the raw MST torque BEFORE combining
% the two datasets, mirroring the fix already used for Fx
% (analyzeFxEnergyMethod_smoothed.m): naive point-by-point combination of
% differently-conditioned sweeps lets per-point mesh noise dominate the
% joint regression. Fix: within each fixed-rz slice, fit torque vs tilt
% with a smooth (odd-only, LOO order-selected) polynomial, evaluate the
% SMOOTHED curve at the sweep's own tilt points, and regress the
% (rz,thx)->a1,a2 model against the smoothed values instead of raw MST.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'centerElectrodeModel.mat'), 'k_fit','a2','a4','a5','g0_z');
k_c = k_fit; a2_fz = a2; a4_fz = a4; a5_fz = a5;   % Fz model coefficients (renamed to avoid clashing with the a2 we're fitting here)

V = 30;
Geff_fun = @(rz,thxr,thyr) g0_z - a5_fz*rz - a2_fz*thxr.^2 - a4_fz*thyr.^2;
Fz_model = @(rz,thxr,thyr) (k_c*V^2) ./ max(Geff_fun(rz,thxr,thyr), 1e-9).^2;

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

%% --- Load original rz=0 sweep (thx-only, thy-only blocks) ---
raw1 = readmatrix(fullfile(scriptDir, 'comsol_centerElectrodeSweep_rzThxThy150_results.txt'), 'NumHeaderLines', 5);
rz1=raw1(:,3); thx1=raw1(:,4); thy1=raw1(:,5); Taux1=raw1(:,10); Tauy1=raw1(:,11);
blkThx0 = 51:100; blkThy0 = 101:150;

%% --- Load rz=+-7um bracket sweep ---
raw2 = readmatrix(fullfile(scriptDir, 'comsol_centerElectrode_rzCrossTerm200_results.txt'), 'NumHeaderLines', 5);
rz2=raw2(:,3); thx2=raw2(:,4); thy2=raw2(:,5); Taux2=raw2(:,10); Tauy2=raw2(:,11);
blkThxBr = 1:100; blkThyBr = 101:200;

candidateOrders = {[1], [1 3]};   % odd-only basis: pure linear, or linear+cubic

%% ==================== Taux(thx), per-rz-slice smoothing ====================
% Slice 1: rz=0 (from raw1)
thxS1 = thx1(blkThx0)*pi/180; TauxS1 = Taux1(blkThx0); rzS1 = zeros(size(thxS1));
% Slice 2/3: rz=-7um, rz=+7um (from raw2, split by sign)
thxBr = thx2(blkThxBr)*pi/180; TauxBr = Taux2(blkThxBr); rzBr = rz2(blkThxBr);
mMinus = rzBr < 0; mPlus = rzBr > 0;
thxS2 = thxBr(mMinus); TauxS2 = TauxBr(mMinus); rzS2 = rzBr(mMinus);
thxS3 = thxBr(mPlus);  TauxS3 = TauxBr(mPlus);  rzS3 = rzBr(mPlus);

fprintf('=== Taux(thx) per-slice smoothing (odd-only poly, LOO order) ===\n');
[TauxS1_sm, ord1] = smoothOdd(thxS1, TauxS1, candidateOrders);
fprintf('  rz=0 slice: LOO-best exponents = %s\n', mat2str(candidateOrders{ord1}));
[TauxS2_sm, ord2] = smoothOdd(thxS2, TauxS2, candidateOrders);
fprintf('  rz=-7um slice: LOO-best exponents = %s\n', mat2str(candidateOrders{ord2}));
[TauxS3_sm, ord3] = smoothOdd(thxS3, TauxS3, candidateOrders);
fprintf('  rz=+7um slice: LOO-best exponents = %s\n', mat2str(candidateOrders{ord3}));

rz_all_x  = [rzS1;  rzS2;  rzS3];
thx_all   = [thxS1; thxS2; thxS3];
Taux_raw_all    = [TauxS1;    TauxS2;    TauxS3];
Taux_smooth_all = [TauxS1_sm; TauxS2_sm; TauxS3_sm];
thy_all_forx = zeros(size(thx_all));   % thy=0 throughout both source sweeps

Fz_at_x = Fz_model(rz_all_x, thx_all, thy_all_forx);
designX = [Fz_at_x.*thx_all, Fz_at_x.*rz_all_x.*thx_all];
cX = designX \ Taux_smooth_all;
a2_leq = cX(1); a1_leq = cX(2);

Taux_model = Fz_at_x .* thx_all .* (a2_leq + a1_leq*rz_all_x);
fprintf('\n=== Taux model: Fz*thx*(a2+a1*rz), fit against SMOOTHED target ===\n');
fprintf('a2=%.6e, a1=%.6e\n', a2_leq, a1_leq);
fprintf('R^2 vs smoothed target = %.4f\n', r2fun(Taux_smooth_all, Taux_model));
fprintf('R^2 vs RAW MST (honest check)  = %.4f\n', r2fun(Taux_raw_all, Taux_model));

%% ==================== Tauy(thy), per-rz-slice smoothing ====================
thyS1 = thy1(blkThy0)*pi/180; TauyS1 = Tauy1(blkThy0); rzS1y = zeros(size(thyS1));
thyBr = thy2(blkThyBr)*pi/180; TauyBr = Tauy2(blkThyBr); rzBry = rz2(blkThyBr);
mMinus2 = rzBry < 0; mPlus2 = rzBry > 0;
thyS2 = thyBr(mMinus2); TauyS2 = TauyBr(mMinus2); rzS2y = rzBry(mMinus2);
thyS3 = thyBr(mPlus2);  TauyS3 = TauyBr(mPlus2);  rzS3y = rzBry(mPlus2);

fprintf('\n=== Tauy(thy) per-slice smoothing (odd-only poly, LOO order) ===\n');
[TauyS1_sm, ordy1] = smoothOdd(thyS1, TauyS1, candidateOrders);
fprintf('  rz=0 slice: LOO-best exponents = %s\n', mat2str(candidateOrders{ordy1}));
[TauyS2_sm, ordy2] = smoothOdd(thyS2, TauyS2, candidateOrders);
fprintf('  rz=-7um slice: LOO-best exponents = %s\n', mat2str(candidateOrders{ordy2}));
[TauyS3_sm, ordy3] = smoothOdd(thyS3, TauyS3, candidateOrders);
fprintf('  rz=+7um slice: LOO-best exponents = %s\n', mat2str(candidateOrders{ordy3}));

rz_all_y  = [rzS1y;  rzS2y;  rzS3y];
thy_all   = [thyS1; thyS2; thyS3];
Tauy_raw_all    = [TauyS1;    TauyS2;    TauyS3];
Tauy_smooth_all = [TauyS1_sm; TauyS2_sm; TauyS3_sm];
thx_all_fory = zeros(size(thy_all));

Fz_at_y = Fz_model(rz_all_y, thx_all_fory, thy_all);
designY = [Fz_at_y.*thy_all, Fz_at_y.*rz_all_y.*thy_all];
cY = designY \ Tauy_smooth_all;
a2_leqy = cY(1); a1_leqy = cY(2);

Tauy_model_own = Fz_at_y .* thy_all .* (a2_leqy + a1_leqy*rz_all_y);
fprintf('\n=== Tauy model (own fit): Fz*thy*(a2+a1*rz), fit against SMOOTHED target ===\n');
fprintf('a2=%.6e, a1=%.6e\n', a2_leqy, a1_leqy);
fprintf('R^2 vs smoothed target = %.4f\n', r2fun(Tauy_smooth_all, Tauy_model_own));
fprintf('R^2 vs RAW MST (honest check)  = %.4f\n', r2fun(Tauy_raw_all, Tauy_model_own));

Tauy_model_shared = Fz_at_y .* thy_all .* (a2_leq + a1_leq*rz_all_y);
fprintf('\nUsing Taux''s shared (a1,a2) on Tauy data:\n');
fprintf('R^2 vs smoothed target = %.4f\nR^2 vs RAW MST = %.4f  (check for x<->y symmetry)\n', ...
    r2fun(Tauy_smooth_all, Tauy_model_shared), r2fun(Tauy_raw_all, Tauy_model_shared));

save(fullfile(scriptDir, 'centerElectrodeTauxTauy_smoothed.mat'), 'a1_leq','a2_leq','Fz_model');

%% ==================== local function ====================
function [ySmooth, bestIdx] = smoothOdd(x, y, candidateOrders)
    % LOO-selects among candidateOrders (each a vector of odd exponents,
    % e.g. [1] or [1 3]), fits x.^exponents \ y with NO intercept (torque
    % is odd/zero-baseline by the symmetric-electrode argument already
    % confirmed for this electrode), returns the smoothed curve
    % evaluated at the input x's own points.
    n = numel(x);
    xScale = max(abs(x));
    xn = x / xScale;
    r2loo = zeros(numel(candidateOrders),1);
    for ci = 1:numel(candidateOrders)
        exps = candidateOrders{ci};
        yLOO = zeros(n,1);
        for i = 1:n
            keep = true(n,1); keep(i) = false;
            Xb = xn(keep).^exps;
            c = Xb \ y(keep);
            yLOO(i) = (xn(i).^exps) * c;
        end
        r2loo(ci) = 1 - sum((y-yLOO).^2)/sum((y-mean(y)).^2);
    end
    [~, bestIdx] = max(r2loo);
    exps = candidateOrders{bestIdx};
    Xb = xn.^exps;
    c = Xb \ y;
    ySmooth = Xb * c;
end
