%% loadFull6DEnergyResults.m
% Loads the overnight full-6D LHS sweep (2000 pts), QA's the data, and
% fits a global quadratic (all pairwise interactions) regression model
% for total electric energy W(rx,ry,rz,thx,thy,thz), validated via LOO.
% Differentiating a quadratic analytically gives affine (linear+const)
% expressions for Fx,Fy,Fz,Taux,Tauy,Tauz -- a first-pass global
% coupling model, to check overall data quality/sanity before refining
% with any physics-informed nonlinear terms (e.g. known rz inverse-gap
% behavior).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_overnightFull6DSweep_results.txt'), 'NumHeaderLines', 5);
fprintf('Loaded %d rows x %d cols\n', size(raw,1), size(raw,2));

rx = raw(:,1); ry = raw(:,2); rz = raw(:,3);
thx = raw(:,4); thy = raw(:,5); thz = raw(:,6);
Fx_MST = raw(:,7); Fy_MST = raw(:,8); Fz_MST = raw(:,9);
Taux_MST = raw(:,11); Tauy_MST = raw(:,12); Tauz_MST = raw(:,13);
W = raw(:,16);

%% QA
fprintf('\n--- QA ---\n');
allCols = [rx ry rz thx thy thz Fx_MST Fy_MST Fz_MST Taux_MST Tauy_MST Tauz_MST W];
colNames = {'rx','ry','rz','thx','thy','thz','Fx','Fy','Fz','Taux','Tauy','Tauz','W'};
nBad = 0;
for k = 1:numel(colNames)
    col = allCols(:,k);
    nNan = sum(isnan(col)); nInf = sum(isinf(col));
    if nNan>0 || nInf>0
        fprintf('%6s: %d NaN, %d Inf\n', colNames{k}, nNan, nInf);
        nBad = nBad + nNan + nInf;
    end
end
if nBad==0
    fprintf('No NaN/Inf found in state vars, MST probes, or W.\n');
end

fprintf('\n%6s  %14s  %14s\n', 'var', 'min', 'max');
for k = 1:numel(colNames)
    fprintf('%6s  %14.4g  %14.4g\n', colNames{k}, min(allCols(:,k)), max(allCols(:,k)));
end

% Flag outliers in W via robust z-score (median/MAD)
medW = median(W); madW = mad(W,1);
zW = abs(W - medW) / (1.4826*madW);
nOutlier = sum(zW > 6);
fprintf('\nW outliers (|robust z|>6): %d / %d\n', nOutlier, numel(W));
if nOutlier>0 && nOutlier<20
    idx = find(zW>6);
    fprintf('Outlier rows: %s\n', mat2str(idx'));
    for i = idx'
        fprintf('  row %d: rx=%.3g ry=%.3g rz=%.3g thx=%.3g thy=%.3g thz=%.3g  W=%.4g\n', ...
            i, rx(i), ry(i), rz(i), thx(i), thy(i), thz(i), W(i));
    end
end

save(fullfile(scriptDir, 'full6D_energyData.mat'), 'rx','ry','rz','thx','thy','thz', ...
    'Fx_MST','Fy_MST','Fz_MST','Taux_MST','Tauy_MST','Tauz_MST','W');
fprintf('\nSaved parsed data to full6D_energyData.mat\n');
