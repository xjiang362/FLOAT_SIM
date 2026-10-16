%% extractForceFromEnergy_full6D.m
% Extracts force/torque data from the 2000-pt full-6D LHS energy sweep
% (full6D_energyData.mat) two ways:
%   1. Direct MST: Fx_MST...Tauz_MST, the COMSOL force/torque probes
%      already in the dataset.
%   2. Energy method, MODEL-FREE: since these 2000 points are scattered
%      LHS (not a grid), a simple central-difference like the dedicated
%      single-axis sweeps used isn't possible -- there's no clean
%      "next point along rz" neighbor. Instead, for each point, a LOCAL
%      LINEAR regression (W ~ W0 + grad.dq) is fit using its k=20
%      nearest neighbors (found in normalized 6D pose space, so all 6
%      axes contribute comparably despite very different physical
%      scales -- meters vs radians). The regression's fitted gradient
%      IS the local force/torque estimate: Fx=dW/drx, ..., Tauz=dW/dthz.
%      This uses ONLY the raw simulated W data, not any fitted physics
%      model (fringing or otherwise) -- a genuine independent,
%      model-free "energy method" reference for this scattered dataset.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));

n = numel(rz);
thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;
Q = [rx, ry, rz, thx_rad, thy_rad, thz_rad];   % raw physical units (m, m, m, rad, rad, rad)

% Normalize each axis by its own std for NEIGHBOR SELECTION ONLY (keeps
% all 6 axes comparable despite scale differences); the actual local
% regression below uses RAW Q so the fitted gradient comes out directly
% in N / N*m.
Qnorm = Q ./ std(Q, 0, 1);

k = 20;   % neighbors per local fit (>6 needed to solve for 6 gradient components + intercept, with margin for noise)
[nnIdx, ~] = knnsearch(Qnorm, Qnorm, 'K', k+1);   % includes self as neighbor 1

Fx_energy = zeros(n,1); Fy_energy = zeros(n,1); Fz_energy = zeros(n,1);
Taux_energy = zeros(n,1); Tauy_energy = zeros(n,1); Tauz_energy = zeros(n,1);

for i = 1:n
    idx = nnIdx(i, 2:end);   % exclude self
    dQ = Q(idx,:) - Q(i,:);
    dW = W(idx) - W(i);
    designMat = [ones(k,1), dQ];   % local affine model: dW = c0 + grad.dQ
    coeffs = designMat \ dW;
    grad = coeffs(2:end);
    Fx_energy(i) = grad(1); Fy_energy(i) = grad(2); Fz_energy(i) = grad(3);
    Taux_energy(i) = grad(4); Tauy_energy(i) = grad(5); Tauz_energy(i) = grad(6);
end

%% Validate against MST
r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
fprintf('Energy-method (local-regression, model-free) vs MST, R^2:\n');
fprintf('  Fx:   %.4f\n', r2fun(Fx_MST, Fx_energy));
fprintf('  Fy:   %.4f\n', r2fun(Fy_MST, Fy_energy));
fprintf('  Fz:   %.4f\n', r2fun(Fz_MST, Fz_energy));
fprintf('  Taux: %.4f\n', r2fun(Taux_MST, Taux_energy));
fprintf('  Tauy: %.4f\n', r2fun(Tauy_MST, Tauy_energy));
fprintf('  Tauz: %.4f\n', r2fun(Tauz_MST, Tauz_energy));

%% Export
T = table(rx, ry, rz, thx, thy, thz, W, ...
    Fx_MST, Fy_MST, Fz_MST, Taux_MST, Tauy_MST, Tauz_MST, ...
    Fx_energy, Fy_energy, Fz_energy, Taux_energy, Tauy_energy, Tauz_energy);
outFile = fullfile(scriptDir, 'full6D_forceTorque_MSTvsEnergy.csv');
writetable(T, outFile);
fprintf('\nSaved %d rows to %s\n', n, outFile);
