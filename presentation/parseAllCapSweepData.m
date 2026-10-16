%% parseAllCapSweepData.m
% Parses all three 300-pose capacitance matrix sweeps (mutual, Maxwell,
% elastance/Cinv) into 21x21x300 arrays, and cross-validates them
% against each other across all 300 poses (Maxwell = row-sum-converted
% mutual; elastance = inv(Maxwell)) as a thorough consistency check.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

files = struct( ...
    'mutual', 'comsol_capacitanceSweep300_mutual_cap_results.txt', ...
    'maxwell', 'comsol_capacitanceSweep300_maxwell_cap_results.txt', ...
    'elastance', 'comsol_capacitanceSweep300_cinv_results.txt');

raw_mutual = readmatrix(fullfile(scriptDir, files.mutual), 'NumHeaderLines', 5);
raw_maxwell = readmatrix(fullfile(scriptDir, files.maxwell), 'NumHeaderLines', 5);
raw_elastance = readmatrix(fullfile(scriptDir, files.elastance), 'NumHeaderLines', 5);

nPose = size(raw_mutual,1)/21;
fprintf('Poses: %d (mutual), %d (maxwell), %d (elastance)\n', ...
    size(raw_mutual,1)/21, size(raw_maxwell,1)/21, size(raw_elastance,1)/21);

Cm_all = zeros(21,21,nPose);
C_all = zeros(21,21,nPose);
S_all = zeros(21,21,nPose);
for p = 1:nPose
    idx = (p-1)*21+1 : p*21;
    Cm_all(:,:,p) = raw_mutual(idx,:);
    C_all(:,:,p) = raw_maxwell(idx,:);
    S_all(:,:,p) = raw_elastance(idx,:);
end

%% Cross-validation 1: Maxwell = row-sum-converted mutual, across all 300 poses
maxErr1 = zeros(nPose,1);
for p = 1:nPose
    Cm = Cm_all(:,:,p);
    Cm_off = Cm - diag(diag(Cm));
    C_from_Cm = diag(Cm*ones(21,1)) - Cm_off;
    Cp = C_all(:,:,p);
    maxErr1(p) = max(abs(C_from_Cm(:) - Cp(:))) / max(abs(Cp(:)));
end
fprintf('\nMaxwell vs row-sum-converted mutual: max relative error across 300 poses = %.4e\n', max(maxErr1));

%% Cross-validation 2: elastance = inv(Maxwell), across all 300 poses
maxErr2 = zeros(nPose,1);
for p = 1:nPose
    S_from_C = inv(C_all(:,:,p));
    Sp = S_all(:,:,p);
    maxErr2(p) = max(abs(S_from_C(:) - Sp(:))) / max(abs(Sp(:)));
end
fprintf('Elastance vs inv(Maxwell): max relative error across 300 poses = %.4e\n', max(maxErr2));

save(fullfile(scriptDir, 'capSweep300_allMatrices.mat'), 'Cm_all', 'C_all', 'S_all', 'nPose');
fprintf('\nSaved all three parsed matrix sets to capSweep300_allMatrices.mat\n');
fprintf('(Cm_all = mutual, C_all = Maxwell, S_all = elastance, each 21x21x%d)\n', nPose);
