%% parseCapSweep300_withGND.m
% Parses the corrected 300-pose, 22-terminal capacitance sweep (with
% bottom Si properly included). Format: 28 columns = 6 pose values
% (only populated on the first row of each 22-row pose block) + 22
% matrix columns (populated on every row). Parses all three matrix
% types (mutual, elastance, Maxwell) and cross-validates them.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

files = struct( ...
    'mutual', 'comsol_capacitanceSweep300_withGND_term_mutual_cap_results.txt', ...
    'elastance', 'comsol_capacitanceSweep300_withGND_term_inv_cap_results.txt', ...
    'maxwell', 'comsol_capacitanceSweep300_withGND_term_maxwell_cap_results.txt');

raw_mutual = readmatrix(fullfile(scriptDir, files.mutual), 'NumHeaderLines', 5);
raw_elastance = readmatrix(fullfile(scriptDir, files.elastance), 'NumHeaderLines', 5);
raw_maxwell = readmatrix(fullfile(scriptDir, files.maxwell), 'NumHeaderLines', 5);

nTerm = 22;
nPose = size(raw_mutual,1) / nTerm;
fprintf('Poses: %d, terminals: %d\n', nPose, nTerm);

% Pose values: first row of each block, columns 1-6
poseRows = 1:nTerm:size(raw_mutual,1);
rx = raw_mutual(poseRows, 1); ry = raw_mutual(poseRows, 2); rz = raw_mutual(poseRows, 3);
thx = raw_mutual(poseRows, 4); thy = raw_mutual(poseRows, 5); thz = raw_mutual(poseRows, 6);
fprintf('Pose ranges: rz [%.3g, %.3g], thx [%.3g, %.3g]\n', min(rz), max(rz), min(thx), max(thx));

Cm_all = zeros(nTerm, nTerm, nPose);
S_all = zeros(nTerm, nTerm, nPose);
C_all = zeros(nTerm, nTerm, nPose);
for p = 1:nPose
    idx = (p-1)*nTerm+1 : p*nTerm;
    Cm_all(:,:,p) = raw_mutual(idx, 7:28);
    S_all(:,:,p) = raw_elastance(idx, 7:28);
    C_all(:,:,p) = raw_maxwell(idx, 7:28);
end

%% Sanity checks
maxAsym = zeros(nPose,1);
for p = 1:nPose
    M = Cm_all(:,:,p);
    maxAsym(p) = max(abs(M(:) - reshape(M',[],1))) / max(abs(M(:)));
end
fprintf('Max relative asymmetry (mutual): %.4e\n', max(maxAsym));

maxErrMaxwell = zeros(nPose,1);
for p = 1:nPose
    Cm = Cm_all(:,:,p);
    Cm_off = Cm - diag(diag(Cm));
    C_from_Cm = diag(Cm*ones(nTerm,1)) - Cm_off;
    Cp = C_all(:,:,p);
    maxErrMaxwell(p) = max(abs(C_from_Cm(:) - Cp(:))) / max(abs(Cp(:)));
end
fprintf('Max relative error (Maxwell vs row-sum-converted mutual): %.4e\n', max(maxErrMaxwell));

save(fullfile(scriptDir, 'capSweep300_withGND_allMatrices.mat'), ...
    'Cm_all', 'S_all', 'C_all', 'rx','ry','rz','thx','thy','thz', 'nPose', 'nTerm');
fprintf('\nSaved to capSweep300_withGND_allMatrices.mat\n');
