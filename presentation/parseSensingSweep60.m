%% parseSensingSweep60.m
% Parses the dedicated 60-pose, 22-terminal rz/thx/thy-only sensing
% sweep (comsol_sensingSweep_rzThxThy60_mutual_cap.txt). Same 28-column
% format as the 300-pt withGND sweep (6 pose values on the first row of
% each 22-row block + 22 matrix columns on every row), but only the
% mutual capacitance matrix was exported this time (no elastance/
% Maxwell cross-check available).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sensingSweep_rzThxThy60_mutual_cap.txt'), 'NumHeaderLines', 5);

nTerm = 22;
nPose = size(raw,1) / nTerm;
fprintf('Poses: %d, terminals: %d\n', nPose, nTerm);

poseRows = 1:nTerm:size(raw,1);
rx = raw(poseRows, 1); ry = raw(poseRows, 2); rz = raw(poseRows, 3);
thx = raw(poseRows, 4); thy = raw(poseRows, 5); thz = raw(poseRows, 6);

Cm_all = zeros(nTerm, nTerm, nPose);
for p = 1:nPose
    idx = (p-1)*nTerm+1 : p*nTerm;
    Cm_all(:,:,p) = raw(idx, 7:28);
end

maxAsym = zeros(nPose,1);
for p = 1:nPose
    M = Cm_all(:,:,p);
    maxAsym(p) = max(abs(M(:) - reshape(M',[],1))) / max(abs(M(:)));
end
fprintf('Max relative asymmetry (mutual): %.4e\n', max(maxAsym));

% Confirm the 3 blocks are rz-only (1-20), thx-only (21-40), thy-only (41-60)
fprintf('Block 1 (rows 1-20):  rz range [%.3g,%.3g], thx/thy all zero: %d\n', ...
    min(rz(1:20)), max(rz(1:20)), all(thx(1:20)==0) && all(thy(1:20)==0));
fprintf('Block 2 (rows 21-40): thx range [%.3g,%.3g], rz/thy all zero: %d\n', ...
    min(thx(21:40)), max(thx(21:40)), all(rz(21:40)==0) && all(thy(21:40)==0));
fprintf('Block 3 (rows 41-60): thy range [%.3g,%.3g], rz/thx all zero: %d\n', ...
    min(thy(41:60)), max(thy(41:60)), all(rz(41:60)==0) && all(thx(41:60)==0));

save(fullfile(scriptDir, 'sensingSweep60_allMatrices.mat'), ...
    'Cm_all', 'rx','ry','rz','thx','thy','thz', 'nPose', 'nTerm');
fprintf('\nSaved to sensingSweep60_allMatrices.mat\n');
