%% checkCapSweepMutualData.m
% Sanity checks on the 300-pose mutual capacitance matrix sweep
% (6300x21 = 300 blocks of 21x21 matrices, one per pose point).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_capacitanceSweep300_mutual_cap_results.txt'), 'NumHeaderLines', 5);
fprintf('Loaded %d rows x %d cols\n', size(raw,1), size(raw,2));

nPose = size(raw,1)/21;
fprintf('Implied number of poses: %.2f\n', nPose);

% Reshape into 21x21x300
Cm_all = zeros(21,21,nPose);
for p = 1:nPose
    Cm_all(:,:,p) = raw((p-1)*21+1 : p*21, :);
end

%% Symmetry check across all poses
maxAsym = zeros(nPose,1);
for p = 1:nPose
    M = Cm_all(:,:,p);
    maxAsym(p) = max(abs(M(:) - reshape(M',[],1))) / max(abs(M(:)));
end
fprintf('\nSymmetry: max relative asymmetry across all 300 poses = %.4e\n', max(maxAsym));

%% Resonator row check (row/col 21 should be much larger than the rest)
resonatorDiag = squeeze(Cm_all(21,21,:));
otherDiag = squeeze(Cm_all(1,1,:));  % terminal 1's diagonal, for comparison
fprintf('\nResonator diagonal (terminal 21): min=%.4e max=%.4e median=%.4e\n', ...
    min(resonatorDiag), max(resonatorDiag), median(resonatorDiag));
fprintf('Terminal 1 diagonal (comparison):  min=%.4e max=%.4e median=%.4e\n', ...
    min(otherDiag), max(otherDiag), median(otherDiag));
fprintf('Ratio (median): %.2fx\n', median(resonatorDiag)/median(otherDiag));

%% NaN/Inf/negative-diagonal checks
fprintf('\nNaN count: %d, Inf count: %d\n', sum(isnan(raw(:))), sum(isinf(raw(:))));
diagVals = zeros(21, nPose);
for p = 1:nPose
    diagVals(:,p) = diag(Cm_all(:,:,p));
end
fprintf('Negative diagonal entries (should be 0, self-capacitance must be positive): %d\n', sum(diagVals(:) < 0));

%% Key quantities for the effective-capacitance formula, across all 300 poses
Cm12 = squeeze(Cm_all(1,2,:));
Cm1_21 = squeeze(Cm_all(1,21,:));
Cm2_21 = squeeze(Cm_all(2,21,:));
S_R = squeeze(sum(Cm_all(21,:,:),2));

fprintf('\nCm(1,2):    min=%.4e max=%.4e\n', min(Cm12), max(Cm12));
fprintf('Cm(1,21):   min=%.4e max=%.4e\n', min(Cm1_21), max(Cm1_21));
fprintf('Cm(2,21):   min=%.4e max=%.4e\n', min(Cm2_21), max(Cm2_21));
fprintf('S_R (row sum, resonator): min=%.4e max=%.4e, range/median = %.2f%%\n', ...
    min(S_R), max(S_R), 100*(max(S_R)-min(S_R))/median(S_R));

save(fullfile(scriptDir, 'capSweep300_parsed.mat'), 'Cm_all', 'nPose');
fprintf('\nSaved parsed data to capSweep300_parsed.mat\n');
