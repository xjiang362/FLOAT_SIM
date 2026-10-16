%% plotSR_vs_rz.m
% Plots S_R = sum_j Cm(21,j) (resonator's mutual row sum) vs rz, using
% the 300-pose capacitance sweep. Pose values come from the
% "Specified combinations" input file (one row per parameter, values
% space-separated in quotes), matched by row order to the capacitance
% matrix blocks (comsol_capacitanceSweep300_mutual_cap_results.txt).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

% Parse the "Specified combinations" pose file to get rz
fid = fopen(fullfile(scriptDir, 'comsol_capacitanceSweep300.txt'));
poseVals = struct();
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);
    if isempty(line), continue; end
    tok = regexp(line, '^(\w+)\s+"([^"]+)"', 'tokens', 'once');
    if isempty(tok), continue; end
    name = tok{1};
    vals = str2double(strsplit(strtrim(tok{2})));
    poseVals.(name) = vals(:);
end
fclose(fid);
rz = poseVals.rz;
fprintf('Parsed %d rz values\n', numel(rz));

% Load the mutual capacitance matrix sweep and compute S_R (resonator row sum)
load(fullfile(scriptDir, 'capSweep300_allMatrices.mat'), 'Cm_all', 'nPose');
fprintf('Loaded %d capacitance matrices\n', nPose);

if numel(rz) ~= nPose
    error('Mismatch: %d pose rows vs %d capacitance matrices', numel(rz), nPose);
end

S_R = squeeze(sum(Cm_all(21,:,:), 2));

figure('Color','w','Name','S_R vs rz');
scatter(rz*1e6, S_R, 20, 'filled', 'MarkerFaceAlpha', 0.5);
grid on;
xlabel('rz [\mum]'); ylabel('S_R = \Sigma_j C_{m,pj} [F]');
title('Resonator row sum (S_R) vs rz');
exportgraphics(gcf, fullfile(scriptDir, 'SR_vs_rz.png'));
fprintf('Saved plot to SR_vs_rz.png\n');

fprintf('\ncorr(S_R, rz) = %.4f\n', corr(S_R, rz));
fprintf('corr(S_R, abs(rz)) = %.4f\n', corr(S_R, abs(rz)));
