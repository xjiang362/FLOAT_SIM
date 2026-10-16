%% checkResonatorDiag_vs_rz.m
% Checks Cm(21,21) (resonator's own mutual self/to-ground term) vs rz
% separately from S_R (full row sum), to see if the resonator's own
% diagonal shows divergence at BOTH rz extremes (as expected if it
% forms a parallel-plate relationship with a ground reference on both
% sides) or just one.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

fid = fopen(fullfile(scriptDir, 'comsol_capacitanceSweep300.txt'));
poseVals = struct();
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);
    if isempty(line), continue; end
    tok = regexp(line, '^(\w+)\s+"([^"]+)"', 'tokens', 'once');
    if isempty(tok), continue; end
    vals = str2double(strsplit(strtrim(tok{2})));
    poseVals.(tok{1}) = vals(:);
end
fclose(fid);
rz = poseVals.rz;

load(fullfile(scriptDir, 'capSweep300_allMatrices.mat'), 'Cm_all', 'nPose');
Cmpp = squeeze(Cm_all(21,21,:));   % resonator's own diagonal
S_R = squeeze(sum(Cm_all(21,:,:), 2));  % full row sum (for comparison)

figure('Color','w','Name','Resonator diagonal vs S_R, both vs rz');
subplot(1,2,1);
scatter(rz*1e6, Cmpp, 20, 'filled', 'MarkerFaceAlpha', 0.5);
grid on; xlabel('rz [\mum]'); ylabel('C_{m,pp} [F]');
title('Resonator own diagonal (Cm_{21,21}) vs rz');

subplot(1,2,2);
scatter(rz*1e6, S_R, 20, 'filled', 'MarkerFaceAlpha', 0.5);
grid on; xlabel('rz [\mum]'); ylabel('S_R [F]');
title('Full row sum (S_R) vs rz');

exportgraphics(gcf, fullfile(scriptDir, 'resonatorDiag_vs_SR_vs_rz.png'));
fprintf('Saved plot to resonatorDiag_vs_SR_vs_rz.png\n');

fprintf('\nCmpp range: [%.4e, %.4e]\n', min(Cmpp), max(Cmpp));
fprintf('Cmpp at most negative rz (5 smallest rz): %s\n', mat2str(Cmpp(rz==min(rz)), 4));
[~, idx] = sort(rz);
fprintf('Cmpp for the 5 most negative rz points: %s\n', mat2str(Cmpp(idx(1:5))', 4));
fprintf('Cmpp for the 5 most positive rz points: %s\n', mat2str(Cmpp(idx(end-4:end))', 4));
