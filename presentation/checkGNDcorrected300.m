%% checkGNDcorrected300.m
% Re-checks S_R (resonator row sum) vs rz and the Cm_pp+Cm_p,bottom
% approximation across the full corrected 300-pt, 22-terminal dataset
% (bottom Si now included), to see if the double-divergence pattern
% found in the 10-pt rz-only check holds up across the full 6D sweep.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'capSweep300_withGND_allMatrices.mat'));

Cmpp = squeeze(Cm_all(21,21,:));
Cm_p_bottom = squeeze(Cm_all(21,22,:));
S_R_full = squeeze(sum(Cm_all(21,:,:), 2));
S_R_approx = Cmpp + Cm_p_bottom;

figure('Color','w','Name','S_R full vs approx, 300pt dataset');
subplot(1,2,1); hold on; grid on;
scatter(rz*1e6, S_R_full, 15, 'b', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'S_R full row sum');
scatter(rz*1e6, S_R_approx, 15, 'r', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'Cm_{pp}+Cm_{p,bottom}');
xlabel('rz [\mum]'); ylabel('[F]');
title('S_R: full vs bottom-plate-only approx');
legend('Location','best');

relErr = 100*(S_R_approx - S_R_full)./S_R_full;
subplot(1,2,2); hold on; grid on;
scatter(rz*1e6, relErr, 15, 'filled', 'MarkerFaceAlpha', 0.4);
yline(0,'k--');
xlabel('rz [\mum]'); ylabel('relative error [%]');
title('Approximation error vs rz');

exportgraphics(gcf, fullfile(scriptDir, 'SR_full_vs_approx_300pt.png'));
fprintf('Saved plot to SR_full_vs_approx_300pt.png\n');
fprintf('\ncorr(S_R_full, rz) = %.4f\n', corr(S_R_full, rz));
