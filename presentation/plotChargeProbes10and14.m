%% plotChargeProbes10and14.m
% Same style as plotChargeProbes.m, for the newly found Probe10/Probe14
% negative-correlation pair (corr=-0.9434).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
raw = readmatrix(fullfile(scriptDir, 'comsol_overnightFull6DSweep_results.txt'), 'NumHeaderLines', 5);
rz = raw(:,3);
p10 = raw(:,15); p14 = raw(:,19);

figure('Color','w','Name','Charge probes 10 and 14', 'Position',[100 100 1400 500]);

subplot(1,2,1); hold on; grid on;
scatter(rz*1e6, p10, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Probe10');
scatter(rz*1e6, p14, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Probe14');
xlabel('rz [\mum]'); ylabel('Charge [C]');
title('Charge vs rz');
legend('Location','best');

subplot(1,2,2); hold on; grid on;
scatter(p10, p14, 10, 'filled', 'MarkerFaceAlpha', 0.3);
lims = [min(p10) max(p10)];
plot(lims, -lims*max(abs(p14))/max(abs(p10)), 'k--', 'LineWidth', 1.2, 'DisplayName', 'reference slope');
xlabel('Probe10 [C]'); ylabel('Probe14 [C]');
title(sprintf('Probe10 vs Probe14 (corr = %.4f)', corr(p10,p14)));
legend('Location','best');

sgtitle('Charge probes: Probe10 and Probe14');
exportgraphics(gcf, fullfile(scriptDir, 'chargeProbes_10and14.png'));
fprintf('Saved plot to chargeProbes_10and14.png\n');
