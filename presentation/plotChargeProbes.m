%% plotChargeProbes.m
% Plots Probe4 (plate terminal charge) and Probe8 (driven-electrode
% terminal charge): both vs rz (to see the physical gap-dependence
% trend, analogous to W~1/gap), and Probe4 vs Probe8 directly (to show
% the near-anti-correlation / charge-conservation relationship).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
raw = readmatrix(fullfile(scriptDir, 'comsol_overnightFull6DSweep_results.txt'), 'NumHeaderLines', 5);
rz = raw(:,3);
p4 = raw(:,10); p8 = raw(:,14);

figure('Color','w','Name','Charge probes 4 and 8', 'Position',[100 100 1400 500]);

subplot(1,2,1); hold on; grid on;
scatter(rz*1e6, p4, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Probe4 (plate)');
scatter(rz*1e6, p8, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Probe8 (electrode)');
xlabel('rz [\mum]'); ylabel('Charge [C]');
title('Charge vs rz');
legend('Location','best');

subplot(1,2,2); hold on; grid on;
scatter(p4, p8, 10, 'filled', 'MarkerFaceAlpha', 0.3);
lims = [min(p4) max(p4)];
plot(lims, -lims, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Probe8 = -Probe4');
xlabel('Probe4 (plate) [C]'); ylabel('Probe8 (electrode) [C]');
title(sprintf('Probe4 vs Probe8 (corr = %.4f)', corr(p4,p8)));
legend('Location','best');

sgtitle('Charge probes: Probe4 (plate) and Probe8 (driven electrode)');
exportgraphics(gcf, fullfile(scriptDir, 'chargeProbes_4and8.png'));
fprintf('Saved plot to chargeProbes_4and8.png\n');
