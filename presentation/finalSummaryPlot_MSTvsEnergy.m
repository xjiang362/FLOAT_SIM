%% finalSummaryPlot_MSTvsEnergy.m
% Final summary: 6 subplots (Fx,Fy,Fz,Taux,Tauy,Tauz), each showing MST
% and energy-derived (fringing model) values, independently sorted
% low-to-high within each series -- x-axis is just data point rank, not
% a physical variable, so this compares the two distributions' overall
% range/shape rather than point-by-point correspondence.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));
load(fullfile(scriptDir, 'fringingModelForces.mat'));

names = {'Fx','Fy','Fz','Taux','Tauy','Tauz'};
mstData = {Fx_MST, Fy_MST, Fz_MST, Taux_MST, Tauy_MST, Tauz_MST};
energyData = {Fx_energy, Fy_energy, Fz_energy, Taux_energy, Tauy_energy, Tauz_energy};
units = {'N','N','N','N*m','N*m','N*m'};

% Sort by MST order, apply the SAME index order to the paired
% energy-derived values -- preserves point-to-point correspondence,
% unlike sorting each series independently (which only compares overall
% range/shape, not whether energy-derived tracks MST's actual ranking).
figure('Color','w','Name','Final summary: MST vs energy-derived, paired sort by MST', 'Position',[50 50 1700 900]);
for i = 1:6
    [mstSorted, idx] = sort(mstData{i});
    energyPaired = energyData{i}(idx);
    subplot(2,3,i); hold on; grid on;
    plot(mstSorted, 'b-', 'LineWidth', 1.2, 'DisplayName', 'MST');
    plot(energyPaired, 'r.', 'MarkerSize', 4, 'DisplayName', 'energy-derived (paired)');
    xlabel('data point (sorted by MST)'); ylabel(sprintf('%s [%s]', names{i}, units{i}));
    title(names{i});
    legend('Location','best','FontSize',8);
    if any(strcmp(names{i}, {'Fx','Fy','Tauz'}))
        ylim([-1e-7, 1e-7]);
    end
end
sgtitle('MST vs energy-derived, paired points sorted by MST order');
exportgraphics(gcf, fullfile(scriptDir, 'finalSummary_MSTvsEnergy.png'));
fprintf('Saved plot to finalSummary_MSTvsEnergy.png\n');
