%% checkNewGndSetup.m
% Checks the updated 10-pt rz-only sweep with the new 22-terminal setup
% (bottom Si now terminal 22, outer boundary Zero Charge instead of
% Ground). Compares the resonator's full row sum (S_R) against the
% simplified "Cm_pp + Cm_p,22" (bottom-plate-only) approximation, and
% checks how the resonator's coupling to the bottom plate specifically
% behaves vs rz.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_capacitanceSweepRz10_cmutual_results.txt'), 'NumHeaderLines', 5);
nPose = size(raw,1)/22;
fprintf('Poses: %d, terminals: %d\n', nPose, size(raw,2));

Cm_all = zeros(22,22,nPose);
for p = 1:nPose
    Cm_all(:,:,p) = raw((p-1)*22+1 : p*22, :);
end

% Symmetry check
maxAsym = zeros(nPose,1);
for p = 1:nPose
    M = Cm_all(:,:,p);
    maxAsym(p) = max(abs(M(:) - reshape(M',[],1))) / max(abs(M(:)));
end
fprintf('Max relative asymmetry: %.4e\n', max(maxAsym));

rz = linspace(-14e-6, 14e-6, 10)';  % matches rzOnlySweep10_gndCheck.m

Cmpp = squeeze(Cm_all(21,21,:));       % resonator's own diagonal
Cm_p_bottom = squeeze(Cm_all(21,22,:)); % resonator <-> bottom Si (terminal 22)
S_R_full = squeeze(sum(Cm_all(21,:,:), 2)); % full row sum (all 22 terminals)
S_R_approx = Cmpp + Cm_p_bottom;        % naive bottom-plate-only approximation

fprintf('\n%8s  %12s  %12s  %12s  %12s  %10s\n', 'rz(um)', 'Cm_pp', 'Cm_p,bot', 'S_R(full)', 'approx', 'err%%');
for i = 1:nPose
    err = 100*(S_R_approx(i)-S_R_full(i))/S_R_full(i);
    fprintf('%8.2f  %12.4e  %12.4e  %12.4e  %12.4e  %10.3f\n', ...
        rz(i)*1e6, Cmpp(i), Cm_p_bottom(i), S_R_full(i), S_R_approx(i), err);
end

figure('Color','w','Name','Resonator row sum: full vs bottom-plate-only approx, vs rz');
hold on; grid on;
plot(rz*1e6, S_R_full, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'S_R full row sum');
plot(rz*1e6, S_R_approx, 'r--x', 'LineWidth', 1.5, 'DisplayName', 'Cm_{pp}+Cm_{p,bottom} approx');
xlabel('rz [\mum]'); ylabel('[F]');
title('Resonator row sum: full vs bottom-plate-only approximation');
legend('Location','best');
exportgraphics(gcf, fullfile(scriptDir, 'SR_full_vs_approx.png'));
fprintf('\nSaved plot to SR_full_vs_approx.png\n');
