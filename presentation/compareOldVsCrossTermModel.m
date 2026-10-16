%% compareOldVsCrossTermModel.m
% Compares three fringing-model variants (old: same-axis quadratics;
% cross-term: rz*thx/rz*thy/rx*thz/ry*thz cross terms; coupled-WL:
% cross-term G_eff + shared-coefficient coupled W_eff/L_eff) against
% MST simulation (the direct COMSOL force/torque probes -- Global
% Variable Probes 1,2,3,5,6,7, "Electromagnetic force"/"Torque") on all
% 6 dedicated single-axis sweeps (comsol_sweep_combined300_results.txt),
% each showing that axis's own natural output (Fx-rx, Fy-ry, Fz-rz,
% Taux-thx, Tauy-thy, Tauz-thz).
%
% NOTE: an "energy method" simulation reference (derived from the
% exported capacitance_pF field via W=0.5*C*V^2) was attempted and
% abandoned -- that field (es.Q0_5/es.V0_5) is negative and ~200x too
% small to be electrode 1's own capacitance, so it's a Maxwell mutual
% (off-diagonal) term to some other terminal, not the driven electrode's
% stored-energy-relevant capacitance. This dedicated sweep does not
% carry a usable energy probe (unlike full6D_energyData.mat's 2000-pt
% sweep, which had a genuine dedicated stored-energy probe) -- confirmed
% with the user. MST is the only simulation reference available here.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

data = parseComsolResults(fullfile(scriptDir, 'comsol_sweep_combined300_results.txt'));
V = 30;

axisBlocks = {
    1:50,    'rx',  'Fx',   1e6, 'rx [\mum]'
    51:100,  'ry',  'Fy',   1e6, 'ry [\mum]'
    101:150, 'rz',  'Fz',   1e6, 'rz [\mum]'
    151:200, 'thx', 'Taux', 1,   '\theta_x [deg]'
    201:250, 'thy', 'Tauy', 1,   '\theta_y [deg]'
    251:300, 'thz', 'Tauz', 1,   '\theta_z [deg]'
};

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

figure('Color','w','Name','Old vs cross-term vs coupled-WL fringing model: dedicated single-axis sweeps', 'Position',[50 50 1700 900]);
for i = 1:6
    idx = axisBlocks{i,1};
    dofName = axisBlocks{i,2};
    outName = axisBlocks{i,3};
    xScale = axisBlocks{i,4};
    xLabelStr = axisBlocks{i,5};

    b = data(idx,:);
    n = height(b);
    U = V*ones(n,1);

    [Fx_o,Fy_o,Fz_o,Taux_o,Tauy_o,Tauz_o] = trapElectrode1Model(b.rx,b.ry,b.rz,b.thx,b.thy,b.thz,U);
    [Fx_n,Fy_n,Fz_n,Taux_n,Tauy_n,Tauz_n] = trapElectrode1Model_crossTerms(b.rx,b.ry,b.rz,b.thx,b.thy,b.thz,U);
    [Fx_c,Fy_c,Fz_c,Taux_c,Tauy_c,Tauz_c] = trapElectrode1Model_coupledWL(b.rx,b.ry,b.rz,b.thx,b.thy,b.thz,U);

    switch outName
        case 'Fx',   yOld=Fx_o; yNew=Fx_n; yCpl=Fx_c; ySim_MST=b.Fx;
        case 'Fy',   yOld=Fy_o; yNew=Fy_n; yCpl=Fy_c; ySim_MST=b.Fy;
        case 'Fz',   yOld=Fz_o; yNew=Fz_n; yCpl=Fz_c; ySim_MST=b.Fz;
        case 'Taux', yOld=Taux_o; yNew=Taux_n; yCpl=Taux_c; ySim_MST=b.Taux;
        case 'Tauy', yOld=Tauy_o; yNew=Tauy_n; yCpl=Tauy_c; ySim_MST=b.Tauy;
        case 'Tauz', yOld=Tauz_o; yNew=Tauz_n; yCpl=Tauz_c; ySim_MST=b.Tauz;
    end

    r2_old = r2fun(ySim_MST, yOld);
    r2_new = r2fun(ySim_MST, yNew);
    r2_cpl = r2fun(ySim_MST, yCpl);
    fprintf('%4s vs %3s-only sweep: old R^2=%.4f, cross-term R^2=%.4f, coupled-WL R^2=%.4f\n', outName, dofName, r2_old, r2_new, r2_cpl);

    x = b.(dofName) * xScale;
    subplot(2,3,i); hold on; grid on;
    scatter(x, ySim_MST, 25, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
    scatter(x, yOld, 18, 'b', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', sprintf('Old model (R^2=%.4f)', r2_old));
    scatter(x, yNew, 18, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', sprintf('Cross-term model (R^2=%.4f)', r2_new));
    scatter(x, yCpl, 18, 'r', 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', sprintf('Coupled-WL model (R^2=%.4f)', r2_cpl));
    xlabel(xLabelStr); ylabel([outName ' [SI]']);
    title([outName ' vs ' dofName ' (dedicated sweep)']);
    legend('Location','best');
end

exportgraphics(gcf, fullfile(scriptDir, 'OldVsCrossTerm_dedicatedSweeps.png'));
fprintf('\nSaved plot to OldVsCrossTerm_dedicatedSweeps.png\n');
