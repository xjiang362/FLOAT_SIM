%% compareModelsVsMST_full6D.m
% Compares all three fringing-model variants (old, cross-term,
% coupled-WL) directly against MST simulation on the full 2000-pt
% scattered 6D LHS dataset (full6D_energyData.mat) -- the dataset all
% three were actually FIT on (via the energy method), for all 6
% force/torque components. No energy-derived reference here (per the
% prior finding that local-regression differentiation on this scattered
% dataset is too noisy to trust) -- MST only.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'full6D_energyData.mat'));

n = numel(rz);
U = 30*ones(n,1);

[Fx_o,Fy_o,Fz_o,Taux_o,Tauy_o,Tauz_o] = trapElectrode1Model(rx,ry,rz,thx,thy,thz,U);
[Fx_n,Fy_n,Fz_n,Taux_n,Tauy_n,Tauz_n] = trapElectrode1Model_crossTerms(rx,ry,rz,thx,thy,thz,U);
[Fx_c,Fy_c,Fz_c,Taux_c,Tauy_c,Tauz_c] = trapElectrode1Model_coupledWL(rx,ry,rz,thx,thy,thz,U);

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

outputs = {
    'Fx',   Fx_MST,   Fx_o,   Fx_n,   Fx_c
    'Fy',   Fy_MST,   Fy_o,   Fy_n,   Fy_c
    'Fz',   Fz_MST,   Fz_o,   Fz_n,   Fz_c
    'Taux', Taux_MST, Taux_o, Taux_n, Taux_c
    'Tauy', Tauy_MST, Tauy_o, Tauy_n, Tauy_c
    'Tauz', Tauz_MST, Tauz_o, Tauz_n, Tauz_c
};

fprintf('%-6s %12s %12s %12s\n', 'output', 'old', 'cross-term', 'coupled-WL');
for i = 1:6
    name = outputs{i,1};
    ySim = outputs{i,2};
    r2_old = r2fun(ySim, outputs{i,3});
    r2_new = r2fun(ySim, outputs{i,4});
    r2_cpl = r2fun(ySim, outputs{i,5});
    fprintf('%-6s %12.4f %12.4f %12.4f\n', name, r2_old, r2_new, r2_cpl);
end

%% Plot: predicted vs actual scatter for all 6 outputs, coupled-WL model
figure('Color','w','Name','Fringing models vs MST: full 2000-pt LHS dataset', 'Position',[50 50 1700 900]);
for i = 1:6
    name = outputs{i,1};
    ySim = outputs{i,2};
    yOld = outputs{i,3}; yNew = outputs{i,4}; yCpl = outputs{i,5};
    r2_old = r2fun(ySim, yOld); r2_new = r2fun(ySim, yNew); r2_cpl = r2fun(ySim, yCpl);

    subplot(2,3,i); hold on; grid on;
    scatter(ySim, yOld, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', sprintf('Old (R^2=%.4f)', r2_old));
    scatter(ySim, yNew, 10, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', sprintf('Cross-term (R^2=%.4f)', r2_new));
    scatter(ySim, yCpl, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', sprintf('Coupled-WL (R^2=%.4f)', r2_cpl));
    lims = [min(ySim) max(ySim)];
    plot(lims, lims, 'k--', 'LineWidth', 1);
    xlabel([name ' (MST)']); ylabel([name ' (model)']);
    title(name);
    legend('Location','best', 'FontSize', 7);
end

exportgraphics(gcf, fullfile(scriptDir, 'ModelsVsMST_full6D.png'));
fprintf('\nSaved plot to ModelsVsMST_full6D.png\n');
