%% analyzeRxSweep_ryThzPairs_smoothed.m
% Same as analyzeRxSweep_ryThzPairs.m, but smooths EACH of the four
% energy tracks (W at ry=-1um, ry=+1um, thz=-thzB(rx), thz=+thzB(rx))
% separately across the rx grid BEFORE taking any difference/gradient --
% rather than differencing/differentiating the raw noisy values
% directly, which amplifies noise (the same issue found on the plain
% rx-only sweep). Each track is fit with a polynomial (order chosen via
% LOO CV, on centered/scaled rx for conditioning), then Fx, Fy, Tauz are
% computed from the smooth fitted curves:
%   Fx(rx)   = d/drx of the mean of the two smoothed ry-bracket curves
%   Fy(rx)   = [smoothed W(ry=+1um) - smoothed W(ry=-1um)] / (2*1um)
%   Tauz(rx) = [smoothed W(thz=+B) - smoothed W(thz=-B)] / (2*B_rad(rx))

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_ryThzPairs_results.txt'), 'NumHeaderLines', 5);
rx_all = raw(:,1); thz_all = raw(:,6);
Fx_all = raw(:,7); Fy_all = raw(:,8); Tauz_all = raw(:,12);
W_all = raw(:,13);

nPtsRx = numel(rx_all)/4;
rxBase = rx_all(1:4:end);
W_ryMinus  = W_all(1:4:end);
W_ryPlus   = W_all(2:4:end);
W_thzMinus = W_all(3:4:end);
W_thzPlus  = W_all(4:4:end);

Fx_MST_avg = (Fx_all(1:4:end) + Fx_all(2:4:end)) / 2;
Fy_MST_avg = (Fy_all(1:4:end) + Fy_all(2:4:end)) / 2;
Tauz_MST_avg = (Tauz_all(3:4:end) + Tauz_all(4:4:end)) / 2;

ry_bracket = 1e-6;
thz_bracket_rad = (thz_all(4:4:end) - thz_all(3:4:end)) / 2 * pi/180;

%% Helper: fit each track vs rx with LOO-CV-selected polynomial order
muRx = mean(rxBase); sigmaRx = std(rxBase);
xhat = (rxBase - muRx) / sigmaRx;
r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

maxOrder = 6;
[p_ryMinus, ord_ryMinus, r2_ryMinus] = fitBestPoly(xhat, W_ryMinus, maxOrder);
[p_ryPlus, ord_ryPlus, r2_ryPlus] = fitBestPoly(xhat, W_ryPlus, maxOrder);
[p_thzMinus, ord_thzMinus, r2_thzMinus] = fitBestPoly(xhat, W_thzMinus, maxOrder);
[p_thzPlus, ord_thzPlus, r2_thzPlus] = fitBestPoly(xhat, W_thzPlus, maxOrder);

fprintf('Track fits (order, LOO R^2): ry-minus (%d, %.4f), ry-plus (%d, %.4f), thz-minus (%d, %.4f), thz-plus (%d, %.4f)\n', ...
    ord_ryMinus, r2_ryMinus, ord_ryPlus, r2_ryPlus, ord_thzMinus, r2_thzMinus, ord_thzPlus, r2_thzPlus);

W_ryMinus_s = polyval(p_ryMinus, xhat);
W_ryPlus_s  = polyval(p_ryPlus, xhat);
W_thzMinus_s = polyval(p_thzMinus, xhat);
W_thzPlus_s  = polyval(p_thzPlus, xhat);

%% Fx: derivative of the mean of the two smoothed ry curves
% REVERTED to +dW/dq: the user confirmed every conductor in this
% energy-method model is held at fixed potential (no floating bodies),
% so the standard coenergy relation F=+dW/dq is theoretically exact for
% all 6 DOF, not just the Fz/Taux/Tauy subset originally validated. The
% earlier -dW/drx "sign correction" (which appeared to help on this
% weak/noisy dataset) is treated as a noise artifact, not a real sign
% law, and is no longer applied.
p_avg = polyfit(xhat, (W_ryMinus+W_ryPlus)/2, ord_ryMinus);  % refit mean directly for a clean single derivative
dp_avg = polyder(p_avg);
Fx_energy = polyval(dp_avg, xhat) / sigmaRx;

%% Fy: difference of smoothed ry curves
Fy_energy = (W_ryPlus_s - W_ryMinus_s) / (2*ry_bracket);

%% Tauz: difference of smoothed thz curves (skip zero-bracket endpoints)
% +dW/dq (see note above -- reverted from the earlier sign flip)
validThz = thz_bracket_rad > 0;
Tauz_energy = nan(nPtsRx,1);
Tauz_energy(validThz) = (W_thzPlus_s(validThz) - W_thzMinus_s(validThz)) ./ (2*thz_bracket_rad(validThz));

fprintf('\nR^2 (Fx_energy, smoothed) vs Fx_MST_avg = %.4f\n', r2fun(Fx_MST_avg, Fx_energy));
fprintf('R^2 (Fy_energy, smoothed) vs Fy_MST_avg = %.4f\n', r2fun(Fy_MST_avg, Fy_energy));
fprintf('R^2 (Tauz_energy, smoothed) vs Tauz_MST_avg = %.4f\n', r2fun(Tauz_MST_avg(validThz), Tauz_energy(validThz)));

%% Plots
figure('Color','w','Name','Fx vs rx (smoothed energy method)');
hold on; grid on;
scatter(rxBase*1e6, Fx_MST_avg, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(rxBase*1e6, Fx_energy, 30, [0 0.6 0], 'filled', 'DisplayName', 'Simulation (smoothed energy method)');
xlabel('rx [\mum]'); ylabel('Fx [N]');
title(sprintf('Fx vs rx: smoothed energy method vs MST (R^2=%.4f)', r2fun(Fx_MST_avg, Fx_energy)));
legend('Location','best');
exportgraphics(gcf, fullfile(scriptDir, 'FxVsRx_smoothedEnergy.png'));
fprintf('Saved plot to FxVsRx_smoothedEnergy.png\n');

figure('Color','w','Name','Fy vs rx (smoothed energy method)');
hold on; grid on;
scatter(rxBase*1e6, Fy_MST_avg, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(rxBase*1e6, Fy_energy, 30, [0 0.6 0], 'filled', 'DisplayName', 'Simulation (smoothed energy method)');
xlabel('rx [\mum]'); ylabel('Fy [N]');
title(sprintf('Fy vs rx: smoothed energy method vs MST (R^2=%.4f)', r2fun(Fy_MST_avg, Fy_energy)));
legend('Location','best');
exportgraphics(gcf, fullfile(scriptDir, 'FyVsRx_smoothedEnergy.png'));
fprintf('Saved plot to FyVsRx_smoothedEnergy.png\n');

figure('Color','w','Name','Tauz vs rx (smoothed energy method)');
hold on; grid on;
scatter(rxBase(validThz)*1e6, Tauz_MST_avg(validThz), 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(rxBase(validThz)*1e6, Tauz_energy(validThz), 30, [0 0.6 0], 'filled', 'DisplayName', 'Simulation (smoothed energy method)');
xlabel('rx [\mum]'); ylabel('Tauz [N*m]');
title(sprintf('Tauz vs rx: smoothed energy method vs MST (R^2=%.4f)', r2fun(Tauz_MST_avg(validThz), Tauz_energy(validThz))));
legend('Location','best');
exportgraphics(gcf, fullfile(scriptDir, 'TauzVsRx_smoothedEnergy.png'));
fprintf('Saved plot to TauzVsRx_smoothedEnergy.png\n');

function [pBest, orderBest, r2loo] = fitBestPoly(xhat, y, maxOrder)
    best_r2loo = -Inf; pBest = []; orderBest = 1;
    for ord = 1:maxOrder
        designMat = xhat .^ (ord:-1:0);
        coeffs = designMat \ y;
        yhat = designMat*coeffs;
        H = designMat * ((designMat'*designMat) \ designMat');
        h = diag(H);
        looResid = (y - yhat) ./ (1 - h);
        r2loo = 1 - sum(looResid.^2)/sum((y-mean(y)).^2);
        if r2loo > best_r2loo
            best_r2loo = r2loo; pBest = polyfit(xhat, y, ord); orderBest = ord;
        end
    end
    r2loo = best_r2loo;
end
