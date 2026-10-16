%% plotThzEnergyFits.m
% Plots the raw W(thz=-B(rx)) and W(thz=+B(rx)) energy tracks vs rx,
% along with their individual polynomial fits (order chosen via LOO CV,
% same as analyzeRxSweep_ryThzPairs_smoothed.m), for visual inspection.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_ryThzPairs_results.txt'), 'NumHeaderLines', 5);
rx_all = raw(:,1);
W_all = raw(:,13);

rxBase = rx_all(1:4:end);
W_thzMinus = W_all(3:4:end);
W_thzPlus  = W_all(4:4:end);

muRx = mean(rxBase); sigmaRx = std(rxBase);
xhat = (rxBase - muRx) / sigmaRx;

maxOrder = 6;
[p_thzMinus, ord_thzMinus, r2_thzMinus] = fitBestPoly(xhat, W_thzMinus, maxOrder);
[p_thzPlus, ord_thzPlus, r2_thzPlus] = fitBestPoly(xhat, W_thzPlus, maxOrder);

xhatFine = linspace(min(xhat), max(xhat), 300)';
rxFine = xhatFine*sigmaRx + muRx;

figure('Color','w','Name','thz-bracket energy fits vs rx', 'Position',[50 50 900 600]);
hold on; grid on;
scatter(rxBase*1e6, W_thzMinus, 30, 'b', 'filled', 'DisplayName', 'W(thz=-B(rx)), raw');
plot(rxFine*1e6, polyval(p_thzMinus, xhatFine), 'b-', 'LineWidth', 2, ...
    'DisplayName', sprintf('fit, order %d (LOO R^2=%.4f)', ord_thzMinus, r2_thzMinus));
scatter(rxBase*1e6, W_thzPlus, 30, 'r', 'filled', 'DisplayName', 'W(thz=+B(rx)), raw');
plot(rxFine*1e6, polyval(p_thzPlus, xhatFine), 'r-', 'LineWidth', 2, ...
    'DisplayName', sprintf('fit, order %d (LOO R^2=%.4f)', ord_thzPlus, r2_thzPlus));
xlabel('rx [\mum]'); ylabel('W [J]');
title('Energy vs rx at thz=-B(rx) and thz=+B(rx), with polynomial fits');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'ThzEnergyFits_vsRx.png'));
fprintf('Saved plot to ThzEnergyFits_vsRx.png\n');

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
