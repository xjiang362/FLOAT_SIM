%% analyzeThzSweep_rxPair.m
% Derives Fx(thz) via the energy method from the thz-primary sweep with
% an rx bracket (comsol_sweep_thzOnly_rxPair_energy_results.txt, 100
% rows = 50 thz values x 2 rows each: rx=-1um, rx=+1um).
%
% Each of the two rx-bracket energy tracks (W at rx=-1um, W at rx=+1um)
% is smoothed with a polynomial (order chosen via LOO CV, on
% centered/scaled thz) BEFORE differencing -- same approach validated
% for the rx-primary sweep, which fixed the raw-gradient noise
% amplification problem there.
%
% SIGN: +dW/dq is used (not the earlier, now-reverted -dW/dq "fix") --
% the user confirmed every conductor in this energy-method model is
% held at fixed potential (no floating bodies), so the standard
% coenergy relation F=+dW/dq is theoretically exact for all 6 DOF:
%   Fx(thz) = +[W(rx=+1um) - W(rx=-1um)] / (2*1um)

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_thzOnly_rxPair_energy_results.txt'), 'NumHeaderLines', 5);
thz_all = raw(:,6);
Fx_all = raw(:,7);
W_all = raw(:,13);

thzBase = thz_all(1:2:end);
W_rxMinus = W_all(1:2:end);   % rx=-1um
W_rxPlus  = W_all(2:2:end);   % rx=+1um
Fx_MST_avg = (Fx_all(1:2:end) + Fx_all(2:2:end)) / 2;

rx_bracket = 1e-6;

muThz = mean(thzBase); sigmaThz = std(thzBase);
xhat = (thzBase - muThz) / sigmaThz;

maxOrder = 6;
[p_rxMinus, ord_rxMinus, r2_rxMinus] = fitBestPoly(xhat, W_rxMinus, maxOrder);
[p_rxPlus, ord_rxPlus, r2_rxPlus] = fitBestPoly(xhat, W_rxPlus, maxOrder);
fprintf('Track fits (order, LOO R^2): rx-minus (%d, %.4f), rx-plus (%d, %.4f)\n', ...
    ord_rxMinus, r2_rxMinus, ord_rxPlus, r2_rxPlus);

W_rxMinus_s = polyval(p_rxMinus, xhat);
W_rxPlus_s  = polyval(p_rxPlus, xhat);

Fx_energy_raw = (W_rxPlus - W_rxMinus) / (2*rx_bracket);           % +dW/dq, unsmoothed
Fx_energy = (W_rxPlus_s - W_rxMinus_s) / (2*rx_bracket);            % +dW/dq, smoothed

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
fprintf('\nR^2 (Fx_energy, raw) vs Fx_MST = %.4f\n', r2fun(Fx_MST_avg, Fx_energy_raw));
fprintf('R^2 (Fx_energy, smoothed) vs Fx_MST = %.4f\n', r2fun(Fx_MST_avg, Fx_energy));
fprintf('corr(Fx_energy_smoothed, Fx_MST) = %.4f\n', corr(Fx_energy, Fx_MST_avg));

%% Plot
figure('Color','w','Name','Fx vs thz (energy method)', 'Position',[50 50 900 600]);
hold on; grid on;
scatter(thzBase, Fx_MST_avg, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(thzBase, Fx_energy_raw, 20, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Energy method (raw)');
scatter(thzBase, Fx_energy, 30, [0 0.6 0], 'filled', 'DisplayName', 'Energy method (smoothed, sign-corrected)');
yline(0, 'k:', 'HandleVisibility', 'off');
xlabel('\theta_z [deg]'); ylabel('Fx [N]');
title(sprintf('Fx vs \\theta_z: energy method vs MST (R^2=%.4f)', r2fun(Fx_MST_avg, Fx_energy)));
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'FxVsThz_energyMethod.png'));
fprintf('\nSaved plot to FxVsThz_energyMethod.png\n');

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
