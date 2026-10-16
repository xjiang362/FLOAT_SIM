%% analyzeTauzVsThz_thzSweep.m
% Derives Tauz(thz) via the energy method from the thz-primary sweep
% (comsol_sweep_thzOnly_rxPair_energy_results.txt): averages the two
% rx-bracket energy tracks (W at rx=-1um, +1um) to approximate
% W(thz,rx=0) to O(rx^2), fits a smooth polynomial to that averaged
% curve vs thz (order via LOO CV), and differentiates analytically:
%   Tauz(thz) = d/dthz of W_avg(thz)
% Checked both sign conventions against MST, given the recurring
% +-dW/dq sign ambiguity found for Fx/Tauz on the rx-primary sweep.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_thzOnly_rxPair_energy_results.txt'), 'NumHeaderLines', 5);
thz_all = raw(:,6);
Tauz_all = raw(:,12);
W_all = raw(:,13);

thzBase = thz_all(1:2:end);
thzBase_rad = thzBase * pi/180;
W_rxMinus = W_all(1:2:end);
W_rxPlus  = W_all(2:2:end);
W_avg = (W_rxMinus + W_rxPlus) / 2;
Tauz_MST_avg = (Tauz_all(1:2:end) + Tauz_all(2:2:end)) / 2;

muThz = mean(thzBase_rad); sigmaThz = std(thzBase_rad);
xhat = (thzBase_rad - muThz) / sigmaThz;

maxOrder = 6;
[p_avg, ord_avg, r2loo_avg] = fitBestPoly(xhat, W_avg, maxOrder);
fprintf('W_avg(thz) fit: order %d, LOO R^2 = %.4f\n', ord_avg, r2loo_avg);

dp_avg = polyder(p_avg);
Tauz_energy_plus = polyval(dp_avg, xhat) / sigmaThz;    % +dW/dthz
Tauz_energy_minus = -Tauz_energy_plus;                  % -dW/dthz

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
fprintf('\nR^2 (+dW/dthz) vs Tauz_MST = %.4f\n', r2fun(Tauz_MST_avg, Tauz_energy_plus));
fprintf('R^2 (-dW/dthz) vs Tauz_MST = %.4f\n', r2fun(Tauz_MST_avg, Tauz_energy_minus));
fprintf('corr(+dW/dthz, Tauz_MST) = %.4f\n', corr(Tauz_energy_plus, Tauz_MST_avg));

r2_plus = r2fun(Tauz_MST_avg, Tauz_energy_plus);
r2_minus = r2fun(Tauz_MST_avg, Tauz_energy_minus);
if r2_plus >= r2_minus
    Tauz_best = Tauz_energy_plus; bestLabel = '+dW/d\theta_z';
else
    Tauz_best = Tauz_energy_minus; bestLabel = '-dW/d\theta_z';
end

%% Plot
figure('Color','w','Name','Tauz vs thz (energy method)', 'Position',[50 50 900 600]);
hold on; grid on;
scatter(thzBase, Tauz_MST_avg, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(thzBase, Tauz_energy_plus, 25, 'b', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', sprintf('+dW/d\\theta_z (R^2=%.4f)', r2_plus));
scatter(thzBase, Tauz_energy_minus, 25, 'r', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', sprintf('-dW/d\\theta_z (R^2=%.4f)', r2_minus));
yline(0, 'k:', 'HandleVisibility', 'off');
xlabel('\theta_z [deg]'); ylabel('Tauz [N*m]');
title('Tauz vs \theta_z: energy method (both signs) vs MST');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'TauzVsThz_energyMethod.png'));
fprintf('\nSaved plot to TauzVsThz_energyMethod.png (better sign: %s)\n', bestLabel);

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
