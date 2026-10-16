%% fitFxModel_kV2GeffSum.m
% Builds Fx(pose) = [g0z/G_eff(rz,thx,thy)] * (h_own(rx) + h_cross(ry) + h_thz(thz))
% analogous to the Tauz model, using three pieces:
%   h_own(rx)   -- Fx's own-axis dependence, from the rx-primary energy
%                  sweep (already-validated +dW/dq energy method)
%   h_cross(.)  -- the CROSS coupling (how the orthogonal lateral axis
%                  feeds into this force). No dedicated "Fx vs ry" sweep
%                  exists, so by the device's x<->y symmetry, this is
%                  borrowed from the DIRECT MST Fy reading vs rx (same
%                  rx-primary sweep, no derivative needed -- Fy is
%                  already a probed quantity) -- Fy's cross-dependence
%                  on rx equals Fx's cross-dependence on ry.
%   h_thz(thz)  -- Fx's thz dependence, from the thz-primary energy
%                  sweep (already-validated +dW/dq energy method)
%
% Fy is then obtained from Fx by the same x<->y reflection symmetry
% already validated for Taux/Tauy (rx<->ry swap, thz->-thz, no extra
% output sign for Fx/Fy since force is a true vector, unlike
% torque/pseudovector Taux/Tauy which DID need a sign flip):
%   Fy(rx,ry,thz) = Fx(ry,rx,-thz) = h_own(ry) + h_cross(rx) - h_thz(thz)
% (using h_thz odd, h_thz(-thz)=-h_thz(thz))

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'simpleFzModel.mat'), 'pfit');
g0_z = 15e-6;
a1_fz = pfit(2); a3_fz = pfit(3); a5_fz = pfit(4);
Geff_fun = @(rz,thx_rad,thy_rad) g0_z - a5_fz*rz - a1_fz*thx_rad - a3_fz*thy_rad;

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

%% --- Load rx-primary sweep ---
raw1 = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_ryThzPairs_results.txt'), 'NumHeaderLines', 5);
rx1_all = raw1(:,1);
Fx1_all = raw1(:,7); Fy1_all = raw1(:,8);
W1_all = raw1(:,13);
rxBase = rx1_all(1:4:end);
W_ryMinus = W1_all(1:4:end); W_ryPlus = W1_all(2:4:end);
Fx_MST_avg = (Fx1_all(1:4:end) + Fx1_all(2:4:end)) / 2;
Fy_MST_avg = (Fy1_all(1:4:end) + Fy1_all(2:4:end)) / 2;

muRx = mean(rxBase); sigmaRx = std(rxBase);
rxhat = rxBase / sigmaRx;

%% h_own(rx): smoothed energy method (+dW/dq), same as before
xhatRx = (rxBase - muRx) / sigmaRx;
[p_ryMinus, ~, ~] = fitBestPoly(xhatRx, W_ryMinus, 6);
[p_ryPlus, ~, ~] = fitBestPoly(xhatRx, W_ryPlus, 6);
p_avg = polyfit(xhatRx, (W_ryMinus+W_ryPlus)/2, 4);
dp_avg = polyder(p_avg);
Fx_energy = polyval(dp_avg, xhatRx) / sigmaRx;   % +dW/drx, smoothed

%% --- Load thz-primary sweep ---
raw2 = readmatrix(fullfile(scriptDir, 'comsol_sweep_thzOnly_rxPair_energy_results.txt'), 'NumHeaderLines', 5);
thz2_all = raw2(:,6);
W2_all = raw2(:,13);
thzBase = thz2_all(1:2:end) * pi/180;
W_rxMinus = W2_all(1:2:end); W_rxPlus = W2_all(2:2:end);
rx_bracket = 1e-6;

muThz = mean(thzBase); sigmaThz = std(thzBase);
xhatThz = (thzBase - muThz) / sigmaThz;
[p_rxMinus, ~, ~] = fitBestPoly(xhatThz, W_rxMinus, 6);
[p_rxPlus, ~, ~] = fitBestPoly(xhatThz, W_rxPlus, 6);
W_rxMinus_s = polyval(p_rxMinus, xhatThz);
W_rxPlus_s  = polyval(p_rxPlus, xhatThz);
Fx_thz_energy = (W_rxPlus_s - W_rxMinus_s) / (2*rx_bracket);   % +dW/drx, smoothed, vs thz

%% --- JOINT fit: Fx = C0 + h_own(rx) + h_thz(thz), ONE shared baseline C0 ---
% NOT odd -- the raw data shows a real, nonzero baseline Fx at rx=0/
% thz=0 (~1.0-1.4e-7 N from each sweep independently), so forcing
% h_own(0)=h_thz(0)=0 was wrong. A single C0 is fit jointly across BOTH
% sweeps stacked together (not two separate constants, which would
% double-count the baseline), while h_own/h_thz remain pure
% "deviation-from-baseline" terms (no constant of their own).
rxhat = rxBase / sigmaRx;
thzhat = thzBase / sigmaThz;
maxOrder = 5;
bestR2loo = -Inf;
for ordRx = 1:maxOrder
    for ordThz = 1:maxOrder
        designRxPart = rxhat .^ (1:ordRx);
        designThzPart = thzhat .^ (1:ordThz);
        nRx = numel(rxhat); nThz = numel(thzhat);
        designJoint = [ones(nRx+nThz,1), ...
            [designRxPart; zeros(nThz,ordRx)], ...
            [zeros(nRx,ordThz); designThzPart]];
        target = [Fx_energy; Fx_thz_energy];
        coeffsJoint = designJoint \ target;
        yhat = designJoint * coeffsJoint;
        H = designJoint * ((designJoint'*designJoint) \ designJoint');
        h = diag(H);
        looResid = (target - yhat) ./ (1 - h);
        r2loo = 1 - sum(looResid.^2)/sum((target-mean(target)).^2);
        if r2loo > bestR2loo
            bestR2loo = r2loo;
            bestCoeffs = coeffsJoint; bestOrdRx = ordRx; bestOrdThz = ordThz;
        end
    end
end
C0 = bestCoeffs(1);
coeffsOwnHat = bestCoeffs(2:1+bestOrdRx);
coeffsThzHat = bestCoeffs(2+bestOrdRx:end);
coeffsOwn = coeffsOwnHat ./ (sigmaRx.^(1:bestOrdRx))';
coeffsThz = coeffsThzHat ./ (sigmaThz.^(1:bestOrdThz))';
h_own = @(x) polyvalOrders(x, coeffsOwn, 1:bestOrdRx);
h_thz = @(x) polyvalOrders(x, coeffsThz, 1:bestOrdThz);

r2_rx = r2fun(Fx_energy, C0 + h_own(rxBase));
r2_thz = r2fun(Fx_thz_energy, C0 + h_thz(thzBase));
fprintf('Joint fit: C0=%.6e N, rx order=%d, thz order=%d, joint LOO R^2=%.4f\n', C0, bestOrdRx, bestOrdThz, bestR2loo);
fprintf('  h_own(rx) coeffs (x^1..x^%d): %s, R^2 (incl. C0) vs rx data=%.4f\n', bestOrdRx, mat2str(coeffsOwn',4), r2_rx);
fprintf('  h_thz(thz) coeffs (x^1..x^%d): %s, R^2 (incl. C0) vs thz data=%.4f\n', bestOrdThz, mat2str(coeffsThz',4), r2_thz);

%% --- Assemble Fx and Fy models ---
% Fx assumed NOT a function of ry (per user instruction) -- only rx,thz:
%   Fx(rx,thz) = [g0z/G_eff] * (C0 + h_own(rx) + h_thz(thz))
% By x<->y reflection symmetry (rx<->ry, thz->-thz; C0 unchanged since
% Fy(0,0,0)=Fx(0,0,0) by the same symmetry; no extra output sign for
% Fx/Fy since force is a true vector, unlike Taux/Tauy):
%   Fy(ry,thz) = [g0z/G_eff] * (C0 + h_own(ry) - h_thz(thz))
FxModel = @(rx,ry,rz,thx_rad,thy_rad,thz_rad) ...
    (g0_z ./ Geff_fun(rz,thx_rad,thy_rad)) .* (C0 + h_own(rx) + h_thz(thz_rad));
FyModel = @(rx,ry,rz,thx_rad,thy_rad,thz_rad) ...
    (g0_z ./ Geff_fun(rz,thx_rad,thy_rad)) .* (C0 + h_own(ry) - h_thz(thz_rad));

save(fullfile(scriptDir, 'FxFyModel_kV2GeffSum.mat'), 'C0', 'coeffsOwn', 'coeffsThz', ...
    'h_own', 'h_thz', 'FxModel', 'FyModel', 'g0_z', 'a1_fz', 'a3_fz', 'a5_fz');
fprintf('\nSaved model to FxFyModel_kV2GeffSum.mat\n');

%% Plot
figure('Color','w','Name','Fx shape functions: h_own(rx), h_thz', 'Position',[50 50 1000 500]);
subplot(1,2,1); hold on; grid on;
scatter(rxBase*1e6, Fx_energy, 30, 'k', 'filled', 'DisplayName', 'Fx energy method (+dW/drx)');
rxFine = linspace(min(rxBase), max(rxBase), 200)';
plot(rxFine*1e6, C0 + h_own(rxFine), 'r-', 'LineWidth', 2, 'DisplayName', sprintf('C0+h_{own}(rx) fit (R^2=%.4f)', r2_rx));
yline(C0, 'b--', 'C0 (shared baseline)');
xlabel('rx [\mum]'); ylabel('Fx contribution [N]');
title('C0 + h_{own}(rx)'); legend('Location','best');

subplot(1,2,2); hold on; grid on;
scatter(thzBase*180/pi, Fx_thz_energy, 30, 'k', 'filled', 'DisplayName', 'Fx energy method (+dW/drx, vs thz)');
thzFine = linspace(min(thzBase), max(thzBase), 200)';
plot(thzFine*180/pi, C0 + h_thz(thzFine), 'r-', 'LineWidth', 2, 'DisplayName', sprintf('C0+h_{thz}(\\theta_z) fit (R^2=%.4f)', r2_thz));
yline(C0, 'b--', 'C0 (shared baseline)');
xlabel('\theta_z [deg]'); ylabel('Fx contribution [N]');
title('C0 + h_{thz}(\theta_z)'); legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'FxShapeFunctions_h.png'));
fprintf('Saved plot to FxShapeFunctions_h.png\n');

function y = polyvalOrders(x, coeffs, orders)
    y = zeros(size(x));
    for i = 1:numel(orders)
        y = y + coeffs(i) * x.^orders(i);
    end
end

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
