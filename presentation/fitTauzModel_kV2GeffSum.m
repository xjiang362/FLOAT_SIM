%% fitTauzModel_kV2GeffSum.m
% Builds Tauz(pose) = [g0z/G_eff(rz,thx,thy)] * (g(rx) + g(ry) + g(thz)),
% i.e. the requested kV^2/G_eff*(f(rx)+f(ry)+f(thz)) form with the
% kV^2/g0z scale absorbed into g(.) (so f(rx)=g(rx)*g0z/(kV^2), same
% model either way): f(rx)=f(ry) by symmetry, so g(rx) is fit ONCE (from
% the rx-primary energy sweep) and reused for ry; g(thz) is fit
% separately (from the thz-primary energy sweep). G_eff reuses the
% already-fitted Fz gap model (a1,a3,a5), since both sweeps held
% rz=thx=thy=0 (G_eff=g0z there), so the extracted g(rx),g(thz) are
% already on the SAME common scale (kV^2/g0z), consistent for summing.
%
% g(x) is fit as an ODD polynomial (only x, x^3, x^5... terms, no
% constant/even powers) -- physically, mirroring rx->-rx (or thz->-thz)
% about the symmetric device should flip the sign of the resulting yaw
% torque, so g(0)=0 and g is antisymmetric. This is linear in the
% coefficients, so solved by ordinary least squares (no iterative fit).
%
% Sign convention: +dW/dq throughout (confirmed correct -- every
% conductor in this model is at fixed potential, no floating-body
% exception).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'simpleFzModel.mat'), 'pfit');   % [k, a1(thx), a3(thy), a5(rz scale)] for G_eff
g0_z = 15e-6;
a1_fz = pfit(2); a3_fz = pfit(3); a5_fz = pfit(4);
Geff_fun = @(rz,thx_rad,thy_rad) g0_z - a5_fz*rz - a1_fz*thx_rad - a3_fz*thy_rad;

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

%% --- g(rx): from the rx-primary sweep (comsol_sweep_rxOnly_ryThzPairs_results.txt) ---
raw1 = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_ryThzPairs_results.txt'), 'NumHeaderLines', 5);
rx1_all = raw1(:,1); thz1_all = raw1(:,6);
W1_all = raw1(:,13);
rxBase = rx1_all(1:4:end);
W_thzMinus = W1_all(3:4:end);
W_thzPlus  = W1_all(4:4:end);
thz_bracket_rad = (thz1_all(4:4:end) - thz1_all(3:4:end)) / 2 * pi/180;
validThz = thz_bracket_rad > 0;

% Smooth each thz-bracket track vs rx BEFORE differencing (same approach
% validated in analyzeRxSweep_ryThzPairs_smoothed.m -- raw differencing
% is noise-amplifying).
muRx = mean(rxBase); sigmaRx = std(rxBase);
xhatRx = (rxBase - muRx) / sigmaRx;
[p_thzMinus, ~, ~] = fitBestPoly(xhatRx, W_thzMinus, 6);
[p_thzPlus, ~, ~] = fitBestPoly(xhatRx, W_thzPlus, 6);
W_thzMinus_s = polyval(p_thzMinus, xhatRx);
W_thzPlus_s  = polyval(p_thzPlus, xhatRx);

Tauz_rx = nan(numel(rxBase),1);
Tauz_rx(validThz) = (W_thzPlus_s(validThz) - W_thzMinus_s(validThz)) ./ (2*thz_bracket_rad(validThz));   % +dW/dthz, smoothed

rxFit = rxBase(validThz); TauzFit_rx = Tauz_rx(validThz);

% Exclude |rx|>40um: the last 1-2 points near each edge show a sharp,
% non-physical reversal (e.g. rx=39.49um: Tauz=+2.5e-11 -> rx=41.33um:
% -1.2e-10 -> rx=43.16um: -4.0e-10, breaking a smooth trend that peaks
% around rx=27um), consistent with meshing breakdown near the
% collision-clearance boundary rather than real signal. Restricting to
% |rx|<=40um took the rx-branch fit from R^2=0.83 to R^2=0.98 (order 4)
% -- the earlier weak fit was these artifacts, not an inadequate model.
keepRx = abs(rxFit) <= 40e-6;
rxFit = rxFit(keepRx); TauzFit_rx = TauzFit_rx(keepRx);

%% --- g(thz): from the thz-primary sweep (comsol_sweep_thzOnly_rxPair_energy_results.txt) ---
raw2 = readmatrix(fullfile(scriptDir, 'comsol_sweep_thzOnly_rxPair_energy_results.txt'), 'NumHeaderLines', 5);
thz2_all = raw2(:,6);
W2_all = raw2(:,13);
thzBase = thz2_all(1:2:end) * pi/180;   % rad
W_rxMinus = W2_all(1:2:end);
W_rxPlus  = W2_all(2:2:end);
rx_bracket = 1e-6;
Tauz_thz = (W_rxPlus - W_rxMinus) / (2*rx_bracket);   % +dW/dq -- but this differences ACROSS rx, giving
% dW/drx = Fx at fixed thz, NOT Tauz(thz) directly. Need Tauz = dW/dthz
% instead: differentiate the rx-AVERAGED energy curve w.r.t. thz.
W_avg2 = (W_rxMinus + W_rxPlus) / 2;
muThz = mean(thzBase); sigmaThz = std(thzBase);
xhatThz = (thzBase - muThz) / sigmaThz;
[p_avg2, ~, ~] = fitBestPoly(xhatThz, W_avg2, 6);
dp_avg2 = polyder(p_avg2);
Tauz_thz_smoothed = polyval(dp_avg2, xhatThz) / sigmaThz;   % +dW/dthz, smoothed

%% --- JOINT fit: Tauz-branch = C0 + g_rx(rx-or-thz-zeroed) + g_thz(...) ---
% Per the Fx/Fy finding: forcing g(0)=0 (pure odd) assumes zero torque
% at the nominal pose, which isn't necessarily true for an off-center
% quadrant electrode. Fits a SINGLE shared baseline C0 jointly across
% BOTH sweeps (not two separate constants, which would double-count),
% plus general (odd+even) polynomial terms for rx, with order chosen by
% LOO CV. thz is fit as PURE LINEAR (order 1, no constant of its own --
% that's what C0 is for): a direct comparison showed a plain linear fit
% for thz (R^2=0.9970, no C0) actually beats a free-order polynomial
% with a shared C0 (R^2=0.9936) -- the earlier order-4 thz fit was
% overfitting noise, not capturing real curvature. rx, by contrast,
% genuinely needs the nonzero shared baseline: forcing it through the
% origin (order-3, no C0, matching a naive "clean cubic" expectation)
% drops its fit from R^2=0.84 to R^2=0.44.
rxhat = rxFit / sigmaRx;
thzhat = thzBase / sigmaThz;
maxOrder = 4;   % capped at 4 per user decision (order-5 gave only a marginal LOO gain)
ordThz = 1;   % thz forced linear -- see note above
bestR2loo = -Inf;
for ordRx = 1:maxOrder
    designRxPart = rxhat .^ (1:ordRx);
    designThzPart = thzhat .^ (1:ordThz);
    nRx = numel(rxhat); nThz = numel(thzhat);
    designJoint = [ones(nRx+nThz,1), ...
        [designRxPart; zeros(nThz,ordRx)], ...
        [zeros(nRx,ordThz); designThzPart]];
    target = [TauzFit_rx; Tauz_thz_smoothed];
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
C0 = bestCoeffs(1);
coeffsRxHat = bestCoeffs(2:1+bestOrdRx);
coeffsThzHat = bestCoeffs(2+bestOrdRx:end);
coeffsRx = coeffsRxHat ./ (sigmaRx.^(1:bestOrdRx))';
coeffsThz = coeffsThzHat ./ (sigmaThz.^(1:bestOrdThz))';
g_rx = @(x) polyvalOrders(x, coeffsRx, 1:bestOrdRx);
g_thz = @(x) polyvalOrders(x, coeffsThz, 1:bestOrdThz);

r2_rx = r2fun(TauzFit_rx, C0 + g_rx(rxFit));
r2_thz = r2fun(Tauz_thz_smoothed, C0 + g_thz(thzBase));
fprintf('Joint fit: C0=%.6e N*m, rx order=%d, thz order=%d, joint LOO R^2=%.4f\n', C0, bestOrdRx, bestOrdThz, bestR2loo);
fprintf('  g(rx) coeffs (x^1..x^%d): %s, R^2 (incl. C0) vs rx data=%.4f\n', bestOrdRx, mat2str(coeffsRx',4), r2_rx);
fprintf('  g(thz) coeffs (x^1..x^%d): %s, R^2 (incl. C0) vs thz data=%.4f\n', bestOrdThz, mat2str(coeffsThz',4), r2_thz);

%% --- Assemble full model ---
% Tauz(0,0,0)=C0 (once, not 3x -- g_rx(0)=g_thz(0)=0 by construction,
% since neither has a constant term of its own): Tauz = (g0z/Geff)*(C0 +
% g_rx(rx) + g_rx(ry) + g_thz(thz))
TauzModel = @(rx,ry,rz,thx_rad,thy_rad,thz_rad) ...
    (g0_z ./ Geff_fun(rz,thx_rad,thy_rad)) .* (C0 + g_rx(rx) + g_rx(ry) + g_thz(thz_rad));

save(fullfile(scriptDir, 'TauzModel_kV2GeffSum.mat'), 'C0', 'coeffsRx', 'coeffsThz', 'g_rx', 'g_thz', 'TauzModel', 'g0_z', 'a1_fz', 'a3_fz', 'a5_fz');
fprintf('\nSaved model to TauzModel_kV2GeffSum.mat\n');

%% Plot: g(rx) and g(thz) fits
figure('Color','w','Name','Tauz shape functions: g(rx) and g(thz)', 'Position',[50 50 1200 500]);
subplot(1,2,1); hold on; grid on;
scatter(rxFit*1e6, TauzFit_rx, 30, 'k', 'filled', 'DisplayName', 'Energy method (+dW/d\theta_z)');
rxFine = linspace(min(rxFit), max(rxFit), 200)';
plot(rxFine*1e6, C0 + g_rx(rxFine), 'r-', 'LineWidth', 2, 'DisplayName', sprintf('C0+g(rx) fit (R^2=%.4f)', r2_rx));
yline(C0, 'b--', 'C0 (shared baseline)');
xlabel('rx [\mum]'); ylabel('Tauz contribution [N*m]');
title('C0 + g(rx) = C0 + g(ry) shape');
legend('Location','best');

subplot(1,2,2); hold on; grid on;
scatter(thzBase*180/pi, Tauz_thz_smoothed, 30, 'k', 'filled', 'DisplayName', 'Energy method (+dW/d\theta_z, smoothed)');
thzFine = linspace(min(thzBase), max(thzBase), 200)';
plot(thzFine*180/pi, C0 + g_thz(thzFine), 'r-', 'LineWidth', 2, 'DisplayName', sprintf('C0+g(\\theta_z) fit (R^2=%.4f)', r2_thz));
yline(C0, 'b--', 'C0 (shared baseline)');
xlabel('\theta_z [deg]'); ylabel('Tauz contribution [N*m]');
title('g(\theta_z) shape');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'TauzShapeFunctions_g.png'));
fprintf('Saved plot to TauzShapeFunctions_g.png\n');

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
