%% plotTauzOrder5Fit.m
% Plots the full 5th-order g(rx) fit (all terms x^1..x^5, shared C0)
% against the order-4 fit and the raw energy-method data. Self-contained
% (reloads the raw sweeps directly) to avoid relying on local functions
% scoped to fitTauzModel_kV2GeffSum.m.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

%% --- rx-primary sweep: Tauz(rx) via smoothed thz-bracket differencing ---
raw1 = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_ryThzPairs_results.txt'), 'NumHeaderLines', 5);
rx1_all = raw1(:,1); thz1_all = raw1(:,6);
W1_all = raw1(:,13);
rxBase = rx1_all(1:4:end);
W_thzMinus = W1_all(3:4:end); W_thzPlus = W1_all(4:4:end);
thz_bracket_rad = (thz1_all(4:4:end) - thz1_all(3:4:end)) / 2 * pi/180;
validThz = thz_bracket_rad > 0;

muRx = mean(rxBase); sigmaRx = std(rxBase);
xhatRx = (rxBase - muRx) / sigmaRx;
p_thzMinus = polyfit(xhatRx, W_thzMinus, 6);
p_thzPlus  = polyfit(xhatRx, W_thzPlus, 6);
W_thzMinus_s = polyval(p_thzMinus, xhatRx);
W_thzPlus_s  = polyval(p_thzPlus, xhatRx);
Tauz_rx = nan(numel(rxBase),1);
Tauz_rx(validThz) = (W_thzPlus_s(validThz) - W_thzMinus_s(validThz)) ./ (2*thz_bracket_rad(validThz));
rxFit = rxBase(validThz); TauzFit_rx = Tauz_rx(validThz);

% Exclude |rx|>40um -- the last 1-2 points near each edge show a sharp,
% non-physical reversal (e.g. rx=39.49um: +2.5e-11 -> rx=41.33um:
% -1.2e-10 -> rx=43.16um: -4.0e-10), consistent with meshing breakdown
% near the collision-clearance boundary rather than real signal.
keepRx = abs(rxFit) <= 40e-6;
rxFit = rxFit(keepRx); TauzFit_rx = TauzFit_rx(keepRx);
rxhat = rxFit / sigmaRx;

%% --- thz-primary sweep: Tauz(thz) via smoothed rx-averaged energy derivative ---
raw2 = readmatrix(fullfile(scriptDir, 'comsol_sweep_thzOnly_rxPair_energy_results.txt'), 'NumHeaderLines', 5);
thz2_all = raw2(:,6);
W2_all = raw2(:,13);
thzBase = thz2_all(1:2:end) * pi/180;
W_rxMinus = W2_all(1:2:end); W_rxPlus = W2_all(2:2:end);
W_avg2 = (W_rxMinus + W_rxPlus) / 2;
muThz = mean(thzBase); sigmaThz = std(thzBase);
xhatThz = (thzBase - muThz) / sigmaThz;
p_avg2 = polyfit(xhatThz, W_avg2, 6);
dp_avg2 = polyder(p_avg2);
Tauz_thz_smoothed = polyval(dp_avg2, xhatThz) / sigmaThz;
thzhat = thzBase / sigmaThz;

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
nRx = numel(rxhat); nThz = numel(thzhat);
target = [TauzFit_rx; Tauz_thz_smoothed];

%% --- Order-4 rx + linear thz, shared C0 (the current model) ---
designRx4 = rxhat .^ (1:4);
designJoint4 = [ones(nRx+nThz,1), [designRx4; zeros(nThz,4)], [zeros(nRx,1); thzhat]];
c4 = designJoint4 \ target;
C0_4 = c4(1);
coeffsRx4 = c4(2:5) ./ (sigmaRx.^(1:4))';
g_rx4 = @(x) coeffsRx4(1)*x + coeffsRx4(2)*x.^2 + coeffsRx4(3)*x.^3 + coeffsRx4(4)*x.^4;
r2_rx4 = r2fun(TauzFit_rx, C0_4 + g_rx4(rxFit));
H4 = designJoint4 * ((designJoint4'*designJoint4) \ designJoint4');
loo4 = (target - designJoint4*c4) ./ (1 - diag(H4));
r2loo4 = 1 - sum(loo4.^2)/sum((target-mean(target)).^2);

%% --- Order-5 rx (all terms) + linear thz, shared C0 ---
designRx5 = rxhat .^ (1:5);
designJoint5 = [ones(nRx+nThz,1), [designRx5; zeros(nThz,5)], [zeros(nRx,1); thzhat]];
c5 = designJoint5 \ target;
C0_5 = c5(1);
coeffsRx5 = c5(2:6) ./ (sigmaRx.^(1:5))';
g_rx5 = @(x) coeffsRx5(1)*x + coeffsRx5(2)*x.^2 + coeffsRx5(3)*x.^3 + coeffsRx5(4)*x.^4 + coeffsRx5(5)*x.^5;
r2_rx5 = r2fun(TauzFit_rx, C0_5 + g_rx5(rxFit));
H5 = designJoint5 * ((designJoint5'*designJoint5) \ designJoint5');
loo5 = (target - designJoint5*c5) ./ (1 - diag(H5));
r2loo5 = 1 - sum(loo5.^2)/sum((target-mean(target)).^2);

fprintf('Order-4: C0=%.4e, R^2(rx)=%.4f, joint LOO=%.4f\n', C0_4, r2_rx4, r2loo4);
fprintf('Order-5: C0=%.4e, R^2(rx)=%.4f, joint LOO=%.4f\n', C0_5, r2_rx5, r2loo5);

%% Plot
figure('Color','w','Name','Tauz g(rx): order-4 vs order-5', 'Position',[50 50 800 600]);
hold on; grid on;
scatter(rxFit*1e6, TauzFit_rx, 30, 'k', 'filled', 'DisplayName', 'Energy method (+dW/d\theta_z)');
rxFine = linspace(min(rxFit), max(rxFit), 300)';
plot(rxFine*1e6, C0_4 + g_rx4(rxFine), 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Order-4 (R^2=%.4f, LOO=%.4f)', r2_rx4, r2loo4));
plot(rxFine*1e6, C0_5 + g_rx5(rxFine), 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Order-5, all terms (R^2=%.4f, LOO=%.4f)', r2_rx5, r2loo5));
xlabel('rx [\mum]'); ylabel('Tauz contribution [N*m]');
title('g(rx): order-4 vs full order-5 fit');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'TauzGrx_order4vs5.png'));
fprintf('Saved plot to TauzGrx_order4vs5.png\n');
