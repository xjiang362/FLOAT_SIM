%% fitSymmetricFxModel_rxOnly.m
% Fits W(rx) with an EVEN-only polynomial (only rx^0, rx^2, rx^4 terms --
% no odd powers), enforcing the physical symmetry that flipping rx -> -rx
% should leave the energy unchanged (electrode/system symmetric about
% rx=0). Differentiating an even polynomial gives a pure ODD polynomial
% for Fx: Fx = dW/drx = 2*w2*rx + 4*w4*rx^3 -- a symmetric (odd),
% 3rd-order model that automatically satisfies Fx(0)=0 and
% Fx(-rx)=-Fx(rx), with no free constant or rx^2 term to introduce a
% spurious offset (which is exactly what made the earlier coupled-WL
% model's Fx prediction flat/wrong: it had no mechanism to force
% Fx(rx=0)=0).
%
% Fit against the raw simulated W(rx) (the validated energy-method
% data), on the dedicated uniform rx-only sweep.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_uniform_energy_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1);
Fx_MST = raw(:,7);
W = raw(:,13);

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

% Center/scale for numerical conditioning (same reasoning as before)
muRx = mean(rx); sigmaRx = std(rx);
xhat = (rx - muRx) / sigmaRx;

%% Even-only fit: W = w0 + w2*xhat^2 + w4*xhat^4
designMat = [ones(size(xhat)), xhat.^2, xhat.^4];
w = designMat \ W;
W_sym = designMat * w;
r2_W = r2fun(W, W_sym);
fprintf('Even-only W fit: w0=%.6e, w2=%.6e, w4=%.6e\n', w(1), w(2), w(3));
fprintf('R^2 (W fit) = %.6f\n', r2_W);

% Fx = dW/drx = dW/dxhat * (1/sigmaRx) = (2*w2*xhat + 4*w4*xhat.^3)/sigmaRx
Fx_sym = (2*w(2)*xhat + 4*w(3)*xhat.^3) / sigmaRx;

%% Comparison: plain unconstrained cubic fit (all 4 terms), no symmetry assumed
pCubic = polyfit(xhat, W, 3);
W_cubic = polyval(pCubic, xhat);
r2_Wcubic = r2fun(W, W_cubic);
dpCubic = polyder(pCubic);
Fx_cubic = polyval(dpCubic, xhat) / sigmaRx;
fprintf('\nUnconstrained cubic W fit: R^2 = %.6f\n', r2_Wcubic);
fprintf('Unconstrained-cubic-derived Fx vs MST, R^2 = %.4f\n', r2fun(Fx_MST, Fx_cubic));

%% Reference: energy method smoothed via movmean (previously validated as
% the more trustworthy estimate) and MST, for comparison
W_movmean = movmean(W, 3);
Fx_movmean = gradient(W_movmean, rx);

fprintf('\nSymmetric (odd, 3rd order) Fx model vs MST, R^2 = %.4f\n', r2fun(Fx_MST, Fx_sym));
fprintf('Symmetric (odd, 3rd order) Fx model vs movmean energy method, R^2 = %.4f\n', r2fun(Fx_movmean, Fx_sym));

%% Plot
figure('Color','w','Name','Symmetric (odd, 3rd-order) Fx model vs rx', 'Position',[50 50 1400 500]);

subplot(1,2,1); hold on; grid on;
scatter(rx*1e6, W, 25, 'k', 'filled', 'DisplayName', 'Simulated W (raw)');
rxFine = linspace(min(rx), max(rx), 200)';
xhatFine = (rxFine - muRx)/sigmaRx;
plot(rxFine*1e6, designMatFine(xhatFine)*w, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Even-only fit (R^2=%.4f)', r2_W));
xlabel('rx [\mum]'); ylabel('W [J]');
title('W(rx): raw simulation vs even-only (symmetric) fit');
legend('Location','best');

subplot(1,2,2); hold on; grid on;
scatter(rx*1e6, Fx_MST, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(rx*1e6, Fx_movmean, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', '3-pt movmean energy method');
plot(rxFine*1e6, (2*w(2)*xhatFine + 4*w(3)*xhatFine.^3)/sigmaRx, 'Color', [0 0.6 0], 'LineWidth', 2, ...
    'DisplayName', sprintf('Symmetric odd-cubic Fx model (R^2 vs MST=%.4f)', r2fun(Fx_MST, Fx_sym)));
yline(0, 'k--', 'HandleVisibility', 'off');
xlabel('rx [\mum]'); ylabel('Fx [N]');
title('Fx vs rx: symmetric odd-cubic model vs MST vs energy method');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'SymmetricFxModel_rxOnly.png'));
fprintf('\nSaved plot to SymmetricFxModel_rxOnly.png\n');

function D = designMatFine(xhat)
D = [ones(size(xhat)), xhat.^2, xhat.^4];
end
