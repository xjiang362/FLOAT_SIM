%% polyFitEnergyFx_rxOnly.m
% Fits a polynomial to the raw simulated W(rx) (smoothing out simulation
% noise), then differentiates the polynomial ANALYTICALLY (not a raw
% numerical gradient on noisy data) to get a smoother, more trustworthy
% energy-method Fx(rx) estimate on the dedicated uniform rx-only sweep.
% Picks the polynomial order via leave-one-out cross-validation on W
% itself, so the order is chosen by predictive fit quality, not by eye.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_sweep_rxOnly_uniform_energy_results.txt'), 'NumHeaderLines', 5);
rx = raw(:,1);
Fx_MST = raw(:,7);
W = raw(:,13);
n = numel(rx);

Fx_energy_raw = gradient(W, rx);   % previous noisy raw-gradient estimate, for comparison

%% Pick polynomial order via LOO CV on W
% rx ~1e-5 with powers up to 6 spans ~30 orders of magnitude, which
% wrecks plain polyfit's conditioning -- center/scale rx first (xhat =
% (rx-mu)/sigma), fit in that normalized variable, and convert the
% derivative back to d/drx via the chain rule (divide by sigma).
maxOrder = 6;
muRx = mean(rx); sigmaRx = std(rx);
xhat = (rx - muRx) / sigmaRx;

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
fprintf('%6s %12s %12s\n', 'order', 'in-sample R2', 'LOO R2');
r2_loo_all = zeros(maxOrder,1);
for ord = 1:maxOrder
    designMat = xhat .^ (ord:-1:0);
    coeffs = designMat \ W;
    What = designMat * coeffs;
    r2 = r2fun(W, What);

    H = designMat * ((designMat'*designMat) \ designMat');
    h = diag(H);
    looResid = (W - What) ./ (1 - h);
    r2_loo = 1 - sum(looResid.^2)/sum((W-mean(W)).^2);
    r2_loo_all(ord) = r2_loo;
    fprintf('%6d %12.6f %12.6f\n', ord, r2, r2_loo);
end
[~, bestOrder] = max(r2_loo_all);
fprintf('\nBest order by LOO R^2: %d\n', bestOrder);

%% Fit at best order (in normalized xhat), differentiate analytically
p = polyfit(xhat, W, bestOrder);
dp = polyder(p);
Fx_poly = polyval(dp, xhat) / sigmaRx;   % chain rule: d/drx = (d/dxhat)*(1/sigmaRx)

r2fun2 = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
fprintf('\nFx: poly-derivative energy method vs MST, R^2 = %.4f\n', r2fun2(Fx_MST, Fx_poly));
fprintf('Fx: raw-gradient energy method vs MST, R^2 = %.4f\n', r2fun2(Fx_MST, Fx_energy_raw));

%% Plot
figure('Color','w','Name','Polynomial-smoothed energy-method Fx vs rx', 'Position',[50 50 1400 500]);

subplot(1,2,1); hold on; grid on;
scatter(rx*1e6, W, 30, 'k', 'filled', 'DisplayName', 'Simulated W (raw)');
rxFine = linspace(min(rx), max(rx), 200)';
xhatFine = (rxFine - muRx) / sigmaRx;
plot(rxFine*1e6, polyval(p, xhatFine), 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Poly fit, order %d', bestOrder));
xlabel('rx [\mum]'); ylabel('W [J]');
title('W(rx): raw simulation vs polynomial fit');
legend('Location','best');

subplot(1,2,2); hold on; grid on;
scatter(rx*1e6, Fx_MST, 30, 'k', 'filled', 'DisplayName', 'Simulation (MST)');
scatter(rx*1e6, Fx_energy_raw, 20, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Energy method (raw gradient)');
plot(rxFine*1e6, polyval(dp, xhatFine)/sigmaRx, 'Color', [0 0.6 0], 'LineWidth', 2, 'DisplayName', sprintf('Energy method (poly-derivative, order %d)', bestOrder));
yline(0, 'k--', 'HandleVisibility', 'off');
xlabel('rx [\mum]'); ylabel('Fx [N]');
title('Fx vs rx: MST vs energy-method (raw vs poly-smoothed)');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'polyFitEnergyFx_rxOnly.png'));
fprintf('\nSaved plot to polyFitEnergyFx_rxOnly.png\n');
