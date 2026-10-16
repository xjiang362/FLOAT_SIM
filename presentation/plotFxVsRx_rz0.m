%% plotFxVsRx_rz0.m
% Fx vs rx at rz = 0 only, energy method.
%
% Fx = +dW/drx at fixed potential. Following the established project
% recipe (polyFitEnergyFx_rxOnly.m), W is FITTED WITH A POLYNOMIAL and the
% polynomial differentiated ANALYTICALLY -- a raw point-to-point gradient
% of W amplifies COMSOL's mesh-rebuild noise until the scatter swamps the
% signal.
%   1. centre/scale rx  (xhat = (rx-mu)/sigma) for conditioning
%   2. pick the order by leave-one-out CV on W itself
%   3. polyder, then chain rule: dW/drx = polyval(dp,xhat)/sigma
%
% Data : comsol_rxFineSweepEnergy_results.txt, rz=0 slice (41 pts,
%        rx = -40..+40 um, 2 um step, ry/thx/thy/thz exactly zero).
%
% SYMBOLS: V [V] drive; g0 [m] gap at rz=0; dY [m] y-extent of the
%   sub-electrode that truncates in x; w [m] fitted edge smear length;
%   W [J]; Fx [N], plotted in nN.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
eps0 = 8.8541878128e-12; V = 30; g0 = 15e-6; um = 1e-6;
dY_a = 1450*um;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

raw = readmatrix(fullfile(scriptDir,'comsol_rxFineSweepEnergy_results.txt'),'NumHeaderLines',5);
m  = raw(:,3)==0;
rx = raw(m,1); W = raw(m,15);
[rx,o] = sort(rx); W = W(o);
g = g0;

%% ---- 1-2. order selection by LOO CV on W, in normalised xhat ----
mu = mean(rx); sg = std(rx); xhat = (rx-mu)/sg;
maxOrder = 6; r2loo = zeros(maxOrder,1);
fprintf('%6s %14s %12s\n','order','in-sample R2','LOO R2');
for ord = 1:maxOrder
    X = xhat.^(ord:-1:0);
    c = X\W; Wh = X*c;
    h = diag(X*((X'*X)\X'));
    looRes = (W-Wh)./(1-h);
    r2loo(ord) = 1 - sum(looRes.^2)/sum((W-mean(W)).^2);
    fprintf('%6d %14.6f %12.6f\n', ord, r2fun(W,Wh), r2loo(ord));
end
[~,looOrder] = max(r2loo);
% LOO is near-flat from order 2 up (differences in the 4th decimal), so it
% does not meaningfully discriminate. Use order 4, matching the order
% already established for the energy method elsewhere in this project
% (polyFitEnergyFx_rxOnly.m / compareMovmeanVsPolyfitEnergyFx.m).
bestOrder = 4;
fprintf('\nLOO would pick order %d (R^2 %.6f); using order %d (project convention, LOO R^2 %.6f)\n', ...
    looOrder, r2loo(looOrder), bestOrder, r2loo(bestOrder));

%% ---- 3. fit and differentiate analytically ----
p  = polyfit(xhat, W, bestOrder);
dp = polyder(p);
rg    = linspace(min(rx),max(rx),400)';
xg    = (rg-mu)/sg;
Fx_poly = polyval(dp, xg)/sg;                 % [N] on the fine grid
Fx_at_data = polyval(dp, xhat)/sg;

% raw gradient kept ONLY to show why it is not used
Fx_raw = gradient(W, rx);

%% ---- closed-form references ----
Fx_sharp = 0.5*V^2*eps0*dY_a/g * (rg < 0);
C = 2*W/V^2;
S = @(r,ww) 0.5*(r - ww*log(2*cosh(r/ww)));
wGrid = logspace(log10(0.5*um), log10(60*um), 600);
best = struct('sse',Inf);
for ww = wGrid
    X = [ones(numel(rx),1), S(rx,ww)];
    q = X\C; sse = sum((C-X*q).^2);
    if sse < best.sse, best = struct('sse',sse,'w',ww,'q',q); end
end
k = best.q(2); w = best.w;
Fx_smear = 0.5*V^2*k*0.5*(1 - tanh(rg/w));

fprintf('\npoly-derivative Fx at rx=-30um : %.1f nN\n', polyval(dp,(-30*um-mu)/sg)/sg*1e9);
fprintf('poly-derivative Fx at rx=+30um : %.1f nN\n', polyval(dp,( 30*um-mu)/sg)/sg*1e9);
fprintf('closed form (sharp) plateau    : %.1f nN\n', 0.5*V^2*eps0*dY_a/g*1e9);
fprintf('edge-smeared fit plateau       : %.1f nN  (w = %.1f um)\n', 0.5*V^2*k*1e9, w*1e6);
fprintf('poly vs edge-smeared model, R^2 = %+.4f\n', ...
    r2fun(Fx_at_data, 0.5*V^2*k*0.5*(1-tanh(rx/w))));

%% ---- plot ----
figure('Color','w','Position',[80 80 780 550]); hold on; grid on; box on;
plot(rx*1e6, Fx_raw*1e9,'.','MarkerSize',9,'Color',[0.78 0.80 0.86], ...
    'DisplayName','raw point-to-point gradient (not used)');
plot(rg*1e6, Fx_sharp*1e9,'--','LineWidth',1.7,'Color',[0.85 0.33 0.10], ...
    'DisplayName','closed form, sharp edge');
plot(rg*1e6, Fx_smear*1e9,'-','LineWidth',1.7,'Color',[0.10 0.55 0.35], ...
    'DisplayName',sprintf('closed form + edge smear (w=%.0f\\mum)', w*1e6));
plot(rg*1e6, Fx_poly*1e9,'-','LineWidth',2.4,'Color',[0.12 0.25 0.60], ...
    'DisplayName',sprintf('energy method: poly%d fit of W, analytic d/dr_x',bestOrder));
yline(0,'k:','HandleVisibility','off'); xline(0,'k:','HandleVisibility','off');
xlabel('r_x  [\mum]'); ylabel('F_x = +dW/dr_x  [nN]');
title({'Quadrant pair: lateral force vs lateral displacement', ...
       'r_z = 0 (gap 15 \mum), V = 30 V, all other DOF zero'});
legend('Location','northeast','FontSize',8.5);
xlim([-42 42]);
outPng = fullfile(scriptDir,'FxVsRx_rz0.png');
exportgraphics(gcf,outPng,'Resolution',150);
fprintf('\nSaved plot to %s\n', outPng);
