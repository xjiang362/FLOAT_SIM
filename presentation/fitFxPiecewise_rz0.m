%% fitFxPiecewise_rz0.m
% Piecewise energy-method Fx at rz=0: separate polynomials either side of
% rx=0, so the fit can represent the CORNER that the overlap-truncation
% physics actually has, instead of smoothing across it the way a single
% global polynomial must.
%
% Energy W must be CONTINUOUS at rx=0 (it is a state function); only its
% DERIVATIVE may jump. So the two branches are tied to a common value at
% the origin:
%
%   W(x) = W0 + sum_k a_k x^k  (x<0)
%        = W0 + sum_k b_k x^k  (x>0),      x = rx/x0, x0 = 40um
%
% One shared W0 enforces continuity automatically. Then
%   Fx = dW/drx = (1/x0) * sum_k k*a_k x^(k-1)   on the left, b on the right
% and Fx(0-) = a_1/x0, Fx(0+) = b_1/x0, whose DIFFERENCE is the force step
% the sharp-edge model predicts.
%
% Order per side chosen by leave-one-out CV on W, as per the project's
% energy-method convention (polyFitEnergyFx_rxOnly.m).
%
% SYMBOLS: V [V] drive; g0 [m] gap at rz=0; dY [m] y-extent of the
%   sub-electrode truncating in x; x0 [m] normalisation half-range;
%   W [J]; Fx [N], plotted in nN.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
eps0 = 8.8541878128e-12; V = 30; g0 = 15e-6; um = 1e-6;
dY_a = 1450*um; x0 = 40*um;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

raw = readmatrix(fullfile(scriptDir,'comsol_rxFineSweepEnergy_results.txt'),'NumHeaderLines',5);
m  = raw(:,3)==0;
rx = raw(m,1); W = raw(m,15);
[rx,o] = sort(rx); W = W(o);
x = rx/x0;
L = x < 0; R = x > 0;                 % x==0 contributes only to W0

design = @(n) [ones(numel(x),1), (x.^(1:n)).*L, (x.^(1:n)).*R];

%% ---- order selection by LOO CV on W ----
fprintf('%6s %8s %14s %12s\n','order','nparam','in-sample R2','LOO R2');
maxN = 4; r2loo = zeros(maxN,1);
for n = 1:maxN
    X = design(n); c = X\W; Wh = X*c;
    h = diag(X*((X'*X)\X'));
    looRes = (W-Wh)./(1-h);
    r2loo(n) = 1 - sum(looRes.^2)/sum((W-mean(W)).^2);
    fprintf('%6d %8d %14.6f %12.6f\n', n, 1+2*n, r2fun(W,Wh), r2loo(n));
end
[~,nBest] = max(r2loo);
fprintf('\nBest per-side order by LOO R^2: %d\n', nBest);

%% ---- fit and differentiate analytically ----
X = design(nBest); c = X\W;
W0 = c(1); a = c(1+(1:nBest)); b = c(1+nBest+(1:nBest));
dpoly = @(co,xx) (1/x0) * sum((1:numel(co)).*co(:)' .* xx.^((1:numel(co))-1), 2);

xg = linspace(-1,1,600)';
FxL = dpoly(a, xg(xg<0)); FxR = dpoly(b, xg(xg>0));

Fx0m = a(1)/x0; Fx0p = b(1)/x0;
Fx_sharp = 0.5*V^2*eps0*dY_a/g0;

fprintf('\n--- piecewise result ---\n');
fprintf('  Fx(0-)  = %7.1f nN      Fx(0+) = %7.1f nN\n', Fx0m*1e9, Fx0p*1e9);
fprintf('  step at rx=0 = %.1f nN   (sharp-edge model predicts %.1f nN)\n', ...
    (Fx0m-Fx0p)*1e9, Fx_sharp*1e9);
fprintf('  Fx at rx=-30um = %7.1f nN,  at rx=+30um = %7.1f nN\n', ...
    dpoly(a,-30*um/x0)*1e9, dpoly(b,30*um/x0)*1e9);
fprintf('  Fx at rx=-40um = %7.1f nN,  at rx=+40um = %7.1f nN\n', ...
    dpoly(a,-1)*1e9, dpoly(b,1)*1e9);
fprintf('  step / sharp-edge prediction = %.2f\n', (Fx0m-Fx0p)/Fx_sharp);

% global poly4 for comparison
mu=mean(rx); sg=std(rx); xh=(rx-mu)/sg;
pg = polyfit(xh,W,4); dpg = polyder(pg);
fprintf('\n  R^2 on W : piecewise(order %d) %.6f   vs   global poly4 %.6f\n', ...
    nBest, r2fun(W, X*c), r2fun(W, polyval(pg,xh)));

%% ---- edge-smeared closed form, for reference ----
C = 2*W/V^2; S = @(r,ww) 0.5*(r - ww*log(2*cosh(r/ww)));
best = struct('sse',Inf);
for ww = logspace(log10(0.5*um),log10(60*um),600)
    Xs = [ones(numel(rx),1), S(rx,ww)]; q = Xs\C;
    sse = sum((C-Xs*q).^2);
    if sse<best.sse, best=struct('sse',sse,'w',ww,'q',q); end
end
k = best.q(2); w = best.w;

%% ---- plot: the ENERGY fit itself, its residuals, and the force ----
% W varies by only ~0.8% of its level across the sweep, so the mean must be
% removed or every candidate fit looks identical.
Wbar = mean(W);
pw_W   = X*c;                                   % piecewise fit of W
pg3    = polyfit(xh,W,3);
glob4_W = polyval(pg,xh); glob3_W = polyval(pg3,xh);

figure('Color','w','Position',[40 40 1400 450]);

subplot(1,3,1); hold on; grid on; box on;
plot(rx*1e6,(W-Wbar)*1e12,'o','MarkerSize',5,'Color',[0.25 0.35 0.55], ...
    'MarkerFaceColor',[0.80 0.85 0.93],'DisplayName','COMSOL W');
plot(rx*1e6,(pw_W-Wbar)*1e12,'-','LineWidth',2.2,'Color',[0.12 0.25 0.60], ...
    'DisplayName',sprintf('piecewise poly%d',nBest));
plot(rx*1e6,(glob4_W-Wbar)*1e12,':','LineWidth',2,'Color',[0.45 0.45 0.45], ...
    'DisplayName','global poly4');
plot(rx*1e6,(glob3_W-Wbar)*1e12,'--','LineWidth',1.4,'Color',[0.85 0.33 0.10], ...
    'DisplayName','global poly3');
xline(0,'k:','HandleVisibility','off');
xlabel('r_x [\mum]'); ylabel('W - mean(W)   [pJ]');
title(sprintf('energy fit  (W \\approx %.4f nJ, swing %.1f pJ)', Wbar*1e9, (max(W)-min(W))*1e12));
legend('Location','southeast','FontSize',8);

subplot(1,3,2); hold on; grid on; box on;
plot(rx*1e6,(W-pw_W)*1e12,'o-','MarkerSize',4,'LineWidth',1.1,'Color',[0.12 0.25 0.60], ...
    'DisplayName',sprintf('piecewise poly%d  (RMS %.3f pJ)',nBest,rms(W-pw_W)*1e12));
plot(rx*1e6,(W-glob4_W)*1e12,':','LineWidth',1.8,'Color',[0.45 0.45 0.45], ...
    'DisplayName',sprintf('global poly4  (RMS %.3f pJ)',rms(W-glob4_W)*1e12));
plot(rx*1e6,(W-glob3_W)*1e12,'--','LineWidth',1.3,'Color',[0.85 0.33 0.10], ...
    'DisplayName',sprintf('global poly3  (RMS %.3f pJ)',rms(W-glob3_W)*1e12));
yline(0,'k-','HandleVisibility','off'); xline(0,'k:','HandleVisibility','off');
xlabel('r_x [\mum]'); ylabel('residual W - fit  [pJ]');
title('energy-fit residuals'); legend('Location','best','FontSize',7.5);

subplot(1,3,3); hold on; grid on; box on;
plot(rx*1e6, gradient(W,rx)*1e9,'.','MarkerSize',9,'Color',[0.80 0.82 0.88], ...
    'DisplayName','raw point-to-point gradient (not used)');
plot(xg*x0*1e6, 0.5*V^2*eps0*dY_a/g0*(xg<0)*1e9,'--','LineWidth',1.6, ...
    'Color',[0.85 0.33 0.10],'DisplayName','closed form, sharp edge');
plot(xg*x0*1e6, 0.5*V^2*k*0.5*(1-tanh(xg*x0/w))*1e9,'-','LineWidth',1.6, ...
    'Color',[0.10 0.55 0.35],'DisplayName',sprintf('closed form + smear (w=%.0f\\mum)',w*1e6));
plot(rx*1e6, polyval(dpg,xh)/sg*1e9,':','LineWidth',2,'Color',[0.45 0.45 0.45], ...
    'DisplayName','global poly4 (for comparison)');
plot(xg(xg<0)*x0*1e6, FxL*1e9,'-','LineWidth',2.6,'Color',[0.12 0.25 0.60], ...
    'DisplayName',sprintf('piecewise poly%d, analytic d/dr_x',nBest));
plot(xg(xg>0)*x0*1e6, FxR*1e9,'-','LineWidth',2.6,'Color',[0.12 0.25 0.60], ...
    'HandleVisibility','off');
plot([0 0],[Fx0p Fx0m]*1e9,'-','LineWidth',2.6,'Color',[0.12 0.25 0.60],'HandleVisibility','off');
yline(0,'k:','HandleVisibility','off'); xline(0,'k:','HandleVisibility','off');
xlabel('r_x  [\mum]'); ylabel('F_x = +dW/dr_x  [nN]');
title('force: analytic derivative of the fit');
legend('Location','northeast','FontSize',7.5); xlim([-42 42]); ylim([-600 800]);

sgtitle(sprintf(['Energy fit and the F_x it implies   (r_z = 0, gap 15 \\mum, ' ...
    'V = 30 V, all other DOF zero)']));
for a = findobj(gcf,'Type','axes')', a.Toolbar.Visible = 'off'; end
outPng = fullfile(scriptDir,'FxVsRx_rz0_piecewise.png');
exportgraphics(gcf,outPng,'Resolution',150);
fprintf('\nSaved plot to %s\n', outPng);
