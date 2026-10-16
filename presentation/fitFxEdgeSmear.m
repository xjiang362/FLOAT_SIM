%% fitFxEdgeSmear.m
% Closed-form Fx (overlap truncation) vs the ENERGY-METHOD Fx from COMSOL.
%
% The sharp-edge closed form predicts a kink in C at rx=0:
%       dX(rx) = dX0 + min(rx,0)   ->   Fx = (1/2)V^2 * eps0*dY/g  for rx<0
%                                             0                    for rx>0
% Real fringing rounds that corner over a length of order the gap. Model
% that with one extra parameter w [m], the EDGE SMEAR LENGTH:
%       S(rx;w) = 0.5*( rx - w*log(2*cosh(rx/w)) )      -> min(rx,0) as w->0
%       C(rx)   = C0 + k*S(rx;w)
%       Fx(rx)  = (1/2)V^2 * k * 0.5*(1 - tanh(rx/w))
% so k is the asymptotic slope the sharp model predicts as eps0*dY/g, and
% w says how far the fringing field smears the plate edge.
%
% If the sharp-edge picture is right, the fit should return k ~ eps0*dY/g
% and w of order g -- and the apparent "model is too big at large gaps"
% effect is then just the smear eating into a +/-40um sweep.
%
% Energy-method Fx only; the MST Fx probe in this file is non-convergent.
%
% SYMBOLS: eps0 [F/m]; V [V]; g0 [m]; a5 [-]; rx,rz [m]; g [m] gap;
%   W [J]; C [F]; k [F/m] slope of C vs rx; w [m] edge smear length;
%   Fx [N], reported in nN.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
eps0 = 8.8541878128e-12; V = 30; g0 = 15e-6; a5 = 0.976; um = 1e-6;
dY_a = 1450*um;          % sub-electrode a's y-extent (it truncates in x)

raw = readmatrix(fullfile(scriptDir,'comsol_rxFineSweepEnergy_results.txt'),'NumHeaderLines',5);
rx = raw(:,1); rz = raw(:,3); C = 2*raw(:,15)/V^2;
rzVals = unique(rz);

S = @(r,w) 0.5*(r - w*log(2*cosh(r/w)));

fprintf('========= closed-form Fx vs energy-method Fx =========\n');
fprintf('sharp-edge model: Fx_inf = (1/2)V^2 * eps0*dY/g,  dY = %.0f um\n\n', dY_a*1e6);
fprintf('%7s %8s | %11s %11s %7s | %9s %7s\n', ...
    'rz[um]','gap[um]','Fx_sim[nN]','Fx_mdl[nN]','sim/mdl','w[um]','w/g');

res = zeros(numel(rzVals),5);
for i = 1:numel(rzVals)
    m = rz==rzVals(i); r = rx(m); c = C(m); [r,o]=sort(r); c=c(o);
    g = g0 - a5*rzVals(i);

    % grid-search w; A0 and k are linear given w
    wGrid = logspace(log10(0.5*um), log10(60*um), 400);
    best = struct('sse',Inf);
    for w = wGrid
        X = [ones(numel(r),1), S(r,w)];
        p = X\c;
        sse = sum((c - X*p).^2);
        if sse < best.sse, best = struct('sse',sse,'w',w,'p',p); end
    end
    k_fit   = best.p(2);
    Fx_sim  = 0.5*V^2*k_fit;              % asymptotic force, rx -> -inf
    Fx_mdl  = 0.5*V^2*eps0*dY_a/g;
    fprintf('%7.1f %8.2f | %11.1f %11.1f %7.2f | %9.2f %7.2f\n', ...
        rzVals(i)*1e6, g*1e6, Fx_sim*1e9, Fx_mdl*1e9, Fx_sim/Fx_mdl, best.w*1e6, best.w/g);
    res(i,:) = [g, Fx_sim, Fx_mdl, best.w, k_fit];
end

fprintf('\nIf the sharp-edge picture is right: sim/mdl -> 1 and w/g = O(1).\n');
fprintf('correlation of w with g across slices: %.3f\n', corr(res(:,1), res(:,4)));

%% ---- plot ----
figure('Color','w','Position',[40 40 1340 430]);
cols = lines(numel(rzVals));

subplot(1,3,1); hold on; grid on;
for i = 1:numel(rzVals)
    m = rz==rzVals(i); r = rx(m); c = C(m); [r,o]=sort(r); c=c(o);
    g = g0 - a5*rzVals(i);
    plot(r*1e6, (c-max(c))*1e12,'o','MarkerSize',3.5,'Color',cols(i,:), ...
        'DisplayName',sprintf('r_z=%+.0f\\mum',rzVals(i)*1e6));
    % sharp-edge closed form, anchored at the same max
    sharp = eps0*dY_a/g*min(r,0);
    plot(r*1e6, sharp*1e12,'--','LineWidth',1.1,'Color',cols(i,:),'HandleVisibility','off');
end
xlabel('r_x [\mum]'); ylabel('C - C_{max} [pF]');
title({'dots = energy-method C','dashed = sharp-edge closed form'});
legend('Location','southeast','FontSize',7);

subplot(1,3,2); hold on; grid on;
plot(res(:,1)*1e6, res(:,2)*1e9,'o','MarkerSize',8,'MarkerFaceColor',[0.2 0.3 0.75], ...
    'Color',[0.2 0.3 0.75],'DisplayName','energy method');
gg = linspace(2,30,200)'*um;
plot(gg*1e6, 0.5*V^2*eps0*dY_a./gg*1e9,'k-','LineWidth',1.6,'DisplayName','closed form');
set(gca,'XScale','log','YScale','log');
xlabel('gap [\mum]'); ylabel('F_x asymptote [nN]');
title('F_x vs gap'); legend('Location','best','FontSize',8);

subplot(1,3,3); hold on; grid on;
plot(res(:,1)*1e6, res(:,4)*1e6,'s','MarkerSize',9,'MarkerFaceColor',[0.85 0.33 0.10], ...
    'Color',[0.85 0.33 0.10],'DisplayName','fitted smear w');
plot([0 30],[0 30],'k--','LineWidth',1.2,'DisplayName','w = g');
xlabel('gap g [\mum]'); ylabel('edge smear w [\mum]');
title('is the smear length set by the gap?'); legend('Location','best','FontSize',8);

sgtitle('Closed-form F_x vs energy-method F_x, and the edge-smear correction');
outPng = fullfile(scriptDir,'FxEdgeSmear.png');
exportgraphics(gcf,outPng);
fprintf('\nSaved plot to %s\n', outPng);
