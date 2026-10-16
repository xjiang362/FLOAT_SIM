%% testOverlapModelFx.m
% Does making dX,dY the PLATE-ELECTRODE OVERLAP (rather than the
% electrode's own extent) produce the lateral force in closed form?
%
%   C(rx,ry) ~= eps0 * dX(rx)*dY(ry) / g,     Fx = +(1/2)V^2 dC/drx
%
% The quadrant-1 pair is (per trapGeometry.m):
%   sub a: x in [1525,2975], y in [  25,1475]   <- reaches the plate edge in X
%   sub b: x in [  25,1475], y in [1525,2975]   <- reaches the plate edge in Y
% and the proof mass half-width is 2975 um, so the electrodes sit FLUSH
% with the plate edge. Displacing the plate by rx<0 therefore uncovers
% sub a's outer edge and shrinks its overlap, giving
%       dX_a(rx) = 1450um + min(rx,0),   dX_b = 1450um (never truncated)
% => a ONE-SIDED, piecewise-linear C(rx): flat for rx>0, falling for rx<0,
%    hence a STEP in Fx at rx=0.
%
% Reference data: comsol_rxFineSweepEnergy_results.txt -- 41 uniform rx
% points (2 um step) at 5 fixed rz slices, with a converged total electric
% energy probe. MST Fx is NOT used: it is non-convergent under mesh
% refinement for this component (sign flips), per rxFineSweepEnergy.m.
%
% SYMBOLS: eps0 [F/m]; V [V]; g0 [m]; a5 [-]; rx,ry,rz [m];
%   W [J] total electric energy; C [F]; Fx [N]; dX,dY [m] overlap extents.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
eps0 = 8.8541878128e-12; V = 30; g0 = 15e-6; a5 = 0.976; um = 1e-6;

raw = readmatrix(fullfile(scriptDir,'comsol_rxFineSweepEnergy_results.txt'),'NumHeaderLines',5);
% 15 cols here: rx,ry,rz,thx,thy,thz, Fx,Fy,Fz, probe4, Taux,Tauy,Tauz,
% probe10, TotalElectricEnergy  (the full-6D export has one extra probe)
rx = raw(:,1); rz = raw(:,3); W = raw(:,15);
C  = 2*W/V^2;
rzVals = unique(rz);
fprintf('Loaded %d rows, %d rz slices: %s um\n\n', numel(rx), numel(rzVals), mat2str(rzVals'*1e6,3));

%% ---- overlap model prediction ----
plateHalf = 2975*um;
subA = struct('x',[1525 2975]*um,'y',[  25 1475]*um);
subB = struct('x',[  25 1475]*um,'y',[1525 2975]*um);
ov = @(e,r) max(0, min(e(2), plateHalf+r) - max(e(1), -plateHalf+r));   % 1-D overlap
Cmodel = @(rxv,g) eps0*( ov(subA.x,rxv).*ov(subA.y,0) + ov(subB.x,rxv).*ov(subB.y,0) )./g;

fprintf('================ shape test: is C(rx) one-sided? ================\n');
fprintf('%8s | %11s %11s %11s | %9s\n','rz[um]','C(-40um)','C(0)','C(+40um)','asym');
for i = 1:numel(rzVals)
    m = rz==rzVals(i); r = rx(m); c = C(m);
    [~,iL] = min(abs(r+40*um)); [~,i0] = min(abs(r)); [~,iR] = min(abs(r-40*um));
    asym = (c(iL)-c(iR))/(c(i0)-0.5*(c(iL)+c(iR)));
    fprintf('%8.1f | %11.5f %11.5f %11.5f | %9.3f\n', rzVals(i)*1e6, ...
        c(iL)*1e12, c(i0)*1e12, c(iR)*1e12, asym);
end
fprintf('\n(one-sided truncation predicts C(-40) << C(0) == C(+40), i.e. strongly\n');
fprintf(' asymmetric; a symmetric dip means the mechanism is NOT edge truncation)\n');

%% ---- energy-derived Fx vs model ----
% Point-to-point differencing of W is noise-dominated here, so instead fit
% a STRAIGHT LINE to C(rx) separately on each side of rx=0 and take the
% slope. The model's claim is sharp: slope = eps0*dY/g for rx<0, and
% slope = 0 for rx>0. Uncertainties are the usual OLS standard errors.
fprintf('\n================ Fx = +(1/2)V^2 dC/drx  vs  overlap model ================\n');
fprintf('%7s | %18s %18s | %12s %7s\n', ...
    'rz[um]','Fx(rx<0) [nN]','Fx(rx>0) [nN]','model [nN]','meas/mdl');
for i = 1:numel(rzVals)
    m = rz==rzVals(i); r = rx(m); c = C(m); [r,o] = sort(r); c = c(o);
    g  = g0 - a5*rzVals(i);
    Fmodel = 0.5*V^2 * eps0*ov(subA.y,0)/g;
    [fn,sn] = slopeFit(r(r<-5*um), c(r<-5*um), V);
    [fp,sp] = slopeFit(r(r> 5*um), c(r> 5*um), V);
    fprintf('%7.1f | %9.1f +/- %5.1f %9.1f +/- %5.1f | %12.1f %7.2f\n', ...
        rzVals(i)*1e6, fn*1e9, sn*1e9, fp*1e9, sp*1e9, Fmodel*1e9, fn/Fmodel);
end
fprintf('\nmodel predicts the right column pair as (nonzero, ZERO).\n');

%% ---- plot ----
figure('Color','w','Position',[40 40 1340 430]);
cols = lines(numel(rzVals));

% The C variation is only ~1-2% of the level, so plot it NORMALISED or the
% shape is invisible.
subplot(1,3,1); hold on; grid on;
for i = 1:numel(rzVals)
    m = rz==rzVals(i); [r,o] = sort(rx(m)); c = C(m); c = c(o);
    plot(r*1e6, 100*(c/max(c)-1),'o-','MarkerSize',3,'Color',cols(i,:), ...
        'DisplayName',sprintf('r_z=%.0f\\mum',rzVals(i)*1e6));
end
xline(0,'k:'); xlabel('r_x [\mum]'); ylabel('C / C_{max} - 1  [%]');
title('measured C vs r_x (normalised)'); legend('Location','south','FontSize',7);

subplot(1,3,2); hold on; grid on;
rg = linspace(-40,40,400)'*um;
for i = 1:numel(rzVals)
    g = g0 - a5*rzVals(i);
    cm = Cmodel(rg,g);
    plot(rg*1e6, 100*(cm/max(cm)-1),'-','LineWidth',1.4,'Color',cols(i,:));
end
xline(0,'k:'); xlabel('r_x [\mum]'); ylabel('C / C_{max} - 1  [%]');
title({'overlap-model prediction','(flat for r_x>0, falling for r_x<0)'});

subplot(1,3,3); hold on; grid on;
for i = 1:numel(rzVals)
    m = rz==rzVals(i); [r,o] = sort(rx(m)); c = C(m); c = c(o);
    g = g0 - a5*rzVals(i);
    Fmodel = 0.5*V^2 * eps0*ov(subA.y,0)/g;
    fn = slopeFit(r(r<-5*um), c(r<-5*um), V);
    fp = slopeFit(r(r> 5*um), c(r> 5*um), V);
    plot(i+[-0.18 0.18], [fn fp]*1e9, 'o','MarkerSize',8,'MarkerFaceColor',cols(i,:), ...
        'Color',cols(i,:));
    plot(i-0.18, Fmodel*1e9, 'k_','MarkerSize',14,'LineWidth',2);
end
set(gca,'XTick',1:numel(rzVals),'XTickLabel',compose('%.0f',rzVals*1e6));
yline(0,'k:'); set(gca,'YScale','log');
xlabel('r_z [\mum]'); ylabel('|F_x| [nN]');
title({'left dot: r_x<0   right dot: r_x>0','black dash = model'});
sgtitle('Does overlap truncation explain F_x? (quadrant pair, V=30 V)');
outPng = fullfile(scriptDir,'OverlapModelFx.png');
exportgraphics(gcf,outPng);
fprintf('\nSaved plot to %s\n', outPng);

function [F, se] = slopeFit(r, c, V)
% OLS slope of C vs rx, converted to force F = 0.5*V^2*dC/drx, with its
% standard error propagated the same way.
X = [r(:), ones(numel(r),1)];
p = X\c(:);
res = c(:) - X*p;
s2 = sum(res.^2)/max(numel(r)-2,1);
Cov = s2 * inv(X'*X);
F  = 0.5*V^2*p(1);
se = 0.5*V^2*sqrt(Cov(1,1));
end
