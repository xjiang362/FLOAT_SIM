%% plotC1RvsModel_3DOF.m
% C1R = Cm(1,21), the sensing port 1 <-> resonator capacitance, vs rx, rz
% and thx on the 63-pose dedicated sweep, compared with the centroid
% analytic model
%
%       C ~= eps0 * A / (g0 - rz +/- thx*Yc -/+ thy*Xc)
%
% Terminal 1 is sense1, a RIGHT TRIANGLE with vertices (1560,1525),
% (2975,1525), (2975,2940) um. Two things matter for the comparison:
%
%   * AREA. The bounding box is dX = dY = 1415 um, but the area is
%     A = dX*dY/2 = 1.0011e6 um^2, NOT dX*dY. Using dX*dY would be 2x
%     too large.
%   * CENTROID. For a triangle the area centroid (mean of the vertices)
%     is (2503.3, 1996.7) um, which is NOT the bounding-box centre
%     (2267.5, 2232.5). The centroid is what the first-order gap
%     expansion needs.
%
% thy = 0 throughout this sweep, so only the thx sign is testable here.
% Both conventions are fitted and compared:
%       S_MINUS : g = g0 - a5*rz - thx*Yc     (capTrapz / quadrant-fit)
%       S_PLUS  : g = g0 - a5*rz + thx*Yc
%
% The exact integral over the true triangular footprint (capTrapz via
% shapeCap) is shown as a third curve.
%
% SYMBOLS: eps0 [F/m]; A [m^2] footprint area; Xc,Yc [m] area centroid;
%   g0 [m] nominal gap; a5 [-] rz->gap scale; rx,rz [m]; thx [rad];
%   C1R [F] (plotted in fF).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
eps0 = 8.8541878128e-12; um = 1e-6; g0 = 15e-6; a5 = 0.976;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

L = load(fullfile(scriptDir,'capMatrixSweep63_parsed.mat'));
Cm = L.Cm; blk = L.blk; nPose = L.nPose;
C1R = squeeze(Cm(1,21,:));
rx = L.rx; rz = L.rz; thx = L.thx*pi/180;

% geometry of sense1
v = [1560 1525; 2975 1525; 2975 2940]*um;
A  = polyarea(v(:,1),v(:,2));
Xc = mean(v(:,1)); Yc = mean(v(:,2));
dX = max(v(:,1))-min(v(:,1)); dY = max(v(:,2))-min(v(:,2));
fprintf('sense1: A = %.4e m^2 (dX*dY = %.4e, ratio %.2f)\n', A, dX*dY, dX*dY/A);
fprintf('        centroid (%.1f, %.1f) um ; bbox centre (%.1f, %.1f) um\n', ...
    Xc/um, Yc/um, (max(v(:,1))+min(v(:,1)))/2/um, (max(v(:,2))+min(v(:,2)))/2/um);

%% ---- the two sign conventions ----
gMinus = g0 - a5*rz - thx*Yc;
gPlus  = g0 - a5*rz + thx*Yc;
Cminus = eps0*A./gMinus;
Cplus  = eps0*A./gPlus;

% exact integral over the true triangle
g = trapGeometry();
Cexact = shapeCap(g.sense1, g0 - a5*rz, thx, 0*thx, eps0);

%% ---- which thx sign? decide on the thx block alone ----
b = blk{3};
fitgain = @(m,d) (m\d);
gm = fitgain(Cminus(b), C1R(b));  gp = fitgain(Cplus(b), C1R(b));
fprintf('\n=== thx block: which sign? (gain fitted, no offset) ===\n');
fprintf('  g = g0 - a5*rz - thx*Yc : gain %.4f  R^2 = %+.5f\n', gm, r2fun(C1R(b), Cminus(b)*gm));
fprintf('  g = g0 - a5*rz + thx*Yc : gain %.4f  R^2 = %+.5f\n', gp, r2fun(C1R(b), Cplus(b)*gp));
fprintf('  exact triangle integral  : gain %.4f  R^2 = %+.5f\n', ...
    fitgain(Cexact(b),C1R(b)), r2fun(C1R(b), Cexact(b)*fitgain(Cexact(b),C1R(b))));

if r2fun(C1R(b),Cminus(b)*gm) >= r2fun(C1R(b),Cplus(b)*gp)
    Cmod = Cminus; signStr = '-\theta_xY_c'; best='MINUS';
else
    Cmod = Cplus;  signStr = '+\theta_xY_c'; best='PLUS';
end
fprintf('  -> %s wins\n', best);

%% ---- per-axis accuracy with the winning sign ----
nm = {'rx','rz','thx'};
fprintf('\n=== accuracy per axis (gain fitted on that block) ===\n');
fprintf('%-5s | %10s %10s | %10s %10s | %10s\n', ...
    'axis','cent gain','cent R^2','exact gain','exact R^2','data swing');
for i=1:3
    bb = blk{i};
    ga = fitgain(Cmod(bb),C1R(bb)); ge = fitgain(Cexact(bb),C1R(bb));
    fprintf('%-5s | %10.4f %+10.5f | %10.4f %+10.5f | %9.1f%%\n', nm{i}, ...
        ga, r2fun(C1R(bb),Cmod(bb)*ga), ge, r2fun(C1R(bb),Cexact(bb)*ge), ...
        100*(max(C1R(bb))-min(C1R(bb)))/mean(C1R(bb)));
end

k = 11;
fprintf('\n=== nominal pose ===\n');
fprintf('  measured C1R        %8.2f fF\n', C1R(k)*1e15);
fprintf('  centroid model      %8.2f fF  (ratio %.4f -> fringing)\n', ...
    Cmod(k)*1e15, C1R(k)/Cmod(k));
fprintf('  exact triangle      %8.2f fF  (ratio %.4f)\n', Cexact(k)*1e15, C1R(k)/Cexact(k));
fprintf('  if dX*dY were used  %8.2f fF  (%.1fx too big)\n', ...
    eps0*dX*dY/g0*1e15, (dX*dY)/A);

%% ---- rx: the model has no rx dependence at all ----
b = blk{1};
fprintf('\n=== rx block: the analytic model is FLAT here ===\n');
fprintf('  measured C1R swings %.1f%% (%.1f .. %.1f fF)\n', ...
    100*(max(C1R(b))-min(C1R(b)))/mean(C1R(b)), min(C1R(b))*1e15, max(C1R(b))*1e15);
fprintf('  model swings        %.3f%%\n', 100*(max(Cmod(b))-min(Cmod(b)))/mean(Cmod(b)));
fprintf('  -> the missing term is plate-electrode OVERLAP truncation:\n');
fprintf('     sense1 reaches x = 2975 um, flush with the plate edge.\n');

%% ---- plot ----
order = [2 1 3]; axVar={rx,rz,thx}; axName={'r_x','r_z','\theta_x'};
axUnit={'\mum','\mum','deg'}; axScale=[1e6 1e6 180/pi];
figure('Color','w','Position',[30 30 1380 460]);
for c=1:3
    i=order(c); bb=blk{i};
    x = axVar{i}(bb)*axScale(i); [x,o]=sort(x); bo=bb(o);
    ga = fitgain(Cmod(bo),C1R(bo));
    subplot(1,3,c); hold on; grid on; box on;
    plot(x, C1R(bo)*1e15,'o','MarkerSize',5,'Color',[0.2 0.3 0.75],'DisplayName','COMSOL C_{1R}');
    plot(x, Cmod(bo)*1e15,'-','LineWidth',1.7,'Color',[0.85 0.33 0.10], ...
        'DisplayName',sprintf('centroid \\epsilon_0A/g, %s',signStr));
    plot(x, Cexact(bo)*1e15,'k--','LineWidth',1.5,'DisplayName','exact \int dA/g');
    if i==2, set(gca,'YScale','log'); end
    xlabel(sprintf('%s [%s]',axName{i},axUnit{i})); ylabel('C_{1R} [fF]');
    title(sprintf('%s   (R^2=%+.4f, gain %.3f)', axName{i}, r2fun(C1R(bo),Cmod(bo)*ga), ga));
    legend('Location','best','FontSize',7.5);
end
sgtitle('C_{1R} (sense1 \leftrightarrow resonator): COMSOL vs centroid analytic model');
for a = findobj(gcf,'Type','axes')', a.Toolbar.Visible='off'; end
outPng = fullfile(scriptDir,'C1R_vs_model_3DOF.png');
exportgraphics(gcf,outPng,'Resolution',150);
fprintf('\nSaved plot to %s\n', outPng);
