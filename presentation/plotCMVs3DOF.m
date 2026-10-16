%% plotCMVs3DOF.m
% C_M -- the total DIRECT branch between sensing ports 1 and 2 after the
% floating resonator has been eliminated -- vs rx, rz and thx, split into
% its two contributions:
%
%       C_M = C12 + C_delta,     C_delta = C1R*C2R/S_R
%
%   C12      the genuine coplanar coupling between the two electrodes,
%            which falls as the plate approaches (shielding).
%   C_delta  the branch TRANSFERRED onto 1<->2 by eliminating the
%            resonator (the Y-delta step), which rises as the plate
%            approaches because C1R*C2R grows as 1/g^2 while S_R grows
%            only as 1/g.
%
% The two therefore run in OPPOSITE directions with gap, and cross over.
% Top row: the three curves. Bottom row: each term as a share of C_M.
%
% Poses with negative off-diagonals in the extracted matrix are marked;
% there C12 is swamped by mesh error and C_M is not trustworthy.
%
% SYMBOLS: C12 [F] direct port-to-port; C1R,C2R [F] port-to-resonator;
%   S_R [F] resonator node row sum; C_delta [F] transferred branch;
%   C_M [F] total direct branch (all plotted in fF); rx,rz [m]; thx [deg].

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
L = load(fullfile(scriptDir,'capMatrixSweep63_parsed.mat'));
Cm = L.Cm; nPose = L.nPose; nTerm = L.nTerm; blk = L.blk;
RES = 21;

C1R = squeeze(Cm(1,RES,:)); C2R = squeeze(Cm(2,RES,:));
C12 = squeeze(Cm(1,2,:));
SR  = C1R + C2R + L.Crg;
Cdel = C1R.*C2R./SR;
CM   = C12 + Cdel;

bad = false(nPose,1);
for p=1:nPose
    M = Cm(:,:,p); off = M - diag(diag(M)); bad(p) = any(off(:)<0);
end
good = ~bad;

axVar={L.rx,L.rz,L.thx}; axName={'r_x','r_z','\theta_x'};
axUnit={'\mum','\mum','deg'}; axScale=[1e6 1e6 1]; order=[2 1 3];

fprintf('%-6s | %-26s | %-26s | %-26s\n','axis','C_M [fF]','C12 [fF]','C_delta [fF]');
for c=1:3
    i=order(c); b=blk{i}; bg=b(good(b));
    f=@(v) sprintf('%7.2f ..%7.2f %5.0f%%', min(v)*1e15, max(v)*1e15, ...
        100*(max(v)-min(v))/abs(mean(v)));
    fprintf('%-6s | %s | %s | %s\n', axName{i}, f(CM(bg)), f(C12(bg)), f(Cdel(bg)));
end
fprintf('(trustworthy poses only)\n');

% crossover on the rz block
b = blk{2}; bg = b(good(b)); [rzs,o] = sort(L.rz(bg)); bo = bg(o);
ix = find(Cdel(bo) > C12(bo), 1, 'first');
if ~isempty(ix)
    fprintf('\nC_delta overtakes C12 between rz = %+.1f and %+.1f um\n', ...
        L.rz(bo(ix-1))*1e6, L.rz(bo(ix))*1e6);
    fprintf('  i.e. at a gap of about %.1f um\n', (15e-6-0.976*L.rz(bo(ix)))*1e6);
end
k=11;
fprintf('\nnominal: C12 %.2f + C_delta %.2f = C_M %.2f fF  (C12 is %.1f%% of C_M)\n', ...
    C12(k)*1e15, Cdel(k)*1e15, CM(k)*1e15, 100*C12(k)/CM(k));

%% ---- plot ----
figure('Color','w','Position',[30 30 1380 720]);
for c=1:3
    i=order(c); b=blk{i};
    x = axVar{i}(b)*axScale(i); [x,o]=sort(x); bo=b(o); bd=bad(bo);

    subplot(2,3,c); hold on; grid on; box on;
    plot(x, CM(bo)*1e15,'o-','MarkerSize',4,'LineWidth',1.7,'Color',[0.12 0.25 0.60], ...
        'DisplayName','C_M = C_{12}+C_\delta');
    plot(x, C12(bo)*1e15,'s-','MarkerSize',3.5,'LineWidth',1.2,'Color',[0.85 0.33 0.10], ...
        'DisplayName','C_{12} (coplanar, shielded)');
    plot(x, Cdel(bo)*1e15,'^-','MarkerSize',3.5,'LineWidth',1.2,'Color',[0.10 0.55 0.35], ...
        'DisplayName','C_\delta = C_{1R}C_{2R}/S_R');
    if any(bd)
        plot(x(bd), CM(bo(bd))*1e15,'rx','MarkerSize',11,'LineWidth',1.8,'DisplayName','untrustworthy');
    end
    xline(0,'k:','HandleVisibility','off'); yline(0,'k-','HandleVisibility','off');
    xlabel(sprintf('%s [%s]',axName{i},axUnit{i})); ylabel('[fF]');
    title(sprintf('C_M vs %s',axName{i})); legend('Location','best','FontSize',7);

    subplot(2,3,3+c); hold on; grid on; box on;
    gd = ~bd;
    plot(x(gd), 100*C12(bo(gd))./CM(bo(gd)),'s-','MarkerSize',3.5,'LineWidth',1.3, ...
        'Color',[0.85 0.33 0.10],'DisplayName','C_{12} share');
    plot(x(gd), 100*Cdel(bo(gd))./CM(bo(gd)),'^-','MarkerSize',3.5,'LineWidth',1.3, ...
        'Color',[0.10 0.55 0.35],'DisplayName','C_\delta share');
    yline(50,'k--','HandleVisibility','off');
    xline(0,'k:','HandleVisibility','off');
    xlabel(sprintf('%s [%s]',axName{i},axUnit{i})); ylabel('share of C_M [%]');
    ylim([0 110]); title('contribution split'); legend('Location','best','FontSize',7);
end
sgtitle({'C_M and its two contributions', ...
    'C_{12} falls with gap (shielding) while C_\delta rises -- they cross'});
for a = findobj(gcf,'Type','axes')', a.Toolbar.Visible='off'; end
outPng = fullfile(scriptDir,'CM_vs_3DOF.png');
exportgraphics(gcf,outPng,'Resolution',150);
fprintf('\nSaved plot to %s\n', outPng);
