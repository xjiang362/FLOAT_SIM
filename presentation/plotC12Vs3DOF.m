%% plotC12Vs3DOF.m
% C12 = Cm(1,2), the DIRECT capacitance between the two quadrant-1
% sensing electrodes (sense1 and sense2), vs rx, rz and thx.
%
% Unlike C1R, this is a COPLANAR coupling: both electrodes sit on the
% stator, facing each other across the x = y diagonal, and neither moves.
% So C12 depends on the plate pose only through SHIELDING -- the proof
% mass hovering above them intercepts the fringing field that would
% otherwise link the two triangles. Bringing the plate closer (increasing
% rz) should therefore REDUCE C12, which is the opposite sense to C1R.
%
% No parallel-plate formula applies here; the comparison made is against
% the gap itself, to see whether the shielding follows a simple law.
%
% Poses with |rz| > ~10 um are marked: the extracted matrix develops
% negative off-diagonals there (mesh under-resolution in the thin gap),
% and C12 is one of the small terms that gets swamped -- it reaches
% -223 fF at rz = +14 um, which is unphysical.
%
% SYMBOLS: Cm [F] mutual capacitance matrix; C12 [F] = Cm(1,2), plotted
%   in fF; G [m] centroid gap = g0 - a5*rz; rx, rz [m]; thx [deg].

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
g0 = 15e-6; a5 = 0.976;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

L = load(fullfile(scriptDir,'capMatrixSweep63_parsed.mat'));
Cm = L.Cm; nPose = L.nPose; nTerm = L.nTerm; blk = L.blk;
C12 = squeeze(Cm(1,2,:));
G = g0 - a5*L.rz;

bad = false(nPose,1);
for p=1:nPose
    M = Cm(:,:,p); off = M - diag(diag(M));
    bad(p) = any(off(:) < 0);
end
good = ~bad;
fprintf('untrustworthy poses: %d (|rz| >= %.1f um)\n', sum(bad), min(abs(L.rz(bad)))*1e6);

axVar = {L.rx, L.rz, L.thx};
axName = {'r_x','r_z','\theta_x'};
axUnit = {'\mum','\mum','deg'};
axScale = [1e6 1e6 1];
order = [2 1 3];

fprintf('\n%-6s %12s %12s %10s\n','axis','C12 min','C12 max','swing');
for c=1:3
    i=order(c); b=blk{i}; v=C12(b(good(b)));
    fprintf('%-6s %9.2f fF %9.2f fF %9.1f%%\n', axName{i}, min(v)*1e15, max(v)*1e15, ...
        100*(max(v)-min(v))/mean(v));
end
fprintf('  (trustworthy poses only)\n');

%% does the rz dependence follow a simple law in the gap?
b = blk{2}; bg = b(good(b));
gg = G(bg); cc = C12(bg);
forms = {'constant', ones(numel(gg),1); ...
         'a + b*G',  [ones(numel(gg),1), gg]; ...
         'a + b/G',  [ones(numel(gg),1), 1./gg]; ...
         'a + b*log(G)', [ones(numel(gg),1), log(gg)]};
fprintf('\n=== C12 vs gap on the rz block (trustworthy poses) ===\n');
for k=1:size(forms,1)
    X = forms{k,2}; p = X\cc;
    fprintf('  %-16s R^2 = %+.5f\n', forms{k,1}, r2fun(cc, X*p));
end
fprintf('  C12 falls from %.2f fF at gap %.1f um to %.2f fF at gap %.1f um\n', ...
    cc(1)*1e15, gg(1)*1e6, cc(end)*1e15, gg(end)*1e6);

k = 11;
fprintf('\nnominal pose: C12 = %.2f fF  (C1R = %.2f fF, so C12/C1R = %.3f)\n', ...
    C12(k)*1e15, Cm(1,21,k)*1e15, C12(k)/Cm(1,21,k));

%% ---- plot ----
figure('Color','w','Position',[30 30 1380 440]);
for c = 1:3
    i = order(c); b = blk{i};
    x = axVar{i}(b)*axScale(i); [x,o] = sort(x); bo = b(o); bd = bad(bo);
    subplot(1,3,c); hold on; grid on; box on;
    plot(x, C12(bo)*1e15,'o-','MarkerSize',4.5,'LineWidth',1.4, ...
        'Color',[0.12 0.25 0.60],'DisplayName','C_{12} = C_m(1,2)');
    if any(bd)
        plot(x(bd), C12(bo(bd))*1e15,'rx','MarkerSize',12,'LineWidth',2, ...
            'DisplayName','negative off-diagonals');
    end
    xline(0,'k:','HandleVisibility','off'); yline(0,'k-','HandleVisibility','off');
    xlabel(sprintf('%s [%s]',axName{i},axUnit{i})); ylabel('C_{12} [fF]');
    title(sprintf('C_{12} vs %s',axName{i}));
    legend('Location','best','FontSize',7.5);
end
sgtitle({'C_{12}: direct coupling between the two quadrant-1 sensing electrodes', ...
    'both electrodes are fixed -- the pose enters only through shielding by the proof mass'});
for a = findobj(gcf,'Type','axes')', a.Toolbar.Visible='off'; end
outPng = fullfile(scriptDir,'C12_vs_3DOF.png');
exportgraphics(gcf,outPng,'Resolution',150);
fprintf('\nSaved plot to %s\n', outPng);
