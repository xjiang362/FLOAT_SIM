%% plotC1gVs3DOF.m
% C1g -- sensing port 1's DIRECT capacitance to ground -- as a function of
% rx, rz and thx, on the 63-pose dedicated capacitance-matrix sweep.
%
% "Ground" here is the reference node plus every fixed-potential terminal,
% i.e. everything except port 1, port 2 and the floating resonator (21).
% From the mutual matrix that is the row sum minus the two live branches:
%
%       C1g = sum_k Cm(1,k) - Cm(1,2) - Cm(1,21)
%           = C_Maxwell(1,1) - C12 - C1R
%
% The lower row breaks C1g into its contributors. At the nominal pose it
% is dominated by the frame (terminal 22) and by the adjacent actuation
% electrode in the same quadrant (terminal 3, C_act1+).
%
% Poses with |rz| > 10 um are marked: there the extracted mutual matrix
% develops NEGATIVE off-diagonal entries, which is unphysical, so C1g is
% not trustworthy in that region (it even goes negative at rz = +14 um).
%
% SYMBOLS: Cm [F] mutual capacitance matrix; C1g [F] port-1-to-ground
%   direct capacitance (plotted in fF); C12 [F] port-to-port direct;
%   C1R [F] port 1 to resonator; rx, rz [m]; thx [deg].

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
L = load(fullfile(scriptDir,'capMatrixSweep63_parsed.mat'));
Cm = L.Cm; nPose = L.nPose; nTerm = L.nTerm; blk = L.blk;
RES = 21; P1 = 1; P2 = 2;

rowsum = squeeze(sum(Cm(P1,:,:),2));
C1g = rowsum - squeeze(Cm(P1,P2,:)) - squeeze(Cm(P1,RES,:));

% contributors
Cfr   = squeeze(Cm(P1,22,:));            % frame
Cact  = squeeze(Cm(P1,3,:));             % C_act1+, adjacent actuation electrode
Cself = squeeze(Cm(P1,P1,:));            % self / to-reference
others = C1g - Cfr - Cact - Cself;

% flag poses whose matrix has negative off-diagonals
bad = false(nPose,1);
for p=1:nPose
    M = Cm(:,:,p); off = M - diag(diag(M));
    bad(p) = any(off(:) < 0);
end
fprintf('poses with negative off-diagonals (untrustworthy): %d -> |rz| >= %.1f um\n', ...
    sum(bad), min(abs(L.rz(bad)))*1e6);

axVar = {L.rx, L.rz, L.thx};
axName = {'r_x','r_z','\theta_x'};
axUnit = {'\mum','\mum','deg'};
axScale = [1e6 1e6 1];
order = [2 1 3];     % display rz, rx, thx

fprintf('\n%-6s %12s %12s %10s\n','axis','C1g min','C1g max','swing');
for c=1:3
    i=order(c); b=blk{i}; v=C1g(b);
    fprintf('%-6s %9.2f fF %9.2f fF %9.1f%%\n', axName{i}, min(v)*1e15, max(v)*1e15, ...
        100*(max(v)-min(v))/abs(mean(v)));
end
g = ~bad;
fprintf('\nrestricted to trustworthy poses: C1g = %.2f .. %.2f fF (median %.2f)\n', ...
    min(C1g(g))*1e15, max(C1g(g))*1e15, median(C1g(g))*1e15);

%% ---- plot ----
figure('Color','w','Position',[30 30 1380 700]);
for c = 1:3
    i = order(c); b = blk{i};
    x = axVar{i}(b)*axScale(i); [x,o] = sort(x); bo = b(o);
    bd = bad(bo);

    subplot(2,3,c); hold on; grid on; box on;
    plot(x, C1g(bo)*1e15,'o-','MarkerSize',4,'LineWidth',1.4,'Color',[0.12 0.25 0.60], ...
        'DisplayName','C_{1g}');
    if any(bd)
        plot(x(bd), C1g(bo(bd))*1e15,'rx','MarkerSize',11,'LineWidth',1.8, ...
            'DisplayName','negative off-diagonals');
    end
    xline(0,'k:','HandleVisibility','off'); yline(0,'k-','HandleVisibility','off');
    xlabel(sprintf('%s [%s]', axName{i}, axUnit{i})); ylabel('C_{1g} [fF]');
    title(sprintf('C_{1g} vs %s', axName{i}));
    if any(bd), legend('Location','best','FontSize',7.5); end

    subplot(2,3,3+c); hold on; grid on; box on;
    plot(x, Cfr(bo)*1e15,'o-','MarkerSize',3.5,'LineWidth',1.2,'DisplayName','frame (22)');
    plot(x, Cact(bo)*1e15,'s-','MarkerSize',3.5,'LineWidth',1.2,'DisplayName','C_{act1+} (3)');
    plot(x, Cself(bo)*1e15,'^-','MarkerSize',3.5,'LineWidth',1.2,'DisplayName','self Cm(1,1)');
    plot(x, others(bo)*1e15,'v-','MarkerSize',3.5,'LineWidth',1.2,'DisplayName','all others');
    xline(0,'k:','HandleVisibility','off'); yline(0,'k-','HandleVisibility','off');
    xlabel(sprintf('%s [%s]', axName{i}, axUnit{i})); ylabel('contribution [fF]');
    title('breakdown of C_{1g}'); legend('Location','best','FontSize',7);
end
sgtitle({'C_{1g}: sensing port 1 direct capacitance to ground', ...
    'ports 1,2 live; resonator (21) floating; all other terminals grounded'});
for a = findobj(gcf,'Type','axes')', a.Toolbar.Visible='off'; end
outPng = fullfile(scriptDir,'C1g_vs_3DOF.png');
exportgraphics(gcf,outPng,'Resolution',150);
fprintf('\nSaved plot to %s\n', outPng);
