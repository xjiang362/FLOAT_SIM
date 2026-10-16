%% deriveCeffFromC1RC2R.m
% Closed-form C_eff(1,2) in terms of the two port-to-resonator
% capacitances, and a quantitative check of each approximation level
% against the exact Schur reduction, on the 63-pose dedicated sweep.
%
% NETWORK (ports 1,2 live; resonator R floating; everything else ground):
%   C1R, C2R  port <-> resonator          C_rg  resonator <-> ground
%   C1g, C2g  port <-> ground (direct)    C12   port <-> port (direct)
%
% LEVEL 1 -- exact for this 6-branch lumped network.
% Eliminating the floating node R (the Y-delta step) with
% S_R = C1R + C2R + C_rg turns the star into three delta branches:
%     1<->2    : C1R*C2R/S_R
%     1<->gnd  : C1R*C_rg/S_R
%     2<->gnd  : C2R*C_rg/S_R
% so, with the direct branches added in parallel,
%     C_M = C12 + C1R*C2R/S_R
%     C1  = C1g + C1R*C_rg/S_R
%     C2  = C2g + C2R*C_rg/S_R
%     C_eff = C_M + C1*C2/(C1+C2)
%
% LEVEL 2 -- use C_rg >> C1R, C2R (measured ratio ~66), so S_R -> C_rg:
%     C_eff ~= C12 + C1R*C2R/C_rg + (C1g+C1R)*(C2g+C2R)/(C1g+C2g+C1R+C2R)
% i.e. the resonator acts as an AC ground and each port simply sees C1R
% (resp. C2R) added to its own direct to-ground capacitance.
%
% LEVEL 3 -- additionally drop C12, C1g, C2g and the C1R*C2R/C_rg term:
%     C_eff ~= C1R*C2R/(C1R + C2R)
% the plain SERIES combination of the two port-to-resonator capacitances.
%
% SYMBOLS: all capacitances in F (reported in fF); S_R [F] resonator node
%   row sum; C_eff [F] measured two-terminal capacitance between ports 1,2.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
L = load(fullfile(scriptDir,'capMatrixSweep63_parsed.mat'));
Cm = L.Cm; nPose = L.nPose; nTerm = L.nTerm; blk = L.blk;
RES = 21; P1 = 1; P2 = 2;
gnd = setdiff(1:nTerm, [P1 P2 RES]);

%% ---- branch values straight from the matrix ----
C1R = squeeze(Cm(P1,RES,:));
C2R = squeeze(Cm(P2,RES,:));
C12 = squeeze(Cm(P1,P2,:));
Crg = L.Crg;
C1g = squeeze(sum(Cm(P1,gnd,:),2)) + squeeze(Cm(P1,P1,:));
C2g = squeeze(sum(Cm(P2,gnd,:),2)) + squeeze(Cm(P2,P2,:));
SR  = C1R + C2R + Crg;

%% ---- exact reference (full 22-terminal Schur reduction) ----
Cex = zeros(nPose,1);
for p = 1:nPose, Cex(p) = effectiveMutualCapacitance(Cm(:,:,p),P1,P2,RES); end

%% ---- the three levels ----
% Level 1
CM1 = C12 + C1R.*C2R./SR;
A1  = C1g + C1R.*Crg./SR;
B1  = C2g + C2R.*Crg./SR;
L1  = CM1 + A1.*B1./(A1+B1);
% Level 2
L2 = C12 + C1R.*C2R./Crg + (C1g+C1R).*(C2g+C2R)./(C1g+C2g+C1R+C2R);
% Level 3
L3 = C1R.*C2R./(C1R+C2R);

err = @(x) 100*(x-Cex)./Cex;
fprintf('=== accuracy vs the exact 22-terminal Schur reduction ===\n');
fprintf('%-42s %10s %10s\n','model','median err','max |err|');
fprintf('%-42s %9.3f%% %9.3f%%\n','L1  6-branch exact (Y-delta)',  median(err(L1)), max(abs(err(L1))));
fprintf('%-42s %9.3f%% %9.3f%%\n','L2  C_rg >> C1R,C2R',            median(err(L2)), max(abs(err(L2))));
fprintf('%-42s %9.3f%% %9.3f%%\n','L3  series C1R*C2R/(C1R+C2R)',   median(err(L3)), max(abs(err(L3))));

k = 11;
fprintf('\n=== at the nominal pose (fF) ===\n');
fprintf('  C1R %7.2f | C2R %7.2f | C1g %7.2f | C2g %7.2f | C12 %6.2f | C_rg %9.1f\n', ...
    C1R(k)*1e15, C2R(k)*1e15, C1g(k)*1e15, C2g(k)*1e15, C12(k)*1e15, Crg(k)*1e15);
fprintf('  exact %7.2f | L1 %7.2f | L2 %7.2f | L3 %7.2f\n', ...
    Cex(k)*1e15, L1(k)*1e15, L2(k)*1e15, L3(k)*1e15);

nm = {'rx','rz','thx'};
fprintf('\n=== max |error| per axis ===\n');
fprintf('%-5s %10s %10s %10s\n','axis','L1','L2','L3');
for i=1:3
    b = blk{i};
    e = @(x) 100*(x-Cex(b))./Cex(b);
    fprintf('%-5s %9.3f%% %9.3f%% %9.3f%%\n', nm{i}, ...
        max(abs(e(L1(b)))), max(abs(e(L2(b)))), max(abs(e(L3(b)))));
end

%% ---- are C1g, C2g usable as constants? ----
fprintf('\n=== how constant are the direct-to-ground terms? (fF) ===\n');
fprintf('  C1g: %7.1f .. %7.1f (median %7.1f)\n', min(C1g)*1e15, max(C1g)*1e15, median(C1g)*1e15);
fprintf('  C2g: %7.1f .. %7.1f (median %7.1f)\n', min(C2g)*1e15, max(C2g)*1e15, median(C2g)*1e15);
fprintf('  C12: %7.2f .. %7.2f (median %7.2f)\n', min(C12)*1e15, max(C12)*1e15, median(C12)*1e15);

% Level 2 with C1g, C2g, C12 frozen at their nominal values -> a model in
% C1R, C2R (and C_rg) alone
c1g = C1g(k); c2g = C2g(k); c12 = C12(k);
L2f = c12 + C1R.*C2R./Crg + (c1g+C1R).*(c2g+C2R)./(c1g+c2g+C1R+C2R);
fprintf('\n  L2 with C1g,C2g,C12 FROZEN at nominal: median err %.3f%%, max |err| %.3f%%\n', ...
    median(err(L2f)), max(abs(err(L2f))));

%% ---- plot ----
figure('Color','w','Position',[40 40 1340 430]);
subplot(1,3,1); hold on; grid on; box on;
plot(Cex*1e15, L1*1e15,'o','MarkerSize',4,'DisplayName','L1 6-branch');
plot(Cex*1e15, L2*1e15,'s','MarkerSize',4,'DisplayName','L2 C_{rg}\gg C_{1R}');
plot(Cex*1e15, L3*1e15,'^','MarkerSize',4,'DisplayName','L3 series only');
lims=[min(Cex) max(Cex)]*1e15; plot(lims,lims,'k--','HandleVisibility','off');
set(gca,'XScale','log','YScale','log');
xlabel('exact C_{eff} [fF]'); ylabel('model [fF]'); title('model vs exact');
legend('Location','northwest','FontSize',8);

subplot(1,3,2); hold on; grid on; box on;
plot(Cex*1e15, err(L1),'o','MarkerSize',4,'DisplayName','L1');
plot(Cex*1e15, err(L2),'s','MarkerSize',4,'DisplayName','L2');
plot(Cex*1e15, err(L3),'^','MarkerSize',4,'DisplayName','L3');
yline(0,'k-','HandleVisibility','off'); set(gca,'XScale','log');
xlabel('exact C_{eff} [fF]'); ylabel('error [%]'); title('relative error');
legend('Location','best','FontSize',8);

subplot(1,3,3); hold on; grid on; box on;
b = blk{2};  % rz block
[rzs,o] = sort(L.rz(b));
plot(rzs*1e6, Cex(b(o))*1e15,'ko','MarkerSize',5,'DisplayName','exact');
plot(rzs*1e6, L1(b(o))*1e15,'-','LineWidth',1.6,'DisplayName','L1');
plot(rzs*1e6, L2(b(o))*1e15,'--','LineWidth',1.5,'DisplayName','L2');
plot(rzs*1e6, L3(b(o))*1e15,':','LineWidth',1.8,'DisplayName','L3');
set(gca,'YScale','log'); xlabel('r_z [\mum]'); ylabel('C_{eff} [fF]');
title('along the r_z block'); legend('Location','northwest','FontSize',8);

sgtitle('C_{eff}(1,2) as a function of C_{1R}, C_{2R}');
for a = findobj(gcf,'Type','axes')', a.Toolbar.Visible='off'; end
outPng = fullfile(scriptDir,'CeffFromC1RC2R.png');
exportgraphics(gcf,outPng,'Resolution',150);
fprintf('\nSaved plot to %s\n', outPng);
