%% evalUserCeffModel.m
% FINAL compact model for the sensed capacitance between sensing ports 1
% and 2, with everything frozen except the pose:
%
%   C_eff,12 ~= C12 + C1g/2
%               + (eps0/2) * dX*dY / (g0 - rz - thx*Yc + thy*Xc)
%               + C1f
%
% All of C12, C1g and C1f are CONSTANTS fixed at the equilibrium pose, so
% the only pose dependence is the single gap term. The three constants
% collapse into one offset C0 = C12 + C1g/2 + C1f, giving the one-line form
%
%   C_eff,12 ~= C0 + (eps0/2) * A / g(pose)
%
% Tilt sign convention: -thx*Yc + thy*Xc, i.e. +thx closes the gap at +y
% and +thy opens it at +x. This is the convention established by the
% quadrant force fit and confirmed directly against C1R (R^2 = +0.9987
% with this sign versus -2.886 with the opposite). thy = 0 throughout this
% sweep, so the thy term is carried on the quadrant evidence, not tested
% here.
%
% Footprint convention (see makeElectrodeGeometryTable.m): dY is the full
% y-extent and dX the width at mid-height, so dX*dY = A exactly even for
% the triangular sensing electrodes. Xc,Yc are the AREA centroid, which is
% a different height from the one dX is measured at.
%
% SYMBOLS: eps0 [F/m]; A = dX*dY [m^2]; Xc,Yc [m] area centroid; g0 [m]
%   nominal gap; rz [m]; thx,thy [rad]; all capacitances [F], shown in fF.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
eps0 = 8.8541878128e-12; um = 1e-6; g0 = 15e-6;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

L = load(fullfile(scriptDir,'capMatrixSweep63_parsed.mat'));
Cm = L.Cm; nPose = L.nPose; blk = L.blk; RES = 21;

% --- frozen footprint constants (sense1) ---
dX = 707.5*um; dY = 1415*um; A = dX*dY;
Xc = 2503.3333*um; Yc = 1996.6667*um;

% --- exact reference, and the trustworthy mask ---
Cex = zeros(nPose,1); bad = false(nPose,1);
for p = 1:nPose
    Cex(p) = effectiveMutualCapacitance(Cm(:,:,p),1,2,RES);
    M = Cm(:,:,p); off = M - diag(diag(M)); bad(p) = any(off(:) < 0);
end
good = ~bad;

% --- frozen electrical constants, read at the equilibrium pose ---
k = 11;
C12  = Cm(1,2,k);
C1g  = sum(Cm(1,:,k)) - Cm(1,2,k) - Cm(1,RES,k);

thx = L.thx*pi/180; thy = L.thy*pi/180; rz = L.rz;
gap = g0 - rz - thx*Yc + thy*Xc;
gapTerm = (eps0/2)*A./gap;

C1f = Cex(k) - (C12 + C1g/2 + gapTerm(k));     % exact at equilibrium
C0  = C12 + C1g/2 + C1f;

fprintf('=== frozen constants (all at the equilibrium pose) ===\n');
fprintf('  C12            %9.4f fF\n', C12*1e15);
fprintf('  C1g            %9.4f fF   -> C1g/2 = %.4f fF\n', C1g*1e15, C1g/2*1e15);
fprintf('  C1f            %9.4f fF\n', C1f*1e15);
fprintf('  ------------------------------\n');
fprintf('  C0 = C12 + C1g/2 + C1f = %.4f fF\n', C0*1e15);
fprintf('  gap term at equilibrium  = %.4f fF\n', gapTerm(k)*1e15);
fprintf('  sum                      = %.4f fF   (exact: %.4f fF)\n', ...
    (C0+gapTerm(k))*1e15, Cex(k)*1e15);
fprintf('\n  geometry: dX = %.1f um, dY = %.1f um, A = %.4e m^2\n', dX/um, dY/um, A);
fprintf('            Xc = %.1f um, Yc = %.1f um, g0 = %.1f um\n', Xc/um, Yc/um, g0/um);

%% --- model and errors ---
Mdl = C0 + gapTerm;
err = 100*(Mdl-Cex)./Cex;
nm = {'rx','rz','thx'};
fprintf('\n=== model vs exact reduction (trustworthy poses) ===\n');
fprintf('  %-6s %6s %12s %12s %12s %10s\n','axis','n','median err','max |err|','RMS err','R^2');
for i=1:3
    b = blk{i}; bg = b(good(b));
    fprintf('  %-6s %6d %11.3f%% %11.3f%% %11.3f%% %10.5f\n', nm{i}, numel(bg), ...
        median(err(bg)), max(abs(err(bg))), rms(err(bg)), r2fun(Cex(bg),Mdl(bg)));
end
bg = good;
fprintf('  %-6s %6d %11.3f%% %11.3f%% %11.3f%% %10.5f\n','ALL', sum(bg), ...
    median(err(bg)), max(abs(err(bg))), rms(err(bg)), r2fun(Cex(bg),Mdl(bg)));

% within the tighter operating band
tight = good & abs(rz) <= 10e-6;
fprintf('  %-6s %6d %11.3f%% %11.3f%% %11.3f%% %10.5f   <- |rz| <= 10 um\n','band', sum(tight), ...
    median(err(tight)), max(abs(err(tight))), rms(err(tight)), r2fun(Cex(tight),Mdl(tight)));

%% --- would a least-squares C1f do better than the equilibrium one? ---
C1f_ls = mean(Cex(good) - (C12 + C1g/2 + gapTerm(good)));
Mls = C12 + C1g/2 + C1f_ls + gapTerm;
els = 100*(Mls-Cex)./Cex;
fprintf('\n=== alternative: C1f fitted over the whole sweep ===\n');
fprintf('  C1f(equilibrium) = %7.4f fF -> median %+.3f%%, max %.3f%%, RMS %.3f%%\n', ...
    C1f*1e15, median(err(good)), max(abs(err(good))), rms(err(good)));
fprintf('  C1f(least sq)    = %7.4f fF -> median %+.3f%%, max %.3f%%, RMS %.3f%%\n', ...
    C1f_ls*1e15, median(els(good)), max(abs(els(good))), rms(els(good)));
fprintf('  (equilibrium calibration keeps the model EXACT at the operating point)\n');

%% --- plot ---
order=[2 1 3]; axVar={L.rx,L.rz,L.thx}; axName={'r_x','r_z','\theta_x'};
axUnit={'\mum','\mum','deg'}; axScale=[1e6 1e6 1];
figure('Color','w','Position',[30 30 1380 700]);
for c=1:3
    i=order(c); b=blk{i};
    x=axVar{i}(b)*axScale(i); [x,o]=sort(x); bo=b(o); bd=bad(bo);
    subplot(2,3,c); hold on; grid on; box on;
    plot(x,Cex(bo)*1e15,'ko','MarkerSize',5,'DisplayName','exact (Schur)');
    plot(x,Mdl(bo)*1e15,'-','LineWidth',1.9,'Color',[0.12 0.25 0.60],'DisplayName','frozen model');
    if any(bd), plot(x(bd),Cex(bo(bd))*1e15,'rx','MarkerSize',11,'LineWidth',1.8, ...
            'DisplayName','untrustworthy'); end
    if i==2, set(gca,'YScale','log'); end
    xlabel(sprintf('%s [%s]',axName{i},axUnit{i})); ylabel('C_{eff,12} [fF]');
    title(sprintf('C_{eff} vs %s',axName{i})); legend('Location','best','FontSize',7.5);

    subplot(2,3,3+c); hold on; grid on; box on;
    gd=~bd;
    plot(x(gd),err(bo(gd)),'o-','MarkerSize',4,'LineWidth',1.4,'Color',[0.12 0.25 0.60]);
    yline(0,'k-'); xline(0,'k:');
    xlabel(sprintf('%s [%s]',axName{i},axUnit{i})); ylabel('model error [%]');
    title('relative error'); ylim([-3 9]);
end
sgtitle(sprintf(['Frozen compact model:  C_{eff,12} = %.2f fF + (\\epsilon_0/2)A/' ...
    '(g_0-r_z-\\theta_xY_c+\\theta_yX_c)'], C0*1e15));
for a=findobj(gcf,'Type','axes')', a.Toolbar.Visible='off'; end
outPng = fullfile(scriptDir,'UserCeffModel_vs_sim.png');
exportgraphics(gcf,outPng,'Resolution',150);
fprintf('\nSaved plot to %s\n', outPng);
