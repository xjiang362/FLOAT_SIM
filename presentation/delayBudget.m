%% delayBudget.m
% How much acquisition delay can the loop actually tolerate?
%
% The textbook bound for an unstable pole p with loop delay td is
% p*td < 1, which here (p = 1144 rad/s) would allow 874 us. That bound is
% NECESSARY, not sufficient: it says no controller can stabilise beyond it,
% not that this controller works below it. Measured here, the real budget
% is far tighter.
%
% The delay has two parts, and they are swept separately first to show that
% only their SUM matters:
%   skew  the four channels are converted dtChannel apart, so a moving
%         plate is sampled at four different poses and the plane fit infers
%         a tilt that is not there
%   lag   tI2C, the transfer after the last conversion, common to all four
%
% The oldest channel is  age = tI2C + dtChannel*(n - 1 + 0.5).
%
% See also TRAPSAMPLERINIT, SIMULATETRAPLEVITATION.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p0 = trapDeviceGeometry();
wu = 1143.68;
x0 = zeros(12,1);
x0(1)=-40e-6; x0(2)=-40e-6; x0(3)=-13e-6; x0(4)=60e-6; x0(5)=-40e-6;
tEnd = 0.6; Vmax = 120; lamR = 50;

run1 = @(T,dc,ti) doRun(p0,T,dc,ti,x0,tEnd,Vmax,lamR);

%% ---- 1. budget vs sample period, simultaneous sampling (zero skew) ----
fprintf('=== budget with SIMULTANEOUS sampling (dtChannel = 0) ===\n');
fprintf('%10s %14s %12s %12s\n','Tisr[us]','max tI2C[us]','p*td','p*(td+T)');
for T = [500 200 100]*1e-6
    lo = 0; hi = 1000e-6;
    if ~isempty(run1(T,0,lo)), lo = NaN; end
    if isnan(lo)
        fprintf('%10.0f %14s\n', T*1e6, 'fails at 0');
    else
        for it = 1:30
            md = (lo+hi)/2;
            if isempty(run1(T,0,md)), lo = md; else, hi = md; end
        end
        fprintf('%10.0f %14.0f %12.3f %12.3f\n', T*1e6, lo*1e6, wu*lo, wu*(lo+T));
    end
end

%% ---- 2. with the realistic sequential schedule ----
fprintf('\n=== budget with SEQUENTIAL sampling (dtChannel = Tisr/4) ===\n');
fprintf('%10s %14s %14s %12s\n','Tisr[us]','max tI2C[us]','oldest[us]','p*td');
for T = [500 200 100]*1e-6
    dc = T/4;
    lo = 0; hi = 1000e-6;
    if ~isempty(run1(T,dc,0))
        fprintf('%10.0f %14s %14.0f %12.3f\n', T*1e6, 'fails at 0', dc*3.5*1e6, wu*dc*3.5);
    else
        for it = 1:30
            md = (lo+hi)/2;
            if isempty(run1(T,dc,md)), lo = md; else, hi = md; end
        end
        old = lo + dc*3.5;
        fprintf('%10.0f %14.0f %14.0f %12.3f\n', T*1e6, lo*1e6, old*1e6, wu*old);
    end
end

%% ---- 3. plot: heave for a few total delays, Tisr = 100 us ----
T = 100e-6; dc = 0;
tis = [0 100 200 300 500]*1e-6;
col = parula(numel(tis)+1);
f = figure('Color','w','Position',[60 60 980 580]);
hold on; grid on; box on;
fprintf('\n=== traces plotted (Tisr = %.0f us, simultaneous sampling) ===\n', T*1e6);
for i = 1:numel(tis)
    [ct, L] = doRun(p0,T,dc,tis(i),x0,tEnd,Vmax,lamR);
    m = L.t*1e3 <= 20;
    if isempty(ct)
        sel = L.t > tEnd/2;
        amp = (max(L.x(sel,3))-min(L.x(sel,3)))/2*1e6;
        nm  = sprintf('t_d = %3.0f \\mus  (\\pm%.3f \\mum)', tis(i)*1e6, amp);
        fprintf('  t_d = %3.0f us (p*td = %.2f): survives, amp %.3f um\n', ...
            tis(i)*1e6, wu*tis(i), amp);
    else
        nm = sprintf('t_d = %3.0f \\mus  (contact %.1f ms)', tis(i)*1e6, ct*1e3);
        fprintf('  t_d = %3.0f us (p*td = %.2f): contact %.1f ms\n', ...
            tis(i)*1e6, wu*tis(i), ct*1e3);
    end
    plot(L.t(m)*1e3, L.x(m,3)*1e6,'-','LineWidth',2.0,'Color',col(i,:),'DisplayName',nm);
    if ~isempty(ct) && ct*1e3 <= 20
        plot(ct*1e3, L.x(end,3)*1e6,'x','Color',col(i,:),'MarkerSize',13, ...
             'LineWidth',2.4,'HandleVisibility','off');
    end
end
yline(0,'k:','HandleVisibility','off');
yline(15,'k--','LineWidth',1.2,'HandleVisibility','off');
yline(-15,'k--','LineWidth',1.2,'HandleVisibility','off');
text(16, 15.9,'wall','FontSize',8);
xlabel('t [ms]'); ylabel('r_z [\mum]'); xlim([0 20]); ylim([-15.5 17]);
title(sprintf(['Acquisition delay budget  \\bullet  T_{ISR} = %.0f \\mus, ' ...
    '\\lambda_r = %d Hz, simultaneous sampling'], T*1e6, lamR),'FontSize',11);
legend('Location','southeast','FontSize',9);
set(gca,'Toolbar',[]);
out = fullfile(root,'presentation','DelayBudget_rz.png');
exportgraphics(f,out,'Resolution',140);
fprintf('\nSaved %s\n', out);

% =====================================================================
function [ct, L] = doRun(p0,T,dc,ti,x0,tEnd,Vmax,lamR)
p = p0;
p.readout.Tisr      = T;
p.readout.dtChannel = dc;
p.readout.tI2C      = ti;
K = trapControlDesign(p,'lambdar',2*pi*lamR);
% With zero total staleness the sampler has nothing to buffer, so run it as
% the undelayed case rather than asking it to interpolate an empty history.
useDelay = (ti + dc) > 0;
L = simulateTrapLevitation(p,K,x0,tEnd,'delay',useDelay,'Vtia',0,'Vmax',Vmax);
if isempty(L.contact), ct = []; else, ct = L.contactTime; end
end
