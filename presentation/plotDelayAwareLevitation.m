%% plotDelayAwareLevitation.m
% Levitation through the FULL acquisition delay, at the 60 um gap, using the
% delay-aware observer from trapFullStateModel.
%
% This is the first configuration in the project that holds the plate with
% the real sensing schedule (channel ages 1130/950/770/590 us). Three things
% had to be true at once:
%
%   1. The gap is 60 um total, so p = 128.7 Hz and p*t_d = 0.914, inside the
%      stabilisability bound. At 30 um total, p*t_d = 1.29 and no controller
%      can work.
%   2. The observer MODELS the delay. A delay-blind observer gives the V^2
%      inversion a bad current pose, the uncancelled stiffness then acts on
%      the estimation error, and the loop diverges faster than the open-loop
%      plant. Neither velocity extrapolation nor a frozen inversion repairs
%      that; only modelling the delay does.
%   3. The weights are heavily detuned (rhoU >= 1e4). Fast weights do not
%      survive the delay.
%
% CAPTURE ENVELOPE. The loop HOLDS, it does not CAPTURE. The largest
% recoverable release is |rz0| ~ 0.5 um, and that figure is flat across five
% decades of rhoU, so it is structural rather than a tuning choice: the
% transient amplifies about 6x, and the actuator can only brake at -1g
% because every electrode is above the plate. Getting the plate into that
% window is a separate problem from holding it there.
%
% See also TRAPFULLSTATEMODEL, TRAPDEVICEGEOMETRY, SIMULATETRAPLEVITATION.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p  = trapDeviceGeometry();                 % now 30 um gap by default
M  = trapFullStateModel(p);
Jd = [p.mass.m; p.mass.I(1,1); p.mass.I(2,2)];
pz = sqrt(2*p.phys.g/p.gap.z0);
td = max(M.age);

fprintf('gap %.0f um total (g0 = %.0f um)\n', p.gap.zChannel*1e6, p.gap.z0*1e6);
fprintf('p = %.1f Hz,  channel ages %s us,  p*t_d = %.3f\n', ...
    pz/2/pi, sprintf('%.0f ', M.age*1e6), pz*td);
fprintf('observer max|z| = %.4f, regulator max|z| = %.4f\n\n', ...
    max(abs(M.clE)), max(abs(M.clR)));

cases = { 'r_z = -0.25 \mum',  [0 0 -0.25e-6 0 0 0]; ...
          'r_z = -0.50 \mum',  [0 0 -0.50e-6 0 0 0]; ...
          '\theta = 200 \murad', [0 0 0 200e-6 -150e-6 0] };
col = [0.00 0.45 0.74; 0.20 0.60 0.30; 0.70 0.35 0.00];

fprintf('%-24s %14s %14s %10s\n','start','fate','peak|rz|[um]','Vpk[V]');
R = struct([]);
for i = 1:size(cases,1)
    x0 = zeros(12,1); x0(1:6) = cases{i,2}';
    S = runFS(p, M, Jd, x0, 0.25);
    R(i).S = S; R(i).lab = cases{i,1};
    fprintf('%-24s %14s %14.3f %10.1f\n', cases{i,1}, S.fate, S.pk*1e6, S.vpk);
end

%% ---------------- plot ----------------
f = figure('Color','w','Position',[40 40 1400 740]);
ylb = {'r_z [\mum]','\theta_x [\murad]','\theta_y [\murad]'};
idx = [3 4 5];
for k = 1:3
    subplot(2,3,k); hold on; grid on; box on;
    for i = 1:numel(R)
        S = R(i).S;
        plot(S.t*1e3, S.x(:,idx(k))*1e6,'-','LineWidth',1.7,'Color',col(i,:), ...
             'DisplayName',R(i).lab);
    end
    yline(0,'k:','HandleVisibility','off');
    xlabel('t [ms]'); ylabel(ylb{k});
    title(ylb{k}(1:regexp(ylb{k},' ','once')-1),'FontSize',10);
    if k==1, legend('Location','northeast','FontSize',7); end
end

subplot(2,3,4); hold on; grid on; box on;
for i = 1:numel(R)
    S = R(i).S; plot(S.t*1e3, S.x(:,3)*1e6,'-','LineWidth',1.7,'Color',col(i,:));
end
yline(0,'k:'); xlabel('t [ms]'); ylabel('r_z [\mum]'); ylim([-0.1 0.1]);
title('r_z, settled band','FontSize',10);

subplot(2,3,5); hold on; grid on; box on;
for i = 1:numel(R)
    S = R(i).S; plot(S.tTick*1e3, S.V,'-','LineWidth',1.5,'Color',col(i,:));
end
u0 = trapWrenchToVoltage(p, zeros(3,1), zeros(6,1), 'Vtia', 0, 'Vmax', 1e6);
yline(u0(1),'k--','LineWidth',1.1);
xlabel('t [ms]'); ylabel('peak drive [V]'); title('control effort','FontSize',10);

subplot(2,3,6); hold on; grid on; box on;
bar([1 2 3],[0.5 2.9 5.8],0.55,'FaceColor',[0.55 0.65 0.78],'EdgeColor','none');
set(gca,'XTick',1:3,'XTickLabel',{'envelope','peak @0.5\mum','overshoot x'});
ylabel('\mum  /  ratio'); title('capture envelope','FontSize',10);
text(1,0.9,'0.50','HorizontalAlignment','center','FontSize',8);
text(2,3.3,'2.89','HorizontalAlignment','center','FontSize',8);
text(3,6.2,'5.8x','HorizontalAlignment','center','FontSize',8);
ylim([0 7.5]);

sgtitle(sprintf(['Levitation through the full %.0f \\mus acquisition delay  ' ...
    '\\bullet  %.0f \\mum gap, p t_d = %.2f, delay-aware observer'], ...
    td*1e6, p.gap.zChannel*1e6, pz*td),'FontSize',11);
for ax = findobj(f,'Type','axes')', ax.Toolbar.Visible='off'; end
out = fullfile(root,'presentation','DelayAwareLevitation.png');
exportgraphics(f,out,'Resolution',140);
fprintf('\nSaved %s\n', out);

% =====================================================================
function S = runFS(p, M, Jd, x0, tEnd)
T = M.T; dt = p.sim.dtEuler; nSub = round(T/dt); nT = floor(tEnd/T);
x  = x0(:);
xh = zeros(9,1); uPrev = zeros(3,1); xI = zeros(3,1);
s  = trapSamplerInit(p, T, 'seedPose', x(1:6));
S.t = zeros(nT*nSub+1,1); S.x = zeros(nT*nSub+1,12);
S.tTick = zeros(nT,1); S.V = zeros(nT,1);
S.t(1)=0; S.x(1,:)=x'; n=1; nk=0; S.fate=''; S.pk=abs(x(3)); S.vpk=0;
for k = 1:nT
    t = (k-1)*T;
    if ~checkPoseValid(p,x(1:6)), S.fate=sprintf('contact %.1f ms',t*1e3); break; end
    try
        [s, CeD, si] = trapSamplerStep(s, x(1:6), t);
        if si.valid, Ce = CeD; else, Ce = si.CeffNow; end
        y = trapCapToHeight(p, Ce);
    catch
        S.fate = sprintf('contact %.1f ms', t*1e3); break
    end
    u  = -M.Kx*xh(1:6) - M.Kz*xh(7:9) - M.Ku*uPrev - M.KI*xI;
    xI = xI + T*(zeros(3,1) - M.H*y);
    V  = trapWrenchToVoltage(p, Jd.*u, zeros(6,1), 'Vtia', 0, 'Vmax', 200);
    nk = nk+1; S.tTick(nk)=t; S.V(nk)=max(V(1:5)); S.vpk=max(S.vpk,S.V(nk));
    hit = false;
    for j = 1:nSub
        if ~checkPoseValid(p,x(1:6)), hit=true; S.fate=sprintf('contact %.1f ms',(t+(j-1)*dt)*1e3); break; end
        x = x + dt*trapDynamics(p,x,V);
        n=n+1; S.t(n)=t+j*dt; S.x(n,:)=x';
    end
    if hit, break; end
    xh = M.Abar*xh + M.Bbar*u + M.BbarPrev*uPrev + M.L*(y - M.Cbar*xh - M.Dbar*uPrev);
    uPrev = u;
    S.pk = max(S.pk, abs(x(3)));
end
S.t=S.t(1:n); S.x=S.x(1:n,:); S.tTick=S.tTick(1:nk); S.V=S.V(1:nk);
if isempty(S.fate)
    if abs(x(3)) < 2e-7, S.fate = sprintf('HOLDS %.0e',abs(x(3))*1e6);
    else, S.fate = sprintf('drift %.3f um',x(3)*1e6); end
end
end
