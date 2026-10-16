%% compareTorqueCentralVsQuadrant.m
% Compares the TORQUE produced by the central electrode against that of a
% quadrant actuation pair, at matched drive (V=30 V), matched tilt range
% (+/-0.129 deg) and matched height (rz=0).
%
% The two are structurally different, and the comparison has to respect
% that:
%   * The quadrant pair's centroid sits at (Xc,Yc)=(1500,1500) um, so a
%     purely vertical force Fz already makes a torque Fz*Yc at ZERO tilt
%     ("moment-arm torque"). This is the term that lets the quadrants
%     actuate tilt.
%   * Central electrode A is centred at the origin, so its moment arm is
%     zero by symmetry. It makes torque ONLY by tilting, through the
%     redistribution of electrostatic pressure across its own footprint
%     ("tilt torque"). NOTE this torque grows WITH the tilt angle
%     (d(Taux)/d(thx) > 0), i.e. it is DESTABILIZING -- the usual
%     fixed-voltage electrostatic negative stiffness, not a restoring
%     spring. Same sign for the quadrants.
%
% Data (both at 30 V):
%   quadrant  : comsol_sweep_combined300_results.txt  (dedicated single-
%               axis sweep, 6 axes x 50 pts, MST force/torque probes)
%   central A : comsol_centerElectrodeSweep_rzThxThy150_results.txt
%
% SYMBOLS:
%   thx,thy  [rad / deg]  plate tilts
%   rz       [m]    plate vertical displacement (+rz closes the gap)
%   Fz       [N]    vertical electrostatic force on the plate
%   Taux     [N*m]  torque about the x axis
%   arm      [m]    effective moment arm, Taux/Fz at zero tilt
%   kTilt    [N*m/rad]  tilt torque stiffness, d(Taux)/d(thx) at thx=0

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
d2r = pi/180;

%% ---------------- quadrant: dedicated single-axis sweep ----------------
Q = readmatrix(fullfile(scriptDir,'comsol_sweep_combined300_results.txt'),'NumHeaderLines',5);
% column order of THIS export: rz,rx,ry,thx,thy,thz,Fx,Fy,Fz,probe4,Taux,Tauy,Tauz,Vf8,Vf9
q.rz=Q(:,1); q.rx=Q(:,2); q.ry=Q(:,3); q.thx=Q(:,4); q.thy=Q(:,5); q.thz=Q(:,6);
q.Fz=Q(:,9); q.Taux=Q(:,11); q.Tauy=Q(:,12);
qThx = find(abs(q.thx)>0 & abs(q.rx)==0 & abs(q.ry)==0 & abs(q.rz)==0 & abs(q.thy)==0 & abs(q.thz)==0);
qZero = find(all(abs(Q(:,1:6))==0,2));
fprintf('quadrant: %d rows, %d in the thx-only block, %d at the exact origin pose\n', ...
    size(Q,1), numel(qThx), numel(qZero));

%% ---------------- central electrode A ----------------
Cc = readmatrix(fullfile(scriptDir,'comsol_centerElectrodeSweep_rzThxThy150_results.txt'),'NumHeaderLines',5);
% column order: rx,ry,rz,thx,thy,thz,Fx,Fy,Fz,Taux,Tauy,Tauz,W
c.rx=Cc(:,1); c.ry=Cc(:,2); c.rz=Cc(:,3); c.thx=Cc(:,4); c.thy=Cc(:,5); c.thz=Cc(:,6);
c.Fz=Cc(:,9); c.Taux=Cc(:,10); c.Tauy=Cc(:,11);
cThx = find(abs(c.thx)>0 & abs(c.rz)==0 & abs(c.thy)==0);
cZero = find(all(abs(Cc(:,1:6))==0,2));
fprintf('central A: %d rows, %d in the thx-only block, %d at the exact origin pose\n\n', ...
    size(Cc,1), numel(cThx), numel(cZero));

%% ---------------- 1. zero-tilt moment arm ----------------
fprintf('================ 1. MOMENT ARM AT ZERO TILT (rz=0) ================\n');
% nearest-to-zero tilt point in each thx block
[~,iq] = min(abs(q.thx(qThx)));  iq = qThx(iq);
[~,ic] = min(abs(c.thx(cThx)));  ic = cThx(ic);
fprintf('%-12s %10s %12s %12s %10s\n','electrode','thx[deg]','Fz[N]','Taux[N*m]','arm[um]');
fprintf('%-12s %10.4f %12.4e %12.4e %10.1f\n','quadrant', q.thx(iq), q.Fz(iq), q.Taux(iq), 1e6*q.Taux(iq)/q.Fz(iq));
fprintf('%-12s %10.4f %12.4e %12.4e %10.1f\n','central A', c.thx(ic), c.Fz(ic), c.Taux(ic), 1e6*c.Taux(ic)/c.Fz(ic));
fprintf('\n -> quadrant arm matches its geometric centroid Yc = 1500 um.\n');
fprintf(' -> central A arm is ~0: centred electrode, no moment arm by symmetry.\n');

%% ---------------- 2. tilt-restoring stiffness ----------------
fprintf('\n================ 2. TILT STIFFNESS d(Taux)/d(thx) at thx=0 ================\n');
kQ = localSlope(q.thx(qThx)*d2r, q.Taux(qThx));
kC = localSlope(c.thx(cThx)*d2r, c.Taux(cThx));
fprintf('%-12s %16s\n','electrode','kTilt [N*m/rad]');
fprintf('%-12s %16.4e\n','quadrant',  kQ);
fprintf('%-12s %16.4e\n','central A', kC);
fprintf('ratio central/quadrant = %.3f\n', kC/kQ);

%% ---------------- 3. peak torque over the tilt range ----------------
fprintf('\n================ 3. TORQUE OVER THE FULL TILT RANGE ================\n');
fprintf('%-12s %10s %12s %12s %12s\n','electrode','|thx|max','Taux(0)','Taux(max)','swing');
fprintf('%-12s %10.4f %12.4e %12.4e %12.4e\n','quadrant', max(abs(q.thx(qThx))), ...
    q.Taux(iq), q.Taux(qThx(end)), max(q.Taux(qThx))-min(q.Taux(qThx)));
fprintf('%-12s %10.4f %12.4e %12.4e %12.4e\n','central A', max(abs(c.thx(cThx))), ...
    c.Taux(ic), c.Taux(cThx(end)), max(c.Taux(cThx))-min(c.Taux(cThx)));

mq = max(abs(q.Taux(qThx)));  mc = max(abs(c.Taux(cThx)));
fprintf('\nmax |Taux| : quadrant %.4e, central A %.4e  -> quadrant is %.1fx larger\n', mq, mc, mq/mc);

%% ---------------- 4. what if all four quadrants are driven equally ----------------
fprintf('\n============ 4. FOUR QUADRANTS DRIVEN EQUALLY (common mode) ============\n');
fprintf('The four pair centroids are at (+-1500,+-1500) um, so the Fz*arm terms\n');
fprintf('CANCEL in the sum and only the tilt-restoring part survives.\n');
fprintf('  single quadrant, moment-arm torque at zero tilt : %.4e N*m\n', abs(q.Taux(iq)));
fprintf('  4 quadrants common mode, net moment-arm torque  : 0 (by symmetry)\n');
fprintf('  4 x quadrant tilt stiffness                     : %.4e N*m/rad\n', 4*kQ);
fprintf('  central A tilt stiffness                        : %.4e N*m/rad  (%.1f%% of that)\n', ...
    kC, 100*kC/(4*kQ));

%% ---------------- plot ----------------
figure('Color','w','Position',[40 40 1340 440]);

subplot(1,3,1); hold on; grid on;
[xs,o]=sort(q.thx(qThx)); plot(xs, q.Taux(qThx(o))*1e9,'o-','MarkerSize',4, ...
    'Color',[0.2 0.3 0.75],'LineWidth',1.2,'DisplayName','quadrant pair');
[xs,o]=sort(c.thx(cThx)); plot(xs, c.Taux(cThx(o))*1e9,'s-','MarkerSize',4, ...
    'Color',[0.85 0.33 0.10],'LineWidth',1.2,'DisplayName','central A');
yline(0,'k:'); xlabel('\theta_x [deg]'); ylabel('Taux [nN\cdotm]');
title('Torque vs tilt, r_z=0, V=30 V'); legend('Location','best','FontSize',8);

subplot(1,3,2); hold on; grid on;
[xs,o]=sort(c.thx(cThx)); plot(xs, c.Taux(cThx(o))*1e9,'s-','MarkerSize',4, ...
    'Color',[0.85 0.33 0.10],'LineWidth',1.3,'DisplayName','central A');
yline(0,'k:'); xlabel('\theta_x [deg]'); ylabel('Taux [nN\cdotm]');
title('central A alone (odd in \theta_x, slope>0 -> DEstabilizing)'); legend('Location','best','FontSize',8);

subplot(1,3,3); hold on; grid on;
[xs,o]=sort(q.thx(qThx));
plot(xs, (q.Taux(qThx(o))-q.Taux(iq))*1e9,'o-','MarkerSize',4,'Color',[0.2 0.3 0.75], ...
    'LineWidth',1.2,'DisplayName','quadrant, arm term removed');
[xs,o]=sort(c.thx(cThx)); plot(xs, c.Taux(cThx(o))*1e9,'s-','MarkerSize',4, ...
    'Color',[0.85 0.33 0.10],'LineWidth',1.2,'DisplayName','central A');
yline(0,'k:'); xlabel('\theta_x [deg]'); ylabel('\DeltaTaux [nN\cdotm]');
title('tilt-dependent part only'); legend('Location','best','FontSize',8);

sgtitle('Torque: central electrode A vs quadrant actuation pair');
outPng = fullfile(scriptDir,'TorqueCentralVsQuadrant.png');
exportgraphics(gcf,outPng);
fprintf('\nSaved plot to %s\n', outPng);

%% ---------------- helper ----------------
function k = localSlope(th, tau)
% slope d(tau)/d(th) at th=0, from the inner third of the sweep range
m = abs(th) <= 0.35*max(abs(th));
p = [th(m), ones(sum(m),1)] \ tau(m);
k = p(1);
end
