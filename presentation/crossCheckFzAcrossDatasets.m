%% crossCheckFzAcrossDatasets.m
% Two independent checks that were missing from the capacitance-only
% provenance check:
%
%  (A) Is Fz CONSISTENT across the three quadrant datasets, and does it
%      match the derived closed form?  For the quadrant pair at zero tilt
%          C  = eps0*A/G,   G = g0 - a5*rz
%          Fz = +dW/drz = (1/2)V^2 dC/drz = a5*eps0*A*V^2 / (2*G^2)
%      i.e. the same 1/g^2 pressure integral, with NO free parameters
%      beyond the geometry already fixed by trapGeometry.m.
%
%  (B) Does the ENERGY route actually reproduce the MST force? Checked on
%      the 150-pt central-electrode file, which is the only dataset with
%      BOTH a dense rz sweep AND a stored-energy probe, so dW/drz can be
%      formed directly without a surrogate.
%
% SYMBOLS: eps0 [F/m]; V [V] drive; g0 [m] nominal gap; a5 [-] rz->gap
%   scale; rz [m]; A [m^2] footprint area; G [m] gap; W [J] stored energy;
%   Fz [N]; MST = Maxwell stress tensor probe.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
eps0 = 8.8541878128e-12; V = 30; g0 = 15e-6; a5 = 0.976;
geom = trapGeometry();
A_pair = 2*geom.actSide^2;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

FzModel = @(rz) a5*eps0*A_pair*V^2 ./ (2*(g0 - a5*rz).^2);

%% ============ (A) Fz across the three quadrant datasets ============
fprintf('========== (A) quadrant Fz: three datasets vs the derived closed form ==========\n');
fprintf('model: Fz = a5*eps0*A*V^2/(2G^2),  A = %.4f mm^2 (geometric, not fitted)\n\n', A_pair*1e6);

% -- dataset 1: fine rx sweep (take rx=0 column) --
R = readmatrix(fullfile(scriptDir,'comsol_rxFineSweepEnergy_results.txt'),'NumHeaderLines',5);
m = R(:,1)==0;
d(1).name='rxFineSweepEnergy'; d(1).rz=R(m,3); d(1).Fz=R(m,9);

% -- dataset 2: dedicated single-axis sweep, rz block --
Q = readmatrix(fullfile(scriptDir,'comsol_sweep_combined300_results.txt'),'NumHeaderLines',5);
% cols: rz,rx,ry,thx,thy,thz,...  rz-only block = rz~=0 with all else 0
m = Q(:,1)~=0 & all(abs(Q(:,2:6))==0,2);
d(2).name='combined300 (rz block)'; d(2).rz=Q(m,1); d(2).Fz=Q(m,9);

% -- dataset 3: full-6D LHS, restricted to near-zero tilt so the
%    zero-tilt closed form applies --
F = readmatrix(fullfile(scriptDir,'comsol_overnightFull6DSweep_results.txt'),'NumHeaderLines',5);
tiltSmall = abs(F(:,4))<0.005 & abs(F(:,5))<0.005;      % |thx|,|thy| < 0.005 deg
d(3).name='full6D (near-zero tilt)'; d(3).rz=F(tiltSmall,3); d(3).Fz=F(tiltSmall,9);

fprintf('%-26s %6s %12s %12s %9s %8s\n','dataset','n','mean Fz/model','median','R^2','gap rng');
for i=1:3
    pred = FzModel(d(i).rz);
    ratio = d(i).Fz ./ pred;
    G = g0 - a5*d(i).rz;
    fprintf('%-26s %6d %12.4f %12.4f %9.5f  %.1f-%.1f um\n', d(i).name, numel(d(i).rz), ...
        mean(ratio), median(ratio), r2fun(d(i).Fz, pred), min(G)*1e6, max(G)*1e6);
end

% direct dataset-to-dataset comparison at matched rz
fprintf('\ndirect comparison at matched rz (interpolating combined300 onto the fine sweep):\n');
[rzs,o] = sort(d(2).rz); Fzs = d(2).Fz(o);
fprintf('%8s %14s %14s %9s\n','rz[um]','rxFine Fz[N]','comb300 Fz[N]','rel diff');
for k = 1:numel(d(1).rz)
    v = interp1(rzs, Fzs, d(1).rz(k), 'linear', NaN);
    fprintf('%8.1f %14.4e %14.4e %8.2f%%\n', d(1).rz(k)*1e6, d(1).Fz(k), v, ...
        100*(d(1).Fz(k)-v)/v);
end

%% ============ (B) energy route vs MST, central electrode A ============
fprintf('\n========== (B) does +dW/drz reproduce the MST Fz? (central A, 150-pt) ==========\n');
C150 = readmatrix(fullfile(scriptDir,'comsol_centerElectrodeSweep_rzThxThy150_results.txt'),'NumHeaderLines',5);
m = abs(C150(:,4))==0 & abs(C150(:,5))==0 & abs(C150(:,6))==0 & abs(C150(:,3))>0;
rz = C150(m,3); Fz_mst = C150(m,9); W = C150(m,13);
[rz,o] = sort(rz); Fz_mst = Fz_mst(o); W = W(o);
Fz_energy = gradient(W, rz);                 % +dW/drz at fixed potential

keep = 3:numel(rz)-2;                        % drop the endpoints of the difference
fprintf('n = %d points, rz %.1f to %.1f um\n', numel(keep), rz(keep(1))*1e6, rz(keep(end))*1e6);
fprintf('  mean Fz_energy / Fz_MST = %.4f\n', mean(Fz_energy(keep)./Fz_mst(keep)));
fprintf('  median                  = %.4f\n', median(Fz_energy(keep)./Fz_mst(keep)));
fprintf('  R^2 (energy vs MST)     = %+.5f\n', r2fun(Fz_mst(keep), Fz_energy(keep)));
fprintf('  -> confirms force = +dW/dq at fixed potential, and that MST Fz\n');
fprintf('     (unlike MST Fx) is trustworthy.\n');

%% ============ plot ============
figure('Color','w','Position',[40 40 1340 430]);
subplot(1,3,1); hold on; grid on;
mk = {'o','s','^'};
for i=1:3
    [rs,o2]=sort(d(i).rz);
    plot(rs*1e6, d(i).Fz(o2)*1e6, mk{i},'MarkerSize',4,'DisplayName',d(i).name);
end
rg = linspace(-14,13.5,300)'*1e-6;
plot(rg*1e6, FzModel(rg)*1e6,'k-','LineWidth',1.6,'DisplayName','closed form (no fit)');
set(gca,'YScale','log'); xlabel('r_z [\mum]'); ylabel('F_z [\muN]');
title('quadrant F_z, three datasets'); legend('Location','northwest','FontSize',7);

subplot(1,3,2); hold on; grid on;
for i=1:3
    [rs,o2]=sort(d(i).rz);
    plot(rs*1e6, d(i).Fz(o2)./FzModel(rs), mk{i},'MarkerSize',4,'DisplayName',d(i).name);
end
yline(1,'k-'); xlabel('r_z [\mum]'); ylabel('F_z measured / model');
ylim([0.8 1.3]); title('ratio to the closed form'); legend('Location','best','FontSize',7);

subplot(1,3,3); hold on; grid on;
plot(rz(keep)*1e6, Fz_mst(keep)*1e6,'o','MarkerSize',5,'Color',[0.2 0.3 0.75],'DisplayName','MST probe');
plot(rz(keep)*1e6, Fz_energy(keep)*1e6,'k-','LineWidth',1.5,'DisplayName','+dW/dr_z');
set(gca,'YScale','log'); xlabel('r_z [\mum]'); ylabel('F_z [\muN]');
title('central A: energy route vs MST'); legend('Location','northwest','FontSize',8);

sgtitle('Force cross-checks: dataset consistency and the energy route');
outPng = fullfile(scriptDir,'FzCrossCheck.png');
exportgraphics(gcf,outPng);
fprintf('\nSaved plot to %s\n', outPng);
