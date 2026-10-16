%% plotTrapLevitationDesign.m
% Closed-loop response at the CURRENT design point:
%   lambdar = 200 Hz, lambdae = 500 Hz, Tisr = 500 us, no acquisition delay,
%   TIA grounded.
%
% Also shows why the region of attraction ends where it does. The local
% unstable pole depends on the gap,
%       k_e(g) = 2W/g,    w_u(g) = sqrt(2W/(g*m))   ~ g^(-1/2)
% so as the plate approaches the electrodes w_u rises and eventually
% exceeds the regulator bandwidth, at which point no amount of voltage
% helps. Going the other way the gap opens, w_u falls, and the limit is
% instead the VOLTAGE needed to carry the weight at that gap.
%
% See also TRAPCONTROLDESIGN, SIMULATETRAPLEVITATION.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p = trapDeviceGeometry();
K = trapControlDesign(p);
W = p.mass.m*p.phys.g;
fprintf('design point: lamR = %.0f Hz, lamE = %.0f Hz, Tisr = %.0f us\n', ...
    K.lambdar/2/pi, K.lambdae/2/pi, K.Tisr*1e6);

%% ---- local unstable pole and hover voltage vs height ----
rzs = linspace(-14, 14, 200)*1e-6;
wu  = sqrt(2*W./((p.gap.z0 - rzs)*p.mass.m));
Vh  = nan(size(rzs));
for i = 1:numel(rzs)
    q = zeros(6,1); q(3) = rzs(i);
    [~,Jc] = trapCapacitance(p,q); [~,~,Bq] = trapWrench(p,Jc,zeros(6,1));
    Vh(i) = sqrt(W/sum(Bq(3,1:5)));
end
fprintf('\n r_z[um]   gap[um]   w_u[Hz]   hover V\n');
for z = [-14 -8 0 8 12 14]*1e-6
    [~,i] = min(abs(rzs-z));
    fprintf(' %+7.1f %9.1f %9.1f %9.2f\n', z*1e6, (p.gap.z0-z)*1e6, wu(i)/2/pi, Vh(i));
end
fprintf('\n regulator bandwidth is %.0f Hz; it is exceeded above r_z = %+.1f um\n', ...
    K.lambdar/2/pi, interp1(wu, rzs, K.lambdar)*1e6);

%% ---- trajectories ----
z0s = [-13 -8 -4 0 4 8]*1e-6;
col = parula(numel(z0s)+1);
res = struct([]);
for i = 1:numel(z0s)
    x0 = zeros(12,1); x0(3) = z0s(i);
    res(i).L = simulateTrapLevitation(p,K,x0,0.15,'delay',false,'Vtia',0);
    res(i).z0 = z0s(i);
end

f = figure('Color','w','Position',[30 30 1400 760]);

subplot(2,3,1); hold on; grid on; box on;
for i=1:numel(z0s)
    L=res(i).L;
    plot(L.t*1e3, L.x(:,3)*1e6,'-','LineWidth',1.7,'Color',col(i,:), ...
        'DisplayName',sprintf('r_z(0) = %+.0f \\mum',z0s(i)*1e6));
end
yline(0,'k:','HandleVisibility','off');
yline(15,'k--','LineWidth',1.1,'HandleVisibility','off');
yline(-15,'k--','LineWidth',1.1,'HandleVisibility','off');
xlabel('t [ms]'); ylabel('r_z [\mum]'); title('heave'); legend('Location','east','FontSize',7);

subplot(2,3,2); hold on; grid on; box on;
for i=1:numel(z0s), L=res(i).L; plot(L.t*1e3, L.x(:,3)*1e6,'-','LineWidth',1.7,'Color',col(i,:)); end
yline(0,'k:'); xlabel('t [ms]'); ylabel('r_z [\mum]');
ylim([-0.3 0.3]); title('heave, zoom on the setpoint');

subplot(2,3,3); hold on; grid on; box on;
for i=1:numel(z0s), L=res(i).L; plot(L.t*1e3, L.x(:,4)*1e6,'-','LineWidth',1.7,'Color',col(i,:)); end
yline(0,'k:'); xlabel('t [ms]'); ylabel('\theta_x [\murad]'); title('tilt');

subplot(2,3,4); hold on; grid on; box on;
for i=1:numel(z0s)
    L=res(i).L; n=L.nTick;
    plot(L.tTick(1:n)*1e3, max(L.u(1:n,1:5),[],2),'-','LineWidth',1.5,'Color',col(i,:));
end
yline(80,'r--','LineWidth',1.2); text(2,83,'V_{max} = 80 V','Color','r','FontSize',8);
xlabel('t [ms]'); ylabel('peak drive [V]'); title('control effort'); ylim([0 95]);

subplot(2,3,5); hold on; grid on; box on;
plot(rzs*1e6, wu/2/pi,'-','LineWidth',2,'Color',[0.2 0.4 0.7],'DisplayName','local \omega_u');
yline(K.lambdar/2/pi,'r--','LineWidth',1.4,'DisplayName','\lambda_r = 200 Hz');
xlabel('r_z [\mum]'); ylabel('[Hz]'); title('local unstable pole vs height');
legend('Location','northwest','FontSize',8); ylim([0 800]);

subplot(2,3,6); hold on; grid on; box on;
plot(rzs*1e6, Vh,'-','LineWidth',2,'Color',[0.2 0.5 0.3],'DisplayName','hover voltage');
yline(80,'r--','LineWidth',1.4,'DisplayName','V_{max} = 80 V');
xlabel('r_z [\mum]'); ylabel('[V]'); title('voltage needed to carry the weight');
legend('Location','northeast','FontSize',8);

sgtitle(sprintf(['Design point: \\lambda_r = %.0f Hz, \\lambda_e = %.0f Hz, ' ...
    'T_{ISR} = %.0f \\mus, no delay, V_{TIA} = 0'], ...
    K.lambdar/2/pi, K.lambdae/2/pi, K.Tisr*1e6),'FontSize',12);
for ax=findobj(f,'Type','axes')', ax.Toolbar.Visible='off'; end
out = fullfile(root,'presentation','TrapLevitationDesign.png');
exportgraphics(f,out,'Resolution',140);
fprintf('\nSaved %s\n', out);
