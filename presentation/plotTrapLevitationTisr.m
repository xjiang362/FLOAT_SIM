%% plotTrapLevitationTisr.m
% Closed-loop levitation at three sample rates, no acquisition delay.
% Isolates the effect of the zero-order hold: the force-to-voltage
% inversion cancels the negative stiffness at the instant it is computed,
% then the voltage is held for one ISR while the plate keeps moving, so
% the uncancelled stiffness acts for the rest of the tick. The governing
% number is wu*T, with wu = 1144 rad/s the open-loop unstable pole.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);
p0 = trapDeviceGeometry();
wu = 1143.68;

Tlist = [50 200 400]*1e-6;
tEnd  = 0.08;
z0    = -0.5e-6;
col   = [0.20 0.55 0.35; 0.20 0.40 0.70; 0.80 0.30 0.22];

res = struct([]);
for i = 1:numel(Tlist)
    p = p0;
    p.readout.Tisr      = Tlist(i);
    p.readout.dtChannel = Tlist(i)/4;
    p.sim.dtEuler       = Tlist(i)/72;
    K = trapControlDesign(p);
    x0 = zeros(12,1); x0(3) = z0;
    L = simulateTrapLevitation(p, K, x0, tEnd, 'delay', false);
    res(i).T = Tlist(i); res(i).L = L; res(i).K = K;
    res(i).wuT = wu*Tlist(i);
    if isempty(L.contact)
        res(i).fate = sprintf('stable to %.0f ms', tEnd*1e3);
    else
        res(i).fate = sprintf('contact at %.1f ms', L.contactTime*1e3);
    end
    fprintf('Tisr %3.0f us  wu*T = %.3f   %s\n', Tlist(i)*1e6, res(i).wuT, res(i).fate);
end

%% ---- plot ----
f = figure('Color','w','Position',[40 40 1320 760]);

% (a) heave
subplot(2,2,1); hold on; grid on; box on;
for i=1:3
    L=res(i).L;
    plot(L.t*1e3, L.x(:,3)*1e6, '-','LineWidth',1.7,'Color',col(i,:), ...
        'DisplayName',sprintf('T=%.0f us (w_uT=%.2f)',res(i).T*1e6,res(i).wuT));
    if ~isempty(L.contact)
        plot(L.contactTime*1e3, L.x(end,3)*1e6,'x','MarkerSize',13, ...
            'LineWidth',2.5,'Color',col(i,:),'HandleVisibility','off');
    end
end
yline(0,'k:','HandleVisibility','off');
yline(15,'k--','LineWidth',1.2,'DisplayName','contact (gap closed)');
xlabel('t [ms]'); ylabel('r_z [\mum]'); ylim([-2 17]);
title('heave'); legend('Location','northwest','FontSize',8);

% (b) heave, zoomed on the settling band
subplot(2,2,2); hold on; grid on; box on;
for i=1:3
    L=res(i).L;
    plot(L.t*1e3, L.x(:,3)*1e6,'-','LineWidth',1.7,'Color',col(i,:));
end
yline(0,'k:');
xlabel('t [ms]'); ylabel('r_z [\mum]'); ylim([-0.6 0.6]); xlim([0 tEnd*1e3]);
title('heave, zoomed: does it hold the setpoint?');

% (c) tilt
subplot(2,2,3); hold on; grid on; box on;
for i=1:3
    L=res(i).L;
    plot(L.t*1e3, L.x(:,4)*1e6,'-','LineWidth',1.7,'Color',col(i,:));
end
yline(0,'k:');
xlabel('t [ms]'); ylabel('\theta_x [\murad]');
title('tilt: where the slow divergence shows first');

% (d) control effort
subplot(2,2,4); hold on; grid on; box on;
for i=1:3
    L=res(i).L; n=L.nTick;
    plot(L.tTick(1:n)*1e3, L.u(1:n,1),'-','LineWidth',1.5,'Color',col(i,:));
end
yline(29.25,'k--','LineWidth',1.2);
text(tEnd*1e3*0.55, 31,'hover 29.25 V','FontSize',8);
xlabel('t [ms]'); ylabel('quadrant 1 voltage [V]');
title('control effort (TIA excluded: fixed at 40 V)'); ylim([0 70]);

sgtitle(sprintf(['Levitation from r_z = %.1f \\mum, no acquisition delay  ' ...
    '\\bullet  the hold lets the uncancelled stiffness act for one tick'], z0*1e6));
for ax=findobj(f,'Type','axes')', ax.Toolbar.Visible='off'; end
out = fullfile(root,'presentation','TrapLevitation_Tisr.png');
exportgraphics(f,out,'Resolution',140);
fprintf('\nSaved %s\n', out);

