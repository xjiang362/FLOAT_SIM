%% plotLevitation60um.m
% Levitation at the WIDENED gap: total channel 60 um, equilibrium g0 = 30 um.
%
% Opening the gap lowers the unstable rate, p = sqrt(2g/g0):
%     g0 = 15 um -> p = 182.0 Hz,  p*t_d = 1.29   (not stabilisable)
%     g0 = 30 um -> p = 128.7 Hz,  p*t_d = 0.91   (inside the bound)
%
% Without acquisition delay this transforms the loop: at 15 um only
% lambda_r >= 100 Hz held the plate, at 30 um every bandwidth from 10 to
% 300 Hz settles. With the delay it still fails, which is the open issue
% documented at the bottom of this file.
%
% COSTS of the wider gap, for the record:
%   hover voltage   24.9 -> 46.7 V        (scales as g0)
%   dC/dz           x0.25                 (scales as 1/g0^2)
%   Ceff            x0.50                 (scales as 1/g0)
% so the position-noise floor degrades by about 4x.
%
% CAVEAT: plant.ceff (C12, Cg, Cf) was extracted from COMSOL at g0 = 15 um.
% The model is internally consistent at any gap -- trapCapToHeight inverts
% trapCapacitance exactly -- but those constants have NOT been re-validated
% against the device at 30 um, so absolute capacitance levels here are
% provisional.
%
% See also TRAPCONTROLDESIGN, SIMULATETRAPLEVITATION, MARGINSTUDY.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p0 = trapDeviceGeometry();
z0 = 30e-6;                       % equilibrium gap; total channel 2*z0
p  = trapDeviceGeometry('z0', z0);   % derived consistently, not patched
g  = p.phys.g;
pz = sqrt(2*g/z0);
td = p.readout.tI2C + p.readout.dtChannel*3.5;

fprintf('total gap %.0f um (g0 = %.0f um)\n', 2*z0*1e6, z0*1e6);
fprintf('p = %.1f rad/s = %.1f Hz,  t_d = %.0f us,  p*t_d = %.3f\n', ...
    pz, pz/2/pi, td*1e6, pz*td);
u0 = trapWrenchToVoltage(p, zeros(3,1), zeros(6,1), 'Vtia', 0, 'Vmax', 1e6);
fprintf('hover voltage %.1f V (was 24.9 V at g0 = 15 um)\n\n', u0(1));

tEnd = 0.25;  Vmax = 200;
cases = { '\lambda_r = 20 Hz, no delay',  20, false; ...
          '\lambda_r = 100 Hz, no delay',100, false; ...
          '\lambda_r = 20 Hz, DELAY',     20, true;  ...
          '\lambda_r = 100 Hz, DELAY',   100, true };
col = [0.00 0.45 0.74; 0.20 0.60 0.30; 0.85 0.33 0.10; 0.60 0.20 0.60];

fprintf('%-32s %18s %14s %10s\n','case','fate','|rz| end','Vpk');
R = struct([]);
for i = 1:size(cases,1)
    K = trapControlDesign(p, 'lambdar', 2*pi*cases{i,2});
    x0 = zeros(12,1); x0(3) = -0.5e-6;
    L = simulateTrapLevitation(p, K, x0, tEnd, 'delay', cases{i,3}, ...
                               'Vtia', 0, 'Vmax', Vmax);
    R(i).L = L; R(i).lab = cases{i,1};
    if isempty(L.contact)
        fate = sprintf('settles %.0f ms', tEnd*1e3);
        fprintf('%-32s %18s %14.3e %10.1f\n', cases{i,1}, fate, ...
            abs(L.x(end,3))*1e6, max(max(L.u(1:L.nTick,1:5))));
    else
        fprintf('%-32s %18s %14s %10.1f\n', cases{i,1}, ...
            sprintf('contact %.1f ms', L.contactTime*1e3), '-', ...
            max(max(L.u(1:L.nTick,1:5))));
    end
end

%% ---------------- plot ----------------
f = figure('Color','w','Position',[40 40 1400 760]);
idx = [3 4 5];
ylb = {'r_z [\mum]','\theta_x [\murad]','\theta_y [\murad]'};
for k = 1:3
    subplot(2,3,k); hold on; grid on; box on;
    for i = 1:numel(R)
        L = R(i).L;
        plot(L.t*1e3, L.x(:,idx(k))*1e6,'-','LineWidth',1.7,'Color',col(i,:), ...
             'DisplayName',R(i).lab);
        if ~isempty(L.contact)
            plot(L.contactTime*1e3, L.x(end,idx(k))*1e6,'x','Color',col(i,:), ...
                 'MarkerSize',12,'LineWidth',2.2,'HandleVisibility','off');
        end
    end
    yline(0,'k:','HandleVisibility','off');
    if k==1
        yline( z0*1e6,'k--','LineWidth',1.1,'HandleVisibility','off');
        yline(-z0*1e6,'k--','LineWidth',1.1,'HandleVisibility','off');
        legend('Location','southwest','FontSize',7);
    end
    xlabel('t [ms]'); ylabel(ylb{k});
    title(ylb{k}(1:regexp(ylb{k},' ','once')-1),'FontSize',10);
end

subplot(2,3,4); hold on; grid on; box on;
for i = 1:numel(R)
    L = R(i).L; plot(L.t*1e3, L.x(:,3)*1e6,'-','LineWidth',1.7,'Color',col(i,:));
end
yline(0,'k:'); xlabel('t [ms]'); ylabel('r_z [\mum]'); ylim([-0.6 0.6]);
title('r_z, zoomed on the setpoint','FontSize',10);

subplot(2,3,5); hold on; grid on; box on;
for i = 1:numel(R)
    L = R(i).L; n = L.nTick;
    plot(L.tTick(1:n)*1e3, max(L.u(1:n,1:5),[],2),'-','LineWidth',1.5,'Color',col(i,:));
end
yline(u0(1),'k--','LineWidth',1.1);
xlabel('t [ms]'); ylabel('peak drive [V]'); title('control effort','FontSize',10);

subplot(2,3,6); hold on; grid on; box on;
gaps = [15 20 25 30 40 50]*1e-6;
plot(2*gaps*1e6, sqrt(2*g./gaps)*td,'o-','LineWidth',2,'Color',[0.2 0.4 0.7]);
yline(1,'r--','LineWidth',1.4);
plot(2*z0*1e6, pz*td,'p','MarkerSize',14,'MarkerFaceColor',[0.85 0.33 0.10], ...
     'MarkerEdgeColor','k');
text(2*z0*1e6+4, pz*td+0.05, 'this design','FontSize',8);
xlabel('total gap [\mum]'); ylabel('p \cdot t_d');
title('stabilisability bound (needs < 1)','FontSize',10); grid on;

sgtitle(sprintf(['Levitation at a %.0f \\mum total gap  \\bullet  ' ...
    'p = %.1f Hz, p t_d = %.2f, hover %.1f V'], ...
    2*z0*1e6, pz/2/pi, pz*td, u0(1)),'FontSize',11);
for ax = findobj(f,'Type','axes')', ax.Toolbar.Visible='off'; end
out = fullfile(root,'presentation','Levitation60um.png');
exportgraphics(f,out,'Resolution',140);
fprintf('\nSaved %s\n', out);

% =====================================================================
function p = widenGap(p0, z0)
% Re-derive everything that depends on the equilibrium gap.
p = p0;
p.gap.z0       = z0;
p.gap.zChannel = 2*z0;
p.gap.zTravel  = z0;
p.damp.c = 4*p0.mass.L^3*p0.env.pressure/(pi^2*p0.env.vbar*z0);
p.damp.q = p0.mass.L^2/12*p.damp.c;
p.damp.C = diag([p.damp.c p.damp.c p.damp.c p.damp.q p.damp.q p.damp.q]);
p.damp.cShear = p.damp.c*(pi^2*z0)/(4*p0.mass.L);
end
