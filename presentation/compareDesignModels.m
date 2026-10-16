%% compareDesignModels.m
% Same plant, same wrench allocation, same estimator -- only the DISCRETE
% MODEL the regulator gains are placed on differs:
%
%   EXACT   A_d = [1  sinh(pT)/p ; 0  cosh(pT)]
%           the true one-tick map when the voltage is held. The V^2
%           inversion cancels the POSITION term (hence the 0 in the corner)
%           but the velocity is still amplified by cosh(pT) every tick.
%
%   DOUBLE  A_d = [1  T ; 0  1]
%           zero-order hold of a double integrator -- what the design used
%           originally, on the assumption that "the inversion cancels the
%           stiffness, so the plant is a double integrator".
%
% At Tisr = 500 us, cosh(pT) = 1.168 (rz), 1.115 (thx), 1.127 (thy): the
% double-integrator model is missing a 12-17% per-tick velocity growth.
%
% Top row    r_z, theta_x, theta_y : the three CONTROLLED degrees of freedom.
% Bottom row r_x, r_y, theta_z     : the three UNCONTROLLED ones, which are
%            passively trapped (11.6 Hz lateral, 5.1 Hz yaw) and essentially
%            undamped. Nothing drives them to the setpoint; they are shown
%            because lateral motion changes the plate/electrode overlap and
%            so leaks into the controlled axes.
%
% See also TRAPCONTROLDESIGN, TRAPWRENCHTOVOLTAGE.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p = trapDeviceGeometry();

% start pose: far lateral offset, low, with a tilt
rng(7);
x0 = zeros(12,1);
x0(1) = -40e-6;  x0(2) = -40e-6;  x0(3) = -13e-6;
x0(4) =  60e-6;  x0(5) = -40e-6;  x0(6) =  0;
fprintf('start: rx=%.0f ry=%.0f rz=%.0f um, thx=%.0f thy=%.0f urad\n\n', ...
    x0(1)*1e6, x0(2)*1e6, x0(3)*1e6, x0(4)*1e6, x0(5)*1e6);

tEnd = 1.0;        % [s] run length
% lambdar = 50 Hz, not the 200 Hz default: at 200 Hz neither design model
% survives this start (both contact at 1.1 ms), because 200 Hz cannot
% capture from rz = -13 um regardless of which model the gains were placed
% on. Comparing two failures says nothing, so the comparison is run at a
% bandwidth where the loop works.
lamR = 50;
mdl  = {'exact','double'};
nmv  = {'exact sampled map','double integrator'};
col  = [0.15 0.45 0.75; 0.85 0.33 0.10];

K0 = trapControlDesign(p, 'lambdar', 2*pi*lamR);
fprintf('Tisr = %.0f us, lamR = %.0f Hz, lamE = %.0f Hz\n', ...
    K0.Tisr*1e6, K0.lambdar/2/pi, K0.lambdae/2/pi);
fprintf('cosh(pT) = %s   (the factor the double-integrator model omits)\n\n', ...
    sprintf('%.4f ', cosh(K0.pRate*K0.Tisr)));

res = struct([]);
for i = 1:2
    K = trapControlDesign(p, 'designModel', mdl{i}, 'lambdar', 2*pi*lamR);
    try
        L = simulateTrapLevitation(p, K, x0, tEnd, 'delay', false, ...
                                   'Vtia', 0, 'Vmax', 120);
    catch ME
        L = []; fprintf('  %s: ERROR %s\n', nmv{i}, ME.identifier);
    end
    res(i).L = L; res(i).K = K;
    if isempty(L.contact)
        fprintf('  %-20s survived %.0f ms: rz=%+9.2e um  thx=%+9.2e urad  Vpk=%5.1f V\n', ...
            nmv{i}, tEnd*1e3, L.x(end,3)*1e6, L.x(end,4)*1e6, ...
            max(max(L.u(1:L.nTick,1:5))));
    else
        fprintf('  %-20s CONTACT at %.1f ms\n', nmv{i}, L.contactTime*1e3);
    end
end

%% ---------------- save, so later zooms need no re-run ----------------
matOut = fullfile(root,'presentation','CompareDesignModels.mat');
save(matOut, 'res', 'x0', 'tEnd', 'mdl', 'nmv', 'col', '-v7.3');
fprintf('Saved %s\n', matOut);

%% ---------------- plot ----------------
f = figure('Color','w','Position',[30 30 1500 760]);
idx   = [3 4 5 1 2 6];                       % rz thx thy | rx ry thz
ttl   = {'r_z  (controlled)','\theta_x  (controlled)','\theta_y  (controlled)', ...
         'r_x  (uncontrolled)','r_y  (uncontrolled)','\theta_z  (uncontrolled)'};
ylb   = {'r_z [\mum]','\theta_x [\murad]','\theta_y [\murad]', ...
         'r_x [\mum]','r_y [\mum]','\theta_z [\murad]'};
sc    = [1e6 1e6 1e6 1e6 1e6 1e6];

for k = 1:6
    subplot(2,3,k); hold on; grid on; box on;
    for i = 1:2
        L = res(i).L; if isempty(L), continue; end
        plot(L.t*1e3, L.x(:,idx(k))*sc(k), '-', 'LineWidth',1.8, ...
             'Color',col(i,:), 'DisplayName',nmv{i});
        if ~isempty(L.contact)
            plot(L.contactTime*1e3, L.x(end,idx(k))*sc(k), 'x', ...
                 'Color',col(i,:),'MarkerSize',13,'LineWidth',2.5, ...
                 'HandleVisibility','off');
        end
    end
    yline(0,'k:','HandleVisibility','off');
    if k == 1
        yline( 15,'k--','LineWidth',1.1,'HandleVisibility','off');
        yline(-15,'k--','LineWidth',1.1,'HandleVisibility','off');
        legend('Location','southeast','FontSize',8);
    end
    if k == 4 || k == 5
        yline( 50,'k--','LineWidth',1.1,'HandleVisibility','off');
        yline(-50,'k--','LineWidth',1.1,'HandleVisibility','off');
    end
    xlabel('t [ms]'); ylabel(ylb{k}); title(ttl{k},'FontSize',10);
end

sgtitle(sprintf(['Regulator design model: exact sampled map vs double integrator ' ...
    '\\bullet  T_{ISR} = %.0f \\mus, \\lambda_r = %.0f Hz, \\lambda_e = %.0f Hz, ' ...
    'start (%.0f, %.0f, %.0f) \\mum'], K0.Tisr*1e6, K0.lambdar/2/pi, ...
    K0.lambdae/2/pi, x0(1)*1e6, x0(2)*1e6, x0(3)*1e6),'FontSize',11);
for ax = findobj(f,'Type','axes')', ax.Toolbar.Visible = 'off'; end
out = fullfile(root,'presentation','CompareDesignModels.png');
exportgraphics(f, out, 'Resolution',140);
fprintf('\nSaved %s\n', out);
