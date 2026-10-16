%% compareLoopSettings.m
% Three head-to-head comparisons on the same plant and the same start pose:
%
%   (1) regulator bandwidth   lambdar = 50 vs 200 Hz   (Tisr 500 us, no delay)
%   (2) sample period         Tisr = 100 vs 500 us     (lambdar 200 Hz, no delay)
%   (3) acquisition delay     off vs on                (Tisr 500 us, both lambdar)
%
% Layout per figure: top row the three CONTROLLED axes (r_z, theta_x,
% theta_y), bottom row the three UNCONTROLLED ones (r_x, r_y, theta_z).
% The uncontrolled axes are passively trapped and essentially undamped, so
% they ring for the whole run and leak into the controlled axes through the
% overlap term -- that leak, not the regulator, sets the steady error from
% a laterally offset start.
%
% DELAY. plant.readout gives age(k) = tI2C + dtChannel*(n-k+0.5), so at
% Tisr = 500 us the oldest channel is 937 us old, and the loop also carries
% one ISR of hold. Against an unstable pole at w_u = 1144 rad/s that is
% w_u*t_d = 1.07 for sensing alone and 1.64 including the hold. The
% stabilisability limit is w_u*t_d ~ 1, so this is expected to fail; the
% comparison quantifies how fast.
%
% See also SIMULATETRAPLEVITATION, TRAPCONTROLDESIGN, TRAPSAMPLERINIT.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p0 = trapDeviceGeometry();
x0 = zeros(12,1);
x0(1) = -40e-6; x0(2) = -40e-6; x0(3) = -13e-6;
x0(4) =  60e-6; x0(5) = -40e-6;
tEnd  = 0.6;
Vmax  = 120;

wu    = 1143.68;
fprintf('start (%.0f, %.0f, %.0f) um, tilt (%.0f, %.0f) urad, %.0f ms\n\n', ...
    x0(1)*1e6, x0(2)*1e6, x0(3)*1e6, x0(4)*1e6, x0(5)*1e6, tEnd*1e3);

%           label                      Tisr    lamR  delay
% The ISR comparison runs at lambdar = 50 Hz because 200 Hz does not
% survive this start at ANY sample period, so comparing 100 vs 500 us there
% would only compare two failures.
cases = { ...
  '\lambda_r = 50 Hz',                500e-6,   50, false; ...
  '\lambda_r = 200 Hz',               500e-6,  200, false; ...
  'T_{ISR} = 100 \mus',               100e-6,   50, false; ...
  'T_{ISR} = 500 \mus',               500e-6,   50, false; ...
  '50 Hz, no delay',                  500e-6,   50, false; ...
  '50 Hz, DELAY',                     500e-6,   50, true;  ...
  '200 Hz, no delay',                 500e-6,  200, false; ...
  '200 Hz, DELAY',                    500e-6,  200, true};

R = struct([]);
fprintf('%-22s %8s %7s %7s %16s %11s %11s %8s\n', ...
    'case','Tisr[us]','lamR','delay','fate','rz amp[um]','thx amp[ur]','Vpk');
for i = 1:size(cases,1)
    p = p0;
    p.readout.Tisr      = cases{i,2};
    p.readout.dtChannel = cases{i,2}/4;
    K = trapControlDesign(p, 'lambdar', 2*pi*cases{i,3});
    try
        L = simulateTrapLevitation(p, K, x0, tEnd, 'delay', cases{i,4}, ...
                                   'Vtia', 0, 'Vmax', Vmax);
    catch ME
        L = []; fprintf('  %s ERROR %s\n', cases{i,1}, ME.identifier);
    end
    R(i).L = L; R(i).lab = cases{i,1};
    if isempty(L)
        fate = 'ERROR'; az = NaN; at = NaN; vp = NaN;
    elseif isempty(L.contact)
        fate = sprintf('survives %.0fms', tEnd*1e3);
        sel = L.t > tEnd/2;
        az = (max(L.x(sel,3))-min(L.x(sel,3)))/2*1e6;
        at = (max(L.x(sel,4))-min(L.x(sel,4)))/2*1e6;
        vp = max(max(L.u(1:L.nTick,1:5)));
    else
        fate = sprintf('contact %.1f ms', L.contactTime*1e3);
        az = NaN; at = NaN; vp = max(max(L.u(1:L.nTick,1:5)));
    end
    dl = {'no','YES'};
    fprintf('%-22s %8.0f %7.0f %7s %16s %11.3f %11.2f %8.1f\n', ...
        cases{i,1}, cases{i,2}*1e6, cases{i,3}, dl{cases{i,4}+1}, fate, az, at, vp);
end

fprintf('\nsensing delay at Tisr = 500 us: oldest channel %.0f us, w_u*t_d = %.2f\n', ...
    (p0.readout.tI2C + 125e-6*3.5)*1e6, wu*(p0.readout.tI2C + 125e-6*3.5));
fprintf('  including one ISR of hold: w_u*t_d = %.2f  (limit is ~1)\n', ...
    wu*(p0.readout.tI2C + 125e-6*3.5 + 500e-6));

%% ---------------- figures ----------------
grp = {[1 2], 'RegulatorBandwidth', '\lambda_r = 50 vs 200 Hz  (T_{ISR} = 500 \mus, no delay)'; ...
       [3 4], 'SamplePeriod',       'T_{ISR} = 100 vs 500 \mus  (\lambda_r = 200 Hz, no delay)'; ...
       [5 6 7 8], 'Delay',          'acquisition delay off vs on  (T_{ISR} = 500 \mus)'};
for g = 1:3
    makeFig(R, grp{g,1}, grp{g,3}, x0, ...
        fullfile(root,'presentation',['CompareLoop_' grp{g,2} '.png']));
end

% =====================================================================
function makeFig(R, sel, ttlTxt, x0, outFile)
col = lines(numel(sel));
idx = [3 4 5 1 2 6];
ttl = {'r_z  (controlled)','\theta_x  (controlled)','\theta_y  (controlled)', ...
       'r_x  (uncontrolled)','r_y  (uncontrolled)','\theta_z  (uncontrolled)'};
ylb = {'r_z [\mum]','\theta_x [\murad]','\theta_y [\murad]', ...
       'r_x [\mum]','r_y [\mum]','\theta_z [\murad]'};
f = figure('Color','w','Position',[25 25 1500 760]);
for k = 1:6
    subplot(2,3,k); hold on; grid on; box on;
    for j = 1:numel(sel)
        L = R(sel(j)).L; if isempty(L), continue; end
        plot(L.t*1e3, L.x(:,idx(k))*1e6, '-', 'LineWidth',1.7, ...
            'Color',col(j,:), 'DisplayName',R(sel(j)).lab);
        if ~isempty(L.contact)
            plot(L.contactTime*1e3, L.x(end,idx(k))*1e6, 'x','Color',col(j,:), ...
                'MarkerSize',12,'LineWidth',2.2,'HandleVisibility','off');
        end
    end
    yline(0,'k:','HandleVisibility','off');
    if k==1, yline(15,'k--','LineWidth',1.0,'HandleVisibility','off');
             yline(-15,'k--','LineWidth',1.0,'HandleVisibility','off');
             legend('Location','best','FontSize',7); end
    if k==4 || k==5, yline(50,'k--','LineWidth',1.0,'HandleVisibility','off');
                     yline(-50,'k--','LineWidth',1.0,'HandleVisibility','off'); end
    xlabel('t [ms]'); ylabel(ylb{k}); title(ttl{k},'FontSize',10);
end
sgtitle([ttlTxt sprintf('  \\bullet  start (%.0f, %.0f, %.0f) \\mum', ...
    x0(1)*1e6, x0(2)*1e6, x0(3)*1e6)],'FontSize',11);
for ax=findobj(f,'Type','axes')', ax.Toolbar.Visible='off'; end
exportgraphics(f, outFile, 'Resolution',140);
fprintf('Saved %s\n', outFile);
end
