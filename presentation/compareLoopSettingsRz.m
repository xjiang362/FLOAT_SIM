%% compareLoopSettingsRz.m
% One figure, four traces: lambdar = 50 Hz throughout, crossing
%   Tisr  = 100 vs 500 us
%   delay = off vs on
%
% All six DOF are simulated; only heave is plotted, over the first tZoom
% milliseconds, which is the capture transient where the settings differ.
%
% NOTE THE DELAY IS NOT INDEPENDENT OF THE ISR. plant.readout sets
% dtChannel = Tisr/4, so the per-channel staleness
%       age(k) = tI2C + dtChannel*(n - k + 0.5)
% shrinks with the sample period: the oldest channel is 938 us old at
% Tisr = 500 us but only 588 us at Tisr = 100 us, because tI2C = 500 us is
% the only fixed part. Against the unstable pole w_u = 1144 rad/s that is
% w_u*t_d = 1.07 versus 0.67 -- the first is past the stabilisability limit
% and the second is not. So shortening the ISR buys delay margin as well as
% bandwidth, and the four traces are not a clean 2x2 of independent knobs.
%
% See also COMPARELOOPSETTINGS, TRAPSAMPLERINIT, SIMULATETRAPLEVITATION.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p0 = trapDeviceGeometry();
x0 = zeros(12,1);
x0(1) = -40e-6; x0(2) = -40e-6; x0(3) = -13e-6;
x0(4) =  60e-6; x0(5) = -40e-6;
tEnd  = 0.6;  Vmax = 120;  lamR = 50;  tZoom = 20;   % [ms]
wu    = 1143.68;

%          Tisr     delay
cases = { 100e-6, false; ...
          500e-6, false; ...
          100e-6, true;  ...
          500e-6, true };
col = [0.00 0.45 0.74;    % 100 us, no delay
       0.85 0.33 0.10;    % 500 us, no delay
       0.47 0.67 0.19;    % 100 us, delay
       0.49 0.18 0.56];   % 500 us, delay

fprintf('start (%.0f, %.0f, %.0f) um, tilt (%.0f, %.0f) urad, lambda_r = %d Hz\n\n', ...
    x0(1)*1e6,x0(2)*1e6,x0(3)*1e6,x0(4)*1e6,x0(5)*1e6,lamR);
fprintf('%-28s %10s %9s %10s %16s %12s\n', ...
    'case','oldest[us]','w_u*t_d','peak[um]','fate','rz amp[um]');

R = struct([]);
for i = 1:4
    T  = cases{i,1};  dly = cases{i,2};
    p  = p0;
    p.readout.Tisr      = T;
    p.readout.dtChannel = T/4;
    K  = trapControlDesign(p,'lambdar',2*pi*lamR);
    L  = simulateTrapLevitation(p,K,x0,tEnd,'delay',dly,'Vtia',0,'Vmax',Vmax);

    tOld = p.readout.tI2C + (T/4)*3.5;
    dl   = {'no delay','DELAY'};
    lab  = sprintf('T_{ISR} = %3.0f \\mus, %s', T*1e6, dl{dly+1});
    if isempty(L.contact)
        sel  = L.t > tEnd/2;
        amp  = (max(L.x(sel,3))-min(L.x(sel,3)))/2*1e6;
        fate = sprintf('survives %.0f ms', tEnd*1e3);
    else
        amp  = NaN;
        fate = sprintf('contact %.1f ms', L.contactTime*1e3);
    end
    R(i).L = L; R(i).lab = lab; R(i).amp = amp; R(i).fate = fate;

    if dly, wt = sprintf('%.2f', wu*tOld); ot = sprintf('%.0f', tOld*1e6);
    else,   wt = '-';                      ot = '-';  end
    fprintf('%-28s %10s %9s %10.3f %16s %12.3f\n', ...
        sprintf('Tisr=%3.0fus, %s', T*1e6, dl{dly+1}), ot, wt, ...
        max(L.x(:,3))*1e6, fate, amp);
end

fprintf('\n(tI2C = %.0f us is fixed; dtChannel = Tisr/4, so the delay shrinks\n', ...
    p0.readout.tI2C*1e6);
fprintf(' with the sample period. The stabilisability limit is w_u*t_d ~ 1.)\n');

%% ---------------- one figure, four traces ----------------
f = figure('Color','w','Position',[60 60 900 560]);
hold on; grid on; box on;
for i = 1:4
    L = R(i).L;
    m = L.t*1e3 <= tZoom;
    if isempty(L.contact)
        nm = sprintf('%s  (\\pm%.3f \\mum)', R(i).lab, R(i).amp);
    else
        nm = sprintf('%s  (contact %.1f ms)', R(i).lab, L.contactTime*1e3);
    end
    plot(L.t(m)*1e3, L.x(m,3)*1e6, '-', 'LineWidth',2.0, 'Color',col(i,:), ...
         'DisplayName',nm);
    if ~isempty(L.contact) && L.contactTime*1e3 <= tZoom
        plot(L.contactTime*1e3, L.x(end,3)*1e6, 'x', 'Color',col(i,:), ...
             'MarkerSize',14,'LineWidth',2.5,'HandleVisibility','off');
    end
end
yline(0,'k:','HandleVisibility','off');
yline( 15,'k--','LineWidth',1.2,'HandleVisibility','off');
yline(-15,'k--','LineWidth',1.2,'HandleVisibility','off');
text(tZoom*0.72, 15.9, 'wall', 'FontSize',8);
xlabel('t [ms]'); ylabel('r_z [\mum]');
xlim([0 tZoom]); ylim([-15.5 17]);
title(sprintf(['Heave capture, \\lambda_r = %d Hz  \\bullet  ' ...
    'start (%.0f, %.0f, %.0f) \\mum, tilt (%.0f, %.0f) \\murad'], ...
    lamR, x0(1)*1e6,x0(2)*1e6,x0(3)*1e6,x0(4)*1e6,x0(5)*1e6),'FontSize',11);
legend('Location','southeast','FontSize',9);
set(gca,'Toolbar',[]);
out = fullfile(root,'presentation','CompareLoop_rzOnly.png');
exportgraphics(f,out,'Resolution',140);
fprintf('\nSaved %s\n', out);
