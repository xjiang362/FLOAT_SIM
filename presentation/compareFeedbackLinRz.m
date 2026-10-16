%% compareFeedbackLinRz.m
% One zoomed heave plot, three traces, all at Tisr = 500 us, lambdar = 50 Hz:
%
%   1. no delay, feedback linearisation NOT compensated  (double integrator)
%   2. no delay, feedback linearisation compensated      (exact sampled map)
%   3. DELAY,    feedback linearisation compensated
%
% WHAT "COMPENSATED" MEANS HERE. The wrench allocation in
% trapWrenchToVoltage cancels the electrostatic stiffness exactly, so it is
% tempting to design the regulator on a double integrator. That is wrong.
% Within one ISR the voltage is frozen and the plate obeys
% qddot = p^2*(q - q_k), whose exact one-tick map is
%
%       A_d = [1  sinh(pT)/p ; 0  cosh(pT)]
%
% The inversion cancels the POSITION term (the zero in the corner) but the
% VELOCITY is still amplified by cosh(pT) every tick -- 1.168 for heave at
% Tisr = 500 us, i.e. 17% per tick. Trace 1 is designed as though that
% factor were 1; traces 2 and 3 carry it.
%
% All six DOF are simulated; only heave is plotted, over the capture window.
%
% See also TRAPCONTROLDESIGN, TRAPWRENCHTOVOLTAGE, COMPAREDESIGNMODELS.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p = trapDeviceGeometry();
p.readout.Tisr      = 500e-6;
p.readout.dtChannel = 500e-6/4;

x0 = zeros(12,1);
x0(1) = -40e-6; x0(2) = -40e-6; x0(3) = -13e-6;
x0(4) =  60e-6; x0(5) = -40e-6;
tEnd = 0.6;  Vmax = 120;  lamR = 50;  tZoom = 20;   % [ms]

K0 = trapControlDesign(p,'lambdar',2*pi*lamR);
fprintf('Tisr = %.0f us, lambda_r = %d Hz, lambda_e = %.0f Hz\n', ...
    K0.Tisr*1e6, lamR, K0.lambdae/2/pi);
fprintf('cosh(pT) = %s   <- the factor the double-integrator design omits\n', ...
    sprintf('%.4f ', cosh(K0.pRate*K0.Tisr)));
fprintf('sensing delay: oldest channel %.0f us, w_u*t_d = %.2f\n\n', ...
    (p.readout.tI2C + p.readout.dtChannel*3.5)*1e6, ...
    1143.68*(p.readout.tI2C + p.readout.dtChannel*3.5));

%          design model   delay   label
cases = { 'double', false, 'no delay, FL NOT compensated (double integrator)'; ...
          'exact',  false, 'no delay, FL compensated (exact sampled map)'; ...
          'exact',  true,  'DELAY, FL compensated' };
col = [0.85 0.33 0.10;      % uncompensated
       0.00 0.45 0.74;      % compensated
       0.47 0.67 0.19];     % compensated + delay

fprintf('%-52s %10s %18s %12s\n','case','peak[um]','fate','rz amp[um]');
R = struct([]);
for i = 1:3
    K = trapControlDesign(p,'lambdar',2*pi*lamR,'designModel',cases{i,1});
    L = simulateTrapLevitation(p,K,x0,tEnd,'delay',cases{i,2},'Vtia',0,'Vmax',Vmax);
    if isempty(L.contact)
        sel  = L.t > tEnd/2;
        amp  = (max(L.x(sel,3))-min(L.x(sel,3)))/2*1e6;
        fate = sprintf('survives %.0f ms', tEnd*1e3);
    else
        amp  = NaN;
        fate = sprintf('contact %.1f ms', L.contactTime*1e3);
    end
    R(i).L = L; R(i).lab = cases{i,3}; R(i).amp = amp;
    fprintf('%-52s %10.3f %18s %12.3f\n', cases{i,3}, max(L.x(:,3))*1e6, fate, amp);
end

%% ---------------- figure ----------------
f = figure('Color','w','Position',[60 60 960 580]);
hold on; grid on; box on;
for i = 1:3
    L = R(i).L;  m = L.t*1e3 <= tZoom;
    if isempty(L.contact)
        nm = sprintf('%s  (\\pm%.3f \\mum)', R(i).lab, R(i).amp);
    else
        nm = sprintf('%s  (contact %.1f ms)', R(i).lab, L.contactTime*1e3);
    end
    plot(L.t(m)*1e3, L.x(m,3)*1e6, '-', 'LineWidth',2.1, 'Color',col(i,:), ...
         'DisplayName',nm);
    if ~isempty(L.contact) && L.contactTime*1e3 <= tZoom
        plot(L.contactTime*1e3, L.x(end,3)*1e6, 'x','Color',col(i,:), ...
             'MarkerSize',14,'LineWidth',2.5,'HandleVisibility','off');
    end
end
yline(0,'k:','HandleVisibility','off');
yline( 15,'k--','LineWidth',1.2,'HandleVisibility','off');
yline(-15,'k--','LineWidth',1.2,'HandleVisibility','off');
text(tZoom*0.80, 15.9,'wall','FontSize',8);
xlabel('t [ms]'); ylabel('r_z [\mum]');
xlim([0 tZoom]); ylim([-15.5 17]);
title(sprintf(['Feedback-linearisation compensation in the design model  ' ...
    '\\bullet  T_{ISR} = %.0f \\mus, \\lambda_r = %d Hz'], ...
    p.readout.Tisr*1e6, lamR),'FontSize',11);
legend('Location','southeast','FontSize',8.5);
set(gca,'Toolbar',[]);
out = fullfile(root,'presentation','CompareFeedbackLin_rz.png');
exportgraphics(f,out,'Resolution',140);
fprintf('\nSaved %s\n', out);

%% ---------------- second figure: the two no-delay traces only ----------
% Isolates the design-model question from the latency question. Both runs
% here have perfect, instantaneous sensing, so the ONLY difference is
% whether the regulator gains were placed on [1 T; 0 1] or on the exact
% [1 sinh(pT)/p; 0 cosh(pT)].
f2 = figure('Color','w','Position',[80 80 960 580]);
hold on; grid on; box on;
for i = 1:2
    L = R(i).L;  m = L.t*1e3 <= tZoom;
    if isempty(L.contact)
        nm = sprintf('%s  (\\pm%.3f \\mum)', R(i).lab, R(i).amp);
    else
        nm = sprintf('%s  (contact %.1f ms)', R(i).lab, L.contactTime*1e3);
    end
    plot(L.t(m)*1e3, L.x(m,3)*1e6, '-', 'LineWidth',2.2, 'Color',col(i,:), ...
         'DisplayName',nm);
    if ~isempty(L.contact) && L.contactTime*1e3 <= tZoom
        plot(L.contactTime*1e3, L.x(end,3)*1e6, 'x','Color',col(i,:), ...
             'MarkerSize',14,'LineWidth',2.5,'HandleVisibility','off');
    end
end
yline(0,'k:','HandleVisibility','off');
yline( 15,'k--','LineWidth',1.2,'HandleVisibility','off');
yline(-15,'k--','LineWidth',1.2,'HandleVisibility','off');
text(tZoom*0.80, 15.9,'wall','FontSize',8);
xlabel('t [ms]'); ylabel('r_z [\mum]');
xlim([0 tZoom]); ylim([-15.5 17]);
title(sprintf(['Design model only, no delay  \\bullet  T_{ISR} = %.0f \\mus, ' ...
    '\\lambda_r = %d Hz, cosh(pT) = %.3f'], ...
    p.readout.Tisr*1e6, lamR, cosh(K0.pRate(1)*K0.Tisr)),'FontSize',11);
legend('Location','southeast','FontSize',9);
set(gca,'Toolbar',[]);
out2 = fullfile(root,'presentation','CompareFeedbackLin_rz_noDelay.png');
exportgraphics(f2,out2,'Resolution',140);
fprintf('Saved %s\n', out2);
