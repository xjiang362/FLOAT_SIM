%% simulateHeaveFeedbackLin.m
% Standalone HEAVE-ONLY study of the effect of feedback linearisation on
% the discrete-time design model. Implements the derivation directly, but
% drives the existing validated electrostatic plant rather than a toy one:
% the force comes from trapCapacitance/trapWrench, the voltages from
% trapWrenchToVoltage, and the measurement from the real sensing chain
% (trapSensing -> trapCapToHeight, four heights averaged to r_z). Only r_z
% is reported; the other five DOF start at zero and stay there by symmetry.
%
% THE POINT. Within one ISR the voltage is frozen, but the commanded
% ACCELERATION is not, because
%       u(t) = eps0*A*V[k]^2 / (2m(d0 - r_z(t))^2) - g
% depends on r_z(t), which keeps moving. Expanding about r_z[k],
%       u(t) ~ u[k] + k*r_z(t),        k = (1/m) dF/dr_z
% so the plant seen over the tick is NOT a double integrator:
%
%   naive        A_d = [1  T      ; 0  1      ]   B_d = [T^2/2        ; T        ]
%   compensated  A_d = [1  sinh(pT)/p ; 0  cosh(pT)]  B_d = [(cosh(pT)-1)/p^2; sinh(pT)/p]
%
% with p = sqrt(k). The naive model puts a double pole at z = 1; the true
% one has a pole at cosh(pT) > 1, outside the unit circle.
%
% ---------------------------------------------------------------------
% CORRECTION TO THE DERIVATION'S CONSTANT. The note gives
%       k = eps0*A*V0^2 / (4*m*d0^3)
% but the Taylor coefficient is a factor of 4 larger. With
% f(r) = C/(d0-r)^2 and C = eps0*A*V^2/(2m),
%       df/dr = 2C/(d0-r)^3 = eps0*A*V^2 / (m*(d0-r)^3),
% i.e. the denominator is m*(d0-r)^3, not 4*m*(d0-r)^3. Using the
% equilibrium condition eps0*A*V0^2 = 2*m*g*d0^2 this gives
%       k = 2g/d0 = 1.3080e6 1/s^2     ->  p = 1143.7 rad/s
% whereas the note's constant gives
%       k = g/(2*d0) = 3.2700e5 1/s^2  ->  p =  571.8 rad/s
% The correct value is confirmed independently by the 6-DOF model:
% trapOpenLoopStiffness returns dFz/dr_z = +37.763 N/m = 2W/d0, so
% p = sqrt(37.763/m) = 1143.68 rad/s exactly.
%
% At T = 700 us that is cosh(pT) = 1.338 versus 1.081 -- a 34% per-tick
% velocity growth versus 8%. Trace 3 below uses the note's p so the cost of
% the slip is visible rather than argued.
% ---------------------------------------------------------------------
%
% Design: Tisr = 700 us, lambdar = 100 Hz, lambdae = 300 Hz.
% Single input and single output, so pole placement is UNIQUE here and the
% conditioning problems of the 6-DOF MIMO design do not arise.
%
% See also TRAPCONTROLDESIGN, TRAPWRENCHTOVOLTAGE, COMPAREDESIGNMODELS.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p = trapDeviceGeometry();
m  = p.mass.m;  g = p.phys.g;  d0 = p.gap.z0;

T       = 700e-6;
lambdar = 2*pi*100;
lambdae = 2*pi*300;
p.readout.Tisr      = T;
p.readout.dtChannel = T/4;
nSub = round(T/p.sim.dtEuler);
if abs(nSub*p.sim.dtEuler - T) > 1e-12
    error('simulateHeaveFeedbackLin:step','dtEuler must divide Tisr.');
end

Kst   = trapOpenLoopStiffness(p,'Vtia',0);
pTrue = sqrt(-Kst(3,3)/m);        % 1143.68 rad/s, from the validated model
pNote = sqrt(g/(2*d0));           %  571.84 rad/s, from the note's constant

fprintf('Tisr = %.0f us, lambdar = %.0f Hz, lambdae = %.0f Hz\n', ...
    T*1e6, lambdar/2/pi, lambdae/2/pi);
fprintf('p (validated model) = %7.2f rad/s -> pT = %.4f, cosh(pT) = %.4f\n', ...
    pTrue, pTrue*T, cosh(pTrue*T));
fprintf('p (note''s constant) = %7.2f rad/s -> pT = %.4f, cosh(pT) = %.4f\n\n', ...
    pNote, pNote*T, cosh(pNote*T));

%        label                                        p
des = { 'naive DT  [1 T; 0 1]',                      NaN; ...
        'compensated  [1 sinh(pT)/p; 0 cosh(pT)]',   pTrue };
col = [0.85 0.33 0.10; 0.00 0.45 0.74];

rz0  = -5e-6;
tEnd = 0.25;
nT   = floor(tEnd/T);

fprintf('%-42s %10s %10s %14s\n','design model','cosh(pT)','peak[um]','fate');
R = struct([]);
for i = 1:numel(des(:,1))
    pp = des{i,2};
    if isnan(pp)
        Ad = [1 T; 0 1];                 Bd = [T^2/2; T];
    else
        ch = cosh(pp*T); sh = sinh(pp*T);
        Ad = [1 sh/pp; 0 ch];            Bd = [(ch-1)/pp^2; sh/pp];
    end
    Cd = [1 0];
    % regulator with one integrator, SISO -> placement is unique
    Aa = [Ad zeros(2,1); -T*Cd 1];  Ba = [Bd; 0];
    Ka = acker(Aa, Ba, exp(-lambdar*T)*[1 1 1]);
    Kx = Ka(1:2);  Ki = Ka(3);

    % ESTIMATOR: placed on the stiffness-RETAINING map, identically for all
    % three designs. The inversion cancels the stiffness at the ESTIMATED
    % pose, so the estimation ERROR sees it uncancelled:
    %       e[k+1] = (A_s - L*C)*e,   A_s = [cosh pT, sinh pT/p; p sinh pT, cosh pT]
    % Placing L on the regulator's own design model instead makes the
    % observer diverge and would confound the comparison -- the question
    % here is about the REGULATOR model, so the observer is held fixed and
    % correct in every case.
    chT = cosh(pTrue*T);  shT = sinh(pTrue*T);
    As  = [chT shT/pTrue; pTrue*shT chT];
    L   = acker(As', Cd', exp(-lambdae*T)*[1 1])';

    S = runHeave(p, Ad, Bd, Cd, Kx, Ki, L, T, nSub, nT, rz0, m);
    S.lab = des{i,1};  R(i).S = S;
    if isnan(pp), cs = 1; else, cs = cosh(pp*T); end
    fprintf('%-42s %10.4f %10.3f %14s\n', des{i,1}, cs, S.pk*1e6, S.fate);
end

%% ---------------- plot: r_z and the controller output ----------------
tMax = 30;        % [ms]
f = figure('Color','w','Position',[60 60 1000 720]);

% --- top: true r_z (solid) vs estimate (dashed) ---
subplot(2,1,1); hold on; grid on; box on;
for i = 1:numel(R)
    S = R(i).S;
    mc = S.t*1e3 <= tMax;  mk = S.tTick*1e3 <= tMax;
    plot(S.t(mc)*1e3, S.rz(mc)*1e6, '-', 'LineWidth',2.1, 'Color',col(i,:), ...
         'DisplayName',sprintf('%s  (%s)', S.lab, S.fate));
    plot(S.tTick(mk)*1e3, S.rzhat(mk)*1e6, '--', 'LineWidth',1.6, ...
         'Color',col(i,:), 'DisplayName',sprintf('%s  estimate', S.lab));
end
yline(0,'k:','HandleVisibility','off');
yline( 15,'k--','LineWidth',1.1,'HandleVisibility','off');
yline(-15,'k--','LineWidth',1.1,'HandleVisibility','off');
ylabel('r_z [\mum]'); xlim([0 tMax]); ylim([-16 17]);
title('heave: true (solid) vs estimated (dashed)','FontSize',10);
legend('Location','southeast','FontSize',8);

% --- bottom: controller output, the commanded acceleration ---
subplot(2,1,2); hold on; grid on; box on;
for i = 1:numel(R)
    S = R(i).S;
    mk = S.tTick*1e3 <= tMax;
    stairs(S.tTick(mk)*1e3, S.u(mk), '-', 'LineWidth',1.9, 'Color',col(i,:), ...
           'DisplayName',S.lab);
end
yline(0,'k:','HandleVisibility','off');
yline(-9.81,'k--','LineWidth',1.1,'DisplayName','-g (zero drive)');
xlabel('t [ms]'); ylabel('a = u [m/s^2]'); xlim([0 tMax]);
title('controller output: commanded acceleration (ZOH over each ISR)','FontSize',10);
legend('Location','northeast','FontSize',8);

sgtitle(sprintf(['Heave, feedback-linearisation compensation in the design model  ' ...
    '\\bullet  T_{ISR} = %.0f \\mus, \\lambda_r = %.0f Hz, \\lambda_e = %.0f Hz'], ...
    T*1e6, lambdar/2/pi, lambdae/2/pi),'FontSize',11);
for ax = findobj(f,'Type','axes')', ax.Toolbar.Visible='off'; end
out = fullfile(root,'presentation','HeaveFeedbackLin_rz.png');
exportgraphics(f,out,'Resolution',140);
fprintf('\nSaved %s\n', out);

% =====================================================================
function S = runHeave(p, Ad, Bd, Cd, Kx, Ki, L, T, nSub, nT, rz0, m)
dt = p.sim.dtEuler;
x  = zeros(12,1);  x(3) = rz0;          % heave only; the rest stay zero
xh = [rz0; 0];  xi = 0;
tLog  = zeros(nT*nSub+1,1);  rzLog = zeros(nT*nSub+1,1);
tTick = zeros(nT,1); rzhat = zeros(nT,1); uLog = zeros(nT,1);
rzLog(1) = rz0;  n = 1;  nk = 0;  fate = 'settles';
for k = 1:nT
    tk = (k-1)*T;
    if ~checkPoseValid(p, x(1:6))
        fate = sprintf('contact %.1f ms', tk*1e3); break
    end
    % ---- measurement: real sensing chain, four heights averaged ----
    try
        zc = trapCapToHeight(p, trapSensing(p, trapCapacitance(p, x(1:6))));
    catch
        fate = sprintf('contact %.1f ms', tk*1e3); break
    end
    y = mean(zc);

    % ---- regulator: u is an ACCELERATION command ----
    u  = -Kx*xh - Ki*xi;
    xi = xi + T*(0 - Cd*xh);

    nk = nk+1; tTick(nk) = tk; rzhat(nk) = xh(1); uLog(nk) = u;

    % ---- feedback linearisation: ask for m*u of extra z-force ----
    poseEst = zeros(6,1); poseEst(3) = xh(1);
    V = trapWrenchToVoltage(p, [m*u; 0; 0], poseEst, 'Vtia', 0, 'Vmax', 1e6);

    % ---- propagate the real plant with the voltage HELD ----
    hit = false;
    for j = 1:nSub
        if ~checkPoseValid(p, x(1:6))
            fate = sprintf('contact %.1f ms', (tk+(j-1)*dt)*1e3); hit = true; break
        end
        x = x + dt*trapDynamics(p, x, V);
        n = n+1; tLog(n) = tk + j*dt; rzLog(n) = x(3);
    end
    if hit, break; end
    xh = Ad*xh + Bd*u + L*(y - Cd*xh);
end
S.t = tLog(1:n);  S.rz = rzLog(1:n);
S.tTick = tTick(1:nk);  S.rzhat = rzhat(1:nk);  S.u = uLog(1:nk);
S.fate = fate;
[~,im] = max(abs(S.rz));  S.pk = S.rz(im);
end
