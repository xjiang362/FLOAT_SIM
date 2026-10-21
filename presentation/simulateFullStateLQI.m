%% simulateFullStateLQI.m
% Closed-loop simulation of the delay-aware full-state LQI design from
% trapFullStateModel, driving the real nonlinear plant.
%
% ACTUATION. The design model in trapFullStateModel RETAINS the
% electrostatic stiffness k in A_c, so the controller must not also cancel
% it. The commanded acceleration is therefore realised through a FROZEN
% inversion at the nominal pose (poseEst = 0), which is small-signal linear
% actuation about hover. As the plate moves, the delivered wrench deviates
% exactly as the physics dictates, and that deviation is what A_c models.
% Passing the live pose estimate instead would cancel the stiffness and
% leave the design model describing a plant that no longer exists.
%
% SENSING. The real chain is used throughout: trapCapacitance ->
% trapSensing -> trapCapToHeight, sampled through trapSamplerInit/Step so
% each channel carries its true age.
%
% See also TRAPFULLSTATEMODEL, SIMULATETRAPLEVITATION.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p = trapDeviceGeometry();
M = trapFullStateModel(p);
m = p.mass.m; Ix = p.mass.I(1,1); Iy = p.mass.I(2,2);
Jd = [m; Ix; Iy];

fprintf('T = %.0f us,  channel ages = %s us\n', M.T*1e6, sprintf('%.0f ', M.age*1e6));
fprintf('regulator max|z| = %.4f,  observer max|z| = %.4f\n\n', ...
    max(abs(M.clR)), max(abs(M.clE)));

tEnd = 0.30;  Vmax = 120;
cases = { 'rz = -5 um',               [0 0 -5e-6 0 0 0]; ...
          'rz = -10 um',              [0 0 -10e-6 0 0 0]; ...
          'rz = -5 um + 100 urad',    [0 0 -5e-6 100e-6 -60e-6 0]; ...
          'lateral -20 um, rz = -5',  [-20e-6 -20e-6 -5e-6 0 0 0] };

fprintf('%-28s %16s %12s %12s %10s\n','start','fate','|rz| end','|thx| end','Vpk');
R = struct([]);
for ic = 1:size(cases,1)
    x0 = zeros(12,1); x0(1:6) = cases{ic,2}';
    S = runLQI(p, M, Jd, x0, tEnd, Vmax);
    R(ic).S = S; R(ic).lab = cases{ic,1};
    if isempty(S.contact)
        fprintf('%-28s %16s %12.3e %12.3e %10.1f\n', cases{ic,1}, ...
            sprintf('survives %.0f ms',tEnd*1e3), abs(S.x(end,3))*1e6, ...
            abs(S.x(end,4))*1e6, S.Vpk);
    else
        fprintf('%-28s %16s %12s %12s %10.1f\n', cases{ic,1}, ...
            sprintf('contact %.1f ms',S.contactTime*1e3), '-', '-', S.Vpk);
    end
end

%% ---------------- plot ----------------
col = lines(numel(R));
f = figure('Color','w','Position',[50 50 1250 720]);
lbl = {'r_z [\mum]','\theta_x [\murad]','\theta_y [\murad]'};
idx = [3 4 5];
for k = 1:3
    subplot(2,3,k); hold on; grid on; box on;
    for i = 1:numel(R)
        S = R(i).S;
        plot(S.t*1e3, S.x(:,idx(k))*1e6,'-','LineWidth',1.7,'Color',col(i,:), ...
             'DisplayName',R(i).lab);
    end
    yline(0,'k:','HandleVisibility','off');
    xlabel('t [ms]'); ylabel(lbl{k}); title(lbl{k}(1:regexp(lbl{k},' ','once')-1));
    if k==1, legend('Location','southeast','FontSize',7); end
end
subplot(2,3,4); hold on; grid on; box on;
for i = 1:numel(R)
    S = R(i).S;
    plot(S.t*1e3, S.x(:,3)*1e6,'-','LineWidth',1.7,'Color',col(i,:));
end
yline(0,'k:'); xlabel('t [ms]'); ylabel('r_z [\mum]'); ylim([-0.5 0.5]);
title('r_z, zoomed');
subplot(2,3,5); hold on; grid on; box on;
for i = 1:numel(R)
    S = R(i).S; nk = numel(S.tTick);
    plot(S.tTick*1e3, S.u(:,1),'-','LineWidth',1.5,'Color',col(i,:));
end
yline(0,'k:'); xlabel('t [ms]'); ylabel('a_z [m/s^2]'); title('commanded a_z');
subplot(2,3,6); hold on; grid on; box on;
for i = 1:numel(R)
    S = R(i).S;
    plot(S.tTick*1e3, S.Vlog,'-','LineWidth',1.5,'Color',col(i,:));
end
yline(Vmax,'r--'); xlabel('t [ms]'); ylabel('peak drive [V]'); title('control effort');
sgtitle(sprintf(['Full-state LQI with per-channel fractional delays  \\bullet  ' ...
    'T_{ISR} = %.0f \\mus, ages %.0f/%.0f/%.0f/%.0f \\mus'], ...
    M.T*1e6, M.age(1)*1e6, M.age(2)*1e6, M.age(3)*1e6, M.age(4)*1e6),'FontSize',11);
for ax = findobj(f,'Type','axes')', ax.Toolbar.Visible='off'; end
out = fullfile(root,'presentation','FullStateLQI.png');
exportgraphics(f,out,'Resolution',140);
fprintf('\nSaved %s\n', out);

% =====================================================================
function S = runLQI(p, M, Jd, x0, tEnd, Vmax)
T = M.T; dt = p.sim.dtEuler; nSub = round(T/dt); nT = floor(tEnd/T);
x = x0(:);
xh = zeros(9,1);                      % [xhat(6); zhat(3)]
uPrev = zeros(3,1); xI = zeros(3,1);
s = trapSamplerInit(p, T);
S.t = zeros(nT*nSub+1,1); S.x = zeros(nT*nSub+1,12);
S.tTick = zeros(nT,1); S.u = zeros(nT,3); S.Vlog = zeros(nT,1);
S.t(1)=0; S.x(1,:)=x'; n=1; nk=0; S.contact=[]; S.Vpk=0;
for k = 1:nT
    t = (k-1)*T;
    if ~checkPoseValid(p, x(1:6)), S.contact=n; S.contactTime=t; break; end
    try
        [s, CeD, sInfo] = trapSamplerStep(s, x(1:6), t);
        if sInfo.valid, Ce = CeD; else, Ce = sInfo.CeffNow; end
        y = trapCapToHeight(p, Ce);
    catch
        S.contact=n; S.contactTime=t; break
    end

    % ---- control law ----
    u = -M.Kx*xh(1:6) - M.Kz*xh(7:9) - M.Ku*uPrev - M.KI*xI;
    xI = xI + T*(zeros(3,1) - M.H*y);

    nk = nk+1; S.tTick(nk)=t; S.u(nk,:)=u';

    % ---- frozen-pose inversion: u is an acceleration ----
    V = trapWrenchToVoltage(p, Jd.*u, zeros(6,1), 'Vtia', 0, 'Vmax', Vmax);
    S.Vlog(nk) = max(V(1:5));  S.Vpk = max(S.Vpk, S.Vlog(nk));

    % ---- propagate ----
    hit = false;
    for j = 1:nSub
        if ~checkPoseValid(p, x(1:6)), hit=true; S.contact=n; S.contactTime=t+(j-1)*dt; break; end
        x = x + dt*trapDynamics(p, x, V);
        n=n+1; S.t(n)=t+j*dt; S.x(n,:)=x';
    end
    if hit, break; end

    % ---- observer (predictor form) ----
    xh = M.Abar*xh + M.Bbar*u + M.L*(y - M.Cbar*xh - M.Dbar*uPrev);
    uPrev = u;
end
S.t=S.t(1:n); S.x=S.x(1:n,:);
S.tTick=S.tTick(1:nk); S.u=S.u(1:nk,:); S.Vlog=S.Vlog(1:nk);
end
