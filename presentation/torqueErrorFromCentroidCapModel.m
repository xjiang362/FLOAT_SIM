%% torqueErrorFromCentroidCapModel.m
% QUESTION: the closed-form capacitance
%       C ~= eps0 * A / (g0 - rz - thx*Yc + thy*Xc)
% keeps only the CENTROID (moment-arm) tilt term. How much torque is lost
% by dropping the area-spread (second-moment) part?
%
% Energy method, every conductor at fixed potential:  tau = +dW/dtheta,
% W = 0.5*C*V^2, so
%       taux = 0.5 * V^2 * dC/dthx
%
% Two contributions to dC/dthx:
%   (a) CENTROID / moment-arm term. Analytic from the formula above:
%           dC/dthx = sum_k eps0*A_k*Yc_k / g_k^2
%       This is what survives in the centroid model. It is ZERO for any
%       electrode centred at the origin (Yc=0).
%   (b) AREA-SPREAD term. The exact integral C = eps0*INT INT dA/g(x,y)
%       has extra tilt dependence from the electrode's own extent, which
%       the centroid formula cannot represent. For a CENTRED electrode
%       this is the ONLY source of torque.
%
% Compared against COMSOL MST torque probes for both electrode types.
%
% SYMBOLS:
%   eps0 [F/m]; V [V] drive; g0 [m] nominal gap; a5 [-] rz->gap scale;
%   rz [m]; thx,thy [rad]; A_k [m^2] sub-rectangle area;
%   Xc_k,Yc_k [m] sub-rectangle centroid; C [F]; taux [N*m];
%   kTilt [N*m/rad] = d(taux)/d(thx) at thx=0.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
eps0 = 8.8541878128e-12; V = 30; g0 = 15e-6; a5 = 0.976; d2r = pi/180;
geom = trapGeometry();

%% ================= model evaluators =================
% centroid model: list of [A_k, Xc_k, Yc_k] with signs folded into A_k
cenList = @(rows) cell2mat(arrayfun(@(k) ...
    [rows(k,1)*rectArea(rows(k,:)), (rows(k,4)+rows(k,6))/2, (rows(k,2)+rows(k,3))/2], ...
    (1:size(rows,1))','UniformOutput',false));

C_cen   = @(L,G,tx,ty) sum(eps0*L(:,1)' ./ (G - tx*L(:,2)'.*0 - tx*L(:,3)' + ty*L(:,2)'), 2);
dC_cen  = @(L,G,tx,ty) sum(eps0*L(:,1)'.*L(:,3)' ./ (G - tx*L(:,3)' + ty*L(:,2)').^2, 2);

C_ex = @(rows,G,tx,ty) eps0*sum(cell2mat(arrayfun(@(k) ...
    rows(k,1)*capTrapz(rows(k,2),rows(k,3),rows(k,4),rows(k,5),rows(k,6),rows(k,7),G,tx,ty), ...
    (1:size(rows,1))','UniformOutput',false)'),2);

%% ================= electrode definitions =================
E(1).name = 'quadrant 1 pair';  E(1).rows = [geom.quad(1).a; geom.quad(1).b];
E(2).name = 'central A';        E(2).rows = [geom.centralA.px; geom.centralA.mx];

%% ================= COMSOL reference data (rz=0, thx block) =================
Q = readmatrix(fullfile(scriptDir,'comsol_sweep_combined300_results.txt'),'NumHeaderLines',5);
qsel = abs(Q(:,4))>0 & all(abs(Q(:,[1 2 3 5 6]))==0,2);
ref(1).thx = Q(qsel,4)*d2r;  ref(1).tau = Q(qsel,11);

Cc = readmatrix(fullfile(scriptDir,'comsol_centerElectrodeSweep_rzThxThy150_results.txt'),'NumHeaderLines',5);
csel = abs(Cc(:,4))>0 & all(abs(Cc(:,[1 2 3 5 6]))==0,2);
ref(2).thx = Cc(csel,4)*d2r;  ref(2).tau = Cc(csel,10);

%% ================= compare =================
h = 1e-8;   % [rad] central-difference step for the exact model
fprintf('========== torque decomposition at rz = 0, V = %g V ==========\n\n', V);
for e = 1:2
    rows = E(e).rows;  L = cenList(rows);
    th = ref(e).thx;   G = g0*ones(size(th));

    tau_cen = 0.5*V^2 * dC_cen(L, G, th, 0*th);
    tau_ex  = 0.5*V^2 * (C_ex(rows,G,th+h,0*th) - C_ex(rows,G,th-h,0*th))/(2*h);
    tau_mst = ref(e).tau;

    % stiffness at thx = 0
    kc = localSlope(th, tau_cen);
    ke = localSlope(th, tau_ex);
    km = localSlope(th, tau_mst);

    fprintf('--- %s ---\n', E(e).name);
    fprintf('  sum(A_k*Yc_k) = %.4e m^3   (the moment-arm weight)\n', sum(L(:,1).*L(:,3)));
    fprintf('  kTilt [N*m/rad]  centroid-only %.4e | exact %.4e | COMSOL MST %.4e\n', kc, ke, km);
    fprintf('  centroid-only captures %.1f %% of the exact tilt stiffness\n', 100*kc/ke);
    fprintf('  exact vs COMSOL MST: %.1f %% \n', 100*ke/km);
    fprintf('  max |tau| over sweep: centroid %.4e | exact %.4e | MST %.4e\n\n', ...
        max(abs(tau_cen)), max(abs(tau_ex)), max(abs(tau_mst)));

    E(e).kc = kc; E(e).ke = ke; E(e).km = km;
    E(e).th = th; E(e).tau_cen = tau_cen; E(e).tau_ex = tau_ex; E(e).tau_mst = tau_mst;
end

%% ================= device-level: 4 quadrants + central A =================
fprintf('========== DEVICE TOTAL: 4 quadrants (common mode) + central A ==========\n');
th = E(1).th; G = g0*ones(size(th));
kq_cen = 0; kq_ex = 0;
for q = 1:4
    rows = [geom.quad(q).a; geom.quad(q).b];  L = cenList(rows);
    kq_cen = kq_cen + localSlope(th, 0.5*V^2*dC_cen(L,G,th,0*th));
    kq_ex  = kq_ex  + localSlope(th, 0.5*V^2*(C_ex(rows,G,th+h,0*th)-C_ex(rows,G,th-h,0*th))/(2*h));
end
kcA_cen = E(2).kc; kcA_ex = E(2).ke;

fprintf('  4 quadrants   : centroid-only %.4e | exact %.4e  (%.1f %% captured)\n', ...
    kq_cen, kq_ex, 100*kq_cen/kq_ex);
fprintf('  central A     : centroid-only %.4e | exact %.4e  (%.1f %% captured)\n', ...
    kcA_cen, kcA_ex, 100*kcA_cen/max(kcA_ex,eps));
fprintf('  DEVICE TOTAL  : centroid-only %.4e | exact %.4e\n', ...
    kq_cen+kcA_cen, kq_ex+kcA_ex);
fprintf('  -> dropping the area-spread term UNDERSTATES total tilt stiffness by %.1f %%\n', ...
    100*(1 - (kq_cen+kcA_cen)/(kq_ex+kcA_ex)));
fprintf('  -> of the missing amount, %.0f %% is the central electrode and %.0f %% the quadrants\n', ...
    100*(kcA_ex-kcA_cen)/((kq_ex-kq_cen)+(kcA_ex-kcA_cen)), ...
    100*(kq_ex-kq_cen)/((kq_ex-kq_cen)+(kcA_ex-kcA_cen)));

%% ================= plot =================
figure('Color','w','Position',[40 40 1340 440]);
for e = 1:2
    subplot(1,3,e); hold on; grid on;
    [xs,o] = sort(E(e).th);
    plot(xs/d2r, E(e).tau_mst(o)*1e9,'o','MarkerSize',4,'Color',[0.2 0.3 0.75],'DisplayName','COMSOL MST');
    plot(xs/d2r, E(e).tau_ex(o)*1e9,'k-','LineWidth',1.6,'DisplayName','exact \intdA/g');
    plot(xs/d2r, E(e).tau_cen(o)*1e9,'-','LineWidth',1.8,'Color',[0.85 0.33 0.10],'DisplayName','centroid only');
    xlabel('\theta_x [deg]'); ylabel('Taux [nN\cdotm]'); title(E(e).name);
    legend('Location','best','FontSize',8);
end
subplot(1,3,3); hold on; grid on;
bar([kq_cen kq_ex; kcA_cen kcA_ex]*1e6);
set(gca,'XTick',1:2,'XTickLabel',{'4 quadrants','central A'});
ylabel('kTilt [\muN\cdotm/rad]'); legend({'centroid only','exact'},'Location','best','FontSize',8);
title('tilt stiffness: what the centroid model keeps');
sgtitle('Effect of dropping the area-spread (non-centroid) term on torque');
outPng = fullfile(scriptDir,'TorqueErrorFromCentroidCapModel.png');
exportgraphics(gcf,outPng);
fprintf('\nSaved plot to %s\n', outPng);

%% ================= helpers =================
function a = rectArea(row)
dp = row(6)-row(4); dq = row(7)-row(5); ya=row(2); yb=row(3);
a = dp*(yb-ya) + dq*(yb^2-ya^2)/2;
end
function k = localSlope(th, tau)
m = abs(th) <= 0.35*max(abs(th));
p = [th(m), ones(sum(m),1)] \ tau(m);
k = p(1);
end
