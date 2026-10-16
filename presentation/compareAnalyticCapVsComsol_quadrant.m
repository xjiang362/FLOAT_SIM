%% compareAnalyticCapVsComsol_quadrant.m
% Tests the TILT terms of the leading-order analytic capacitance
%
%       C ~= eps0 * (dX*dY) / (g0 - rz - thx*Yc -/+ thy*Xc)
%
% on a QUADRANT actuation electrode, whose centroid is at (1500,1500) um
% so the tilt terms are FIRST ORDER (unlike central electrode A, where
% Xc=Yc=0 makes them vanish identically).
%
% Data: comsol_overnightFull6DSweep_results.txt, 2000-pt full-6D LHS,
% which carries a verified Total-Electric-Energy probe (col 16) plus the
% two terminal-charge probes of the differential pair (cols 10 and 14).
%
% The quadrant is a DIFFERENTIAL PAIR of 1450 um squares at (2250,750)
% and (750,2250), driven at +V and -V. Because the mutual capacitance
% between them is negligible (verified here: 2W/V^2 == (|Qa|+|Qb|)/V),
% the pair's effective capacitance is simply Ca + Cb, and the analytic
% formula must be applied PER SUB-ELECTRODE and summed -- each has its
% own centroid, and the pair centroid (1500,1500) is only a first-order
% stand-in for the two.
%
% SYMBOLS (SI unless suffixed):
%   eps0     [F/m]  vacuum permittivity
%   V        [V]    magnitude of the differential drive on each sub-electrode
%   g0       [m]    gap at nominal pose (rz=0)
%   a5       [-]    calibrated scale of rz onto gap change
%   rz       [m]    plate vertical displacement (+rz closes the gap)
%   thx,thy  [rad]  plate tilts
%   Xc,Yc    [m]    sub-electrode centroid, plate-centred coords
%   A_sub    [m^2]  one sub-electrode's area (1450 um square)
%   W        [J]    COMSOL total electric energy
%   Qa,Qb    [C]    terminal charges of the two sub-electrodes
%   C_*      [F]    capacitance (reported in pF)

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

eps0 = 8.8541878128e-12;
V    = 30;
g0   = 15e-6;
a5   = 0.976;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

%% ---------------- data ----------------
raw = readmatrix(fullfile(scriptDir,'comsol_overnightFull6DSweep_results.txt'),'NumHeaderLines',5);
rx = raw(:,1); ry = raw(:,2); rz = raw(:,3);
thx = raw(:,4)*pi/180; thy = raw(:,5)*pi/180; thz = raw(:,6)*pi/180;
W  = raw(:,16); Qa = raw(:,10); Qb = raw(:,14);

C_comsol = 2*W/V^2;
C_charge = (abs(Qa)+abs(Qb))/V;

% keep poses where the plate has not been driven into contact
G = g0 - a5*rz;
ok = G > 1.0e-6 & isfinite(C_comsol) & C_comsol > 0;
fprintf('Loaded %d rows; %d usable (gap > 1 um).\n', numel(rz), sum(ok));
fprintf('Energy vs charge consistency: max rel diff = %.3e  -> mutual C negligible,\n', ...
    max(abs(C_comsol(ok)-C_charge(ok))./C_charge(ok)));
fprintf('  so C_eff = 2W/V^2 = Ca + Cb is the right quantity to model.\n');
fprintf('gap range %.2f to %.2f um;  C range %.2f to %.2f pF\n', ...
    min(G(ok))*1e6, max(G(ok))*1e6, min(C_comsol(ok))*1e12, max(C_comsol(ok))*1e12);

%% ------- scan all 4 quadrants x both thy sign conventions -------
geom = trapGeometry();
A_sub = geom.actSide^2;

fprintf('\n============ identifying quadrant AND thy sign ============\n');
fprintf('model: C = sum_sub eps0*A_sub/(g0 - a5*rz - thx*Yc + s*thy*Xc), gain+offset fitted\n');
fprintf('%5s %6s | %9s %10s %10s\n','quad','s(thy)','R^2','gain','offset[pF]');
best = struct('r2',-Inf);
for q = 1:4
    cen = [centroidOf(geom.quad(q).a); centroidOf(geom.quad(q).b)];
    for s = [-1 +1]
        Cm = zeros(numel(rz),1);
        for k = 1:2
            gk = G - thx*cen(k,2) + s*thy*cen(k,1);
            Cm = Cm + eps0*A_sub./max(gk,1e-9);
        end
        M = [Cm(ok), ones(sum(ok),1)];  p = M\C_comsol(ok);
        r2 = r2fun(C_comsol(ok), M*p);
        fprintf('%5d %6+d | %+9.5f %10.4f %10.4f\n', q, s, r2, p(1), p(2)*1e12);
        if r2 > best.r2
            best = struct('r2',r2,'q',q,'s',s,'gain',p(1),'off',p(2),'C',Cm,'cen',cen);
        end
    end
end
fprintf('\nBEST: quadrant %d, thy sign %+d  (R^2 = %.5f)\n', best.q, best.s, best.r2);

%% ------- how much do the tilt terms actually matter? -------
Cm_notilt = 2*eps0*A_sub./G;      % drop BOTH tilt terms
M = [Cm_notilt(ok), ones(sum(ok),1)];  p0 = M\C_comsol(ok);
r2_notilt = r2fun(C_comsol(ok), M*p0);

cen = best.cen; s = best.s;
% pair-centroid version: treat the pair as ONE block at (1500,1500)
Xcp = mean(cen(:,1)); Ycp = mean(cen(:,2));
Cm_pair = 2*eps0*A_sub./max(G - thx*Ycp + s*thy*Xcp, 1e-9);
M = [Cm_pair(ok), ones(sum(ok),1)];  pp = M\C_comsol(ok);
r2_pair = r2fun(C_comsol(ok), M*pp);

% exact integral over the two squares
rowsQ = [geom.quad(best.q).a; geom.quad(best.q).b];
I = zeros(numel(rz),1);
for k = 1:size(rowsQ,1)
    r = rowsQ(k,:);
    I = I + r(1)*capTrapz(r(2),r(3),r(4),r(5),r(6),r(7), G, thx, s*thy);
end
C_exact = eps0*I;
M = [C_exact(ok), ones(sum(ok),1)];  pe = M\C_comsol(ok);
r2_exact = r2fun(C_comsol(ok), M*pe);

fprintf('\n================ model comparison (gain+offset fitted) ================\n');
fprintf('%-34s %10s %9s %11s\n','model','R^2','gain','offset[pF]');
fprintf('%-34s %+10.5f %9.4f %11.4f\n','no tilt terms at all',    r2_notilt, p0(1), p0(2)*1e12);
fprintf('%-34s %+10.5f %9.4f %11.4f\n','pair centroid (1500,1500)',r2_pair,  pp(1), pp(2)*1e12);
fprintf('%-34s %+10.5f %9.4f %11.4f\n','per-sub-electrode centroid',best.r2, best.gain, best.off*1e12);
fprintf('%-34s %+10.5f %9.4f %11.4f\n','exact INT dA/g',           r2_exact, pe(1), pe(2)*1e12);
fprintf('\nimplied effective area from gain (per sub-electrode model): %.4f mm^2 vs geometric %.4f\n', ...
    best.gain*A_sub*1e6, A_sub*1e6);

%% ------- does C depend on rx, ry, thz as the formula says it should not? -------
resid = C_comsol(ok) - [best.C(ok) ones(sum(ok),1)]*[best.gain; best.off];
fprintf('\n--- formula predicts NO rx, ry, thz dependence; test on residual ---\n');
nm = {'rx','ry','thz'}; vv = {rx(ok), ry(ok), thz(ok)};
for i=1:3
    cc = corr(vv{i}, resid);
    fprintf('  corr(resid, %4s) = %+.4f\n', nm{i}, cc);
end
fprintf('  residual RMS = %.4f pF  (%.2f %% of mean C)\n', ...
    rms(resid)*1e12, 100*rms(resid)/mean(C_comsol(ok)));

%% ---------------- plots ----------------
figure('Color','w','Position',[40 40 1400 430]);

subplot(1,3,1); hold on; grid on;
Cfit = [best.C(ok) ones(sum(ok),1)]*[best.gain; best.off];
plot(Cfit*1e12, C_comsol(ok)*1e12, '.','MarkerSize',6,'Color',[0.2 0.3 0.75]);
lims = [0 max(C_comsol(ok))*1.05e12]; plot(lims,lims,'k--','LineWidth',1.2);
xlabel('analytic C [pF]'); ylabel('COMSOL C = 2W/V^2 [pF]');
set(gca,'XScale','log','YScale','log'); axis tight;
title(sprintf('per-sub-electrode model  (R^2=%.4f)', best.r2));

subplot(1,3,2); hold on; grid on;
% isolate the tilt effect: residual of the NO-TILT model vs the tilt term
tiltTerm = -(thx(ok)*Ycp - s*thy(ok)*Xcp)./G(ok);      % fractional gap change
resid0 = C_comsol(ok) - [Cm_notilt(ok) ones(sum(ok),1)]*p0;
scatter(tiltTerm, resid0*1e12, 7, G(ok)*1e6, 'filled');
cb=colorbar; cb.Label.String='gap [\mum]';
xlabel('tilt term  -(\theta_xY_c - \theta_yX_c)/g'); ylabel('no-tilt model residual [pF]');
title('what the tilt terms have to explain');

subplot(1,3,3); hold on; grid on;
plot(G(ok)*1e6, resid*1e12, '.','MarkerSize',6,'Color',[0.85 0.33 0.10]);
yline(0,'k-'); xlabel('gap g_0-a_5r_z [\mum]'); ylabel('residual [pF]');
title(sprintf('residual of best model (RMS %.3f pF)', rms(resid)*1e12));

sgtitle(sprintf('Quadrant %d actuation pair: analytic capacitance vs COMSOL (V=%g V, 2000-pt full-6D)', best.q, V));
outPng = fullfile(scriptDir,'AnalyticCapVsComsol_quadrant.png');
exportgraphics(gcf, outPng);
fprintf('\nSaved plot to %s\n', outPng);

%% ---------------- helper ----------------
function c = centroidOf(row)
% centroid of a single rectangle row [sign ya yb pL qL pR qR]
c = [ (row(4)+row(6))/2 , (row(2)+row(3))/2 ];
end
