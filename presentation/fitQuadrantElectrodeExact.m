%% fitQuadrantElectrodeExact.m
% Applies the derived closed-form method (which fixed the central
% electrode) to the QUADRANT electrode, replacing the empirical
%     Fz   = k V^2 / G_eff^2,  G_eff = g0 - a5 rz - a1 thx - a3 thy
%     Taux = Leqx*Fz,  Leqx = a0 + a1 rz + a2 thx + a3 thy   (4 fitted coeffs)
%
% A quadrant electrode is OFFSET from the rotation axis, so its integral
% differs from the centered case. For a rectangle spanning y in [y1,y2],
% x in [x1,x2], gap g(y) = g - y*th:
%
%   Fz  = (eps V^2/2) dx INT dy/(g-y th)^2
%       = k V^2 / [ (g - y1 th)(g - y2 th) ]
%         -> the product of the NEAR-edge and FAR-edge gaps, not one
%            lumped effective gap squared.
%
%   tau = (eps V^2/2) dx INT y dy/(g-y th)^2
%       = k V^2 / (th*g) * f(u1,u2),   u_i = y_i th/g
%       f(u1,u2) = 1/((1-u1)(1-u2)) + ln((1-u2)/(1-u1))/(u2-u1)
%
% CRITICAL PROPERTY: f depends only on y1,y2 and k -- all three are
% already pinned down by the Fz fit. So the torque model has ZERO free
% parameters. It is a prediction, not a fit, and any agreement with the
% MST torque is genuine validation.
%
% Small-angle limit: f -> (u1+u2)/2, giving tau -> Fz*(y1+y2)/2, i.e. the
% lever arm is the electrode CENTROID -- recovering the a0 of the old
% moment-arm model and explaining what a0 physically was.
%
% Bonus: y1,y2,x1,x2 come out in metres and can be checked directly
% against the CAD electrode footprint.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
V = 30; g0_z = 15e-6;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

data = parseComsolResults(fullfile(scriptDir,'comsol_sweep_combined300_results.txt'));
b = data([101:150, 151:200, 201:250], :);
rz=b.rz; thx=b.thx*pi/180; thy=b.thy*pi/180;
Fz=b.Fz; Taux=b.Taux; Tauy=b.Tauy;
blk = {1:50, 51:100, 101:150}; blkName = {'rz-only','thx-only','thy-only'};

%% ============ EXACT Fz: 1/Fz linear in [g^2, -g*thx, thx^2, -g*thy, thy^2] ============
% The gap scale a5 in g = g0 - a5*rz is fitted, exactly as the empirical
% model does. Physically a5 should be 1 (rz IS the gap change), but the
% empirical model fits it, so leaving it at 1 would be an unfair
% comparison. Only the rz block constrains it (the tilt blocks are all at
% rz=0), so a cheap 1-D grid search is enough; everything else stays
% closed-form linear given a5.
a5Grid = linspace(0.85, 1.25, 401);
bestR2 = -inf; a5_q = 1; c = []; X = []; g = [];
for a5t = a5Grid
    gt = g0_z - a5t*rz;
    Xt = [gt.^2, -gt.*thx, thx.^2, -gt.*thy, thy.^2];
    ct = Xt \ (1./Fz);
    rr = r2fun(Fz, 1./(Xt*ct));
    if rr > bestR2, bestR2 = rr; a5_q = a5t; c = ct; X = Xt; g = gt; end
end
k_q  = 1/(c(1)*V^2);
Sy = c(2)/c(1); Py = c(3)/c(1);      % y1+y2, y1*y2
Sx = c(4)/c(1); Px = c(5)/c(1);
yr = roots([1 -Sy Py]);  xr = roots([1 -Sx Px]);

% Model prediction: 1/Fz is exactly the fitted linear combination.
% (Single-axis blocks: only one tilt is nonzero per row, so this reduces
%  to the correct near-edge x far-edge product for whichever axis is live.)
Fz_hat = 1./(X*c);

fprintf('=============== QUADRANT ELECTRODE: EXACT vs EMPIRICAL ===============\n\n');
fprintf('--- EXACT Fz = k V^2/((g-y1 th)(g-y2 th)) ---\n');
fprintf('  k  = %.6e,  a5 (gap scale) = %.4f  [physically expected 1.0]\n', k_q, a5_q);
fprintf('  RECOVERED GEOMETRY (check these against the CAD footprint!):\n');
fprintf('    y edges: %8.1f um  to %8.1f um   (span %7.1f um, centroid %7.1f um)\n', ...
    min(yr)*1e6, max(yr)*1e6, abs(diff(yr))*1e6, mean(yr)*1e6);
fprintf('    x edges: %8.1f um  to %8.1f um   (span %7.1f um, centroid %7.1f um)\n', ...
    min(xr)*1e6, max(xr)*1e6, abs(diff(xr))*1e6, mean(xr)*1e6);
if isreal(yr) && isreal(xr)
    fprintf('    -> both real: consistent with a physical rectangular footprint\n');
else
    fprintf('    -> COMPLEX roots: the rectangle model is not capturing this geometry\n');
end

%% ============ EMPIRICAL Fz (current model) for comparison ============
S1 = load(fullfile(scriptDir,'simpleFzModel.mat'),'pfit');
p = S1.pfit;   % [k, a1(thx), a3(thy), a5(rz)]
Geff = max(g0_z - p(4)*rz - p(2)*thx - p(3)*thy, 1e-9);
Fz_emp = p(1)*V^2 ./ Geff.^2;

fprintf('\n--- Fz fit quality, per block ---\n');
fprintf('%-12s %14s %14s\n','block','EXACT R^2','EMPIRICAL R^2');
for i=1:3
    j = blk{i};
    fprintf('%-12s %14.5f %14.5f\n', blkName{i}, r2fun(Fz(j),Fz_hat(j)), r2fun(Fz(j),Fz_emp(j)));
end
fprintf('%-12s %14.5f %14.5f\n','ALL', r2fun(Fz,Fz_hat), r2fun(Fz,Fz_emp));

%% ============ TORQUE: ZERO free parameters ============
% Taux from the y-extent; Tauy from the x-extent (with the established
% sign convention Tauy = -Leqy*Fz).
% NOTE on the Tauy sign: the x-side fit absorbs the opposite tilt
% convention (a positive thy raises the gap at positive x), so it returns
% NEGATIVE x roots. The torque tau_y = -INT x P dA then already comes out
% negative through those roots -- applying a further explicit minus sign
% double-counts it. Both axes therefore use the SAME +k V^2 f/(g th)
% form, and the small-angle limit tau -> Fz*centroid reproduces the
% established sign convention (Taux>0, Tauy<0) automatically.
Taux_pred = k_q*V^2 .* ffun(yr(1)*thx./g, yr(2)*thx./g) ./ (g .* thx);
Tauy_pred = k_q*V^2 .* ffun(xr(1)*thy./g, xr(2)*thy./g) ./ (g .* thy);
% theta -> 0 rows: use the centroid limit tau = Fz*(edge1+edge2)/2
z = abs(thx) < 1e-12; Taux_pred(z) = Fz_hat(z)*mean(yr);
z = abs(thy) < 1e-12; Tauy_pred(z) = Fz_hat(z)*mean(xr);

S2 = load(fullfile(scriptDir,'tauxTauyMomentArm.mat'),'a0','a1','a2','a3');
Leqx = S2.a0 + S2.a1*rz + S2.a2*thx + S2.a3*thy;
Leqy = S2.a0 + S2.a1*rz + S2.a2*thy + S2.a3*thx;
Taux_emp = Leqx .* Fz;  Tauy_emp = -Leqy .* Fz;

fprintf('\n--- Torque: EXACT (0 free params) vs EMPIRICAL (4 fitted coeffs) ---\n');
fprintf('%-12s %16s %16s\n','block','EXACT R^2','EMPIRICAL R^2');
for i=1:3
    j = blk{i};
    fprintf('%-12s %16.5f %16.5f    (Taux)\n', blkName{i}, r2fun(Taux(j),Taux_pred(j)), r2fun(Taux(j),Taux_emp(j)));
end
fprintf('%-12s %16.5f %16.5f    (Taux, ALL)\n','ALL', r2fun(Taux,Taux_pred), r2fun(Taux,Taux_emp));
for i=1:3
    j = blk{i};
    fprintf('%-12s %16.5f %16.5f    (Tauy)\n', blkName{i}, r2fun(Tauy(j),Tauy_pred(j)), r2fun(Tauy(j),Tauy_emp(j)));
end
fprintf('%-12s %16.5f %16.5f    (Tauy, ALL)\n','ALL', r2fun(Tauy,Tauy_pred), r2fun(Tauy,Tauy_emp));

fprintf('\nImplied small-angle lever arm (electrode centroid): y=%.1f um, x=%.1f um\n', ...
    mean(yr)*1e6, mean(xr)*1e6);
fprintf('Old empirical a0 (fitted lever arm)              : %.1f um\n', S2.a0*1e6);

save(fullfile(scriptDir,'quadrantExactModel.mat'),'k_q','yr','xr','g0_z','V');

%% ============ Plot ============
figure('Color','w','Position',[40 40 1400 760]);
xv = {rz(blk{1})*1e6, thx(blk{2})*180/pi, thy(blk{3})*180/pi};
xl = {'rz [\mum]','\theta_x [deg]','\theta_y [deg]'};
for i=1:3
    j = blk{i};
    subplot(2,3,i); hold on; grid on;
    plot(xv{i}, Fz(j), 'ko','MarkerSize',5,'DisplayName','MST');
    plot(xv{i}, Fz_hat(j), 'r-','LineWidth',1.6,'DisplayName',sprintf('exact (R^2=%.4f)',r2fun(Fz(j),Fz_hat(j))));
    plot(xv{i}, Fz_emp(j), 'b--','LineWidth',1.3,'DisplayName',sprintf('empirical (R^2=%.4f)',r2fun(Fz(j),Fz_emp(j))));
    xlabel(xl{i}); ylabel('Fz [N]'); title(['Fz, ' blkName{i}]);
    if i==1, set(gca,'YScale','log'); end
    legend('Location','best','FontSize',7);
end
tauD = {Taux, Taux, Tauy}; tauE = {Taux_pred,Taux_pred,Tauy_pred}; tauM = {Taux_emp,Taux_emp,Tauy_emp};
tl = {'Taux, rz-only','Taux, thx-only','Tauy, thy-only'};
for i=1:3
    j = blk{i};
    subplot(2,3,3+i); hold on; grid on;
    plot(xv{i}, tauD{i}(j), 'ko','MarkerSize',5,'DisplayName','MST');
    plot(xv{i}, tauE{i}(j), 'r-','LineWidth',1.6,'DisplayName',sprintf('exact, 0 params (R^2=%.4f)',r2fun(tauD{i}(j),tauE{i}(j))));
    plot(xv{i}, tauM{i}(j), 'b--','LineWidth',1.3,'DisplayName',sprintf('empirical (R^2=%.4f)',r2fun(tauD{i}(j),tauM{i}(j))));
    xlabel(xl{i}); ylabel('Torque [N\cdotm]'); title(tl{i});
    legend('Location','best','FontSize',7);
end
sgtitle('Quadrant electrode: derived closed-form vs current empirical model');
exportgraphics(gcf, fullfile(scriptDir,'QuadrantElectrode_ExactVsEmpirical.png'));
fprintf('\nSaved plot to QuadrantElectrode_ExactVsEmpirical.png\n');

%% ============ local functions ============
function out = ffun(u1,u2)
% f(u1,u2) = 1/((1-u1)(1-u2)) + ln((1-u2)/(1-u1))/(u2-u1)
% Series branch avoids the catastrophic cancellation between the two
% terms as u->0 (they each tend to +-1 while their sum tends to 0).
    small = max(abs(u1),abs(u2)) < 1e-4;
    out = zeros(size(u1));
    out(small) = (u1(small)+u2(small))/2 + ...
                 (2/3)*(u1(small).^2 + u1(small).*u2(small) + u2(small).^2);
    ns = ~small;
    out(ns) = 1./((1-u1(ns)).*(1-u2(ns))) + ...
              (log(1-u2(ns)) - log(1-u1(ns)))./(u2(ns)-u1(ns));
end
