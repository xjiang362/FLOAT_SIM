%% fitCenterElectrodeExact.m
% v2 fixed the cross-run scale bug but still fitted badly (n=2.52,
% R^2=0.63) because it regressed EXTRACTED SLOPES -- which presupposes
% tau is linear in theta. It isn't: the dimensionless tilt u = b*theta/g
% runs 0.3-0.5 in every slice, so the small-angle expansion behind the
% 1/g^3 law is invalid. Here we use the EXACT integral instead.
%
% Geometry: centered electrode, half-widths a (x) and b (y), plate tilted
% by theta about the centroid, local gap g(y) = g - y*theta.
%
%   Fz   = (eps V^2/2) INT dA/(g-y th)^2
%        = eps V^2 a b 2/(g^2 - b^2 th^2)
%        => Fz = k V^2 / (g^2 - b^2 th^2)        [EXACT, one extra param b]
%
%   tau  = (eps V^2/2) INT y dA/(g-y th)^2
%        = eps V^2 a / th^2 * [ 2u/(1-u^2) - 2 atanh(u) ],  u = b th/g
%        => tau = K V^2 / th^2 * [ 2u/(1-u^2) - 2 atanh(u) ]   [EXACT]
%
% Small-u series: [.] = (4/3)u^3 + (8/5)u^5 + ...  -> tau ~ (4/3)K V^2 b^3 th/g^3,
% recovering the 1/g^3 law only in the limit u<<1 (never satisfied here).
%
% This makes b a SHARED, physically meaningful parameter: it is fitted
% from Fz alone (closed form, linear after inverting), then reused in the
% torque model where K is then the single remaining free parameter --
% a genuine cross-validation, not a curve-fit knob.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir,'centerElectrodeModel.mat'),'a5','g0_z');

V = 30; Lplate = 6000e-6;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

raw1 = readmatrix(fullfile(scriptDir,'comsol_centerElectrodeSweep_rzThxThy150_results.txt'),'NumHeaderLines',5);
raw2 = readmatrix(fullfile(scriptDir,'comsol_centerElectrode_rzCrossTerm200_results.txt'),'NumHeaderLines',5);
rz1=raw1(:,3); thx1=raw1(:,4); thy1=raw1(:,5); Fz1=raw1(:,9); Taux1=raw1(:,10); Tauy1=raw1(:,11);
rz2=raw2(:,3); thx2=raw2(:,4); thy2=raw2(:,5); Fz2=raw2(:,9); Taux2=raw2(:,10); Tauy2=raw2(:,11);
scale = 0.2773;                      % measured in fitCenterElectrodeTauxTauy_v2 (spread 1.7%)

gapOf = @(rz) g0_z - a5*rz;

%% ===== STEP 1: exact Fz model -> k and the shared half-width b =====
% 1/Fz = g^2/(k V^2) - b^2 th^2/(k V^2)  : linear in [1/(kV^2), b^2/(kV^2)]
g1 = gapOf(rz1); th1sq = (thx1*pi/180).^2 + (thy1*pi/180).^2;
X = [g1.^2, -th1sq];
c = X \ (1./Fz1);
k_c = 1/(c(1)*V^2);
b_Fz = sqrt(c(2)/c(1));

Fz_exact = @(rz,thsq,k,b) k*V^2 ./ max(gapOf(rz).^2 - b^2*thsq, 1e-18);

fprintf('===== STEP 1: exact Fz = k V^2/(g^2 - b^2 th^2), fitted on run 1 =====\n');
fprintf('k = %.6e,  b = %.2f um   (plate half-width is %.0f um)\n', k_c, b_Fz*1e6, Lplate/2*1e6);
blk = struct('rz',1:50,'thx',51:100,'thy',101:150);
for f = fieldnames(blk)'
    i = blk.(f{1});
    fprintf('  %-4s-only block R^2 = %+.4f\n', f{1}, r2fun(Fz1(i), Fz_exact(rz1(i),th1sq(i),k_c,b_Fz)));
end
% Held-out check on run 2 (rescaled)
g2 = gapOf(rz2); th2sq = (thx2*pi/180).^2 + (thy2*pi/180).^2;
fprintf('  run 2 (held out, rescaled) R^2 = %+.4f\n\n', ...
    r2fun(Fz2/scale, Fz_exact(rz2,th2sq,k_c,b_Fz)));

%% ===== STEP 2: exact torque model, b taken from the Fz fit =====
% tau = K V^2/th^2 * [2u/(1-u^2) - 2 atanh(u)],  u = b*th/g
% Numerically: use the series for tiny u to avoid 0/0.
shapeFun = @(u) (abs(u) < 1e-3) .* ((4/3)*u.^3 + (8/5)*u.^5) + ...
                (abs(u) >= 1e-3) .* (2*u./max(1-u.^2,1e-12) - 2*atanh(min(abs(u),0.999999).*sign(u)));
basisFun = @(rz,thr,b) V^2 .* shapeFun(b*thr./gapOf(rz)) ./ max(thr.^2, 1e-30);

% Assemble every torque point, run 2 rescaled onto run 1's V=30 footing.
rzX  = [rz1(blk.thx); rz2(1:100)];
thrX = [thx1(blk.thx); thx2(1:100)]*pi/180;
tauX = [Taux1(blk.thx); Taux2(1:100)/scale];
rzY  = [rz1(blk.thy); rz2(101:200)];
thrY = [thy1(blk.thy); thy2(101:200)]*pi/180;
tauY = [Tauy1(blk.thy); Tauy2(101:200)/scale];

% Physical validity: the integral diverges as u->1 (plate touching the
% electrode edge). Keep a margin.
uX = b_Fz*abs(thrX)./gapOf(rzX);  okX = uX < 0.85 & abs(thrX) > 1e-9;
uY = b_Fz*abs(thrY)./gapOf(rzY);  okY = uY < 0.85 & abs(thrY) > 1e-9;

KX = basisFun(rzX(okX),thrX(okX),b_Fz) \ tauX(okX);
KY = basisFun(rzY(okY),thrY(okY),b_Fz) \ tauY(okY);

fprintf('===== STEP 2: exact torque, b locked to the Fz-derived %.2f um =====\n', b_Fz*1e6);
fprintf('K_x = %.6e, K_y = %.6e   (x<->y symmetry: %.2f%% apart)\n', KX, KY, 200*abs(KX-KY)/(KX+KY));
fprintf('u range used: %.3f to %.3f  (%d of %d pts kept)\n', min(uX(okX)), max(uX(okX)), sum(okX), numel(okX));

tauX_hat = basisFun(rzX,thrX,b_Fz)*KX;
tauY_hat = basisFun(rzY,thrY,b_Fz)*KY;
fprintf('\n%-30s %10s\n','dataset','R^2');
report = @(nm,m,y,yh) fprintf('%-30s %+10.4f\n', nm, r2fun(y(m), yh(m)));
n1 = numel(blk.thx);
report('run1 thx-only (rz=0)',      [true(n1,1);false(100,1)] & okX, tauX, tauX_hat);
report('run2 thx rz=-7um',  [false(n1,1); rz2(1:100)<0] & okX, tauX, tauX_hat);
report('run2 thx rz=+7um',  [false(n1,1); rz2(1:100)>0] & okX, tauX, tauX_hat);
report('POOLED Taux',                okX, tauX, tauX_hat);
report('POOLED Tauy',                okY, tauY, tauY_hat);

%% ===== STEP 3: also let b float for the torque, as a consistency check =====
bGrid = linspace(0.2*b_Fz, 2.0*b_Fz, 400);
best = -inf; bBest = b_Fz;
for bb = bGrid
    uu = bb*abs(thrX)./gapOf(rzX); m = uu < 0.85 & abs(thrX)>1e-9;
    if sum(m) < 60, continue; end
    B = basisFun(rzX(m),thrX(m),bb); Kt = B \ tauX(m);
    rr = r2fun(tauX(m), B*Kt);
    if rr > best, best = rr; bBest = bb; end
end
fprintf('\n===== STEP 3: torque-only best-fit b =====\n');
fprintf('b from Fz    = %.2f um\n', b_Fz*1e6);
fprintf('b from torque= %.2f um   (R^2 = %+.4f)   -> agreement %.1f%%\n', ...
    bBest*1e6, best, 100*(1-abs(bBest-b_Fz)/b_Fz));

save(fullfile(scriptDir,'centerElectrodeExactModel.mat'), ...
    'k_c','b_Fz','KX','KY','scale','g0_z','a5','V');

%% ===== Plot =====
figure('Color','w','Position',[50 50 1320 540]);
rzPts = [-7e-6 0 7e-6]; cols = [0.80 0.25 0.20; 0.20 0.30 0.75; 0.15 0.60 0.30];
subplot(1,2,1); hold on; grid on;
for i=1:3
    m = abs(rzX-rzPts(i))<1e-12 & okX;
    scatter(thrX(m)*180/pi, tauX(m), 26, cols(i,:), 'filled', 'DisplayName',sprintf('data rz=%+.0f\\mum',rzPts(i)*1e6));
    tg = linspace(min(thrX(m)),max(thrX(m)),200)';
    plot(tg*180/pi, basisFun(rzPts(i)*ones(size(tg)),tg,b_Fz)*KX, '-', ...
        'Color',cols(i,:),'LineWidth',1.8,'HandleVisibility','off');
end
xlabel('\theta_x [deg]'); ylabel('Taux [N\cdotm] (V=30 footing)');
title(sprintf('Exact model, b=%.0f\\mum locked from Fz', b_Fz*1e6));
legend('Location','best','FontSize',8);

subplot(1,2,2); hold on; grid on;
scatter(tauX(okX), tauX_hat(okX), 26, [0.2 0.3 0.75],'filled');
lim = [min(tauX(okX)) max(tauX(okX))];
plot(lim,lim,'k--','LineWidth',1.3);
xlabel('Taux simulated [N\cdotm]'); ylabel('Taux model [N\cdotm]');
title(sprintf('Parity, pooled  (R^2 = %.4f)', r2fun(tauX(okX),tauX_hat(okX))));
axis equal; grid on;
sgtitle('Center electrode: exact closed-form torque model');
exportgraphics(gcf, fullfile(scriptDir,'CenterElectrode_ExactFit.png'));
fprintf('\nSaved plot to CenterElectrode_ExactFit.png\n');
