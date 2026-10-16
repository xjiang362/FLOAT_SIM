%% fitCentralElectrodeA_final.m
% FINAL model for the central electrode captured by the 150-pt file
% ("electrode A" -- the LARGER of the two central electrodes, b~1580um,
% so most likely the OUTER one; confirm against the COMSOL model before
% naming it as such in the report).
%
% Forensics (forensicCenterElectrodeBlocks.m) established:
%   * the 150-pt file is internally consistent -- k agrees to 0.1% across
%     all three blocks, b agrees to 3% between its thx and thy blocks
%   * the 200-pt file is a DIFFERENT electrode (k smaller by 3.567x,
%     consistent with an area ratio at the same 30 V), so it must NOT be
%     combined with this one. It is modeled separately.
%
% Model (exact, derived -- not fitted shapes). For a centered electrode
% of half-widths a (x) and b (y) at gap g = g0 - a5*rz, tilted by theta:
%
%   Fz  = k V^2 / (g^2 - b^2 th^2)
%   tau = K V^2 / th^2 * [ 2u/(1-u^2) - 2 atanh(u) ],   u = b*th/g
%
% obtained by integrating eps V^2/(2 (g - y th)^2) and
% y eps V^2/(2 (g - y th)^2) in closed form over the electrode. Because
% both are EXACT in theta and carry g analytically, the rz=0 tilt blocks
% alone pin down (k, b, K) and the model then extrapolates over rz with
% no further data -- which is why one clean slice suffices here.
%
% Small-u limit: [.] -> (4/3)u^3, i.e. tau -> (4/3) K V^2 b^3 th / g^3,
% the familiar 1/g^3 law. u reaches 0.24 here, so the exact form matters.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
V = 30; g0_z = 15e-6; a5 = 0.976;
gapOf = @(rz) g0_z - a5*rz;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);
shapeFun = @(u) (abs(u)<1e-3).*((4/3)*u.^3 + (8/5)*u.^5) + ...
                (abs(u)>=1e-3).*(2*u./max(1-u.^2,1e-12) - 2*atanh(min(abs(u),0.999999).*sign(u)));

raw = readmatrix(fullfile(scriptDir,'comsol_centerElectrodeSweep_rzThxThy150_results.txt'),'NumHeaderLines',5);
rz=raw(:,3); thx=raw(:,4); thy=raw(:,5); Fz=raw(:,9); Taux=raw(:,10); Tauy=raw(:,11);
blkRz=1:50; blkThx=51:100; blkThy=101:150;

%% ===== Fz: k and b, closed form on all 150 points =====
g = gapOf(rz); thsq = (thx*pi/180).^2 + (thy*pi/180).^2;
c = [g.^2, -thsq] \ (1./Fz);
k_A = 1/(c(1)*V^2);  b_A = sqrt(c(2)/c(1));
Fz_fun = @(rz,thxr,thyr) k_A*V^2 ./ max(gapOf(rz).^2 - b_A^2*(thxr.^2+thyr.^2), 1e-18);
Fz_hat = Fz_fun(rz, thx*pi/180, thy*pi/180);

fprintf('================ CENTRAL ELECTRODE A (150-pt file) ================\n');
fprintf('Fz = k V^2/(g^2 - b^2 th^2)\n');
fprintf('  k = %.6e   b = %.1f um\n', k_A, b_A*1e6);
fprintf('  R^2  rz-only=%.4f  thx-only=%.4f  thy-only=%.4f  ALL=%.4f\n\n', ...
    r2fun(Fz(blkRz),Fz_hat(blkRz)), r2fun(Fz(blkThx),Fz_hat(blkThx)), ...
    r2fun(Fz(blkThy),Fz_hat(blkThy)), r2fun(Fz,Fz_hat));

%% ===== Torque: single coefficient K, b locked from Fz =====
basis = @(rz,thr) V^2 .* shapeFun(b_A*thr./gapOf(rz)) ./ max(thr.^2, 1e-30);

ix = blkThx(abs(thx(blkThx)) > 1e-9);
iy = blkThy(abs(thy(blkThy)) > 1e-9);
Bx = basis(rz(ix), thx(ix)*pi/180);  KX = Bx \ Taux(ix);
By = basis(rz(iy), thy(iy)*pi/180);  KY = By \ Tauy(iy);

fprintf('tau = K V^2/th^2 * [2u/(1-u^2) - 2 atanh(u)],  u = b th/g\n');
fprintf('  K_x = %.6e   R^2 = %+.4f   (%d pts)\n', KX, r2fun(Taux(ix), Bx*KX), numel(ix));
fprintf('  K_y = %.6e   R^2 = %+.4f   (%d pts)\n', KY, r2fun(Tauy(iy), By*KY), numel(iy));
fprintf('  x<->y symmetry: %.2f%% apart  -> shared K = %.6e\n', ...
    200*abs(KX-KY)/(KX+KY), mean([KX KY]));
K_A = mean([KX KY]);
fprintf('  R^2 with the SHARED K:  Taux=%+.4f  Tauy=%+.4f\n', ...
    r2fun(Taux(ix), Bx*K_A), r2fun(Tauy(iy), By*K_A));

% How much does the exact form beat the pure 1/g^3 (small-u) limit?
Bx_lin = V^2*(4/3)*(b_A*(thx(ix)*pi/180)./gapOf(rz(ix))).^3 ./ max((thx(ix)*pi/180).^2,1e-30);
Kl = Bx_lin \ Taux(ix);
fprintf('\n  Small-angle 1/g^3 approximation instead: R^2 = %+.4f (exact form gains %+.4f)\n', ...
    r2fun(Taux(ix), Bx_lin*Kl), r2fun(Taux(ix),Bx*KX) - r2fun(Taux(ix), Bx_lin*Kl));
fprintf('  max u reached = %.3f  (1/g^3 needs u<<1)\n', max(abs(b_A*(thx(ix)*pi/180)./gapOf(rz(ix)))));

save(fullfile(scriptDir,'centralElectrodeA_model.mat'),'k_A','b_A','K_A','g0_z','a5','V');
fprintf('\nSaved coefficients to centralElectrodeA_model.mat\n');

%% ===== Plot =====
figure('Color','w','Position',[50 50 1340 460]);

subplot(1,3,1); hold on; grid on;
[rzs,o] = sort(rz(blkRz));
plot(rzs*1e6, Fz(blkRz(o)), 'o', 'MarkerSize',5, 'Color',[0.2 0.3 0.75], 'DisplayName','COMSOL');
plot(rzs*1e6, Fz_fun(rzs,0*rzs,0*rzs), 'k-','LineWidth',1.7,'DisplayName','model');
set(gca,'YScale','log'); xlabel('rz [\mum]'); ylabel('Fz [N]');
title(sprintf('Fz vs rz   (R^2=%.4f)', r2fun(Fz(blkRz),Fz_hat(blkRz))));
legend('Location','best','FontSize',8);

subplot(1,3,2); hold on; grid on;
scatter(thx(blkThx), Taux(blkThx), 26, [0.80 0.25 0.20],'filled','DisplayName','COMSOL Taux');
tg = linspace(min(thx(blkThx)),max(thx(blkThx)),200)';
plot(tg, basis(0*tg, tg*pi/180)*K_A, 'k-','LineWidth',1.8,'DisplayName','exact model');
xlabel('\theta_x [deg]'); ylabel('Taux [N\cdotm]');
title(sprintf('Taux at rz=0   (R^2=%.4f)', r2fun(Taux(ix),Bx*K_A)));
legend('Location','best','FontSize',8);

subplot(1,3,3); hold on; grid on;
% Model's rz extrapolation for torque -- no data yet, this is the prediction
rzg = linspace(-14e-6,14e-6,200)';
thFix = 0.05*pi/180;
plot(rzg*1e6, basis(rzg, thFix*ones(size(rzg)))*K_A, 'k-','LineWidth',1.8, ...
    'DisplayName','model prediction');
plot(0, basis(0,thFix)*K_A, 'o','MarkerSize',10,'MarkerFaceColor',[0.80 0.25 0.20], ...
    'Color',[0.80 0.25 0.20],'DisplayName','only rz with data');
set(gca,'YScale','log'); xlabel('rz [\mum]'); ylabel(sprintf('Taux at \\theta_x=%.3f deg [N\\cdotm]',thFix*180/pi));
title({'Torque vs rz: analytic extrapolation','(UNVALIDATED - needs one rz sweep to confirm)'});
legend('Location','best','FontSize',8);

sgtitle('Central electrode A (150-pt file, b\approx1580\mum) - exact closed-form model');
exportgraphics(gcf, fullfile(scriptDir,'CentralElectrodeA_FinalFit.png'));
fprintf('Saved plot to CentralElectrodeA_FinalFit.png\n');
