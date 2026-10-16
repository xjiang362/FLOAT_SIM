%% testRun1ExtraElectrodes.m
% Fz is consistent between the two runs after a single uniform factor
% 0.2773, but TORQUE is not -- run 2's two rz slices agree with the exact
% centered-electrode model (R^2 0.91/0.89) while run 1's rz=0 slice sits
% ~6x too high (R^2 0.31). A pure voltage difference cannot do that:
% Fz and tau both scale as V^2, so one factor would fix both.
%
% Hypothesis: run 1 had MORE THAN THE CENTER ELECTRODE energized (center
% + the 4 quadrants), run 2 only the center. That predicts exactly this
% asymmetry: the quadrant electrodes sit at the same gap, so their Fz
% just inflates the total by a near-constant factor (a uniform Fz scale,
% as observed), but they are laterally OFFSET, so their torque has a big
% lever arm and a completely different law -- which no uniform scale can
% absorb.
%
% Decisive test: the 4 quadrant models are already fitted and validated.
% Evaluate them at run 1's own poses and see whether they account for
% run 1's Fz excess AND its torque excess simultaneously.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir,'centerElectrodeExactModel.mat'),'k_c','b_Fz','KX','scale','g0_z','a5','V');

r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);
gapOf = @(rz) g0_z - a5*rz;
shapeFun = @(u) (abs(u)<1e-3).*((4/3)*u.^3+(8/5)*u.^5) + ...
                (abs(u)>=1e-3).*(2*u./max(1-u.^2,1e-12) - 2*atanh(min(abs(u),0.999999).*sign(u)));
tauCenter = @(rz,thr,K,b) K*V^2 .* shapeFun(b*thr./gapOf(rz)) ./ max(thr.^2,1e-30);
FzCenter  = @(rz,thsq,k,b) k*V^2 ./ max(gapOf(rz).^2 - b^2*thsq, 1e-18);

raw1 = readmatrix(fullfile(scriptDir,'comsol_centerElectrodeSweep_rzThxThy150_results.txt'),'NumHeaderLines',5);
rx1=raw1(:,1); ry1=raw1(:,2); rz1=raw1(:,3); thx1=raw1(:,4); thy1=raw1(:,5); thz1=raw1(:,6);
Fz1=raw1(:,9); Taux1=raw1(:,10);
blkThx = 51:100;

%% ---- Sum the 4 quadrant electrode models at run 1's poses (V=30) ----
Fz_q = zeros(150,1); Taux_q = zeros(150,1);
for e = 1:4
    fn = str2func(sprintf('trapElectrode%dModel_v2', e));
    [~,~,fz_e,tx_e,~,~] = fn(rx1,ry1,rz1,thx1,thy1,thz1,V);
    Fz_q = Fz_q + fz_e; Taux_q = Taux_q + tx_e;
end

i = blkThx; thr = thx1(i)*pi/180;

%% ---- Test A: does the quadrant sum explain run 1's Fz excess? ----
% Center-only Fz on run-1 footing would be, from the run-2-validated
% model.  If run 1 = center + quads, then Fz1 ~ Fz_center + Fz_quads.
fprintf('===== TEST A: Fz budget at run 1''s thx-only block (rz=0) =====\n');
Fz_cen = FzCenter(rz1(i), thr.^2, k_c, b_Fz);
fprintf('mean Fz run1 (simulated)      = %.5e N\n', mean(Fz1(i)));
fprintf('mean Fz 4 quadrant models     = %.5e N\n', mean(Fz_q(i)));
fprintf('mean Fz center model (run1 k) = %.5e N\n', mean(Fz_cen));
fprintf('  NOTE: k was fitted ON run 1, so Fz_cen already absorbs any extra electrodes.\n');
fprintf('  Center-only share implied by the measured scale factor = %.4f\n', scale);
fprintf('  Quadrant share implied  = %.4f   vs from quadrant models = %.4f\n', ...
    1-scale, mean(Fz_q(i))/mean(Fz1(i)));

%% ---- Test B: does the quadrant sum explain run 1's TORQUE excess? ----
fprintf('\n===== TEST B: torque budget, same rows =====\n');
% Center-only torque predicted from the run-2-fitted model, put on a
% 30V/center-only footing: KX was fitted on run-2-rescaled data, i.e.
% already on run 1's V=30 footing but representing the CENTER only if
% run 2 was center-only.
Taux_cen = tauCenter(rz1(i), thr, KX, b_Fz);
slope = @(th,t) th \ t;
fprintf('measured   d(Taux)/dthx run1        = %.5e N*m/rad\n', slope(thr, Taux1(i)));
fprintf('center model alone                  = %.5e\n', slope(thr, Taux_cen));
fprintf('4 quadrant models alone             = %.5e\n', slope(thr, Taux_q(i)));
fprintf('center + quadrants                  = %.5e\n', slope(thr, Taux_cen+Taux_q(i)));
fprintf('\nR^2 of run1 Taux vs:\n');
fprintf('  center only          = %+.4f\n', r2fun(Taux1(i), Taux_cen));
fprintf('  quadrants only       = %+.4f\n', r2fun(Taux1(i), Taux_q(i)));
fprintf('  center + quadrants   = %+.4f\n', r2fun(Taux1(i), Taux_cen+Taux_q(i)));

%% ---- Test C: refit the center model on RUN 2 ONLY (self-consistent) ----
fprintf('\n===== TEST C: center model fitted on run 2 ALONE =====\n');
raw2 = readmatrix(fullfile(scriptDir,'comsol_centerElectrode_rzCrossTerm200_results.txt'),'NumHeaderLines',5);
rz2=raw2(:,3); thx2=raw2(:,4); thy2=raw2(:,5); Fz2=raw2(:,9); Taux2=raw2(:,10); Tauy2=raw2(:,11);

% Fz on run 2 alone -> its own k2, b2 (no cross-run assumption at all)
g2 = gapOf(rz2); th2sq=(thx2*pi/180).^2+(thy2*pi/180).^2;
c2 = [g2.^2, -th2sq] \ (1./Fz2);
k2 = 1/(c2(1)*V^2); b2 = sqrt(c2(2)/c2(1));
fprintf('run2-only Fz fit:  k2 = %.6e, b2 = %.2f um   (R^2 = %+.4f)\n', ...
    k2, b2*1e6, r2fun(Fz2, k2*V^2./max(g2.^2-b2^2*th2sq,1e-18)));
fprintf('  (run1-only gave b = %.2f um -- geometry should MATCH if same electrode)\n', b_Fz*1e6);

thrX2 = thx2(1:100)*pi/180; rzX2 = rz2(1:100);
u2 = b2*abs(thrX2)./gapOf(rzX2); ok2 = u2 < 0.85 & abs(thrX2) > 1e-9;
B2 = V^2*shapeFun(b2*thrX2(ok2)./gapOf(rzX2(ok2)))./max(thrX2(ok2).^2,1e-30);
K2 = B2 \ Taux2(ok2);
fprintf('run2-only torque fit: K2 = %.6e,  R^2 = %+.4f  (%d pts)\n', ...
    K2, r2fun(Taux2(ok2), B2*K2), sum(ok2));

thrY2 = thy2(101:200)*pi/180; rzY2 = rz2(101:200);
uy2 = b2*abs(thrY2)./gapOf(rzY2); oky2 = uy2 < 0.85 & abs(thrY2) > 1e-9;
By2 = V^2*shapeFun(b2*thrY2(oky2)./gapOf(rzY2(oky2)))./max(thrY2(oky2).^2,1e-30);
Ky2 = By2 \ Tauy2(oky2);
fprintf('run2-only Tauy fit:   K2y= %.6e,  R^2 = %+.4f  (x<->y %.2f%% apart)\n', ...
    Ky2, r2fun(Tauy2(oky2), By2*Ky2), 200*abs(K2-Ky2)/(K2+Ky2));

save(fullfile(scriptDir,'centerElectrodeRun2Model.mat'),'k2','b2','K2','Ky2','g0_z','a5','V');
