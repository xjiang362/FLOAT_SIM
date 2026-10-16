%% diagnoseCenterElectrodeSlopeShape.m
% The smoothed-target refit barely moved a1/a2 (R^2 still ~0.6) even
% though each individual rz-slice is itself a near-perfect line
% (LOO picked plain linear every time) -- so per-point mesh noise isn't
% the bottleneck. This means the FUNCTIONAL FORM slope(rz) = a2+a1*rz is
% probably wrong. Directly compute the per-slice slope at rz=0,-7,+7um
% (3 known points) and see what shape connects them.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw1 = readmatrix(fullfile(scriptDir, 'comsol_centerElectrodeSweep_rzThxThy150_results.txt'), 'NumHeaderLines', 5);
thx1=raw1(:,4); Taux1=raw1(:,10);
blkThx0 = 51:100;
p0 = polyfit(thx1(blkThx0)*pi/180, Taux1(blkThx0), 1);

raw2 = readmatrix(fullfile(scriptDir, 'comsol_centerElectrode_rzCrossTerm200_results.txt'), 'NumHeaderLines', 5);
rz2=raw2(:,3); thx2=raw2(:,4); Taux2=raw2(:,10);
blkThxBr = 1:100;
thxBr = thx2(blkThxBr)*pi/180; TauxBr = Taux2(blkThxBr); rzBr = rz2(blkThxBr);
mMinus = rzBr < 0; mPlus = rzBr > 0;
pMinus = polyfit(thxBr(mMinus), TauxBr(mMinus), 1);
pPlus  = polyfit(thxBr(mPlus),  TauxBr(mPlus),  1);

rzPts = [-7e-6; 0; 7e-6];
slopes = [pMinus(1); p0(1); pPlus(1)];
fprintf('rz [um]   slope [Nm/rad]\n');
for i=1:3
    fprintf('%8.2f   %.6e\n', rzPts(i)*1e6, slopes(i));
end
fprintf('\nRatio slope(+7)/slope(-7) = %.4f\n', slopes(3)/slopes(1));
fprintf('Ratio slope(0)/slope(-7)  = %.4f\n', slopes(2)/slopes(1));
fprintf('Ratio slope(+7)/slope(0)  = %.4f\n', slopes(3)/slopes(2));

% Candidate shapes using the already-fitted Fz model's Geff (same rz
% dependence source): does slope ~ 1/Geff(rz)^n for some n?
load(fullfile(scriptDir, 'centerElectrodeModel.mat'), 'k_fit','a2','a4','a5','g0_z');
Geff0 = @(rz) g0_z - a5*rz;   % thx=thy=0 along this rz axis
fprintf('\nGeff(rz) at the 3 points (thx=thy=0):\n');
for i=1:3
    fprintf('  rz=%.2fum: Geff=%.4e m\n', rzPts(i)*1e6, Geff0(rzPts(i)));
end
fprintf('\nGeff ratios: Geff(-7)/Geff(+7)=%.4f, Geff(-7)/Geff(0)=%.4f, Geff(0)/Geff(+7)=%.4f\n', ...
    Geff0(rzPts(1))/Geff0(rzPts(3)), Geff0(rzPts(1))/Geff0(rzPts(2)), Geff0(rzPts(2))/Geff0(rzPts(3)));

for n = 1:4
    pred = 1./Geff0(rzPts).^n;
    predRatio = pred/pred(2)*slopes(2);
    fprintf('\n1/Geff^%d prediction (anchored at rz=0): %.4e, %.4e, %.4e  (actual: %.4e, %.4e, %.4e)\n', ...
        n, predRatio(1), predRatio(2), predRatio(3), slopes(1), slopes(2), slopes(3));
end
