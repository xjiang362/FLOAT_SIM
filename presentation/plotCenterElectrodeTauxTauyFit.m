%% plotCenterElectrodeTauxTauyFit.m
% Visualizes the Taux/Tauy = Fz_model*theta*(a2+a1*rz) fit against raw MST
% data at the 3 known rz slices (-7,0,+7um), plus the slope(rz) shape
% diagnostic that explains why a linear-in-rz model tops out at R^2~0.6:
% the true slope peaks at rz=0 and falls off asymmetrically on both sides,
% which a straight line a2+a1*rz cannot represent.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'centerElectrodeTauxTauy_smoothed.mat'), 'a1_leq','a2_leq','Fz_model');

raw1 = readmatrix(fullfile(scriptDir, 'comsol_centerElectrodeSweep_rzThxThy150_results.txt'), 'NumHeaderLines', 5);
thx1=raw1(:,4); thy1=raw1(:,5); Taux1=raw1(:,10); Tauy1=raw1(:,11);
blkThx0 = 51:100; blkThy0 = 101:150;

raw2 = readmatrix(fullfile(scriptDir, 'comsol_centerElectrode_rzCrossTerm200_results.txt'), 'NumHeaderLines', 5);
rz2=raw2(:,3); thx2=raw2(:,4); thy2=raw2(:,5); Taux2=raw2(:,10); Tauy2=raw2(:,11);
blkThxBr = 1:100; blkThyBr = 101:200;
thxBr = thx2(blkThxBr); TauxBr = Taux2(blkThxBr); rzBr = rz2(blkThxBr);
mMinus = rzBr < 0; mPlus = rzBr > 0;
thyBr = thy2(blkThyBr); TauyBr = Tauy2(blkThyBr); rzBr2 = rz2(blkThyBr);
mMinus2 = rzBr2 < 0; mPlus2 = rzBr2 > 0;

fitCurve = @(thDeg, rzVal) Fz_model(rzVal, thDeg*pi/180, zeros(size(thDeg))) ...
    .* (thDeg*pi/180) .* (a2_leq + a1_leq*rzVal);

figure('Color','w','Name','Center electrode Taux/Tauy: data vs fit','Position',[50 50 1300 900]);

%% Taux(thx) at 3 rz slices, data vs fit
subplot(2,2,1); hold on; grid on;
thxGrid = linspace(-max(abs(thx1(blkThx0))), max(abs(thx1(blkThx0))), 100)';
scatter(thx1(blkThx0), Taux1(blkThx0), 24, [0.2 0.2 0.8], 'filled', 'DisplayName','MST, rz=0');
plot(thxGrid, fitCurve(thxGrid, 0), '-', 'Color',[0.2 0.2 0.8], 'LineWidth',1.5, 'HandleVisibility','off');
scatter(thxBr(mMinus), TauxBr(mMinus), 24, [0.8 0.2 0.2], 'filled', 'DisplayName','MST, rz=-7\mum');
plot(thxGrid, fitCurve(thxGrid, -7e-6), '-', 'Color',[0.8 0.2 0.2], 'LineWidth',1.5, 'HandleVisibility','off');
scatter(thxBr(mPlus), TauxBr(mPlus), 24, [0.2 0.7 0.3], 'filled', 'DisplayName','MST, rz=+7\mum');
plot(thxGrid, fitCurve(thxGrid, 7e-6), '-', 'Color',[0.2 0.7 0.3], 'LineWidth',1.5, 'HandleVisibility','off');
xlabel('\theta_x [deg]'); ylabel('Taux [N\cdotm]');
title('Taux(\theta_x): MST data (dots) vs Fz\cdot\theta_x\cdot(a2+a1\cdotrz) fit (lines)');
legend('Location','best');

%% Tauy(thy) at 3 rz slices, data vs fit
subplot(2,2,2); hold on; grid on;
thyGrid = linspace(-max(abs(thy1(blkThy0))), max(abs(thy1(blkThy0))), 100)';
fitCurveY = @(thDeg, rzVal) Fz_model(rzVal, zeros(size(thDeg)), (thDeg*pi/180)) ...
    .* (thDeg*pi/180) .* (a2_leq + a1_leq*rzVal);
scatter(thy1(blkThy0), Tauy1(blkThy0), 24, [0.2 0.2 0.8], 'filled', 'DisplayName','MST, rz=0');
plot(thyGrid, fitCurveY(thyGrid, 0), '-', 'Color',[0.2 0.2 0.8], 'LineWidth',1.5, 'HandleVisibility','off');
scatter(thyBr(mMinus2), TauyBr(mMinus2), 24, [0.8 0.2 0.2], 'filled', 'DisplayName','MST, rz=-7\mum');
plot(thyGrid, fitCurveY(thyGrid, -7e-6), '-', 'Color',[0.8 0.2 0.2], 'LineWidth',1.5, 'HandleVisibility','off');
scatter(thyBr(mPlus2), TauyBr(mPlus2), 24, [0.2 0.7 0.3], 'filled', 'DisplayName','MST, rz=+7\mum');
plot(thyGrid, fitCurveY(thyGrid, 7e-6), '-', 'Color',[0.2 0.7 0.3], 'LineWidth',1.5, 'HandleVisibility','off');
xlabel('\theta_y [deg]'); ylabel('Tauy [N\cdotm]');
title('Tauy(\theta_y): MST data (dots) vs fit (lines)');
legend('Location','best');

%% Slope(rz) shape diagnostic: 3 known points vs the fitted line a2+a1*rz
p0 = polyfit(thx1(blkThx0)*pi/180, Taux1(blkThx0), 1);
pMinus = polyfit(thxBr(mMinus)*pi/180, TauxBr(mMinus), 1);
pPlus  = polyfit(thxBr(mPlus)*pi/180,  TauxBr(mPlus),  1);
rzPts = [-7e-6; 0; 7e-6];
slopeRaw = [pMinus(1); p0(1); pPlus(1)];

subplot(2,2,3); hold on; grid on;
scatter(rzPts*1e6, slopeRaw, 80, 'k', 'filled', 'DisplayName', 'actual d(Taux)/d\theta_x per slice');
rzGrid = linspace(-14,14,100)';
plot(rzGrid, a2_leq + a1_leq*rzGrid*1e-6, 'r-', 'LineWidth',1.5, 'DisplayName', 'fitted a2+a1\cdotrz (this line CANNOT show a peak)');
xlabel('rz [\mum]'); ylabel('Taux slope [N\cdotm/rad]');
title({'Why R^2\approx0.6: true slope(rz) PEAKS at rz=0,','a straight line structurally cannot match that'});
legend('Location','best');

subplot(2,2,4); axis off;
text(0,0.9, sprintf('slope(-7um) = %.3e', slopeRaw(1)), 'FontSize',11);
text(0,0.7, sprintf('slope(0)    = %.3e  <- peak', slopeRaw(2)), 'FontSize',11);
text(0,0.5, sprintf('slope(+7um) = %.3e', slopeRaw(3)), 'FontSize',11);
text(0,0.25, sprintf('Only 3 rz points known so far.'), 'FontSize',10);
text(0,0.1, sprintf('Fine rz sweep (comsol_centerElectrodeTauRzFine50.txt)'), 'FontSize',10);
text(0,-0.05, sprintf('queued to resolve the true shape.'), 'FontSize',10);

exportgraphics(gcf, fullfile(scriptDir, 'CenterElectrode_TauxTauyFit_diagnostic.png'));
fprintf('Saved plot to CenterElectrode_TauxTauyFit_diagnostic.png\n');
