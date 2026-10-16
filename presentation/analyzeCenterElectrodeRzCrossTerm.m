%% analyzeCenterElectrodeRzCrossTerm.m
% Compares Taux(thx) slope at rz=-7um vs rz=+7um (and Tauy(thy)
% similarly) to test/quantify the rz*thx, rz*thy cross-coupling that
% the original rz=0-only sweep couldn't see.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

raw = readmatrix(fullfile(scriptDir, 'comsol_centerElectrode_rzCrossTerm200_results.txt'), 'NumHeaderLines', 5);
rz=raw(:,3); thx=raw(:,4); thy=raw(:,5);
Taux=raw(:,10); Tauy=raw(:,11);

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

%% Taux(thx) at rz=-7um vs +7um
blk1 = 1:100;
thxB = thx(blk1)*pi/180; rzB = rz(blk1); TauxB = Taux(blk1);
mMinus = rzB < 0; mPlus = rzB > 0;

pMinus = polyfit(thxB(mMinus), TauxB(mMinus), 1);
pPlus  = polyfit(thxB(mPlus),  TauxB(mPlus),  1);
fprintf('Taux(thx) slope at rz=-7um: %.6e  (R^2=%.4f)\n', pMinus(1), r2fun(TauxB(mMinus), polyval(pMinus,thxB(mMinus))));
fprintf('Taux(thx) slope at rz=+7um: %.6e  (R^2=%.4f)\n', pPlus(1),  r2fun(TauxB(mPlus),  polyval(pPlus, thxB(mPlus))));
fprintf('Slope change: %.4f%% \n\n', 100*(pPlus(1)-pMinus(1))/pMinus(1));

% Joint linear model: Taux = a1*thx + a2*(rz*thx) -- solve directly
designX = [thxB, rzB.*thxB];
cX = designX \ TauxB;
fprintf('Joint fit Taux = a1*thx + a2*(rz*thx): a1=%.6e, a2=%.6e, R^2=%.4f\n\n', cX(1), cX(2), r2fun(TauxB, designX*cX));

%% Tauy(thy) at rz=-7um vs +7um
blk2 = 101:200;
thyB = thy(blk2)*pi/180; rzB2 = rz(blk2); TauyB = Tauy(blk2);
mMinus2 = rzB2 < 0; mPlus2 = rzB2 > 0;

qMinus = polyfit(thyB(mMinus2), TauyB(mMinus2), 1);
qPlus  = polyfit(thyB(mPlus2),  TauyB(mPlus2),  1);
fprintf('Tauy(thy) slope at rz=-7um: %.6e  (R^2=%.4f)\n', qMinus(1), r2fun(TauyB(mMinus2), polyval(qMinus,thyB(mMinus2))));
fprintf('Tauy(thy) slope at rz=+7um: %.6e  (R^2=%.4f)\n', qPlus(1),  r2fun(TauyB(mPlus2),  polyval(qPlus, thyB(mPlus2))));
fprintf('Slope change: %.4f%% \n\n', 100*(qPlus(1)-qMinus(1))/qMinus(1));

designY = [thyB, rzB2.*thyB];
cY = designY \ TauyB;
fprintf('Joint fit Tauy = b1*thy + b2*(rz*thy): b1=%.6e, b2=%.6e, R^2=%.4f\n', cY(1), cY(2), r2fun(TauyB, designY*cY));

%% Plot
figure('Color','w','Name','Taux(thx), Tauy(thy) at two rz values', 'Position',[50 50 1200 500]);
subplot(1,2,1); hold on; grid on;
scatter(thxB(mMinus)*180/pi, TauxB(mMinus), 30, 'b', 'filled', 'DisplayName', 'rz=-7\mum');
scatter(thxB(mPlus)*180/pi, TauxB(mPlus), 30, 'r', 'filled', 'DisplayName', 'rz=+7\mum');
xlabel('\theta_x [deg]'); ylabel('Taux [N*m]');
title('Taux(\theta_x) at rz=\pm7\mum'); legend('Location','best');

subplot(1,2,2); hold on; grid on;
scatter(thyB(mMinus2)*180/pi, TauyB(mMinus2), 30, 'b', 'filled', 'DisplayName', 'rz=-7\mum');
scatter(thyB(mPlus2)*180/pi, TauyB(mPlus2), 30, 'r', 'filled', 'DisplayName', 'rz=+7\mum');
xlabel('\theta_y [deg]'); ylabel('Tauy [N*m]');
title('Tauy(\theta_y) at rz=\pm7\mum'); legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'CenterElectrode_TauRzCrossTerm.png'));
fprintf('\nSaved plot to CenterElectrode_TauRzCrossTerm.png\n');
