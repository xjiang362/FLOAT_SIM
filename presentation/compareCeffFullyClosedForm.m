%% compareCeffFullyClosedForm.m
% Fully closed-form Ceff(terminals 1,2) model -- NO simulation lookups
% at all. Every input quantity comes from a fitted function of pose:
%   Cm1R(rz,thx,thy), Cm2R(rz,thx,thy), S_R(rx,ry,rz)  -- earlier fits
%   Cm12(rz,thx,thy), S1(rz,thx,thy), S2(rz,thx,thy)   -- new fits
%
% Tested two ways:
%   1. On the dedicated sensingSweep60 sweep (the data Cm12/S1/S2 were
%      fit on -- an honest in-sample check for those three, but still a
%      genuine out-of-sample test for Cm1R/Cm2R/S_R).
%   2. On the independent 300-pt scattered LHS dataset, which was NEVER
%      used to fit Cm12/S1/S2 -- the real out-of-sample test for the
%      complete assembled model.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'SR_physicsModel.mat'));      % pfit, modelFun
SR_pfit = pfit; SR_fun = modelFun;
load(fullfile(scriptDir, 'Cm1p_fringingFit.mat'));     % pfit, modelFun
Cm1p_pfit = pfit; Cm1p_fun = modelFun;
load(fullfile(scriptDir, 'Cm2p_fringingFit.mat'));     % pfit, modelFun
Cm2p_pfit = pfit; Cm2p_fun = modelFun;
load(fullfile(scriptDir, 'Cm12_S1_S2_models.mat'));    % results, modelFun (Cm12/S1/S2 model)
Cm12_coeffs = results.Cm12.coeffs;
S1_coeffs = results.S1.coeffs;
S2_coeffs = results.S2.coeffs;
Q_fun = modelFun;   % Q_fun(coeffs, [rz, thx_rad, thy_rad])

r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
floatingTerms = [5,6,9,10,13,14,21];

ceffModel = @(rx,ry,rz,thx_rad,thy_rad,thz_rad) localCeff(rx,ry,rz,thx_rad,thy_rad,thz_rad, ...
    SR_pfit,SR_fun,Cm1p_pfit,Cm1p_fun,Cm2p_pfit,Cm2p_fun,Q_fun,Cm12_coeffs,S1_coeffs,S2_coeffs);

%% Test 1: dedicated sensingSweep60
S = load(fullfile(scriptDir, 'sensingSweep60_allMatrices.mat'));   % Cm_all, rx,ry,rz,thx,thy,thz, nPose
rx60=S.rx; ry60=S.ry; rz60=S.rz; thx60=S.thx; thy60=S.thy; thz60=S.thz; nPose60=S.nPose; Cm_all60=S.Cm_all;
thx60_rad = thx60*pi/180; thy60_rad = thy60*pi/180; thz60_rad = thz60*pi/180;
Ceff_gt60 = zeros(nPose60,1);
for pp = 1:nPose60
    Ceff_gt60(pp) = effectiveMutualCapacitance(Cm_all60(:,:,pp), 1, 2, floatingTerms);
end
Ceff_model60 = ceffModel(rx60,ry60,rz60,thx60_rad,thy60_rad,thz60_rad);
r2_60 = r2fun(Ceff_gt60, Ceff_model60);
fprintf('Test 1 (dedicated sensingSweep60): overall R^2 = %.4f\n', r2_60);
blocks = {1:20, 21:40, 41:60};
blockNames = {'rz-only','thx-only','thy-only'};
for b = 1:3
    fprintf('  %8s block R^2 = %.4f\n', blockNames{b}, r2fun(Ceff_gt60(blocks{b}), Ceff_model60(blocks{b})));
end

%% Test 2: independent 300-pt scattered LHS dataset (never used to fit Cm12/S1/S2)
T = load(fullfile(scriptDir, 'capSweep300_withGND_allMatrices.mat'));
rx3=T.rx; ry3=T.ry; rz3=T.rz; thx3=T.thx; thy3=T.thy; thz3=T.thz; nPose3=T.nPose; Cm_all3=T.Cm_all;
thx3_rad = thx3*pi/180; thy3_rad = thy3*pi/180; thz3_rad = thz3*pi/180;
Ceff_gt300 = zeros(nPose3,1);
for pp = 1:nPose3
    Ceff_gt300(pp) = effectiveMutualCapacitance(Cm_all3(:,:,pp), 1, 2, floatingTerms);
end
Ceff_model300 = ceffModel(rx3,ry3,rz3,thx3_rad,thy3_rad,thz3_rad);
r2_300 = r2fun(Ceff_gt300, Ceff_model300);
fprintf('\nTest 2 (independent 300-pt scattered LHS): overall R^2 = %.4f\n', r2_300);

%% Plot: Test 1 (dedicated sweep, per-axis) + Test 2 (scattered, vs rz)
figure('Color','w','Name','Fully closed-form Ceff model vs COMSOL', 'Position',[50 50 1600 800]);
blockX = {rz60*1e6, thx60, thy60};
blockXNames = {'rz [\mum]', '\theta_x [deg]', '\theta_y [deg]'};
for b = 1:3
    subplot(2,3,b); hold on; grid on;
    idx = blocks{b};
    r2b = r2fun(Ceff_gt60(idx), Ceff_model60(idx));
    scatter(blockX{b}(idx), Ceff_gt60(idx), 30, 'k', 'filled', 'DisplayName', 'Ground truth');
    scatter(blockX{b}(idx), Ceff_model60(idx), 30, 'r', 'filled', 'DisplayName', sprintf('Fully closed-form model (R^2=%.4f)', r2b));
    xlabel(blockXNames{b}); ylabel('C_{eff} [F]');
    title(['Test 1: dedicated sweep, ' blockNames{b}]);
    legend('Location','best');
end

subplot(2,3,[4 5 6]); hold on; grid on;
scatter(rz3*1e6, Ceff_gt300, 15, 'k', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'Ground truth (COMSOL)');
scatter(rz3*1e6, Ceff_model300, 15, 'r', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', sprintf('Fully closed-form model (R^2=%.4f)', r2_300));
xlabel('rz [\mum]'); ylabel('C_{eff} [F]');
title('Test 2: independent 300-pt scattered LHS (never used to fit Cm12/S1/S2)');
legend('Location','best');

exportgraphics(gcf, fullfile(scriptDir, 'Ceff_fullyClosedForm.png'));
fprintf('\nSaved plot to Ceff_fullyClosedForm.png\n');

function Ceff = localCeff(rx,ry,rz,thx_rad,thy_rad,thz_rad, ...
        SR_pfit,SR_fun,Cm1p_pfit,Cm1p_fun,Cm2p_pfit,Cm2p_fun,Q_fun,Cm12_coeffs,S1_coeffs,S2_coeffs)
    Cm1R = Cm1p_fun(Cm1p_pfit, [rx, ry, rz, thx_rad, thy_rad, thz_rad]);
    Cm2R = Cm2p_fun(Cm2p_pfit, [rx, ry, rz, thx_rad, thy_rad, thz_rad]);
    SR   = SR_fun(SR_pfit, [rx, ry, rz]);
    Cm12 = Q_fun(Cm12_coeffs, [rz, thx_rad, thy_rad]);
    S1   = Q_fun(S1_coeffs,   [rz, thx_rad, thy_rad]);
    S2   = Q_fun(S2_coeffs,   [rz, thx_rad, thy_rad]);

    C_M = Cm12 + Cm1R.*Cm2R./SR;
    C1 = S1 - Cm12 - Cm1R.*(Cm1R+Cm2R)./SR;
    C2 = S2 - Cm12 - Cm2R.*(Cm1R+Cm2R)./SR;
    Ceff = C_M + C1.*C2./(C1+C2);
end
