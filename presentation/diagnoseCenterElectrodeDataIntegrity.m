%% diagnoseCenterElectrodeDataIntegrity.m
% Before changing the fitting equation again, check whether the DATA is
% self-consistent. Two independent worries:
%
%  (1) The rz=0 slice comes from one COMSOL run (rzThxThy150) and the
%      rz=+-7um slices from a different run (rzCrossTerm200). If the two
%      runs differ in setup (voltage, probes, energized electrode), the
%      slope(rz) "shape" is an artifact of stitching, not physics.
%      TEST: the Fz model was fitted ONLY on run 1 and validated to
%      R^2~0.98-0.99 there. Predict Fz on run 2's points. If run 2 is
%      the same setup, the run-1 model must still predict it well.
%
%  (2) Both sweeps sized their max tilt from the rz=0 gap (g0=15um). At
%      rz=+7um the gap is only ~8um, so the same tilt swings the plate
%      edge to a near-collision minimum gap. Parallel-plate pressure
%      ~1/g^2 blows up there, so tau vs thx stops being linear and a
%      linear polyfit slope becomes meaningless.
%      TEST: compute min gap under the plate at every point.
%
% Analytic expectation for a symmetric centered electrode (derived by
% expanding tau = integral of y*eps*V^2/(2*(g-y*th)^2) dA over a
% symmetric span): tau ~ (4/3)*b^3*th/g^3, i.e. slope ~ 1/g^3,
% MONOTONICALLY INCREASING as rz increases (gap shrinks).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
load(fullfile(scriptDir, 'centerElectrodeModel.mat'), 'k_fit','a2','a4','a5','g0_z');

V = 30;
L = 6000e-6;                         % plate side [m]
Geff_fun = @(rz,thxr,thyr) g0_z - a5*rz - a2*thxr.^2 - a4*thyr.^2;
Fz_model = @(rz,thxr,thyr) (k_fit*V^2) ./ max(Geff_fun(rz,thxr,thyr),1e-9).^2;
r2fun = @(y,yhat) 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);

%% ---------- Load both runs ----------
raw1 = readmatrix(fullfile(scriptDir,'comsol_centerElectrodeSweep_rzThxThy150_results.txt'),'NumHeaderLines',5);
raw2 = readmatrix(fullfile(scriptDir,'comsol_centerElectrode_rzCrossTerm200_results.txt'),'NumHeaderLines',5);

get = @(R) deal(R(:,3),R(:,4),R(:,5),R(:,9),R(:,10),R(:,11));
[rz1,thx1,thy1,Fz1,Taux1,Tauy1] = get(raw1);
[rz2,thx2,thy2,Fz2,Taux2,Tauy2] = get(raw2);

%% ---------- TEST 1: is run 2 the same physical setup as run 1? ----------
fprintf('===== TEST 1: cross-run consistency via the Fz model =====\n');
Fz1_pred = Fz_model(rz1, thx1*pi/180, thy1*pi/180);
Fz2_pred = Fz_model(rz2, thx2*pi/180, thy2*pi/180);
fprintf('Run 1 (fitted on this run):   R^2 = %+.4f\n', r2fun(Fz1,Fz1_pred));
fprintf('Run 2 (held out entirely):    R^2 = %+.4f\n', r2fun(Fz2,Fz2_pred));
fprintf('Run 2 median ratio Fz_sim/Fz_pred = %.4f  (1.00 => same setup)\n', median(Fz2./Fz2_pred));
fprintf('Run 2 Fz range: [%.4e, %.4e] N\n', min(Fz2), max(Fz2));
fprintf('Run 1 Fz range: [%.4e, %.4e] N\n\n', min(Fz1), max(Fz1));

%% ---------- TEST 2: geometric clearance at every point ----------
% Tilt about the centroid lifts one plate edge and lowers the other by
% (L/2)*sin(th). Minimum gap anywhere under the plate:
minGap = @(rz,thxDeg,thyDeg) (g0_z - rz) - (L/2)*( abs(sind(thxDeg)) + abs(sind(thyDeg)) );

fprintf('===== TEST 2: minimum gap under the plate =====\n');
fprintf('%-28s %12s %12s %10s\n','block','minGap_min','minGap_med','#<2um');
report = @(name,rz,thxD,thyD) fprintf('%-28s %10.3f um %10.3f um %10d\n', name, ...
    min(minGap(rz,thxD,thyD))*1e6, median(minGap(rz,thxD,thyD))*1e6, ...
    sum(minGap(rz,thxD,thyD) < 2e-6));

report('run1 thx-only (rz=0)',      rz1(51:100), thx1(51:100), thy1(51:100));
report('run1 thy-only (rz=0)',      rz1(101:150),thx1(101:150),thy1(101:150));
mM = rz2(1:100)<0; mP = rz2(1:100)>0;
idxM = find(mM); idxP = find(mP);
report('run2 thx-block rz=-7um',    rz2(idxM),   thx2(idxM),   thy2(idxM));
report('run2 thx-block rz=+7um',    rz2(idxP),   thx2(idxP),   thy2(idxP));
fprintf('\n');

%% ---------- TEST 3: is Taux vs thx actually LINEAR in each slice? ----------
% A meaningful "slope" requires linearity. Compare linear vs cubic fit,
% and re-fit the slope using only geometrically safe points.
fprintf('===== TEST 3: per-slice linearity of Taux(thx) =====\n');
slices = struct('name',{},'rz',{},'th',{},'tau',{});
slices(1) = struct('name','rz = 0   (run 1)','rz',0,     'th',thx1(51:100),'tau',Taux1(51:100));
slices(2) = struct('name','rz = -7um (run 2)','rz',-7e-6,'th',thx2(idxM),  'tau',Taux2(idxM));
slices(3) = struct('name','rz = +7um (run 2)','rz',+7e-6,'th',thx2(idxP),  'tau',Taux2(idxP));

fprintf('%-20s %11s %11s %13s %13s %8s\n','slice','R2_linear','R2_cubic','slope_all','slope_safe','n_safe');
slopeAll = zeros(3,1); slopeSafe = zeros(3,1); gapNom = zeros(3,1);
for i = 1:3
    s = slices(i);
    thr = s.th*pi/180;
    pL = polyfit(thr,s.tau,1);  pC = polyfit(thr,s.tau,3);
    safe = minGap(s.rz, s.th, 0) > 4e-6;      % keep >=4um clearance
    pS = polyfit(thr(safe), s.tau(safe), 1);
    slopeAll(i)=pL(1); slopeSafe(i)=pS(1); gapNom(i)=g0_z-s.rz;
    fprintf('%-20s %11.4f %11.4f %13.4e %13.4e %8d\n', s.name, ...
        r2fun(s.tau,polyval(pL,thr)), r2fun(s.tau,polyval(pC,thr)), pL(1), pS(1), sum(safe));
end

%% ---------- TEST 4: does slope ~ 1/g^3 hold on the SAFE slopes? ----------
fprintf('\n===== TEST 4: slope(rz) vs analytic 1/g^3 law =====\n');
fprintf('%-12s %10s %14s %14s %14s\n','rz[um]','gap[um]','slope_all','slope_safe','1/g^3 pred');
for n = 3
    pred = (1./gapNom.^n);
    pred = pred / pred(1) * slopeSafe(1);      % anchor at the clean rz=0 slice
    for i = 1:3
        fprintf('%-12.1f %10.2f %14.4e %14.4e %14.4e\n', ...
            (g0_z-gapNom(i))*1e6, gapNom(i)*1e6, slopeAll(i), slopeSafe(i), pred(i));
    end
end

%% ---------- Plot: raw Taux vs thx per slice, collision zone shaded ----------
figure('Color','w','Position',[60 60 1250 480]);
cols = [0.20 0.30 0.75; 0.80 0.25 0.20; 0.15 0.60 0.30];
for i = 1:3
    subplot(1,3,i); hold on; grid on;
    s = slices(i); thr = s.th*pi/180;
    safe = minGap(s.rz, s.th, 0) > 4e-6;
    scatter(s.th(safe),  s.tau(safe),  26, cols(i,:), 'filled', 'DisplayName','clearance > 4\mum');
    scatter(s.th(~safe), s.tau(~safe), 40, [0.6 0.6 0.6], 'x', 'LineWidth',1.2, 'DisplayName','near-collision');
    tg = linspace(min(s.th),max(s.th),100)';
    plot(tg, polyval(polyfit(thr(safe),s.tau(safe),1), tg*pi/180), '-', ...
        'Color',cols(i,:), 'LineWidth',1.6, 'DisplayName','linear fit (safe pts)');
    xlabel('\theta_x [deg]'); ylabel('Taux [N\cdotm]');
    title(sprintf('%s\nnominal gap = %.1f\\mum', s.name, gapNom(i)*1e6));
    legend('Location','best','FontSize',7);
end
sgtitle('Center electrode Taux(\theta_x) per rz slice, near-collision points flagged');
exportgraphics(gcf, fullfile(scriptDir,'CenterElectrode_DataIntegrity.png'));
fprintf('\nSaved plot to CenterElectrode_DataIntegrity.png\n');
