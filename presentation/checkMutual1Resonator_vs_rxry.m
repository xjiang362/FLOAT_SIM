%% checkMutual1Resonator_vs_rxry.m
% Does the MUTUAL capacitance between terminal 1 and the resonator (the
% proof mass, terminal 21) depend on the in-plane motion rx/ry?
%
% Physics expectation: the parallel-plate term depends only on the gap, so
% it is blind to rx/ry. BUT the quadrant sub-electrodes sit FLUSH with the
% proof-mass edge (plate half-width 2975 um; quad1-a spans x in
% [1525,2975]), so lateral motion changes the plate-electrode OVERLAP and
% therefore must move Cm(1,21). This is the same mechanism that produces
% the lateral force Fx, so the two must be consistent.
%
% Method: rz and the tilts dominate Cm by orders of magnitude, so fit and
% remove that dependence first, then test the RESIDUAL against rx and ry.
% Testing raw Cm against rx would just measure the LHS sampling.
%
% Data: comsol_capacitanceSweep300_* (300-pose full-6D LHS, 21x21 matrices),
% parsed by parseAllCapSweepData.m into capSweep300_allMatrices.mat.
%
% SYMBOLS: Cm [F] mutual capacitance matrix entry; rx,ry,rz [m] plate
%   translations; thx,thy,thz [deg as stored, rad in use] plate rotations;
%   G [m] centroid gap = g0 - rz; eps0 [F/m]; resid [F] Cm after the
%   gap/tilt model is removed.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
eps0 = 8.8541878128e-12; g0 = 15e-6; um = 1e-6;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

%% ---- poses ----
fid = fopen(fullfile(scriptDir,'comsol_capacitanceSweep300.txt'));
P = struct();
while true
    line = fgetl(fid); if ~ischar(line), break; end
    tok = regexp(strtrim(line),'^(\w+)\s+"([^"]+)"','tokens','once');
    if isempty(tok), continue; end
    P.(tok{1}) = str2double(strsplit(strtrim(tok{2})))';
end
fclose(fid);
rx=P.rx; ry=P.ry; rz=P.rz; thx=P.thx*pi/180; thy=P.thy*pi/180; thz=P.thz*pi/180;

load(fullfile(scriptDir,'capSweep300_allMatrices.mat'),'Cm_all','nPose');
fprintf('Loaded %d poses, matrices %s\n', nPose, mat2str(size(Cm_all,[1 2])));

RES = 21;                                  % resonator / proof mass terminal
Cm1R = squeeze(Cm_all(1,RES,:));
fprintf('Cm(1,%d): mean %.4e F, range %.4e to %.4e\n', RES, mean(Cm1R), min(Cm1R), max(Cm1R));

%% ---- which electrode is terminal 1? ----
g = trapGeometry(); G = g0 - rz;
shapes = struct('name',{},'S',{});
shapes(end+1) = struct('name','sense1','S',g.sense1);
shapes(end+1) = struct('name','sense2','S',g.sense2);
shapes(end+1) = struct('name','centralA','S',[g.centralA.px; g.centralA.mx]);
shapes(end+1) = struct('name','centralB','S',[g.centralB.my; g.centralB.py]);
for q=1:4
    shapes(end+1) = struct('name',sprintf('quad%d-a',q),'S',g.quad(q).a);
    shapes(end+1) = struct('name',sprintf('quad%d-b',q),'S',g.quad(q).b);
end
% Correlation CANNOT discriminate here: every footprint scales as 1/G, so
% they all correlate ~0.99 with each other. Identify instead by which
% shape leaves the SMALLEST residual after a gain+offset fit, and by
% whether that gain is near 1 (a shape of the wrong area needs a gain far
% from unity to compensate).
fprintf('\n%-12s %9s %10s %12s %10s\n','shape','corr','gain','resid RMS','resid %%');
best = struct('rms',Inf);
for j=1:numel(shapes)
    Ca = shapeCap(shapes(j).S, G, thx, thy, eps0);
    Xj = [Ca, ones(nPose,1)];
    bj = Xj\abs(Cm1R);
    rj = abs(Cm1R) - Xj*bj;
    fprintf('%-12s %9.4f %10.3f %12.4e %9.1f%%\n', shapes(j).name, ...
        corr(abs(Cm1R),Ca), bj(1), rms(rj), 100*rms(rj)/mean(abs(Cm1R)));
    if rms(rj) < best.rms
        best = struct('rms',rms(rj),'j',j,'b',bj,'resid',rj,'gain',bj(1));
    end
end
kb = best.j;
fprintf('\n-> terminal 1 identified as %s (lowest residual, gain %.3f)\n', ...
    shapes(kb).name, best.gain);

%% ---- remove the gap/tilt dependence, then test rx, ry ----
Cana = shapeCap(shapes(kb).S, G, thx, thy, eps0);
X = [Cana, ones(nPose,1)];
beta = best.b;
resid = best.resid;
fprintf('\ngap/tilt model R^2 on |Cm(1,%d)| = %.5f ; residual RMS = %.4e F (%.2f%% of mean)\n', ...
    RES, r2fun(abs(Cm1R), X*beta), rms(resid), 100*rms(resid)/mean(abs(Cm1R)));

fprintf('\n--- residual vs in-plane DOF ---\n');
nm = {'rx','ry','thz'}; vv = {rx, ry, thz};
for i=1:3
    [cc2,pv] = corr(vv{i}, resid);
    fprintf('  corr(resid, %4s) = %+.4f   (p = %.3g)\n', nm{i}, cc2, pv);
end

% joint linear fit of the residual on rx, ry (and thz), with t-stats
Z = [rx, ry, thz, ones(nPose,1)];
b = Z\resid;
r = resid - Z*b;
s2 = sum(r.^2)/(nPose-size(Z,2));
se = sqrt(diag(s2*inv(Z'*Z)));
fprintf('\n  joint fit resid = a*rx + b*ry + c*thz + d\n');
fprintf('  %-6s %14s %14s %9s\n','term','coeff','std err','t');
lbl = {'rx [F/m]','ry [F/m]','thz [F/rad]','const [F]'};
for i=1:4
    fprintf('  %-12s %14.4e %14.4e %9.2f\n', lbl{i}, b(i), se(i), b(i)/se(i));
end
fprintf('  R^2 of the residual explained by rx,ry,thz = %.4f\n', r2fun(resid, Z*b));

% predicted sensitivity if the mechanism is overlap truncation:
% dCm/drx = eps0*dY/G  (one truncating edge), evaluated at the mean gap
dY = g.actSide;
fprintf('\n  overlap-truncation prediction |dCm/drx| = eps0*dY/G = %.4e F/m (at G = %.1f um)\n', ...
    eps0*dY/mean(G), mean(G)*1e6);
fprintf('  measured |rx| coefficient                = %.4e F/m  -> ratio %.2f\n', ...
    abs(b(1)), abs(b(1))/(eps0*dY/mean(G)));

%% ---- plot ----
figure('Color','w','Position',[50 50 1340 430]);
subplot(1,3,1); hold on; grid on;
scatter(rz*1e6, abs(Cm1R)*1e12, 18, 'filled','MarkerFaceAlpha',0.5);
xlabel('r_z [\mum]'); ylabel('|C_m(1,21)| [pF]');
title('raw: dominated by the gap');

subplot(1,3,2); hold on; grid on;
scatter(rx*1e6, resid*1e15, 20, 'filled','MarkerFaceAlpha',0.55,'MarkerFaceColor',[0.2 0.3 0.75]);
xg = linspace(min(rx),max(rx),100)';
plot(xg*1e6, (b(1)*xg + b(4))*1e15,'r-','LineWidth',1.8);
xlabel('r_x [\mum]'); ylabel('residual [fF]');
title(sprintf('residual vs r_x  (t = %.1f)', b(1)/se(1)));

subplot(1,3,3); hold on; grid on;
scatter(ry*1e6, resid*1e15, 20, 'filled','MarkerFaceAlpha',0.55,'MarkerFaceColor',[0.85 0.33 0.10]);
yg = linspace(min(ry),max(ry),100)';
plot(yg*1e6, (b(2)*yg + b(4))*1e15,'r-','LineWidth',1.8);
xlabel('r_y [\mum]'); ylabel('residual [fF]');
title(sprintf('residual vs r_y  (t = %.1f)', b(2)/se(2)));

sgtitle('Mutual capacitance terminal 1 <-> resonator: is there an r_x/r_y dependence?');
outPng = fullfile(scriptDir,'Mutual1Resonator_vs_rxry.png');
exportgraphics(gcf,outPng);
fprintf('\nSaved plot to %s\n', outPng);
