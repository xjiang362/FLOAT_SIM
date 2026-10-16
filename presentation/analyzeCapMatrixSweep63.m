%% analyzeCapMatrixSweep63.m
% Dedicated single-axis capacitance-matrix sweep (63 poses: rx, rz, thx,
% 21 points each with the other five DOF at exactly zero).
%
% NETWORK ASSUMPTION (as specified): terminals 1 and 2 are the two live
% sensing ports; the RESONATOR (terminal 21) FLOATS; every other terminal
% -- including the frame, terminal 22 -- is held at ground. A fixed-voltage
% terminal contributes no AC current, so it is simply dropped from the
% reduction; the floating resonator must instead be ELIMINATED (Schur /
% Kron reduction, i.e. the Y-delta conversion).
%
% Two quantities are extracted at each pose:
%   C12eff  the measured two-terminal capacitance between ports 1 and 2,
%           via effectiveMutualCapacitance(Cm,1,2,21)
%   Crg     the resonator node's total capacitance to GROUND, i.e. the
%           row sum of the resonator minus its couplings to the two live
%           ports:   Crg = S_R - C1R - C2R
%           This is the term that sits in the denominator of the Y-delta
%           branch, C_delta = C1R*C2R/S_R, so it modulates the sensed
%           signal without appearing in it directly.
%
% SYMBOLS (all capacitances in F; plotted in fF or pF as marked)
%   Cm       22x22 mutual (lumped) capacitance matrix at one pose     [F]
%   C1R,C2R  mutual capacitance, port 1 or 2 <-> resonator            [F]
%   S_R      resonator node row sum, sum_k Cm(21,k)                   [F]
%   Crg      resonator <-> ground, = S_R - C1R - C2R                  [F]
%   CRF      resonator <-> frame (terminal 22), a component of Crg    [F]
%   C_delta  branch transferred onto 1<->2 by eliminating the resonator[F]
%   C12      direct port-to-port mutual capacitance, Cm(1,2)          [F]
%   C12eff   measured two-terminal capacitance between ports 1 and 2   [F]
%   rx,rz    plate translations                                       [m]
%   thx      plate tilt about x                                      [deg]

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

nTerm = 22; RES = 21; FR = 22; P1 = 1; P2 = 2;
raw = readmatrix(fullfile(scriptDir,'comsol_capMatrixSweep_rxRzThx63_results.txt'),'NumHeaderLines',5);
nPose = size(raw,1)/nTerm;
fprintf('Parsed %d poses x %d terminals (%d columns)\n', nPose, nTerm, size(raw,2));

poseRows = 1:nTerm:size(raw,1);
rx = raw(poseRows,1); ry = raw(poseRows,2); rz = raw(poseRows,3);
thx = raw(poseRows,4); thy = raw(poseRows,5); thz = raw(poseRows,6);

Cm = zeros(nTerm,nTerm,nPose);
for p = 1:nPose
    Cm(:,:,p) = raw((p-1)*nTerm+1 : p*nTerm, 7:28);
end

% symmetry sanity check
asym = zeros(nPose,1);
for p=1:nPose
    M = Cm(:,:,p); asym(p) = max(abs(M(:)-reshape(M',[],1)))/max(abs(M(:)));
end
fprintf('max relative matrix asymmetry: %.3e\n', max(asym));

%% ---- blocks (generator order: rx, rz, thx) ----
blk = {1:21, 22:42, 43:63};
blkVar = {rx, rz, thx};
blkName = {'rx','rz','thx'};
blkUnit = {'\mum','\mum','deg'};
blkScale = [1e6, 1e6, 1];
for i=1:3
    v = blkVar{i}(blk{i});
    fprintf('block %-4s rows %2d-%2d : %.4g .. %.4g\n', blkName{i}, blk{i}(1), blk{i}(end), ...
        min(v)*blkScale(i), max(v)*blkScale(i));
end

%% ---- extract the two quantities ----
C12eff = zeros(nPose,1); CM = zeros(nPose,1);
C1R = squeeze(Cm(RES,P1,:)); C2R = squeeze(Cm(RES,P2,:));
SR  = squeeze(sum(Cm(RES,:,:),2));
CRF = squeeze(Cm(RES,FR,:));
C12 = squeeze(Cm(P1,P2,:));
for p=1:nPose
    [C12eff(p), CM(p)] = effectiveMutualCapacitance(Cm(:,:,p), P1, P2, RES);
end
Crg    = SR - C1R - C2R;
Cdelta = C1R.*C2R./SR;

%% ---- zero-pose repeatability (rows 11, 32, 53 are the same pose) ----
z = [11, 32, 53];
fprintf('\n--- repeatability at the identical all-zero pose (poses %d, %d, %d) ---\n', z);
fprintf('  C12eff : %s fF   spread %.3f %%\n', ...
    sprintf('%.4f ', C12eff(z)*1e15), 100*(max(C12eff(z))-min(C12eff(z)))/mean(C12eff(z)));
fprintf('  Crg    : %s pF   spread %.3f %%\n', ...
    sprintf('%.5f ', Crg(z)*1e12), 100*(max(Crg(z))-min(Crg(z)))/mean(Crg(z)));

%% ---- per-axis summary ----
fprintf('\n--- variation of each quantity along each axis ---\n');
fprintf('%-5s | %12s %12s %9s | %12s %12s %9s\n', ...
    'axis','C12eff min','C12eff max','swing','Crg min','Crg max','swing');
for i=1:3
    b = blk{i};
    fprintf('%-5s | %10.3f fF %10.3f fF %8.1f%% | %10.4f pF %10.4f pF %8.1f%%\n', blkName{i}, ...
        min(C12eff(b))*1e15, max(C12eff(b))*1e15, ...
        100*(max(C12eff(b))-min(C12eff(b)))/mean(C12eff(b)), ...
        min(Crg(b))*1e12, max(Crg(b))*1e12, ...
        100*(max(Crg(b))-min(Crg(b)))/mean(Crg(b)));
end

fprintf('\n--- composition at the nominal pose (pose %d) ---\n', z(1));
k = z(1);
fprintf('  C1R %.4f pF | C2R %.4f pF | Crg %.4f pF | CRF %.4f pF (%.1f%% of Crg)\n', ...
    C1R(k)*1e12, C2R(k)*1e12, Crg(k)*1e12, CRF(k)*1e12, 100*CRF(k)/Crg(k));
fprintf('  S_R %.4f pF | C12 direct %.4f fF | C_delta %.4f fF | C_M %.4f fF | C12eff %.3f fF\n', ...
    SR(k)*1e12, C12(k)*1e15, Cdelta(k)*1e15, CM(k)*1e15, C12eff(k)*1e15);

%% ---- plot: 2 rows (C12eff, Crg) x 3 cols (rz, rx, thx) ----
order = [2 1 3];                       % display as rz, rx, thx
figure('Color','w','Position',[30 30 1380 720]);
for c = 1:3
    i = order(c); b = blk{i};
    x = blkVar{i}(b)*blkScale(i);
    [x,o] = sort(x);

    subplot(2,3,c); hold on; grid on; box on;
    plot(x, C12eff(b(o))*1e15,'o-','MarkerSize',4,'LineWidth',1.4,'Color',[0.12 0.25 0.60]);
    xline(0,'k:'); xlabel(sprintf('%s [%s]', blkName{i}, blkUnit{i}));
    ylabel('C_{12,eff} [fF]');
    title(sprintf('C_{12,eff} vs %s', blkName{i}));

    subplot(2,3,3+c); hold on; grid on; box on;
    plot(x, Crg(b(o))*1e12,'s-','MarkerSize',4,'LineWidth',1.4,'Color',[0.85 0.33 0.10]);
    xline(0,'k:'); xlabel(sprintf('%s [%s]', blkName{i}, blkUnit{i}));
    ylabel('C_{rg} [pF]');
    title(sprintf('C_{rg} (resonator\\leftrightarrowground) vs %s', blkName{i}));
end
sgtitle({'Dedicated 63-pose capacitance-matrix sweep', ...
    'ports 1,2 live; resonator (21) floating; all other terminals incl. frame (22) grounded'});
for a = findobj(gcf,'Type','axes')', a.Toolbar.Visible='off'; end
outPng = fullfile(scriptDir,'CapMatrixSweep63_C12eff_Crg.png');
exportgraphics(gcf,outPng,'Resolution',150);
fprintf('\nSaved plot to %s\n', outPng);

save(fullfile(scriptDir,'capMatrixSweep63_parsed.mat'), ...
    'Cm','rx','ry','rz','thx','thy','thz','nPose','nTerm', ...
    'C12eff','CM','Crg','SR','C1R','C2R','CRF','C12','Cdelta','blk');
fprintf('Saved parsed data to capMatrixSweep63_parsed.mat\n');
