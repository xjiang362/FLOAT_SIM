%% inspectCapMatrices.m
% Identify which Maxwell-matrix terminal corresponds to which physical
% electrode, by correlating each matrix row/column against the ANALYTIC
% capacitance of each decoded footprint over the 60-pose sweep. The
% analytic curves come from geometry alone, so a high correlation is a
% genuine identification rather than a fit.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
g = trapGeometry(); um = g.um;
eps0 = 8.8541878128e-12;

S = load(fullfile(scriptDir,'sensingSweep60_allMatrices.mat'));
fprintf('Variables: %s\n', strjoin(fieldnames(S)', ', '));
fprintf('Cm_all size: %s\n', mat2str(size(S.Cm_all)));
fprintf('nPose = %d\n\n', S.nPose);

rz = S.rz; thx = S.thx*pi/180; thy = S.thy*pi/180;
g0 = 15e-6; G = g0 - rz;

%% Analytic capacitance of every decoded electrode, over the sweep
shapes = struct('name',{},'S',{});
shapes(end+1) = struct('name','sense1 (Cm1p)','S',g.sense1);
shapes(end+1) = struct('name','sense2 (Cm2p)','S',g.sense2);
shapes(end+1) = struct('name','centralA (outer)','S',[g.centralA.px; g.centralA.mx]);
shapes(end+1) = struct('name','centralB (inner)','S',[g.centralB.my; g.centralB.py]);
for q=1:4
    shapes(end+1) = struct('name',sprintf('quad%d-a',q),'S',g.quad(q).a);
    shapes(end+1) = struct('name',sprintf('quad%d-b',q),'S',g.quad(q).b);
end

nS = numel(shapes);
Cana = zeros(S.nPose, nS);
for j = 1:nS
    Cana(:,j) = shapeCap(shapes(j).S, G, thx, thy, eps0);
end

fprintf('Analytic capacitance at the nominal pose (rz=thx=thy=0):\n');
Cnom = zeros(nS,1);
for j=1:nS
    Cnom(j) = shapeCap(shapes(j).S, g0, 0, 0, eps0);
    fprintf('  %-20s %12.5e F\n', shapes(j).name, Cnom(j));
end

%% Correlate each matrix entry against each analytic curve
n = size(S.Cm_all,1);
fprintf('\nMatrix is %dx%d. Best analytic match for each |Cm(i,j)| curve:\n', n, n);
fprintf('%-10s %-22s %9s %12s %12s\n','entry','best-match shape','corr','mean|Cm|','mean C_ana');

rows = {};
for i = 1:n
    for j = i:n
        v = squeeze(S.Cm_all(i,j,:));
        if all(abs(v) < 1e-18) || std(v) == 0, continue; end
        cc = zeros(nS,1);
        for k = 1:nS
            cc(k) = corr(abs(v), Cana(:,k));
        end
        [mx, kb] = max(cc);
        if mx > 0.90
            rows(end+1,:) = {sprintf('(%d,%d)',i,j), shapes(kb).name, mx, mean(abs(v)), mean(Cana(:,kb))}; %#ok<SAGROW>
        end
    end
end
for r = 1:size(rows,1)
    fprintf('%-10s %-22s %9.4f %12.4e %12.4e\n', rows{r,:});
end

%% Which terminal is the resonator (proof mass)?
% The proof mass couples to EVERY electrode, so its row has the most
% large off-diagonal entries.
fprintf('\nOff-diagonal magnitude per terminal (the proof mass should dominate):\n');
for i = 1:n
    v = squeeze(sum(abs(S.Cm_all(i,:,:)),2)) - squeeze(abs(S.Cm_all(i,i,:)));
    fprintf('  term %2d: mean off-diag sum = %.4e\n', i, mean(v));
end
