%% reportNumbers.m
% Collects every number quoted in the capacitance / force-torque report,
% straight from the current model, so the report cannot drift from the code.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'estimator'), fullfile(root,'presentation'));
p = trapDeviceGeometry();
R = struct();

%% ---- geometry table ----
fprintf('=== GEOMETRY ===\n');
for i = [1 3 17 19]
    e = p.term(i);
    fprintf('%2d %-9s A=%.4f mm^2 P=%.3f mm Xc=%+8.1f Yc=%+8.1f dX=%7.1f dY=%7.1f\n', ...
        i, e.name, e.A*1e6, e.P*1e3, e.Xc*1e6, e.Yc*1e6, e.dX*1e6, e.dY*1e6);
end
R.mass_mg = p.mass.m*1e6; R.Ixx = p.mass.I(1,1);
fprintf('mass %.4f mg, Ixx %.4e kg m^2, weight %.2f uN\n', ...
    R.mass_mg, R.Ixx, p.mass.m*p.phys.g*1e6);

%% ---- capacitance accuracy, all 20 terminals ----
D = load(fullfile(root,'presentation','capMatrixSweep63_parsed.mat'));
Cm = D.Cm; nP = D.nPose; RES = 21;
bad = false(nP,1);
for q=1:nP, M=Cm(:,:,q); o=M-diag(diag(M)); bad(q)=any(o(:)<0); end
gd = find(~bad);
allE = []; perT = zeros(20,1);
for i=1:20
    ee = zeros(numel(gd),1);
    for ii=1:numel(gd)
        q = gd(ii);
        pose = [D.rx(q);D.ry(q);D.rz(q);D.thx(q)*pi/180;D.thy(q)*pi/180;D.thz(q)*pi/180];
        C = trapCapacitance(p,pose);
        ee(ii) = 100*(C(i)-Cm(i,RES,q))/Cm(i,RES,q);
    end
    perT(i)=median(ee); allE=[allE; ee]; %#ok<AGROW>
end
R.capMed = median(allE); R.capRMS = rms(allE); R.capMax = max(abs(allE));
R.nPoses = numel(gd);
fprintf('\n=== CAPACITANCE, 20 terminals x %d poses ===\n', R.nPoses);
fprintf('median %+.2f%%  RMS %.2f%%  max %.2f%%\n', R.capMed, R.capRMS, R.capMax);
R.defAct = median(perT(p.idx.actuation));
R.defSns = median(perT(p.idx.sensing));
R.defCp  = median(perT([17 18]));
R.defTIA = median(perT([19 20]));
fprintf('deficit: act %+.2f%%, Cp %+.2f%%, sns %+.2f%%, TIA %+.2f%%\n', ...
    R.defAct, R.defCp, R.defSns, R.defTIA);

%% ---- sensing chain ----
C0 = trapCapacitance(p, zeros(6,1));
Ce0 = trapSensing(p, C0);
R.Ceff0 = Ce0(1); R.CeffMeas = 414.5503e-15;
eS=[];
for ii=1:numel(gd)
    q=gd(ii); pose=[D.rx(q);D.ry(q);D.rz(q);D.thx(q)*pi/180;D.thy(q)*pi/180;D.thz(q)*pi/180];
    Ce = trapSensing(p, trapCapacitance(p,pose));
    eS(end+1)=100*(Ce(1)-effectiveMutualCapacitance(Cm(:,:,q),1,2,RES))/effectiveMutualCapacitance(Cm(:,:,q),1,2,RES); %#ok<AGROW>
end
R.sensMed=median(eS); R.sensRMS=rms(eS); R.sensMax=max(abs(eS));
fprintf('\n=== SENSING chain ===\n');
fprintf('Ceff %.3f fF (meas %.3f), chain median %+.2f%% RMS %.2f%% max %.2f%%\n', ...
    R.Ceff0*1e15, R.CeffMeas*1e15, R.sensMed, R.sensRMS, R.sensMax);

%% ---- hover + allocation ----
pose0 = zeros(6,1);
[~,J0] = trapCapacitance(p,pose0); [~,~,B] = trapWrench(p,J0,zeros(6,1));
sel=[1 1 1 1 0 0]'; W=p.mass.m*p.phys.g;
R.vHover = W/(B(3,:)*sel); R.uHover = sqrt(R.vHover);
R.B = B; R.rankB = rank(B); R.svdB = svd(B);
fprintf('\n=== WRENCH ===\n');
fprintf('hover %.2f V, rank(B) = %d\n', R.uHover, R.rankB);
u = R.uHover*sel;

%% ---- stiffness ----
h=[1e-9 1e-9 1e-9 1e-7 1e-7 1e-7]; K=zeros(6,6);
for j=1:6
    qp=pose0; qp(j)=qp(j)+h(j); qm=pose0; qm(j)=qm(j)-h(j);
    [~,Jp]=trapCapacitance(p,qp); [Fp,Tp]=trapWrench(p,Jp,u);
    [~,Jm]=trapCapacitance(p,qm); [Fm,Tm]=trapWrench(p,Jm,u);
    K(:,j) = -([Fp;Tp]-[Fm;Tm])/(2*h(j));
end
R.K = K;
M=[p.mass.m p.mass.m p.mass.m diag(p.mass.I)'];
R.freq = zeros(6,1); R.stable = false(6,1);
for i=1:6
    R.freq(i)=sqrt(abs(K(i,i))/M(i)); R.stable(i)=K(i,i)>0;
end
fprintf('\n=== STIFFNESS ===\n');
lbl={'rx','ry','rz','thx','thy','thz'};
for i=1:6
    kind = {'UNSTABLE','restoring'};
    fprintf('%5s K=%+11.4g  %-9s  %.1f rad/s (%.1f Hz)\n', lbl{i}, K(i,i), ...
        kind{R.stable(i)+1}, R.freq(i), R.freq(i)/2/pi);
end
R.poleHeave = R.freq(3);

%% ---- delay ----
R.age = p.readout.tI2C + p.readout.dtChannel*(4-(1:4)');
R.maxAge = max(R.age); R.skew = max(R.age)-min(R.age);
R.pTau = R.poleHeave*R.maxAge;
R.tauBudget = 1/R.poleHeave;
fprintf('\n=== DELAY ===\n');
fprintf('ages %s us, skew %.0f us, p*tau %.3f, budget %.0f us\n', ...
    sprintf('%.0f ',R.age*1e6), R.skew*1e6, R.pTau, R.tauBudget*1e6);

save(fullfile(root,'presentation','reportNumbers.mat'),'R');
fprintf('\nSaved to reportNumbers.mat\n');
