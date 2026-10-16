%% forensicCenterElectrodeBlocks.m
% The device has TWO central electrodes (inner + outer), and COMSOL may
% have overwritten one of the exports. So treat every BLOCK of every
% results file as an unknown-provenance dataset and fingerprint it.
%
% Fingerprint: fit the exact centered-electrode force law
%     Fz = k*V^2 / (g^2 - b^2*theta^2)
% which is linear in parameters after inversion:
%     1/Fz = g^2/(k V^2) - b^2*theta^2/(k V^2)
% b is a pure GEOMETRY parameter (the electrode half-width). It cannot
% change between runs of the same electrode, so blocks that disagree on
% b are different electrodes. k is proportional to electrode AREA, so it
% is a second, independent fingerprint.
%
% Blocks fingerprinted:
%   R1-rz   rows   1- 50  (theta=0 -> only k identifiable)
%   R1-thx  rows  51-100
%   R1-thy  rows 101-150
%   R2-thx  rows   1-100
%   R2-thy  rows 101-200

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
V = 30; g0_z = 15e-6; a5 = 0.976;            % a5 from the validated Fz fit
gapOf = @(rz) g0_z - a5*rz;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

raw1 = readmatrix(fullfile(scriptDir,'comsol_centerElectrodeSweep_rzThxThy150_results.txt'),'NumHeaderLines',5);
raw2 = readmatrix(fullfile(scriptDir,'comsol_centerElectrode_rzCrossTerm200_results.txt'),'NumHeaderLines',5);

B = struct('name',{},'rz',{},'thx',{},'thy',{},'Fz',{},'Taux',{},'Tauy',{});
mk = @(nm,R,idx) struct('name',nm,'rz',R(idx,3),'thx',R(idx,4),'thy',R(idx,5), ...
                        'Fz',R(idx,9),'Taux',R(idx,10),'Tauy',R(idx,11));
B(1) = mk('R1-rz   (150file)', raw1,   1:50);
B(2) = mk('R1-thx  (150file)', raw1,  51:100);
B(3) = mk('R1-thy  (150file)', raw1, 101:150);
B(4) = mk('R2-thx  (200file)', raw2,   1:100);
B(5) = mk('R2-thy  (200file)', raw2, 101:200);

fprintf('=========== GEOMETRY FINGERPRINT PER BLOCK ===========\n');
fprintf('%-18s %6s %10s %12s %12s %9s %14s\n', ...
    'block','n','rz vals','b [um]','k','R2(Fz)','Fz@(0,0) [N]');
kAll = nan(5,1); bAll = nan(5,1);
for i = 1:numel(B)
    s = B(i);
    g = gapOf(s.rz); thsq = (s.thx*pi/180).^2 + (s.thy*pi/180).^2;
    nRz = numel(uniquetol(s.rz,1e-9));
    if max(thsq) < 1e-12
        % theta == 0 throughout: only k is identifiable
        kk = median(s.Fz .* g.^2 / V^2);
        bb = NaN;
        rr = r2fun(s.Fz, kk*V^2./g.^2);
    else
        c = [g.^2, -thsq] \ (1./s.Fz);
        kk = 1/(c(1)*V^2);
        bb = sqrt(max(c(2),0)/c(1));
        rr = r2fun(s.Fz, kk*V^2./max(g.^2-bb^2*thsq,1e-18));
    end
    kAll(i) = kk; bAll(i) = bb;
    fprintf('%-18s %6d %10d %12.2f %12.4e %9.4f %14.4e\n', ...
        s.name, numel(s.rz), nRz, bb*1e6, kk, rr, kk*V^2/g0_z^2);
end

fprintf('\nk is proportional to electrode AREA; b is its half-width.\n');
fprintf('Blocks sharing (k,b) are the SAME electrode.\n');

fprintf('\n=========== PAIRWISE k RATIOS ===========\n');
fprintf('%-18s', ''); fprintf('%12s', B.name); fprintf('\n');
for i=1:5
    fprintf('%-18s', B(i).name);
    for j=1:5, fprintf('%12.4f', kAll(i)/kAll(j)); end
    fprintf('\n');
end

fprintf('\n=========== TORQUE SIGN / MAGNITUDE PER BLOCK ===========\n');
% Slope of the PRIMARY torque vs its own tilt, through the origin, split
% by rz so blocks with two rz values report both.
fprintf('%-18s %8s %16s %16s\n','block','rz[um]','dTau/dth [Nm/rad]','sign');
for i = 2:5
    s = B(i);
    isX = contains(s.name,'thx');
    th = (isX .* s.thx + ~isX .* s.thy)*pi/180;
    tau = isX*0; if isX, tau = s.Taux; else, tau = s.Tauy; end
    for rzv = uniquetol(s.rz,1e-9)'
        m = abs(s.rz-rzv) < 1e-9 & abs(th) > 1e-9;
        if sum(m) < 5, continue; end
        sl = th(m) \ tau(m);
        fprintf('%-18s %8.1f %16.5e %16s\n', s.name, rzv*1e6, sl, string(sign(sl)));
    end
end

fprintf('\n=========== CROSS-BLOCK CONSISTENCY WITHIN EACH FILE ===========\n');
fprintf('If a file is internally consistent, its thx and thy blocks must\n');
fprintf('agree on b and k (x<->y symmetry of a centered electrode).\n\n');
pairs = {[2 3],'run 1 (150-pt file): thx vs thy'; [4 5],'run 2 (200-pt file): thx vs thy'};
for p = 1:2
    ii = pairs{p,1};
    fprintf('%-38s b: %8.1f vs %8.1f um  (%6.1f%% apart)   k: %6.3f ratio\n', ...
        pairs{p,2}, bAll(ii(1))*1e6, bAll(ii(2))*1e6, ...
        100*abs(bAll(ii(1))-bAll(ii(2)))/mean(bAll(ii)), kAll(ii(1))/kAll(ii(2)));
end
fprintf('\nAlso: R1-rz''s k vs R1-thx''s k ratio = %.4f (same file, must be 1.00)\n', kAll(1)/kAll(2));
