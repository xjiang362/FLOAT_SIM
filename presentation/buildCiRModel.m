%% buildCiRModel.m
% Capacitance model C_iR (terminal i <-> resonator) for ALL 20 electrode
% terminals, built purely from the footprint geometry in
% electrodeGeometryTable.tex / makeElectrodeGeometryTable.m:
%
%   C_iR ~= k_i * eps0 * A_i / (g0 - rz - thx*Yc_i + thy*Xc_i)
%
% A_i is the NET footprint area (outer minus any subtracted well) and
% (Xc_i, Yc_i) the NET area centroid, so the hollow central electrodes are
% handled correctly rather than by their outer rectangle alone. k_i is a
% per-electrode gain that should represent FRINGING only -- if the model
% is right, k_i must come out slightly above 1 and should scale with the
% perimeter-to-area ratio P_i/A_i, since fringing is a perimeter effect.
%
% This simultaneously validates the terminal<->geometry assignment in the
% table: a terminal whose assumed footprint is wrong cannot fit, because
% the area sets the magnitude and the centroid sets the tilt slope, and
% neither is free.
%
% Reference data: comsol_capMatrixSweep_rxRzThx63_results.txt (63 poses,
% rx / rz / thx one at a time). Poses whose matrix has negative
% off-diagonals are excluded.
%
% SYMBOLS: eps0 [F/m]; A_i [m^2] net footprint area; P_i [m] net
%   perimeter (hole edges included); Xc_i,Yc_i [m] net area centroid;
%   g0 [m] nominal gap; rz [m]; thx,thy [rad]; k_i [-] fitted gain;
%   C_iR [F] (reported in fF).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
eps0 = 8.8541878128e-12; um = 1e-6; g0 = 15e-6;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

%% ================= geometry for terminals 1..20 =================
% quadrant-1 primitives (micrometres), matching the table
v1 = [1560 1525; 2975 1525; 2975 2940];     % sense, x>y half
v2 = [1525 1560; 1525 2975; 2940 2975];     % sense, x<y half
triG = @(v) struct('A',polyarea(v(:,1),v(:,2)), ...
                   'P',sum(sqrt(sum(diff(v([1:end 1],:)).^2,2))), ...
                   'Xc',mean(v(:,1)),'Yc',mean(v(:,2)));
rectG = @(cx,cy,w,h) struct('A',w*h,'P',2*(w+h),'Xc',cx,'Yc',cy);

S1 = triG(v1); S2 = triG(v2);
A1 = rectG(2250,750,1450,1450);             % act, x>y half
A2 = rectG(750,2250,1450,1450);             % act, x<y half
q1 = {S1,S2,A1,A2};
sx = [1 -1 -1 1]; sy = [1 1 -1 -1];         % quadrant sign maps
% Terminals advance by a 90-degree rotation, so quadrants 2 and 4 take the
% pair in the opposite order (see makeElectrodeGeometryTable.m).
pairOrder = {[1 2 3 4],[2 1 4 3],[1 2 3 4],[2 1 4 3]};

% hollow central pairs: net area, net centroid, perimeter = outer + hole
hollow = @(o,h) struct('A',o.A-h.A,'P',o.P+h.P, ...
                       'Xc',(o.A*o.Xc-h.A*h.Xc)/(o.A-h.A), ...
                       'Yc',(o.A*o.Yc-h.A*h.Yc)/(o.A-h.A));
Cp_o  = rectG((75+1475)/2,0,1400,2950);   Cp_h  = rectG((75+737)/2,0,662,1450);
TIA_o = rectG(0,(25+675)/2,1400,650);     TIA_h = rectG(0,(25+100)/2,200,75);
mirX = @(s) setfield(s,'Xc',-s.Xc);  mirY = @(s) setfield(s,'Yc',-s.Yc); %#ok<SFLD>

E = cell(20,1); nameE = cell(20,1);
for q=1:4
    for j=1:4
        e = q1{ pairOrder{q}(j) }; e.Xc = sx(q)*e.Xc; e.Yc = sy(q)*e.Yc;
        E{4*(q-1)+j} = e;
        tag = {'sns','sns','act','act'}; pm = {'+','-','+','-'};
        nameE{4*(q-1)+j} = sprintf('C_%s%d%s', tag{j}, q, pm{j});
    end
end
% 17 Cp+ is the -x half, 18 Cp- the +x half (verified from the coupling test)
E{17} = mirX(hollow(Cp_o,Cp_h));  nameE{17} = 'Cp+';
E{18} =       hollow(Cp_o,Cp_h);  nameE{18} = 'Cp-';
E{19} =       hollow(TIA_o,TIA_h);nameE{19} = 'C_TIA+';
E{20} = mirY(hollow(TIA_o,TIA_h));nameE{20} = 'C_TIA-';

%% ================= reference data =================
L = load(fullfile(scriptDir,'capMatrixSweep63_parsed.mat'));
Cm = L.Cm; nPose = L.nPose; blk = L.blk; RES = 21;
bad = false(nPose,1);
for p=1:nPose, M=Cm(:,:,p); o=M-diag(diag(M)); bad(p)=any(o(:)<0); end
good = ~bad;
rz = L.rz; thx = L.thx*pi/180; thy = L.thy*pi/180;

%% ================= fit one gain per terminal =================
fprintf('%4s %-9s %9s %9s %9s %9s | %7s %9s %9s\n', ...
    'term','label','A[mm^2]','P[mm]','Xc[um]','Yc[um]','k_i','R^2','P/A[1/mm]');
k = zeros(20,1); R2 = zeros(20,1); PA = zeros(20,1);
for i=1:20
    e = E{i};
    Ai = e.A*um^2; Pi = e.P*um; Xci = e.Xc*um; Yci = e.Yc*um;
    gi = g0 - rz - thx*Yci + thy*Xci;
    Cmod = eps0*Ai./gi;
    meas = squeeze(Cm(i,RES,:));
    k(i)  = Cmod(good)\meas(good);
    R2(i) = r2fun(meas(good), Cmod(good)*k(i));
    PA(i) = Pi/Ai*1e-3;
    fprintf('%4d %-9s %9.4f %9.3f %9.1f %9.1f | %7.4f %9.5f %9.4f\n', ...
        i, nameE{i}, Ai*1e6, Pi*1e3, e.Xc, e.Yc, k(i), R2(i), PA(i));
end

fprintf('\ngain k_i: min %.4f, median %.4f, max %.4f\n', min(k), median(k), max(k));
fprintf('all R^2 >= %.5f\n', min(R2));

%% ================= does the gain track perimeter/area? =================
X = [ones(20,1), PA]; b = X\k;
fprintf('\n=== fringing hypothesis: k_i = 1 + c*(P_i/A_i) ===\n');
fprintf('  fit k = %.4f + %.5f * (P/A)   R^2 = %.4f\n', b(1), b(2), r2fun(k, X*b));
fprintf('  corr(k, P/A) = %.4f\n', corr(k,PA));
fprintf('  intercept should be ~1 if the gain is pure fringing: %.4f\n', b(1));

%% ================= save & plot =================
geomTable = struct('term',num2cell((1:20)'),'name',nameE, ...
    'A',cellfun(@(e) e.A*um^2, E,'UniformOutput',false), ...
    'P',cellfun(@(e) e.P*um,  E,'UniformOutput',false), ...
    'Xc',cellfun(@(e) e.Xc*um,E,'UniformOutput',false), ...
    'Yc',cellfun(@(e) e.Yc*um,E,'UniformOutput',false), ...
    'k',num2cell(k),'R2',num2cell(R2));
save(fullfile(scriptDir,'CiRModel.mat'),'geomTable','k','R2','g0');
fprintf('\nSaved model coefficients to CiRModel.mat\n');

figure('Color','w','Position',[40 40 1380 440]);
subplot(1,3,1); hold on; grid on; box on;
for i=1:20
    e=E{i}; Ai=e.A*um^2; gi=g0-rz-thx*e.Yc*um+thy*e.Xc*um;
    plot(eps0*Ai./gi(good)*k(i)*1e15, squeeze(Cm(i,RES,good))*1e15,'.','MarkerSize',7);
end
lims=[20 8000]; plot(lims,lims,'k--'); set(gca,'XScale','log','YScale','log');
xlabel('model k_i\epsilon_0A_i/g_i [fF]'); ylabel('COMSOL C_{iR} [fF]');
title(sprintf('all 20 terminals x %d poses', sum(good)));

subplot(1,3,2); hold on; grid on; box on;
scatter(PA, k, 46, 1:20, 'filled');
pg = linspace(min(PA),max(PA),50); plot(pg, b(1)+b(2)*pg,'k-','LineWidth',1.5);
yline(1,'k:'); colorbar('Ticks',[1 5 9 13 17],'TickLabels',{'1','5','9','13','17'});
xlabel('P_i/A_i [1/mm]'); ylabel('fitted gain k_i');
title(sprintf('gain vs perimeter/area (R^2=%.3f)', r2fun(k,X*b)));

subplot(1,3,3); hold on; grid on; box on;
bar(1:20, R2); ylim([0.99 1.0005]); xlabel('terminal'); ylabel('R^2');
title('per-terminal fit quality'); xticks(1:2:20);

sgtitle('C_{iR} model from footprint geometry alone, all 20 electrode terminals');
for a=findobj(gcf,'Type','axes')', if isprop(a,'Toolbar'), a.Toolbar.Visible='off'; end; end
outPng = fullfile(scriptDir,'CiRModel_allTerminals.png');
exportgraphics(gcf,outPng,'Resolution',150);
fprintf('Saved plot to %s\n', outPng);
