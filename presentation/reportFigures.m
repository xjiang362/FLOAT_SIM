%% reportFigures.m
% Figures for the capacitance / force-torque report, generated from the
% CURRENT parameter-free model so nothing predates it.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'estimator'), fullfile(root,'presentation'));
p = trapDeviceGeometry();
L = load(fullfile(root,'presentation','reportNumbers.mat')); R = L.R;

D = load(fullfile(root,'presentation','capMatrixSweep63_parsed.mat'));
Cm = D.Cm; nP = D.nPose; RES = 21;
bad=false(nP,1); for q=1:nP, M=Cm(:,:,q); o=M-diag(diag(M)); bad(q)=any(o(:)<0); end
gd=find(~bad);

%% ===== Figure 1: capacitance model =====
fig = figure('Color','w','Position',[40 40 1300 430]);

% (a) model vs COMSOL, all terminals
subplot(1,3,1); hold on; grid on; box on;
mdl=[]; mea=[]; pa=[];
for i=1:20
    for ii=1:numel(gd)
        q=gd(ii);
        pose=[D.rx(q);D.ry(q);D.rz(q);D.thx(q)*pi/180;D.thy(q)*pi/180;D.thz(q)*pi/180];
        C=trapCapacitance(p,pose);
        mdl(end+1)=C(i); mea(end+1)=Cm(i,RES,q); pa(end+1)=p.term(i).P/p.term(i).A*1e-3; %#ok<AGROW>
    end
end
scatter(mdl*1e15, mea*1e15, 5, pa, 'filled','MarkerFaceAlpha',0.45);
lims=[200 7000]; plot(lims,lims,'k--','LineWidth',1.2);
set(gca,'XScale','log','YScale','log'); xlim(lims); ylim(lims);
cb=colorbar; cb.Label.String='P/A [1/mm]';
xlabel('model  \epsilon_0A/g  [fF]'); ylabel('COMSOL C_{iR} [fF]');
title(sprintf('20 terminals \\times %d poses',numel(gd)));

% (b) deficit vs perimeter/area
subplot(1,3,2); hold on; grid on; box on;
pav = arrayfun(@(i) p.term(i).P/p.term(i).A*1e-3, 1:20);
defv = zeros(20,1);
for i=1:20
    ee=zeros(numel(gd),1);
    for ii=1:numel(gd)
        q=gd(ii); pose=[D.rx(q);D.ry(q);D.rz(q);D.thx(q)*pi/180;D.thy(q)*pi/180;D.thz(q)*pi/180];
        C=trapCapacitance(p,pose); ee(ii)=100*(C(i)-Cm(i,RES,q))/Cm(i,RES,q);
    end
    defv(i)=median(ee);
end
scatter(pav, defv, 55, 1:20, 'filled');
pf=polyfit(pav,defv',1); xg=linspace(min(pav),max(pav),20);
plot(xg, polyval(pf,xg),'k-','LineWidth',1.4);
yline(0,'k:');
xlabel('P_i/A_i  [1/mm]'); ylabel('median deficit [%]');
title({'deficit tracks perimeter/area','(the un-modelled fringing)'});

% (c) stiffness
subplot(1,3,3); hold on; grid on; box on;
lbl={'r_x','r_y','r_z','\theta_x','\theta_y','\theta_z'};
fr = R.freq/2/pi; st = R.stable;
bh = bar(1:6, fr, 0.6, 'FaceColor','flat');
for i=1:6
    if st(i), bh.CData(i,:) = [0.20 0.55 0.35];
    else,     bh.CData(i,:) = [0.80 0.25 0.20]; end
end
set(gca,'YScale','log','XTick',1:6,'XTickLabel',lbl);
ylabel('rate [Hz]');
title({'green = passively trapped','red = open-loop unstable'});
ylim([1 400]);

sgtitle('Closed-form capacitance model and the stiffness it implies');
for a=findobj(fig,'Type','axes')', if isprop(a,'Toolbar'), a.Toolbar.Visible='off'; end; end
f1 = fullfile(root,'presentation','report_fig1_capacitance.png');
exportgraphics(fig,f1,'Resolution',140);
fprintf('Saved %s\n', f1);

%% ===== Figure 2: wrench + delay =====
fig2 = figure('Color','w','Position',[40 40 1300 430]);

% (a) lateral force: one-sided per quadrant, two-sided for the set
subplot(1,3,1); hold on; grid on; box on;
rxv = linspace(-45,45,181)'*1e-6;
u1 = R.uHover*[1 0 0 0 0 0]';          % quadrant 1 only
uAll = R.uHover*[1 1 1 1 0 0]';        % all four
Fx1 = zeros(size(rxv)); FxA = zeros(size(rxv));
for i=1:numel(rxv)
    pose=zeros(6,1); pose(1)=rxv(i);
    [~,Jl]=trapCapacitance(p,pose);
    f1=trapWrench(p,Jl,u1);  Fx1(i)=f1(1);
    fa=trapWrench(p,Jl,uAll); FxA(i)=fa(1);
end
plot(rxv*1e6, Fx1*1e9,'-','LineWidth',1.6,'Color',[0.85 0.33 0.10], ...
    'DisplayName','quadrant 1 alone');
plot(rxv*1e6, FxA*1e9,'-','LineWidth',2.0,'Color',[0.12 0.25 0.60], ...
    'DisplayName','all four quadrants');
yline(0,'k:','HandleVisibility','off'); xline(0,'k:','HandleVisibility','off');
xlabel('r_x [\mum]'); ylabel('F_x [nN]');
title({'lateral force from overlap truncation','one-sided alone, restoring in the set'});
legend('Location','northeast','FontSize',8);

% (b) singular values
subplot(1,3,2); hold on; grid on; box on;
sv = R.svdB; sv(sv<1e-30)=1e-30;
bar(1:6, sv, 0.6, 'FaceColor',[0.2 0.3 0.75]);
set(gca,'YScale','log'); xlabel('index'); ylabel('singular value');
title({'B is effectively rank 3','three directions are unreachable'});
ylim([1e-28 1e-6]);

% (c) delay budget
subplot(1,3,3); hold on; box on;
dtv=(0:10:200)*1e-6; i2c=(0:40:800)*1e-6;
[DT,IC]=meshgrid(dtv*1e6,i2c*1e6);
PT = R.poleHeave*(IC*1e-6 + DT*1e-6*3);
contourf(DT,IC,PT,0:0.1:2.0,'LineColor','none');
[cc,hh]=contour(DT,IC,PT,[0.5 1.0],'w-','LineWidth',2);
clabel(cc,hh,'Color','w','FontWeight','bold');
plot(p.readout.dtChannel*1e6, p.readout.tI2C*1e6,'rp','MarkerSize',16,'MarkerFaceColor','r');
colorbar; xlabel('dtChannel [\mus]'); ylabel('tI2C [\mus]');
title({sprintf('delay budget p\\tau (now %.2f)',R.pTau),'red star = current setting'});

sgtitle('Actuation authority and the sensing delay budget');
for a=findobj(fig2,'Type','axes')', if isprop(a,'Toolbar'), a.Toolbar.Visible='off'; end; end
f2 = fullfile(root,'presentation','report_fig2_wrench.png');
exportgraphics(fig2,f2,'Resolution',140);
fprintf('Saved %s\n', f2);
fprintf('deficit-vs-P/A slope: %.3f %%/(1/mm), intercept %.3f%%\n', pf(1), pf(2));
