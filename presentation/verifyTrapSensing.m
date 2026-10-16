%% verifyTrapSensing.m
% Validates estimator/trapSensing.m against the exact 22-terminal Schur
% reduction on the 63-pose dedicated capacitance-matrix sweep.
%
% The measured per-terminal C_iR = Cm(i,21) is fed in, so this isolates
% the SENSING stage: any error here belongs to the Y-delta reduction and
% the lumped constants, not to the capacitance model upstream.
%
% Poses whose extracted matrix has negative off-diagonals (mesh error in
% the thin gap, |rz| >= 11.2 um) are excluded throughout.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
root = fileparts(scriptDir);
addpath(scriptDir, fullfile(root,'config'), fullfile(root,'estimator'));
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

plant = trapDeviceGeometry();
D = load(fullfile(scriptDir,'capMatrixSweep63_parsed.mat'));
Cm = D.Cm; nPose = D.nPose; blk = D.blk; RES = plant.idx.resonator;

bad = false(nPose,1);
for p=1:nPose, M=Cm(:,:,p); o=M-diag(diag(M)); bad(p)=any(o(:)<0); end
good = ~bad;

%% ---- run the sensing stage at every pose, quadrant 1 pair ----
Cmodel = zeros(nPose,1); Cexact = zeros(nPose,1);
for p = 1:nPose
    CiR = squeeze(Cm(1:20,RES,p));          % as the capacitance stage would supply
    Ce  = trapSensing(plant, CiR);
    Cmodel(p) = Ce(1);                      % quadrant 1
    Cexact(p) = effectiveMutualCapacitance(Cm(:,:,p),1,2,RES);
end
err = 100*(Cmodel-Cexact)./Cexact;

fprintf('=== trapSensing vs exact Schur reduction (quadrant 1 pair) ===\n');
fprintf('  %-6s %6s %12s %12s %12s %10s\n','axis','n','median','max |err|','RMS','R^2');
nm = {'rx','rz','thx'};
for i=1:3
    b = blk{i}; bg = b(good(b));
    fprintf('  %-6s %6d %11.3f%% %11.3f%% %11.3f%% %10.5f\n', nm{i}, numel(bg), ...
        median(err(bg)), max(abs(err(bg))), rms(err(bg)), r2fun(Cexact(bg),Cmodel(bg)));
end
bg = good;
fprintf('  %-6s %6d %11.3f%% %11.3f%% %11.3f%% %10.5f\n','ALL', sum(bg), ...
    median(err(bg)), max(abs(err(bg))), rms(err(bg)), r2fun(Cexact(bg),Cmodel(bg)));

k0 = 11;
fprintf('\n  at equilibrium: model %.4f fF vs exact %.4f fF (%.4f%%)\n', ...
    Cmodel(k0)*1e15, Cexact(k0)*1e15, err(k0));

%% ---- all four quadrants should agree at a symmetric pose ----
CiR0 = squeeze(Cm(1:20,RES,k0));
[Ce0, info] = trapSensing(plant, CiR0);
fprintf('\n=== all four quadrant pairs at equilibrium ===\n');
for q = 1:4
    fprintf('  quadrant %d (terms %2d,%2d): Ceff = %8.4f fF  (series %8.4f fF)\n', ...
        q, info.pair(q,1), info.pair(q,2), Ce0(q)*1e15, info.series(q)*1e15);
end
fprintf('  spread across quadrants: %.3f%%  (symmetry check)\n', ...
    100*(max(Ce0)-min(Ce0))/mean(Ce0));

%% ---- how much of Ceff is the through-ground path? ----
fprintf('\n=== composition at equilibrium ===\n');
fprintf('  C12 %7.3f fF + Cf %6.3f fF + series %8.3f fF = %8.3f fF\n', ...
    plant.ceff.C12*1e15, plant.ceff.Cf*1e15, info.series(1)*1e15, Ce0(1)*1e15);
fprintf('  series term is %.1f%% of Ceff\n', 100*info.series(1)/Ce0(1));

%% ---- plot ----
order=[2 1 3]; axVar={D.rx,D.rz,D.thx}; axName={'r_x','r_z','\theta_x'};
axUnit={'\mum','\mum','deg'}; axScale=[1e6 1e6 1];
figure('Color','w','Position',[40 40 1380 700]);
for c=1:3
    i=order(c); b=blk{i};
    x=axVar{i}(b)*axScale(i); [x,o]=sort(x); bo=b(o); bd=bad(bo);
    subplot(2,3,c); hold on; grid on; box on;
    plot(x,Cexact(bo)*1e15,'ko','MarkerSize',5,'DisplayName','exact (Schur)');
    plot(x,Cmodel(bo)*1e15,'-','LineWidth',1.8,'Color',[0.12 0.25 0.60],'DisplayName','trapSensing');
    if any(bd), plot(x(bd),Cexact(bo(bd))*1e15,'rx','MarkerSize',11,'LineWidth',1.8, ...
            'DisplayName','excluded'); end
    if i==2, set(gca,'YScale','log'); end
    xlabel(sprintf('%s [%s]',axName{i},axUnit{i})); ylabel('C_{eff} [fF]');
    title(sprintf('C_{eff} vs %s',axName{i})); legend('Location','best','FontSize',7.5);

    subplot(2,3,3+c); hold on; grid on; box on;
    gd2=~bd;
    plot(x(gd2),err(bo(gd2)),'o-','MarkerSize',4,'LineWidth',1.4,'Color',[0.12 0.25 0.60]);
    yline(0,'k-'); xline(0,'k:');
    xlabel(sprintf('%s [%s]',axName{i},axUnit{i})); ylabel('error [%]');
    title('relative error'); ylim([-7 7]);
end
sgtitle('estimator/trapSensing.m: sensed capacitance vs exact reduction');
for a=findobj(gcf,'Type','axes')', a.Toolbar.Visible='off'; end
outPng = fullfile(scriptDir,'TrapSensing_validation.png');
exportgraphics(gcf,outPng,'Resolution',150);
fprintf('\nSaved plot to %s\n', outPng);
