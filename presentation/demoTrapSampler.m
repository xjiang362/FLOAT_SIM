%% demoTrapSampler.m
% Exercises the discrete acquisition model and quantifies what the
% delay costs, on a plant whose unstable poles sit near 180 Hz.
%
% The plate is driven along a known trajectory (it is NOT integrated here
% -- this isolates the sampler), the sampler is stepped once per ISR, and
% the delivered capacitances are compared with the true ones.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), fullfile(root,'estimator'));
p = trapDeviceGeometry();

Tisr = 1/4000;                       % 4 kHz, the rate the unstable pole needs
s = trapSamplerInit(p, Tisr);   % timing defaults come from p.readout

fprintf('=== acquisition schedule ===\n');
fprintf('  ISR %.0f us, channel spacing %.0f us, I2C %.0f us\n', ...
    Tisr*1e6, s.dtChannel*1e6, s.tI2C*1e6);
for k = 1:s.nCh
    pr = p.idx.sensePair(k,:);
    fprintf('  ch%d (terminals %2d,%2d): age %6.0f us\n', k, pr(1), pr(2), s.age(k)*1e6);
end
fprintf('  skew between first and last channel: %.0f us\n', s.skew*1e6);
fprintf('  buffer: %d ISR slots\n', s.nBuf);

%% ---- drive a known rz oscillation and watch the lag ----
f = 50;                              % Hz, well inside the control band
amp = 1e-6;                          % 1 um
T = 0.02;
N = round(T/Tisr);
tv = (0:N-1)'*Tisr;
Cd = nan(N,s.nCh); Cn = nan(N,s.nCh); vld = false(N,1);
for i = 1:N
    pose = zeros(6,1);
    pose(3) = amp*sin(2*pi*f*tv(i));
    [s, Ce, info] = trapSamplerStep(s, pose, tv(i));
    Cd(i,:) = Ce'; Cn(i,:) = info.CeffNow'; vld(i) = info.valid;
end

g = vld;
fprintf('\n=== lag error on a %.0f Hz, %.0f um heave oscillation ===\n', f, amp*1e6);
fprintf('  signal swing (true): %.3f fF peak-to-peak\n', (max(Cn(g,1))-min(Cn(g,1)))*1e15);
fprintf('%6s %14s %14s %12s\n','ch','RMS lag err','as %% of swing','phase lag');
for k = 1:s.nCh
    e = Cd(g,k)-Cn(g,k);
    sw = max(Cn(g,k))-min(Cn(g,k));
    fprintf('%6d %11.4f fF %13.1f%% %10.1f deg\n', k, rms(e)*1e15, ...
        100*rms(e)/sw, 360*f*s.age(k));
end

%% ---- the control-theoretic cost ----
fprintf('\n=== delay vs the unstable pole ===\n');
pose0 = zeros(6,1);
[~,J0] = trapCapacitance(p,pose0); [~,~,B0] = trapWrench(p,J0,zeros(6,1));
sel = [1 1 1 1 0 0]'; vH = p.mass.m*p.phys.g/(B0(3,:)*sel); u = sqrt(vH)*sel;
h = 1e-9; qp=pose0; qp(3)=h; qm=pose0; qm(3)=-h;
[~,Jp]=trapCapacitance(p,qp); [Fp,~]=trapWrench(p,Jp,u);
[~,Jm]=trapCapacitance(p,qm); [Fm,~]=trapWrench(p,Jm,u);
kz = -(Fp(3)-Fm(3))/(2*h);
pole = sqrt(abs(kz)/p.mass.m);
fprintf('  heave unstable pole  p = %.0f rad/s (%.0f Hz)\n', pole, pole/2/pi);
fprintf('  worst sensing delay  tau = %.0f us\n', s.maxAge*1e6);
fprintf('  p*tau = %.3f\n', pole*s.maxAge);
fprintf('\n  A plant with an unstable pole p and loop delay tau cannot be\n');
fprintf('  stabilised with any reasonable margin once p*tau approaches 1\n');
fprintf('  (the delay eats the phase the loop needs above the pole).\n');
if pole*s.maxAge > 1
    fprintf('  -> p*tau = %.2f EXCEEDS 1: not stabilisable as configured.\n', pole*s.maxAge);
    fprintf('     max tolerable delay at this pole: %.0f us\n', 1e6/pole);
    fprintf('     i.e. the acquisition chain must come down from %.0f us.\n', s.maxAge*1e6);
end

%% ---- sweep the two knobs ----
% plant.readout.dtChannel and .tI2C are the levers. Map p*tau over them so
% the feasible corner is visible rather than guessed at.
fprintf('\n=== delay budget: p*tau over the two knobs ===\n');
dtv = (0:25:200)*1e-6;
i2c = (0:100:800)*1e-6;
PT = zeros(numel(i2c), numel(dtv));
for a = 1:numel(i2c)
    for b = 1:numel(dtv)
        PT(a,b) = pole*(i2c(a) + dtv(b)*(s.nCh-1));
    end
end
fprintf('  rows = tI2C [us], cols = dtChannel [us], entries = p*tau\n');
fprintf('%8s','');  fprintf('%8.0f', dtv*1e6); fprintf('\n');
for a = 1:numel(i2c)
    fprintf('%8.0f', i2c(a)*1e6);
    fprintf('%8.2f', PT(a,:));
    fprintf('\n');
end
fprintf('\n  current setting (%.0f, %.0f) gives p*tau = %.2f\n', ...
    s.tI2C*1e6, s.dtChannel*1e6, pole*s.maxAge);
fprintf('  need p*tau well below 1; p*tau <= 0.5 is a workable target.\n');

figure('Color','w','Position',[80 80 620 460]); hold on; box on;
[DT,IC] = meshgrid(dtv*1e6, i2c*1e6);
contourf(DT, IC, PT, 0:0.1:2.5, 'LineColor','none');
[cc,hh] = contour(DT, IC, PT, [0.5 1.0], 'w-','LineWidth',2);
clabel(cc,hh,'Color','w','FontWeight','bold');
plot(s.dtChannel*1e6, s.tI2C*1e6, 'rp','MarkerSize',18,'MarkerFaceColor','r');
text(s.dtChannel*1e6, s.tI2C*1e6-45, sprintf('  current, p\\tau=%.2f',pole*s.maxAge), ...
    'Color','r','FontWeight','bold');
colorbar; colormap(parula);
xlabel('dtChannel [\mus]'); ylabel('tI2C [\mus]');
title({'Delay budget: p\tau for the heave pole', ...
       'white contours at p\tau = 0.5 and 1.0'});
ax=gca; ax.Toolbar.Visible='off';
out2 = fullfile(root,'presentation','TrapSampler_delayBudget.png');
exportgraphics(gcf,out2,'Resolution',150);
fprintf('  Saved budget map to %s\n', out2);

%% ---- plot ----
figure('Color','w','Position',[60 60 1200 440]);
subplot(1,2,1); hold on; grid on; box on;
plot(tv(g)*1e3, Cn(g,1)*1e15,'k-','LineWidth',1.6,'DisplayName','true C_{12}(t)');
plot(tv(g)*1e3, Cd(g,1)*1e15,'-','LineWidth',1.5,'Color',[0.85 0.33 0.10], ...
    'DisplayName',sprintf('delivered ch1 (%.0f us old)',s.age(1)*1e6));
plot(tv(g)*1e3, Cd(g,4)*1e15,'-','LineWidth',1.5,'Color',[0.12 0.25 0.60], ...
    'DisplayName',sprintf('delivered ch4 (%.0f us old)',s.age(4)*1e6));
xlabel('t [ms]'); ylabel('C_{eff} [fF]');
title(sprintf('%.0f Hz heave, %.0f um amplitude',f,amp*1e6));
legend('Location','best','FontSize',8);

subplot(1,2,2); hold on; grid on; box on;
for k=1:s.nCh
    plot(tv(g)*1e3, (Cd(g,k)-Cn(g,k))*1e15,'-','LineWidth',1.3, ...
        'DisplayName',sprintf('ch%d, %.0f us',k,s.age(k)*1e6));
end
yline(0,'k-','HandleVisibility','off');
xlabel('t [ms]'); ylabel('delivered - true [fF]');
title('lag error, per channel'); legend('Location','best','FontSize',8);
sgtitle('Discrete capacitance acquisition: per-channel delay and skew');
for a=findobj(gcf,'Type','axes')', a.Toolbar.Visible='off'; end
out = fullfile(root,'presentation','TrapSampler_demo.png');
exportgraphics(gcf,out,'Resolution',150);
fprintf('\nSaved plot to %s\n', out);

