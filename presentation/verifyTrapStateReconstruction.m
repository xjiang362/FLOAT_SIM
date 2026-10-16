%% verifyTrapStateReconstruction.m
% OPEN-LOOP accuracy of the sensing chain. No controller, no estimator.
%
% The voltages are FROZEN, the plate is released, and the continuous-time
% plant is integrated. At every ISR tick the four effective capacitances
% are sampled through trapSamplerStep, inverted to four heights by
% trapCapToHeight, and combined into a pose by the plane relations
%
%   rz   = (xs1 + xs2 + xs3 + xs4)/4
%   thx  = (xs1 + xs2 - xs3 - xs4)/(4*l)
%   thy  = (xs2 + xs3 - xs1 - xs4)/(4*l),      l = 2250 um
%
% which are exactly pinv(Cs) for the Cs used by the controller. The
% reconstruction is then compared against the true state from the dynamics.
%
% Separating the two error sources matters, so both are reported:
%   NO DELAY  reconstruction from the capacitances at the current instant.
%             Whatever error remains is pure MODEL error in the
%             capacitance -> height inversion.
%   DELAYED   reconstruction from the sampler output, i.e. four channels of
%             different ages (tI2C + dtChannel*(n-k+0.5)). The extra error
%             is pure LAG, and on a moving plate the channel SKEW also
%             fakes a tilt that is not there.
%
% READ THE "NO DELAY" NUMBERS WITH CARE. When the plate is laterally
% centred, trapCapToHeight is the EXACT algebraic inverse of the forward
% model in trapCapacitance/trapSensing: same gap expression, same C0, same
% area. Its error is then zero by construction, and reporting that as an
% accuracy result would be circular. The inversion only becomes
% approximate when the plate moves LATERALLY, because the overlap
% truncation changes the area while the inversion assumes it fixed. Case D
% is therefore the one that measures real model error; in A-C the
% informative quantity is the LAG.
%
% CASES (each released at the voltage that balances the weight AT its own
% start pose, except A which is the 40 V asked for)
%   A  40 V, from rz = -14 um, Tisr = 50 us. The gap to the electrodes is
%      29 um there, so 40 V gives only 195 uN against a 283 uN weight: the
%      plate FALLS, from 1 um off the floor. It lands in 0.75 ms, which is
%      LESS than the 544 us sensing delay -- so the delayed reconstruction
%      never catches up at all. Run at 50 us because at 720 us the whole
%      flight is one tick.
%   B  48.20 V from rz = -14 um, Tisr = 720 us. An UNSTABLE equilibrium, so
%      the plate drifts off it slowly and sweeps a range of heights. This is
%      the case that exercises the inversion over many samples.
%   C  as B but released at rz = -8 um with a 100 urad tilt, to exercise
%      the two tilt channels and expose the sampler SKEW. Released higher
%      than B because at -14 um a 100 urad tilt leaves the corner only
%      0.4 um of clearance and the run ends almost at once.
%   D  as C but with a 20 um lateral offset instead of tilt. Lateral motion
%      is not observable in these four channels, but it LEAKS into the
%      height through the overlap term, so this is the genuine model error.
%
% See also TRAPCAPTOHEIGHT, TRAPSAMPLERSTEP, TRAPSENSING.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p0 = trapDeviceGeometry();
W  = p0.mass.m * p0.phys.g;

% balance voltage at an arbitrary pose, derived from the geometry
vbalAt = @(q) localVbal(p0, q, W);

fprintf('sensing delay: tI2C %.0f us + channel spacing -> oldest channel %.0f us\n', ...
    p0.readout.tI2C*1e6, (p0.readout.tI2C + p0.readout.dtChannel*3.5)*1e6);
fprintf('balance voltage at rz = -14 um : %.2f V\n', vbalAt([0;0;-14e-6;0;0;0]));
fprintf('balance voltage at rz =  -8 um : %.2f V\n\n', vbalAt([0;0;-8e-6;0;0;0]));

%          name                          T        rz0      thx0    rxy0
cases = { ...
  'A: 40 V, falls to floor',           50e-6,  -14e-6,       0,      0,  40; ...
  'B: balance V, slow drift',         720e-6,  -14e-6,       0,      0, NaN; ...
  'C: balance V + 100 urad tilt',     720e-6,   -8e-6,  100e-6,      0, NaN; ...
  'D: balance V + 20 um lateral',     720e-6,   -8e-6,       0,  20e-6, NaN};

res = struct([]);
for ic = 1:size(cases,1)
    T    = cases{ic,2};
    rz0  = cases{ic,3};
    thx0 = cases{ic,4};
    rxy0 = cases{ic,5};
    Vq   = cases{ic,6};
    if isnan(Vq)
        q0 = zeros(6,1); q0(3)=rz0; q0(4)=thx0; q0(1)=rxy0; q0(2)=rxy0;
        Vq = vbalAt(q0);
    end

    p = p0;
    p.readout.Tisr = T;  p.readout.dtChannel = T/4;  p.sim.dtEuler = T/72;
    dt = p.sim.dtEuler;  nSub = round(T/dt);
    ell = abs((p.term(p.idx.sensePair(1,1)).Yc + p.term(p.idx.sensePair(1,2)).Yc)/2);

    u = zeros(6,1);  u(1:5) = Vq;      % quadrants + Cp held, TIA grounded
    s = trapSamplerInit(p, T);

    x = zeros(12,1);  x(3) = rz0;  x(4) = thx0;
    x(1) = rxy0;  x(2) = rxy0;

    nT = 4000;
    rec = nan(nT,3); recND = nan(nT,3); tru = nan(nT,3);
    tt = nan(nT,1);  xsLog = nan(nT,4); rxLog = nan(nT,1); nk = 0;

    for kTick = 1:nT
        t = (kTick-1)*T;
        if ~checkPoseValid(p, x(1:6)), break; end
        try
            [s, CeD, sInfo] = trapSamplerStep(s, x(1:6), t);
            [~, p3nd] = trapCapToHeight(p, sInfo.CeffNow);
            if sInfo.valid
                [xs, p3d] = trapCapToHeight(p, CeD);
            else
                xs = nan(4,1); p3d = nan(3,1);
            end
        catch ME
            if contains(ME.identifier,'trapCapToHeight') || ...
               contains(ME.identifier,'trapCapacitance')
                break
            else
                rethrow(ME)
            end
        end
        nk = nk + 1;  tt(nk) = t;
        tru(nk,:)   = [x(3) x(4) x(5)];
        rxLog(nk,1) = x(1);
        recND(nk,:) = p3nd';
        rec(nk,:)   = p3d';
        xsLog(nk,:) = xs';

        for j = 1:nSub
            if ~checkPoseValid(p, x(1:6)), break; end
            x = x + dt*trapDynamics(p, x, u);
        end
    end

    res(ic).name = cases{ic,1};  res(ic).T = T;  res(ic).V = Vq;
    res(ic).t = tt(1:nk);  res(ic).tru = tru(1:nk,:);
    res(ic).rec = rec(1:nk,:);  res(ic).recND = recND(1:nk,:);
    res(ic).ell = ell;  res(ic).nValid = sum(~isnan(rec(1:nk,1)));

    rms_ = @(v) sqrt(mean(v(~isnan(v)).^2));
    mx_  = @(v) max(abs(v(~isnan(v))));
    eND = recND(1:nk,1) - tru(1:nk,1);
    eD  = rec(1:nk,1)   - tru(1:nk,1);

    fprintf('%-30s V = %.2f, Tisr %3.0f us\n', cases{ic,1}, Vq, T*1e6);
    fprintf('    %3d ticks (%3d with a valid delayed set), %6.2f ms, rz %+6.2f -> %+6.2f um\n', ...
        nk, res(ic).nValid, tt(nk)*1e3, tru(1,1)*1e6, tru(nk,1)*1e6);
    fprintf('    rz err, no delay (MODEL only) : RMS %8.2f nm   max %8.2f nm\n', ...
        1e9*rms_(eND), 1e9*mx_(eND));
    if res(ic).nValid > 0
        fprintf('    rz err, delayed  (model+LAG)  : RMS %8.2f nm   max %8.2f nm\n', ...
            1e9*rms_(eD), 1e9*mx_(eD));
    else
        fprintf('    rz err, delayed               : NEVER VALID -- the run is shorter\n');
        fprintf('                                    than the sensing delay itself\n');
    end
    if thx0 ~= 0
        etx = rec(1:nk,2) - tru(1:nk,2);
        ety = rec(1:nk,3) - tru(1:nk,3);
        fprintf('    thx err, delayed              : RMS %8.2f urad  max %8.2f urad\n', ...
            1e6*rms_(etx), 1e6*mx_(etx));
        fprintf('    thy err  (true thy = 0, so this is PURE channel-skew artefact)\n');
        fprintf('                                  : RMS %8.2f urad  max %8.2f urad\n', ...
            1e6*rms_(ety), 1e6*mx_(ety));
    end
    if rxy0 ~= 0
        fprintf('    lateral: released at rx = ry = %.0f um; rx swings %+.2f to %+.2f um\n', ...
            rxy0*1e6, min(rxLog(1:nk))*1e6, max(rxLog(1:nk))*1e6);
        fprintf('    (lateral is NOT observable in these 4 channels, but it leaks into rz)\n');
    end
    fprintf('\n');
end

%% ---------------- plot ----------------
nC = numel(res);
f = figure('Color','w','Position',[20 20 1420 240*nC]);
cT = [0.15 0.15 0.15]; cN = [0.20 0.50 0.75]; cD = [0.85 0.35 0.15];

for ic = 1:nC
    r = res(ic);
    % --- left: trajectory ---
    subplot(nC,3,3*ic-2); hold on; grid on; box on;
    plot(r.t*1e3, r.tru(:,1)*1e6, '-',  'Color',cT,'LineWidth',2.4,'DisplayName','true r_z (CT dynamics)');
    plot(r.t*1e3, r.recND(:,1)*1e6,'--','Color',cN,'LineWidth',1.6,'DisplayName','reconstructed, no delay');
    plot(r.t*1e3, r.rec(:,1)*1e6, ':', 'Color',cD,'LineWidth',2.0,'DisplayName','reconstructed, delayed');
    xlabel('t [ms]'); ylabel('r_z [\mum]');
    title(sprintf('%s  (%.1f V, T_{ISR} = %.0f \\mus)', r.name, r.V, r.T*1e6),'FontSize',9);
    if ic==1, legend('Location','southwest','FontSize',7); end

    % --- middle: error vs time ---
    subplot(nC,3,3*ic-1); hold on; grid on; box on;
    plot(r.t*1e3, (r.recND(:,1)-r.tru(:,1))*1e9,'-','Color',cN,'LineWidth',1.6,'DisplayName','model error only');
    plot(r.t*1e3, (r.rec(:,1)-r.tru(:,1))*1e9, '-','Color',cD,'LineWidth',1.6,'DisplayName','model + lag');
    yline(0,'k:','HandleVisibility','off');
    xlabel('t [ms]'); ylabel('r_z error [nm]');
    title('reconstruction error','FontSize',9);
    if ic==1, legend('Location','best','FontSize',7); end

    % --- right: tilt channels, or model error vs height ---
    subplot(nC,3,3*ic); hold on; grid on; box on;
    if any(r.tru(:,2) ~= 0) || any(r.tru(:,3) ~= 0)
        plot(r.t*1e3, r.tru(:,2)*1e6,'-','Color',cT,'LineWidth',2.2,'DisplayName','true \theta_x');
        plot(r.t*1e3, r.rec(:,2)*1e6,':','Color',cD,'LineWidth',2.0,'DisplayName','recon \theta_x');
        plot(r.t*1e3, r.tru(:,3)*1e6,'-','Color',[0.4 0.6 0.4],'LineWidth',1.4,'DisplayName','true \theta_y = 0');
        plot(r.t*1e3, r.rec(:,3)*1e6,':','Color',[0.5 0.2 0.6],'LineWidth',1.8,'DisplayName','recon \theta_y (skew)');
        xlabel('t [ms]'); ylabel('tilt [\murad]');
        title('tilt channels','FontSize',9); legend('Location','best','FontSize',6);
    else
        plot(r.tru(:,1)*1e6, (r.recND(:,1)-r.tru(:,1))*1e9,'.-','Color',cN,'LineWidth',1.2,'MarkerSize',8);
        yline(0,'k:');
        xlabel('true r_z [\mum]'); ylabel('error, no delay [nm]');
        title('model error vs height','FontSize',9);
    end
end

sgtitle(['Open-loop state reconstruction: frozen voltages, no controller, no estimator  \bullet  ' ...
         'r_z = (x_{s1}+x_{s2}+x_{s3}+x_{s4})/4'],'FontSize',11);
for ax = findobj(f,'Type','axes')', ax.Toolbar.Visible='off'; end
out = fullfile(root,'presentation','TrapStateReconstruction.png');
exportgraphics(f, out, 'Resolution',140);
fprintf('Saved %s\n', out);

% =====================================================================
function V = localVbal(p, pose, W)
%LOCALVBAL Voltage on quadrants + Cp (TIA grounded) that balances W at pose.
[~, J]   = trapCapacitance(p, pose);
[~,~,B]  = trapWrench(p, J, zeros(6,1));
V = sqrt(W / sum(B(3,1:5)));
end
