%% verifyForceToVoltageRoundTrip.m
% Round trip  f -> trapForceToVoltage -> u -> trapWrench -> wrench,
% at a FIXED pose, comparing what the controller COMMANDED against what the
% device actually delivers.
%
% WHAT IS "COMMANDED". trapForceToVoltage takes the regulator increment f
% (5-by-1, order [fp; f1; f2; f3; f4]) and adds the gravity bias, so the
% total per-pair force it is trying to produce is Ftot = Fbias + f. The
% wrench the controller BELIEVES that produces is the design allocation,
%       Fz   = sum(Ftot)
%       taux = a*( Ftot_q1 + Ftot_q2 - Ftot_q3 - Ftot_q4 )
%       tauy = a*(-Ftot_q1 + Ftot_q2 + Ftot_q3 - Ftot_q4 )
% with a = 1500 um, which is K.B rows 2, 4, 6 scaled by the inertias.
%
% WHAT IS DELIVERED. u = sqrt(Ftot./b(pose)) per pair, then the true wrench
% is trapWrench(plant, dCdq(pose), u), i.e.
%       Q = sum_p 0.5*(dC_p1/dq + dC_p2/dq) * u_p^2.
%
% WHY Fz MUST BE EXACT HERE. The inversion divides by the same b that
% trapWrench multiplies by, so the z-force round trip is an identity as
% long as the pose used for the inversion equals the pose used for the
% wrench. Any Fz mismatch is therefore CLIPPING, not model error.
%
% WHY THE TORQUES ARE EXACT AT ZERO TILT. A pair's effective arm is the
% gap-weighted mean of its two members' coordinates,
%       arm_p = [sum_e A_e Y_e/g_e^2] / [sum_e A_e/g_e^2].
% With no tilt every g_e is equal, so the weights are areas alone and the
% arm collapses to the pair CENTROID, +/-1500 um, exactly the design value
% -- at ANY rz. The arm only departs from the centroid once a tilt makes
% the two members see different gaps. So this fixed-pose sweep tests the
% allocation and the clipping, and is deliberately blind to the tilt
% effect; sweep thx/thy to see that.
%
% CLIPPING. At rz = -14 um the electrode gap is 29 um, so b is about
% (15/29)^2 = 0.27 of its hover value and the voltages are correspondingly
% higher. The sweep is run twice: once at the default 80 V clamp and once
% with the clamp lifted, so saturation is separable from everything else.
%
% See also TRAPFORCETOVOLTAGE, TRAPWRENCH, TRAPCAPACITANCE.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p = trapDeviceGeometry();
W  = p.mass.m * p.phys.g;
a  = 1500e-6;
Vtia = 0;

pose = zeros(6,1);  pose(3) = -14e-6;
[~, dCdq] = trapCapacitance(p, pose);
[~,~,Bp]  = trapWrench(p, dCdq, zeros(6,1));
fprintf('pose: rz = %.1f um -> electrode gap %.1f um\n', pose(3)*1e6, (p.gap.z0-pose(3))*1e6);
fprintf('weight W = %.3f uN;  sweep each f_k from 0 to W\n\n', W*1e6);

% bias split, as trapForceToVoltage computes it
map = [5 1 2 3 4];
A5  = arrayfun(@(i) p.term(p.idx.actPair(i,1)).A + p.term(p.idx.actPair(i,2)).A, map)';
Fbias = (W - Bp(3,6)*Vtia^2) * A5/sum(A5);

chName = {'f_p (Cp)','f_1 (quad1)','f_2 (quad2)','f_3 (quad3)','f_4 (quad4)'};
sgn_x  = [0  1  1 -1 -1];      % taux sign per controller channel
sgn_y  = [0 -1  1  1 -1];      % tauy sign per controller channel

nS = 200;
fs = linspace(0, W, nS);

for clampCase = 1:2
    if clampCase == 1, Vmax = 80;   tag = 'Vmax = 80 V (default clamp)';
    else,              Vmax = 1e6;  tag = 'clamp lifted';
    end

    cmd = nan(nS,3,5); act = nan(nS,3,5); Vpk = nan(nS,5); clip = false(nS,5);
    for k = 1:5
        for i = 1:nS
            f = zeros(5,1); f(k) = fs(i);
            [u, info] = trapForceToVoltage(p, f, pose, 'Vtia', Vtia, 'Vmax', Vmax);
            Ftot = Fbias + f;
            cmd(i,:,k) = [sum(Ftot), a*(sgn_x*Ftot), a*(sgn_y*Ftot)];
            [F, tau] = trapWrench(p, dCdq, u);
            act(i,:,k) = [F(3), tau(1), tau(2)];
            Vpk(i,k) = max(u(1:5));  clip(i,k) = info.clipped;
        end
    end

    fprintf('=== %s ===\n', tag);
    fprintf('%-13s %12s %14s %14s %14s\n','channel','clip at f=','max|dFz|[uN]','max|dtaux|[pNm]','max|dtauy|[pNm]');
    for k = 1:5
        ic = find(clip(:,k),1,'first');
        if isempty(ic), cs = 'never'; sel = true(nS,1);
        else, cs = sprintf('%.1f uN', fs(ic)*1e6); sel = (1:nS)' < ic; end
        d = act(sel,:,k) - cmd(sel,:,k);
        fprintf('%-13s %12s %14.3e %14.3e %14.3e\n', chName{k}, cs, ...
            max(abs(d(:,1)))*1e6, max(abs(d(:,2)))*1e12, max(abs(d(:,3)))*1e12);
    end
    fprintf('\n');

    if clampCase == 1, R1.cmd=cmd; R1.act=act; R1.Vpk=Vpk; R1.clip=clip; end
end

%% ---------------- plot (default clamp) ----------------
cmd = R1.cmd; act = R1.act; Vpk = R1.Vpk; clip = R1.clip;
col = lines(5);
ylabs = {'F_z [\muN]','\tau_x [pN\cdotm]','\tau_y [pN\cdotm]'};
scl   = [1e6 1e12 1e12];

f = figure('Color','w','Position',[40 40 1500 520]);
for j = 1:3
    subplot(1,3,j); hold on; grid on; box on;
    for k = 1:5
        plot(fs*1e6, cmd(:,j,k)*scl(j), '--', 'Color',col(k,:), 'LineWidth',1.3, ...
            'HandleVisibility','off');
        plot(fs*1e6, act(:,j,k)*scl(j), '-',  'Color',col(k,:), 'LineWidth',1.9, ...
            'DisplayName',chName{k});
        ic = find(clip(:,k),1,'first');
        if ~isempty(ic)
            plot(fs(ic)*1e6, act(ic,j,k)*scl(j), 'o', 'Color',col(k,:), ...
                'MarkerFaceColor','w','MarkerSize',6,'HandleVisibility','off');
        end
    end
    xlabel('commanded f_k [\muN]'); ylabel(ylabs{j});
    title(ylabs{j}(1:regexp(ylabs{j},' ','once')-1),'FontSize',11);
    if j==1
        legend('Location','northwest','FontSize',8);
        text(5, 1.02*max(max(cmd(:,1,:)))*scl(1), ...
            'dashed = commanded, solid = delivered, o = clip onset', 'FontSize',8);
    end
end
sgtitle(sprintf(['f \\rightarrow trapForceToVoltage \\rightarrow u \\rightarrow trapWrench, ' ...
    'fixed pose r_z = %.0f \\mum, V_{TIA} = %g V'], pose(3)*1e6, Vtia),'FontSize',12);
for ax = findobj(f,'Type','axes')', ax.Toolbar.Visible='off'; end
out = fullfile(root,'presentation','ForceToVoltageRoundTrip.png');
exportgraphics(f, out, 'Resolution',140);
fprintf('Saved %s\n', out);
