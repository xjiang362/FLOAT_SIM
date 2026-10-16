%% compareFringingVsEdgeSmear.m
% Two different things have been called "fringing" in this analysis:
%
%  (1) FRINGING CAPACITANCE EXCESS. Field lines bulge outside the
%      electrode perimeter, so C exceeds the parallel-plate eps0*A/g at
%      EVERY pose. Classic result: the excess is proportional to the
%      electrode PERIMETER and depends on the gap only LOGARITHMICALLY,
%      so it should look like a near-constant offset, NOT like a 1/g term.
%
%  (2) EDGE SMEAR w. The rounding of the plate-electrode overlap
%      transition as the plate edge slides past the electrode edge --
%      a property of the DERIVATIVE near the crossing point, fitted as
%      w = 25.1 um at g = 15 um from the rx sweep.
%
% Both come from the field not stopping at a conductor boundary, so they
% should share one length scale of order the gap. This script measures (1)
% directly from the gap dependence and compares its length scale to w.
%
% Test: dC(g) = C_measured - eps0*A_geom/g.
%   if fringing   -> dC roughly constant (weak log dependence in g)
%   if effective-area inflation -> dC proportional to 1/g
%
% SYMBOLS: eps0 [F/m]; A [m^2] geometric footprint; P [m] footprint
%   perimeter (hole edges included); g [m] gap; dC [F] excess over
%   parallel plate; lambda [m] implied effective edge extension; w [m]
%   fitted edge smear from the rx sweep.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
eps0 = 8.8541878128e-12; V = 30; g0 = 15e-6; a5 = 0.976; um = 1e-6;
geom = trapGeometry();
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

%% ---------------- geometry: area and perimeter ----------------
% every shape here is a rectangle; a subtracted hole ADDS edge length
perim = @(rows) sum(arrayfun(@(k) 2*((rows(k,6)-rows(k,4)) + (rows(k,3)-rows(k,2))), 1:size(rows,1)));
area  = @(rows) sum(arrayfun(@(k) rows(k,1)*(rows(k,6)-rows(k,4))*(rows(k,3)-rows(k,2)), 1:size(rows,1)));

E(1).name='central A';       E(1).rows=[geom.centralA.px; geom.centralA.mx];
E(2).name='quadrant 1 pair'; E(2).rows=[geom.quad(1).a; geom.quad(1).b];
for i=1:2
    E(i).A = area(E(i).rows); E(i).P = perim(E(i).rows);
    fprintf('%-16s A = %.4f mm^2, P = %.1f mm, P/A = %.4f /mm\n', ...
        E(i).name, E(i).A*1e6, E(i).P*1e3, E(i).P/E(i).A*1e-3);
end

%% ---------------- C vs gap, both electrodes ----------------
% central A: 150-pt file, rz block (zero tilt)
C150 = readmatrix(fullfile(scriptDir,'comsol_centerElectrodeSweep_rzThxThy150_results.txt'),'NumHeaderLines',5);
m = all(abs(C150(:,[1 2 4 5 6]))==0,2);
E(1).g = g0 - a5*C150(m,3);  E(1).C = 2*C150(m,13)/V^2;

% quadrant: full-6D sweep, near-zero tilt. NOTE this is an LHS sweep, so
% rx/ry/thz still vary freely and add scatter to dC -- the central-A
% single-axis sweep is the clean test and the quadrant is shown only for
% magnitude comparison.
F = readmatrix(fullfile(scriptDir,'comsol_overnightFull6DSweep_results.txt'),'NumHeaderLines',5);
m = abs(F(:,4))<0.004 & abs(F(:,5))<0.004;
E(2).g = g0 - a5*F(m,3);  E(2).C = 2*F(m,16)/V^2;

%% ---------------- is the excess constant, or 1/g ? ----------------
fprintf('\n============ excess over parallel plate:  dC = C - eps0*A/g ============\n');
for i=1:2
    g = E(i).g; C = E(i).C;
    dC = C - eps0*E(i).A./g;
    ok = g > 2*um;                         % drop the contact-adjacent end
    g = g(ok); dC = dC(ok);

    % candidate forms, each fitted by least squares
    m1 = [ones(size(g))];              p1 = m1\dC;   % constant
    m2 = [1./g];                       p2 = m2\dC;   % pure 1/g
    m3 = [ones(size(g)), log(g)];      p3 = m3\dC;   % a + b ln g  (classic fringing)
    fprintf('\n--- %s  (n = %d, gap %.1f to %.1f um) ---\n', E(i).name, numel(g), min(g)*1e6, max(g)*1e6);
    fprintf('  dC range: %.4f to %.4f pF\n', min(dC)*1e12, max(dC)*1e12);
    fprintf('  %-22s R^2 = %+.4f\n','constant',      r2fun(dC, m1*p1));
    fprintf('  %-22s R^2 = %+.4f\n','proportional 1/g', r2fun(dC, m2*p2));
    fprintf('  %-22s R^2 = %+.4f\n','a + b*ln(g)',   r2fun(dC, m3*p3));

    % classic fringing predicts the ln(g) slope to be ~ eps0*P/pi
    fprintf('  ln(g) coefficient: fitted %.4f pF, predicted eps0*P/pi = %.4f pF, ratio %.2f\n', ...
        p3(2)*1e12, eps0*E(i).P/pi*1e12, p3(2)/(eps0*E(i).P/pi));
    fprintf('  dC at g=15um: %.3f pF (%.1f%% of C);  at g=3um: %.3f pF (%.1f%% of C)\n', ...
        (p3(1)+p3(2)*log(15e-6))*1e12, ...
        100*(p3(1)+p3(2)*log(15e-6))/(eps0*E(i).A/15e-6 + p3(1)+p3(2)*log(15e-6)), ...
        (p3(1)+p3(2)*log(3e-6))*1e12, ...
        100*(p3(1)+p3(2)*log(3e-6))/(eps0*E(i).A/3e-6 + p3(1)+p3(2)*log(3e-6)));
    E(i).gfit=g; E(i).dC=dC; E(i).p3=p3;
end

%% ---------------- plot ----------------
w_fit = 25.1*um;      % edge smear fitted on the rz=0 rx sweep (fitFxEdgeSmear.m)
figure('Color','w','Position',[50 50 1340 430]);

subplot(1,3,1); hold on; grid on;
cc = [0.85 0.33 0.10; 0.20 0.30 0.75];
for i=1:2
    plot(E(i).gfit*1e6, E(i).dC*1e12,'o','MarkerSize',4,'Color',cc(i,:),'DisplayName',E(i).name);
    gg = linspace(min(E(i).gfit),max(E(i).gfit),200)';
    plot(gg*1e6, (E(i).p3(1)+E(i).p3(2)*log(gg))*1e12,'-','LineWidth',1.6, ...
        'Color',cc(i,:),'HandleVisibility','off');
end
set(gca,'XScale','log'); xlabel('gap [\mum]'); ylabel('\DeltaC = C - \epsilon_0A/g  [pF]');
title({'excess over parallel plate','lines = a + b ln(g) fit'}); legend('Location','best','FontSize',8);

subplot(1,3,2); hold on; grid on;
i=1; g=E(i).gfit;
Cpp = eps0*E(i).A./g;
plot(g*1e6, 100*E(i).dC./(Cpp+E(i).dC),'o','MarkerSize',4,'Color',cc(1,:));
xlabel('gap g [\mum]'); ylabel('fringing share of C  [%]');
title({'central A: fringing as a fraction of C','(big at large gap, small at small gap)'});

subplot(1,3,3); hold on; grid on;
bar([1 2],[w_fit/g0, 1]); set(gca,'XTick',[1 2], ...
    'XTickLabel',{sprintf('edge smear\nw/g = %.2f',w_fit/g0),'fringing scale\n~ g'});
ylabel('length scale / gap'); ylim([0 2]);
title({'both are O(g)','same physics, different observables'});

sgtitle('Fringing capacitance excess vs the overlap edge smear');
outPng = fullfile(scriptDir,'FringingVsEdgeSmear.png');
exportgraphics(gcf,outPng);
fprintf('\nSaved plot to %s\n', outPng);
