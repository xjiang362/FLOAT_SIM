%% refineC1RGapModel.m
% C_eff cannot beat the quality of its inputs. The analytic sensing-to-
% proof-mass capacitance C1R has the right MAGNITUDE (0.72% on the mean,
% with zero fitted parameters) but only R^2=0.90 in its POSE dependence,
% which caps the whole C_eff model. Diagnose that here on C1R alone,
% where the COMSOL ground truth is a single clean matrix entry.
%
% Candidates, in increasing order of freedom:
%   M0  bare geometry, G = g0 - rz                    (0 params)
%   M1  + scale alpha                                 (1 param)
%   M2  + gap offset and scale: G = gz - a5*rz        (3 params)
%   M3  + additive perimeter fringing term            (4 params)
%
% The fringing term is physical: field lines bulge past the electrode
% edge, adding roughly eps0*beta per unit edge length. It is ADDITIVE and
% only weakly pose dependent, which is exactly the signature of a term
% that shifts the mean without tracking the gap.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
g = trapGeometry(); eps0 = 8.8541878128e-12;
r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

T = load(fullfile(scriptDir,'sensingSweep60_allMatrices.mat'));
rz = T.rz; thx = T.thx*pi/180; thy = T.thy*pi/180;
C1R = squeeze(abs(T.Cm_all(1,21,:)));

% Perimeter of the sensing triangle (legs 1415 um + hypotenuse)
Ltri = 1415e-6; P = 2*Ltri + Ltri*sqrt(2);

I = @(gz,a5) shapeCap(g.sense1, gz - a5*rz, thx, thy, 1);

fprintf('================ C1R MODEL LADDER ================\n');
fprintf('%-52s %10s\n','model','R^2');

% M0
m0 = eps0*I(15e-6,1);
fprintf('%-52s %+10.4f\n','M0  bare geometry, G=g0-rz, 0 params', r2fun(C1R,m0));

% M1: scale only
a = I(15e-6,1)\ (C1R/eps0);
fprintf('%-52s %+10.4f   (alpha=%.4f)\n','M1  + scale alpha', r2fun(C1R, eps0*a*I(15e-6,1)), a);

% M2: gap offset + scale
cost2 = @(p) 1 - r2fun(C1R, eps0*p(3)*I(p(1),p(2)));
p2 = fminsearch(cost2,[15e-6,1,1],optimset('Display','off','TolX',1e-14,'TolFun',1e-14,'MaxFunEvals',5e4,'MaxIter',5e4));
fprintf('%-52s %+10.4f   (gz=%.3f um, a5=%.4f, alpha=%.4f)\n','M2  + gap offset & scale', ...
    r2fun(C1R, eps0*p2(3)*I(p2(1),p2(2))), p2(1)*1e6, p2(2), p2(3));

% M3: + additive perimeter fringing
cost3 = @(p) 1 - r2fun(C1R, eps0*(p(3)*I(p(1),p(2)) + p(4)*P));
p3 = fminsearch(cost3,[p2(1),p2(2),p2(3),0.5],optimset('Display','off','TolX',1e-14,'TolFun',1e-14,'MaxFunEvals',5e4,'MaxIter',5e4));
m3 = eps0*(p3(3)*I(p3(1),p3(2)) + p3(4)*P);
fprintf('%-52s %+10.4f   (gz=%.3f um, a5=%.4f, alpha=%.4f, beta=%.3f)\n','M3  + perimeter fringing', ...
    r2fun(C1R,m3), p3(1)*1e6, p3(2), p3(3), p3(4));

fprintf('\nfringing share of C1R at nominal pose: %.1f%%\n', ...
    100*eps0*p3(4)*P / mean(m3));

%% Per-axis breakdown of the best model
blocks = {1:20,21:40,41:60}; bn={'rz-only','thx-only','thy-only'};
fprintf('\nBest model (M3) per axis:\n');
for b=1:3
    fprintf('  %-10s R^2 = %+.4f\n', bn{b}, r2fun(C1R(blocks{b}), m3(blocks{b})));
end

%% Does rx/ry/thz matter? (assumed negligible earlier -- verify)
fprintf('\nResidual correlation with the ignored DOFs:\n');
res = C1R - m3;
fprintf('  corr(residual, rx)  = %+.4f\n', corr(res, T.rx));
fprintf('  corr(residual, ry)  = %+.4f\n', corr(res, T.ry));
fprintf('  corr(residual, thz) = %+.4f\n', corr(res, T.thz));

save(fullfile(scriptDir,'c1rGapModel.mat'),'p3','P','eps0');

figure('Color','w','Position',[60 60 1250 420]);
xs={rz*1e6,thx*180/pi,thy*180/pi}; xl={'rz [\mum]','\theta_x [deg]','\theta_y [deg]'};
for b=1:3
    subplot(1,3,b); hold on; grid on; i=blocks{b};
    plot(xs{b}(i), C1R(i)*1e15,'ko','MarkerSize',6,'DisplayName','COMSOL');
    plot(xs{b}(i), m0(i)*1e15,'b--','LineWidth',1.3,'DisplayName','M0 bare geometry');
    plot(xs{b}(i), m3(i)*1e15,'r-','LineWidth',1.7,'DisplayName','M3 + fringing');
    xlabel(xl{b}); ylabel('C_{1R} [fF]'); title(bn{b}); legend('Location','best','FontSize',7);
end
sgtitle('Sensing-to-proof-mass capacitance: analytic ladder vs COMSOL');
exportgraphics(gcf, fullfile(scriptDir,'C1R_modelLadder.png'));
fprintf('\nSaved plot to C1R_modelLadder.png\n');
