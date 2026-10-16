%% compareAnalyticCapVsComsol.m
% Compares the LEADING-ORDER ("centroid") analytic capacitance
%
%       C ~= eps0 * (dX*dY) / (g0 - rz - thx*Yc - thy*Xc)
%
% against COMSOL, using the 150-pt central-electrode sweep that carries a
% real Total-Electric-Energy probe. Capacitance is recovered from energy as
% C = 2W/V^2 (valid because every conductor is held at fixed potential).
%
% Also computes the EXACT integral C = eps0 * INT INT dA/g(x,y) via
% capTrapz, so the error of the centroid approximation can be separated
% from the error of the parallel-plate (no-fringing) assumption itself.
%
% SYMBOLS (all SI unless a _deg/_um suffix says otherwise):
%   eps0      [F/m]  vacuum permittivity
%   V         [V]    potential applied to the driven central electrode
%   g0        [m]    plate-to-electrode gap at the nominal pose (rz=0)
%   a5        [-]    calibrated scale mapping rz onto gap change
%   rz        [m]    plate vertical displacement (+rz CLOSES the gap)
%   thx,thy   [rad]  plate tilts; +thx closes the gap at +y,
%                    +thy OPENS the gap at +x (capTrapz convention)
%   dX,dY     [m]    electrode footprint extents (x and y)
%   Xc,Yc     [m]    electrode footprint centroid, plate-centred coords
%   A         [m^2]  electrode footprint area
%   W         [J]    COMSOL total electric energy at the pose
%   C_*       [F]    capacitance (plotted in pF)

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

eps0 = 8.8541878128e-12;
V    = 30;
g0   = 15e-6;
a5   = 0.976;

r2fun = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

%% ---------------- COMSOL data ----------------
f   = fullfile(scriptDir,'comsol_centerElectrodeSweep_rzThxThy150_results.txt');
raw = readmatrix(f,'NumHeaderLines',5);
fprintf('Loaded %s\n  %d rows x %d cols\n', 'comsol_centerElectrodeSweep_rzThxThy150_results.txt', ...
    size(raw,1), size(raw,2));

rz      = raw(:,3);
thx     = raw(:,4)*pi/180;
thy     = raw(:,5)*pi/180;
W       = raw(:,13);
C_comsol = 2*W/V^2;

blk    = {1:50, 51:100, 101:150};
blkNam = {'rz','thx','thy'};

%% ---------------- Geometry of central electrode A ----------------
geom = trapGeometry();
rowsA = [geom.centralA.px; geom.centralA.mx];   % signed trapezoid list

[A, Xc, Yc, dX, dY] = trapMoments(rowsA);
fprintf('\nCentral electrode A footprint (signed sum of %d rectangles):\n', size(rowsA,1));
fprintf('  area A   = %.4f mm^2\n', A*1e6);
fprintf('  centroid = (%.3f, %.3f) um   <-- both ~0 by symmetry\n', Xc*1e6, Yc*1e6);
fprintf('  extents  dX = %.1f um, dY = %.1f um  (bounding box)\n', dX*1e6, dY*1e6);
fprintf('  NOTE: because Xc=Yc=0, the centroid formula predicts NO tilt\n');
fprintf('        dependence at all -- the thx/thy blocks test exactly that.\n');

%% ---------------- Models ----------------
G = g0 - a5*rz;                      % [m] centroid gap vs rz

% (1) user's leading-order centroid formula
C_centroid = eps0 * A ./ (G - thx*Yc - thy*Xc);

% (2) exact geometric integral over the true footprint
I = zeros(numel(rz),1);
for k = 1:size(rowsA,1)
    r = rowsA(k,:);
    I = I + r(1)*capTrapz(r(2),r(3),r(4),r(5),r(6),r(7), G, thx, thy);
end
C_exact = eps0 * I;

%% ---------------- Report ----------------
fprintf('\n================ C = 2W/V^2  vs  analytic ================\n');
fprintf('%6s | %10s %10s %10s | %9s %9s\n', ...
    'block','C_comsol','centroid','exact','ratio','exact/c');
for i = 1:3
    b = blk{i};
    fprintf('%6s | %9.4f  %9.4f  %9.4f | %8.3f  %8.3f\n', blkNam{i}, ...
        mean(C_comsol(b))*1e12, mean(C_centroid(b))*1e12, mean(C_exact(b))*1e12, ...
        mean(C_comsol(b))/mean(C_centroid(b)), mean(C_exact(b))/mean(C_centroid(b)));
end

fprintf('\nR^2 of each analytic model against C_comsol:\n');
fprintf('%6s | %9s %9s | %9s %9s\n','block','centroid','exact','cent+fit','exact+fit');
for i = 1:3
    b = blk{i};
    % with a free gain+offset (tests SHAPE alone, absorbing stray C and fringing gain)
    Mc = [C_centroid(b), ones(numel(b),1)];  pc = Mc\C_comsol(b);
    Me = [C_exact(b),    ones(numel(b),1)];  pe = Me\C_comsol(b);
    fprintf('%6s | %+8.4f  %+8.4f | %+8.4f  %+8.4f\n', blkNam{i}, ...
        r2fun(C_comsol(b), C_centroid(b)), r2fun(C_comsol(b), C_exact(b)), ...
        r2fun(C_comsol(b), Mc*pc), r2fun(C_comsol(b), Me*pe));
end

% How much does C actually move across each block?
fprintf('\nFractional swing of C_comsol within each block (max/min - 1):\n');
for i = 1:3
    b = blk{i};
    fprintf('  %4s : %6.2f %%\n', blkNam{i}, 100*(max(C_comsol(b))/min(C_comsol(b)) - 1));
end

% Global gain+offset on the rz block -> implied effective area / stray C
b = blk{1};
M = [C_centroid(b), ones(numel(b),1)];  p = M\C_comsol(b);
fprintf('\nrz block, C_comsol = gain*C_centroid + offset:\n');
fprintf('  gain   = %.4f   -> implied effective area %.4f mm^2 (geometric %.4f)\n', ...
    p(1), p(1)*A*1e6, A*1e6);
fprintf('  offset = %.4f pF  (pose-independent stray capacitance)\n', p(2)*1e12);

%% ---------------- Plot ----------------
figure('Color','w','Position',[40 40 1400 430]);
for i = 1:3
    b = blk{i};
    switch i
        case 1, xv = rz(b)*1e6;        xl = 'r_z [\mum]';
        case 2, xv = thx(b)*180/pi;    xl = '\theta_x [deg]';
        case 3, xv = thy(b)*180/pi;    xl = '\theta_y [deg]';
    end
    [xs,o] = sort(xv);
    subplot(1,3,i); hold on; grid on;
    plot(xs, C_comsol(b(o))*1e12, 'o','MarkerSize',5,'Color',[0.2 0.3 0.75], ...
        'DisplayName','COMSOL  2W/V^2');
    plot(xs, C_centroid(b(o))*1e12,'-','LineWidth',1.8,'Color',[0.85 0.33 0.10], ...
        'DisplayName','centroid \epsilon_0A/g(X_c,Y_c)');
    plot(xs, C_exact(b(o))*1e12,'k--','LineWidth',1.6,'DisplayName','exact \int dA/g');
    xlabel(xl); ylabel('C [pF]');
    title(sprintf('%s block   (R^2_{centroid}=%+.3f)', blkNam{i}, r2fun(C_comsol(b),C_centroid(b))));
    legend('Location','best','FontSize',8);
end
sgtitle('Central electrode A: analytic capacitance vs COMSOL (V=30 V)');
outPng = fullfile(scriptDir,'AnalyticCapVsComsol.png');
exportgraphics(gcf, outPng);
fprintf('\nSaved plot to %s\n', outPng);

%% ---------------- helper ----------------
function [A, Xc, Yc, dX, dY] = trapMoments(rows)
% Area and centroid of a signed trapezoid list, rows = [sign ya yb pL qL pR qR].
A = 0; Mx = 0; My = 0;
xs = []; ys = [];
for k = 1:size(rows,1)
    s=rows(k,1); ya=rows(k,2); yb=rows(k,3);
    pL=rows(k,4); qL=rows(k,5); pR=rows(k,6); qR=rows(k,7);
    dp = pR-pL;  dq = qR-qL;
    % width(y) = dp + dq*y
    a  = dp*(yb-ya) + dq*(yb^2-ya^2)/2;
    my = dp*(yb^2-ya^2)/2 + dq*(yb^3-ya^3)/3;          % INT y*width dy
    % INT xmid*width dy, xmid = (pL+pR)/2 + (qL+qR)/2*y
    c0 = (pL+pR)/2; c1 = (qL+qR)/2;
    mx = c0*dp*(yb-ya) + (c0*dq + c1*dp)*(yb^2-ya^2)/2 + c1*dq*(yb^3-ya^3)/3;
    A = A + s*a;  Mx = Mx + s*mx;  My = My + s*my;
    xs = [xs, pL+qL*ya, pL+qL*yb, pR+qR*ya, pR+qR*yb];
    ys = [ys, ya, yb];
end
Xc = Mx/A;  Yc = My/A;
dX = max(xs)-min(xs);  dY = max(ys)-min(ys);
end
