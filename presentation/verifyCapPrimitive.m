%% verifyCapPrimitive.m
% Nothing downstream is trustworthy unless capTrapz is exact, so check it
% against brute-force 2-D quadrature (integral2) across every branch:
% generic (both tilts), thy->0, thx->0, both->0, and the tilt magnitudes
% actually reached by the device (where y*thx is COMPARABLE to the gap,
% not a small perturbation -- so the exact form genuinely matters).

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rng(7);

g0 = 15e-6;
thMax = 0.1289*pi/180;      % established tilt limit
um = 1e-6;

fprintf('=========== capTrapz vs integral2 ===========\n');
fprintf('%-34s %14s %14s %11s\n','case','closed form','integral2','rel err');

cases = {};
% name, ya, yb, pL, qL, pR, qR, G, thx, thy
cases(end+1,:) = {'rect, no tilt',            25*um, 1475*um, 1525*um,0, 2975*um,0, g0, 0, 0};
cases(end+1,:) = {'rect, thx only',           25*um, 1475*um, 1525*um,0, 2975*um,0, g0, thMax, 0};
cases(end+1,:) = {'rect, thy only',           25*um, 1475*um, 1525*um,0, 2975*um,0, g0, 0, thMax};
cases(end+1,:) = {'rect, both tilts',         25*um, 1475*um, 1525*um,0, 2975*um,0, g0, thMax, -0.7*thMax};
cases(end+1,:) = {'rect, tiny thy (branch)',  25*um, 1475*um, 1525*um,0, 2975*um,0, g0, thMax, 1e-13};
cases(end+1,:) = {'rect, tiny thx (B~0)',     25*um, 1475*um, 1525*um,0, 2975*um,0, g0, 1e-13, thMax};
cases(end+1,:) = {'rect, rz=+10um gap',       25*um, 1475*um, 1525*um,0, 2975*um,0, 5e-6, thMax, 0.5*thMax};
% sensing triangle 1: y 1525->2940, left edge x=y+35, right edge x=2975
cases(end+1,:) = {'triangle, no tilt',       1525*um,2940*um, 35*um,1, 2975*um,0, g0, 0, 0};
cases(end+1,:) = {'triangle, thx only',      1525*um,2940*um, 35*um,1, 2975*um,0, g0, thMax, 0};
cases(end+1,:) = {'triangle, both tilts',    1525*um,2940*um, 35*um,1, 2975*um,0, g0, thMax, -0.7*thMax};
% Gap must stay positive over the whole footprint: the triangle reaches
% y=2940um, so y*thx alone is 6.6um at full tilt -- pick G and the tilts
% so the plate never intersects the electrode.
cases(end+1,:) = {'triangle, tight gap',     1525*um,2940*um, 35*um,1, 2975*um,0, 12e-6, 0.5*thMax, 0.3*thMax};

maxRel = 0;
for i = 1:size(cases,1)
    [nm,ya,yb,pL,qL,pR,qR,G,tx,ty] = cases{i,:};
    Cc = capTrapz(ya,yb,pL,qL,pR,qR,G,tx,ty);
    % integral2(FUN,XMIN,XMAX,YMIN,YMAX) makes the FIRST argument the
    % OUTER variable, with the inner limits given as functions of it. Our
    % region sweeps y outer and x inner, so y must occupy the first slot
    % and the integrand's arguments are swapped accordingly. (Getting
    % this backwards is invisible at zero tilt, where the integrand is
    % symmetric in the two variables -- which is exactly how the first
    % run of this check produced passing no-tilt rows and failing
    % tilted ones.)
    Cn = integral2(@(yy,xx) 1./(G - yy*tx + xx*ty), ya, yb, ...
                   @(yy) pL+qL*yy, @(yy) pR+qR*yy, ...
                   'AbsTol',1e-16,'RelTol',1e-12);
    rel = abs(Cc-Cn)/abs(Cn);
    maxRel = max(maxRel, rel);
    fprintf('%-34s %14.6e %14.6e %11.2e\n', nm, Cc, Cn, rel);
end

fprintf('\nworst relative error across all branches: %.3e\n', maxRel);
if maxRel < 1e-9
    fprintf('PASS -- closed form is exact to quadrature precision.\n');
else
    fprintf('FAIL -- investigate before building anything on top of this.\n');
end

%% How large is the tilt effect? (is the exact form actually needed?)
fprintf('\n=========== is the exact form necessary? ===========\n');
ya=25*um; yb=1475*um; pL=1525*um; pR=2975*um;
C0   = capTrapz(ya,yb,pL,0,pR,0, g0, 0, 0);
Cmax = capTrapz(ya,yb,pL,0,pR,0, g0, thMax, 0);
A    = (yb-ya)*(pR-pL);
fprintf('flat-plate A/g            = %.6e\n', A/g0);
fprintf('exact at max thx          = %.6e  (%.2f%% above flat)\n', Cmax, 100*(Cmax-C0)/C0);
fprintf('max |y*thx| / g           = %.3f   <- perturbation parameter\n', yb*thMax/g0);
fprintf('  (order 1 => the tilt is NOT a small correction; a linearized\n');
fprintf('   effective-gap model cannot represent this faithfully)\n');
