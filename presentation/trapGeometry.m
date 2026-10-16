function geom = trapGeometry()
%TRAPGEOMETRY Electrode footprints of the FLOATS trap, as signed
% trapezoid lists (units: metres, plate-centred coordinates).
%
% Every electrode is stored as an N-by-7 array of rows
%       [sign, ya, yb, pL, qL, pR, qR]
% meaning: add (sign=+1) or subtract (sign=-1) the region swept in y from
% ya to yb between the straight edges xL(y)=pL+qL*y and xR(y)=pR+qR*y.
% This is exactly the argument list of capTrapz, so the capacitance of
% any electrode is the signed sum of that one primitive. Hollow shapes
% (the central electrodes) are a rectangle minus a rectangle; the sensing
% electrodes are right triangles expressed with a sloped edge (q = +-1).
%
% Source: dimensions supplied by the user. Proof mass 5950 um square.
% Actuation electrodes are square DIFFERENTIAL PAIRS of side 1450 um --
% each quadrant drives two sub-electrodes at opposite polarity so the
% proof mass floats near zero potential.

um = 1e-6;
geom.um = um;
geom.plateHalf = 2975*um;          % proof mass half-width (5950 um square)
geom.actSide   = 1450*um;

rect = @(x1,x2,y1,y2) [1, y1, y2, x1, 0, x2, 0];
neg  = @(x1,x2,y1,y2) [-1, y1, y2, x1, 0, x2, 0];

%% ---------- Quadrant actuation electrodes (differential pairs) ----------
% Quadrant 1 pair: squares centred at (2250,750) and (750,2250).
s = geom.actSide/2;
c1 = [2250, 750]*um;      % sub-electrode "a"
c2 = [ 750, 2250]*um;     % sub-electrode "b"
mkSq = @(c) rect(c(1)-s, c(1)+s, c(2)-s, c(2)+s);

geom.quad(1).a = mkSq(c1);    geom.quad(1).b = mkSq(c2);
% Quadrants 2,3,4 by the established mirror / 180-degree symmetry.
mirX = @(R) [R(1), R(2), R(3), -R(6), -R(7), -R(4), -R(5)];             % x -> -x
mirY = @(R) [R(1), -R(3), -R(2), R(4), -R(5), R(6), -R(7)];             % y -> -y
geom.quad(2).a = mirX(geom.quad(1).a);  geom.quad(2).b = mirX(geom.quad(1).b);
geom.quad(3).a = mirY(mirX(geom.quad(1).a)); geom.quad(3).b = mirY(mirX(geom.quad(1).b));
geom.quad(4).a = mirY(geom.quad(1).a);  geom.quad(4).b = mirY(geom.quad(1).b);

%% ---------- Sensing electrodes (right triangles), quadrant 1 ----------
% Cm1p: vertices (1560,1525) (2975,1525) (2975,2940).
%   Right angle at (2975,1525); hypotenuse has slope 1 -> x = y + 35.
%   Sweep y from 1525 to 2940 with xL(y) = y+35, xR = 2975.
geom.sense1 = [1, 1525*um, 2940*um,  35*um, 1, 2975*um, 0];
% Cm2p: the mirror image of Cm1p across the line x = y, i.e. vertices
%   (1525,1560) (1525,2975) (2940,2975). Sweep y from 1560 to 2975 with
%   xL = 1525, xR(y) = y - 35.
geom.sense2 = [1, 1560*um, 2975*um, 1525*um, 0, -35*um, 1];

%% ---------- Central electrode A (OUTER): hollow square pair ----------
% +x half: outer rectangle minus inner rectangle (the hollow).
% NOTE: no space between the function name and "(" -- inside [ ] a space
% is a column separator and would split the call into two elements.
geom.centralA.px = [ rect(75*um, 1475*um, -1475*um, 1475*um);
                     neg(75*um, 737*um, -725*um, 725*um) ];
% -x half: mirror about the y-axis (pair A is y-axis symmetric).
geom.centralA.mx = [ mirX(geom.centralA.px(1,:));
                     mirX(geom.centralA.px(2,:)) ];

%% ---------- Central electrode B (INNER): smaller hollow pair ----------
% -y half: outer rectangle minus a notch at the top centre.
geom.centralB.my = [ rect(-700*um, 700*um, -675*um, -25*um);
                     neg(-100*um, 100*um, -100*um, -25*um) ];
% +y half: mirror about the x-axis (pair B is x-axis symmetric).
geom.centralB.py = [ mirY(geom.centralB.my(1,:));
                     mirY(geom.centralB.my(2,:)) ];
end
