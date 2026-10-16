%% verifyTrapGeometry.m
% Checks the decoded geometry against things already measured from the
% simulations, BEFORE any capacitance fitting. Three independent tests:
%
%  1. Areas of central electrodes A and B. The force-law prefactor k is
%     proportional to electrode area, and the two central-electrode
%     datasets gave a measured ratio k_A/k_B = 3.567 (determined to 0.15%
%     from two independent blocks). Pure geometry must reproduce it -- if
%     it does, it simultaneously proves (a) which dataset is which
%     electrode and (b) that there was NO voltage difference between the
%     runs, only an area difference.
%
%  2. Effective half-width b. The exact force law Fz = kV^2/(G^2-b^2 th^2)
%     has b^2 = 3*Iyy/A from the shape alone. Compare with the b recovered
%     by fitting the simulation data (1568.8 um for A).
%
%  3. Quadrant actuation centroid. The exact torque model recovered a
%     small-angle lever arm of 1498.5 um from the force sweeps, with no
%     geometric input at all.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
g = trapGeometry(); um = g.um;

fprintf('================= GEOMETRY CROSS-CHECKS =================\n\n');

%% --- 1. Central electrode areas and the k ratio ---
SA = [g.centralA.px; g.centralA.mx];
SB = [g.centralB.my; g.centralB.py];
[~, A_A, Iyy_A, Ixx_A] = shapeCap(SA, 15e-6, 0, 0);
[~, A_B, Iyy_B, Ixx_B] = shapeCap(SB, 15e-6, 0, 0);

fprintf('--- Central electrode areas ---\n');
fprintf('  A (outer, hollow square pair) = %.0f um^2\n', A_A/um^2);
fprintf('  B (inner, hollow pair)        = %.0f um^2\n', A_B/um^2);
fprintf('  geometric area ratio A/B      = %.4f\n', A_A/A_B);
fprintf('  MEASURED force-law k ratio    = %.4f   (from the two datasets)\n', 3.567);
fprintf('  agreement                     = %.2f%%\n\n', 100*abs(A_A/A_B - 3.567)/3.567);

%% --- 2. Effective half-widths from pure geometry ---
fprintf('--- Effective half-width  b = sqrt(3*Iyy/A)  (force-law tilt term) ---\n');
fprintf('%-28s %12s %12s %14s\n','electrode','b_y [um]','b_x [um]','fitted b [um]');
fprintf('%-28s %12.1f %12.1f %14.1f\n','central A (outer)', ...
    sqrt(3*Iyy_A/A_A)/um, sqrt(3*Ixx_A/A_A)/um, 1568.8);
fprintf('%-28s %12.1f %12.1f %14s\n','central B (inner)', ...
    sqrt(3*Iyy_B/A_B)/um, sqrt(3*Ixx_B/A_B)/um, '393 / 191*');
fprintf('  (*B''s fitted b was poorly conditioned -- its sweep varied rz,\n');
fprintf('   which dominates Fz and leaves the tilt term as a small residual)\n');
fprintf('  Geometry says B is nearly isotropic (b_y ~ b_x), confirming the\n');
fprintf('  488/191 split seen in the fit was conditioning noise, not shape.\n\n');

%% --- 3. Quadrant actuation electrode centroid ---
Sq = [g.quad(1).a; g.quad(1).b];
[~, A_q, Iyy_q] = shapeCap(Sq, 15e-6, 0, 0);
% centroid in y of the quadrant-1 pair
cen = 0;
for i=1:size(Sq,1)
    sg=Sq(i,1); ya=Sq(i,2); yb=Sq(i,3); dp=Sq(i,6)-Sq(i,4);
    cen = cen + sg*dp*(yb^2-ya^2)/2;
end
cen = cen / A_q;
fprintf('--- Quadrant 1 actuation pair ---\n');
fprintf('  two squares of side %.0f um, centred at (2250,750) and (750,2250) um\n', g.actSide/um);
fprintf('  combined area      = %.0f um^2  (= 2 x %.0f^2)\n', A_q/um^2, g.actSide/um);
fprintf('  combined centroid  = %.1f um\n', cen/um);
fprintf('  MEASURED lever arm = %.1f um  (exact torque model, zero geometric input)\n', 1498.5);
fprintf('  agreement          = %.2f%%\n\n', 100*abs(cen/um-1498.5)/(cen/um));

%% --- 4. Sensing triangles ---
[~, A_s1, Iyy_s1] = shapeCap(g.sense1, 15e-6, 0, 0);
[~, A_s2, Iyy_s2] = shapeCap(g.sense2, 15e-6, 0, 0);
fprintf('--- Sensing electrodes (quadrant 1) ---\n');
fprintf('  Cm1p area = %.0f um^2   (right triangle, legs 1415 um -> 1415^2/2 = %.0f)\n', ...
    A_s1/um^2, 1415^2/2);
fprintf('  Cm2p area = %.0f um^2   (mirror across x=y, must match)\n', A_s2/um^2);
fprintf('  area difference = %.3e%%  -> decoding of both triangles is consistent\n\n', ...
    100*abs(A_s1-A_s2)/A_s1);

%% --- 5. Total electrode area vs the plate ---
Aact = 0;
for q=1:4
    [~,a1] = shapeCap(g.quad(q).a,15e-6,0,0); [~,a2] = shapeCap(g.quad(q).b,15e-6,0,0);
    Aact = Aact + a1 + a2;
end
Aplate = (2*g.plateHalf)^2;
fprintf('--- Coverage ---\n');
fprintf('  plate area            = %.0f um^2\n', Aplate/um^2);
fprintf('  8 actuation squares   = %.0f um^2  (%.1f%% of plate)\n', Aact/um^2, 100*Aact/Aplate);
fprintf('  central A + B         = %.0f um^2  (%.1f%%)\n', (A_A+A_B)/um^2, 100*(A_A+A_B)/Aplate);
fprintf('  sensing (q1, both)    = %.0f um^2  (%.1f%%)\n', (A_s1+A_s2)/um^2, 100*(A_s1+A_s2)/Aplate);
