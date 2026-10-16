%% makeElectrodeGeometryTable.m
% Emits a LaTeX table of the electrode footprint geometry (Xc, Yc, dX, dY)
% using the TERMINAL NUMBERING of the source/terminal sweep.
%
% Numbering convention (per quadrant q = 1..4, four terminals each):
%   4q-3  C_sns q+   sensing   electrode, outer-x half of the quadrant
%   4q-2  C_sns q-   sensing   electrode, outer-y half of the quadrant
%   4q-1  C_act q+   actuation electrode, outer-x half of the quadrant
%   4q    C_act q-   actuation electrode, outer-y half of the quadrant
% "+" is the member of the pair with |Xc| > |Yc| (further out along x),
% "-" the one with |Yc| > |Xc|. In quadrant 1 this reduces exactly to the
% x>y / x<y split.
%   17,18  Cp+,  Cp-     larger inner pair  (central electrode A)
%   19,20  C_TIA+, C_TIA- smaller inner pair (central electrode B)
%   21     resonator (proof mass)      22  frame
%
% Hollow shapes are given as TWO rows: the solid outer rectangle (marked
% with a star) and the subtracted inner rectangle. The net footprint is
% outer MINUS inner. For Cp the subtracted region is the well that houses
% the TIA pair.
%
% Frame coupling Cm(21,22) is decomposed as five parallel plate capacitors:
% four side faces Crs1..4 (plate edge to side wall, lateral gap) and one
% bottom face Crb1 (plate underside to floor, vertical gap).
%
% All lengths in micrometres. Source of truth: trapGeometry.m for the
% electrode footprints; sampleTrapPoseLHS.m for the gaps; trapDeviceParams.m
% for the plate thickness.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
g = trapGeometry(); um = 1e-6;

plateHalf = g.plateHalf/um;        % 2975 um
gz  = 15;                          % [um] vertical gap, plate centred in a 2*gz channel
gxy = 50;                          % [um] lateral gap
tPlate = 350;                      % [um] plate thickness

%% ---------- quadrant-1 primitives, in micrometres ----------
% Sensing triangles. Xc,Yc are the AREA centroid (mean of the vertices) --
% that is the moment the first-order gap expansion needs, and it is what
% the thx fit validated at R^2 = 0.9987.
%
% dY is the full y-extent. dX is the width measured at the BOUNDING-BOX
% MID-HEIGHT, i.e. XR(Ymid) - XL(Ymid), which for a right triangle is
% exactly half the base. That choice makes dX*dY equal the true area A for
% triangles just as it does for rectangles, so the area column is uniform
% and C ~= eps0*dX*dY/g needs no shape-dependent factor.
%
% NOTE the width is NOT evaluated at the area centroid: at Yc the width is
% 2/3 of the base (943.33 um here), which would make dX*dY = 1.333*A. The
% two heights must not be conflated -- Yc is for the tilt term, Ymid is
% for the area.
v1 = [1560 1525; 2975 1525; 2975 2940];        % Cm1p, the x>y half
v2 = [1525 1560; 1525 2975; 2940 2975];        % Cm2p, its mirror across x=y
tri = @(v) struct('Xc',mean(v(:,1)),'Yc',mean(v(:,2)), ...
                  'dX',polyarea(v(:,1),v(:,2))/(max(v(:,2))-min(v(:,2))), ...
                  'dY',max(v(:,2))-min(v(:,2)), ...
                  'A',polyarea(v(:,1),v(:,2)),'shape','triangle');
S1 = tri(v1); S2 = tri(v2);

s = (g.actSide/um)/2;                           % 725 um
rct = @(cx,cy,w,h) struct('Xc',cx,'Yc',cy,'dX',w,'dY',h,'A',w*h,'shape','rectangle');
A1 = rct(2250, 750, 2*s, 2*s);                  % C_act1+
A2 = rct( 750,2250, 2*s, 2*s);                  % C_act1-

q1 = {S1, S2, A1, A2};
q1name = {'C_{\mathrm{sns}1+}','C_{\mathrm{sns}1-}','C_{\mathrm{act}1+}','C_{\mathrm{act}1-}'};

% Quadrant sign maps: Q1 (+,+), Q2 (-,+), Q3 (-,-), Q4 (+,-).
sx = [ 1 -1 -1  1];
sy = [ 1  1 -1 -1];

% WITHIN-PAIR ORDER. The terminals advance by a 90-degree ROTATION from
% quadrant to quadrant, not by mirroring, so |Xc| and |Yc| swap roles every
% step. Quadrants 1 and 3 are 180 degrees apart and keep the quadrant-1
% order; quadrants 2 and 4 take the pair in the opposite order.
% Verified against the sweep: fitting C_iR with the quadrant-1 order gives
% R^2 = 0.996 for quadrants 1 and 3 but only 0.927 for 2 and 4, and the
% swapped order reverses that exactly. Equivalently, rotating quadrant 1's
% C_act1+ centroid (2250,750) by +90 degrees gives (-750,2250), which is
% the SMALLER-|Xc| member of the quadrant-2 pair.
pairOrder = {[1 2 3 4], [2 1 4 3], [1 2 3 4], [2 1 4 3]};

%% ---------- central electrodes ----------
% Cp (central A): outer rect minus the well that houses the TIA
Cp_out_p = rct((75+1475)/2, 0, 1475-75, 2*1475);
Cp_in_p  = rct((75+737)/2,  0,  737-75, 2*725);
% C_TIA (central B): outer rect minus the small notch
TIA_out_p = rct(0,  (25+675)/2, 1400, 675-25);
TIA_in_p  = rct(0, (25+100)/2,   200, 100-25);

%% ---------- assemble the row list ----------
rows = {};   % {term, label, shape, Xc, Yc, dX, dY, Area, note}
for q = 1:4
    for k = 1:4
        e = q1{ pairOrder{q}(k) };          % 90-degree rotation, not mirror
        term = 4*(q-1) + k;
        lbl = strrep(q1name{k}, '1', sprintf('%d',q));
        rows(end+1,:) = {term, lbl, e.shape, sx(q)*e.Xc, sy(q)*e.Yc, e.dX, e.dY, e.A, ''}; %#ok<SAGROW>
    end
end
% ORIENTATION OF THE Cp PAIR -- verified against the matrix, not assumed.
% Electrodes at +x (terminals 1,2,3) couple to terminal 18 at 2.39-49.66 fF
% and to terminal 17 at only 0.03-0.20 fF; electrodes at -x (terminals
% 5,6,7) show exactly the mirrored pattern. Coupling falls with distance,
% so terminal 18 is the +x half and terminal 17 is the -x half. The names
% Cp+ / Cp- therefore do NOT follow the sign of Xc. (The TIA pair does:
% terminal 19 is the +y half, 20 the -y half, confirmed the same way.)
rows(end+1,:) = {17,'C_{p+}^{*}',   'rectangle', -Cp_out_p.Xc, 0, Cp_out_p.dX, Cp_out_p.dY, Cp_out_p.A, 'solid outer'};
rows(end+1,:) = {17,'C_{p+}^{\,\mathrm{well}}','rectangle',-Cp_in_p.Xc, 0, Cp_in_p.dX, Cp_in_p.dY, Cp_in_p.A, 'subtracted'};
rows(end+1,:) = {18,'C_{p-}^{*}',   'rectangle',  Cp_out_p.Xc, 0, Cp_out_p.dX, Cp_out_p.dY, Cp_out_p.A, 'solid outer'};
rows(end+1,:) = {18,'C_{p-}^{\,\mathrm{well}}','rectangle', Cp_in_p.Xc, 0, Cp_in_p.dX, Cp_in_p.dY, Cp_in_p.A, 'subtracted'};
rows(end+1,:) = {19,'C_{\mathrm{TIA}+}^{*}','rectangle', 0,  TIA_out_p.Yc, TIA_out_p.dX, TIA_out_p.dY, TIA_out_p.A,'solid outer'};
rows(end+1,:) = {19,'C_{\mathrm{TIA}+}^{\,\mathrm{notch}}','rectangle',0, TIA_in_p.Yc, TIA_in_p.dX, TIA_in_p.dY, TIA_in_p.A,'subtracted'};
rows(end+1,:) = {20,'C_{\mathrm{TIA}-}^{*}','rectangle', 0, -TIA_out_p.Yc, TIA_out_p.dX, TIA_out_p.dY, TIA_out_p.A,'solid outer'};
rows(end+1,:) = {20,'C_{\mathrm{TIA}-}^{\,\mathrm{notch}}','rectangle',0,-TIA_in_p.Yc, TIA_in_p.dX, TIA_in_p.dY, TIA_in_p.A,'subtracted'};

%% ---------- console check ----------
fprintf('%5s %-28s %-10s %9s %9s %8s %8s %12s\n', ...
    'term','label','shape','Xc','Yc','dX','dY','A [um^2]');
for i=1:size(rows,1)
    fprintf('%5d %-28s %-10s %9.1f %9.1f %8.0f %8.0f %12.0f\n', ...
        rows{i,1}, rows{i,2}, rows{i,3}, rows{i,4}, rows{i,5}, rows{i,6}, rows{i,7}, rows{i,8});
end
fprintf('\nnet Cp per half  = %.0f um^2   net TIA per half = %.0f um^2\n', ...
    Cp_out_p.A-Cp_in_p.A, TIA_out_p.A-TIA_in_p.A);
fprintf('total actuation (8 squares) = %.0f um^2, total sensing (8 triangles) = %.0f um^2\n', ...
    8*A1.A, 8*S1.A);

%% ---------- write LaTeX ----------
texFile = fullfile(scriptDir,'electrodeGeometryTable.tex');
fid = fopen(texFile,'wt');
fprintf(fid, '%% Electrode footprint geometry, terminal-sweep numbering.\n');
fprintf(fid, '%% Generated by makeElectrodeGeometryTable.m -- do not edit by hand.\n');
fprintf(fid, '%% Requires \\usepackage{booktabs,siunitx}\n\n');
fprintf(fid, '\\begin{table}[htbp]\n  \\centering\n');
fprintf(fid, '  \\caption{Electrode footprint geometry in plate-centred coordinates.\n');
fprintf(fid, '  $X_c,Y_c$ are the area centroid. $\\Delta Y$ is the full $y$-extent and\n');
fprintf(fid, '  $\\Delta X$ the width at mid-height, $X_R(Y_\\mathrm{mid})-X_L(Y_\\mathrm{mid})$,\n');
fprintf(fid, '  which for the triangular sensing electrodes is half the base. With this\n');
fprintf(fid, '  convention $\\Delta X\\,\\Delta Y = A$ for every footprint, triangles\n');
fprintf(fid, '  included. Note $\\Delta X$ is \\emph{not} the width at $Y_c$: there a\n');
fprintf(fid, '  triangle is $\\tfrac{2}{3}$ of its base. Rows marked $^{*}$ give the solid\n');
fprintf(fid, '  outer rectangle; the row beneath is the region subtracted from it, so the\n');
fprintf(fid, '  net footprint is outer minus inner. All lengths in \\si{\\micro\\metre}.}\n');
fprintf(fid, '  \\label{tab:electrode-geometry}\n');
fprintf(fid, '  \\begin{tabular}{r l l S[table-format=-5.1] S[table-format=-5.1] S[table-format=5.1] S[table-format=4.0] S[table-format=8.1]}\n');
fprintf(fid, '    \\toprule\n');
fprintf(fid, '    Term. & Label & Shape & {$X_c$} & {$Y_c$} & {$\\Delta X$} & {$\\Delta Y$} & {Area} \\\\\n');
fprintf(fid, '    \\midrule\n');
lastq = 0;
for i=1:size(rows,1)
    t = rows{i,1};
    q = ceil(t/4);
    if t<=16 && q~=lastq && i>1, fprintf(fid,'    \\addlinespace\n'); end
    if t==17 && lastq<=4, fprintf(fid,'    \\midrule\n'); end
    lastq = q;
    fprintf(fid, '    %d & $%s$ & %s & %.1f & %.1f & %.1f & %.0f & %.1f \\\\\n', ...
        t, rows{i,2}, rows{i,3}, rows{i,4}, rows{i,5}, rows{i,6}, rows{i,7}, rows{i,8});
end

% ---- frame decomposition ----
fprintf(fid, '    \\midrule\n');
fprintf(fid, '    \\multicolumn{8}{l}{\\emph{Resonator--frame coupling $C_{m}(21,22)$, five parallel plate capacitors}}\\\\\n');
fprintf(fid, '    \\addlinespace[2pt]\n');
sideName = {'C_{rs1}','C_{rs2}','C_{rs3}','C_{rs4}'};
sideFace = {'$+x$ edge','$+y$ edge','$-x$ edge','$-y$ edge'};
sideXc = [ plateHalf, 0, -plateHalf, 0];
sideYc = [ 0, plateHalf, 0, -plateHalf];
for k=1:4
    if mod(k,2)==1   % +-x faces: extent is along y, thickness in z
        dX = tPlate; dY = 2*plateHalf;
    else             % +-y faces: extent along x
        dX = 2*plateHalf; dY = tPlate;
    end
    fprintf(fid, '    22 & $%s$ & side & %.1f & %.1f & %.1f & %.0f & %.1f \\\\\n', ...
        sideName{k}, sideXc(k), sideYc(k), dX, dY, 2*plateHalf*tPlate);
end
fprintf(fid, '    22 & $C_{rb1}$ & bottom & 0.0 & 0.0 & %.1f & %.0f & %.1f \\\\\n', ...
    2*plateHalf, 2*plateHalf, (2*plateHalf)^2);
fprintf(fid, '    \\bottomrule\n  \\end{tabular}\n\n');
fprintf(fid, '  \\vspace{2pt}\n');
fprintf(fid, '  \\footnotesize Nominal gaps: \\SI{%d}{\\micro\\metre} vertical (the plate sits\n', gz);
fprintf(fid, '  centred in a \\SI{%d}{\\micro\\metre} channel, so top and bottom gaps are equal),\n', 2*gz);
fprintf(fid, '  \\SI{%d}{\\micro\\metre} lateral. Plate \\SI{%d}{\\micro\\metre} square,\n', gxy, 2*plateHalf);
fprintf(fid, '  \\SI{%d}{\\micro\\metre} thick. For the four side capacitors $C_{rs1..4}$ the\n', tPlate);
fprintf(fid, '  listed $\\Delta X,\\Delta Y$ are the in-plane edge length and the plate\n');
fprintf(fid, '  thickness respectively, and the gap is lateral rather than vertical.\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('\nWrote %s\n', texFile);
