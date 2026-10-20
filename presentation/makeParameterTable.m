%% makeParameterTable.m
% Emits a LaTeX table of every parameter used by the simulation, EXCEPT the
% per-electrode footprint geometry (Xc, Yc, dX, dY), which has its own table
% from makeElectrodeGeometryTable.m.
%
% Every value is read from trapDeviceGeometry()/trapControlDesign() rather
% than written out by hand, so the table cannot drift from the code. Derived
% quantities (weight, hover voltage, open-loop stiffnesses, the unstable
% rate p and cosh(pT)) are recomputed here from the same source.
%
% Output: presentation/trap_parameters_table.tex
% Requires: \usepackage{booktabs} in the document preamble.
%
% See also MAKEELECTRODEGEOMETRYTABLE, TRAPDEVICEGEOMETRY, TRAPCONTROLDESIGN.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p  = trapDeviceGeometry();
K  = trapControlDesign(p);
W  = p.mass.m*p.phys.g;
Kst = trapOpenLoopStiffness(p, 'Vtia', 0);
u0  = trapWrenchToVoltage(p, zeros(3,1), zeros(6,1), 'Vtia', 0);
Ce0 = mean(trapSensing(p, trapCapacitance(p, zeros(6,1))));
C0  = p.ceff.C12 + p.ceff.Cg/2 + p.ceff.Cf/2;
pr  = K.pRate;

% rows: {symbol, description, value, unit}  -- 'Â§' starts a new section
R = {
'Â§','Physical constants',[],''
'\varepsilon_0','vacuum permittivity',p.phys.eps0,'\si{\farad\per\metre}'
'g','gravitational acceleration',p.phys.g,'\si{\metre\per\second\squared}'

'Â§','Proof mass',[],''
'L','plate length',p.mass.L,'\si{\micro\metre}'
'W_{\mathrm{p}}','plate width',p.mass.W,'\si{\micro\metre}'
't','plate thickness',p.mass.t,'\si{\micro\metre}'
'\rho','silicon density',p.mass.rho,'\si{\kilo\gram\per\metre\cubed}'
'm','plate mass',p.mass.m,'\si{\milli\gram}'
'I_{xx}=I_{yy}','transverse moment of inertia',p.mass.I(1,1),'\si{\kilo\gram\metre\squared}'
'I_{zz}','yaw moment of inertia',p.mass.I(3,3),'\si{\kilo\gram\metre\squared}'

'Â§','Gaps and travel',[],''
'g_0','nominal electrode gap at equilibrium',p.gap.z0,'\si{\micro\metre}'
'2g_0','total vertical channel height',p.gap.zChannel,'\si{\micro\metre}'
'xy_0','nominal lateral clearance',p.gap.xy0,'\si{\micro\metre}'
'2xy_0','total lateral channel width',p.gap.xyChannel,'\si{\micro\metre}'

'Â§','Environment',[],''
'P','chamber pressure ($10^{-4}$\,Torr $=0.1$\,mTorr)',p.env.pressure,'\si{\pascal}'
'\bar{v}','mean molecular speed',p.env.vbar,'\si{\metre\per\second}'
'T_{\mathrm{amb}}','ambient temperature',p.env.temperature,'\si{\kelvin}'

'Â§','Squeeze-film damping',[],''
'c','translational damping, $c=4L^3P/(\pi^2\bar{v}g_0)$',p.damp.c,'\si{\newton\second\per\metre}'
'q','rotational damping, $q=(L^2/12)\,c$',p.damp.q,'\si{\newton\metre\second\per\radian}'
'c_{\mathrm{shear}}','lateral (shear) damping',p.damp.cShear,'\si{\newton\second\per\metre}'

'Â§','Capacitance and sensing model',[],''
'\delta/g','conformal-mapping edge extension, $2/\pi$',p.fringe.deltaOverGap,'--'
'C_{12}','pair-to-pair direct capacitance',p.ceff.C12,'\si{\femto\farad}'
'C_{g}','electrode-to-ground capacitance',p.ceff.Cg,'\si{\femto\farad}'
'C_{f}','unmodelled fringing capacitance',p.ceff.Cf,'\si{\femto\farad}'
'C_{0}','pose-independent floor, $C_{12}+C_g/2+C_f/2$',C0,'\si{\femto\farad}'
'C_{\mathrm{eff}}','effective pair capacitance at equilibrium',Ce0,'\si{\femto\farad}'
'\ell','sensing-pair moment arm',K.ell,'\si{\micro\metre}'
'a','actuation-pair moment arm',K.a,'\si{\micro\metre}'

'Â§','Timing',[],''
'T_{\mathrm{ISR}}','controller interrupt period',p.readout.Tisr,'\si{\micro\second}'
'\Delta t_{\mathrm{ch}}','inter-channel conversion spacing',p.readout.dtChannel,'\si{\micro\second}'
't_{\mathrm{I2C}}','bus transfer after the last conversion',p.readout.tI2C,'\si{\micro\second}'
't_{d}','oldest-channel staleness, $t_{\mathrm{I2C}}+3.5\Delta t_{\mathrm{ch}}$',p.readout.tI2C+3.5*p.readout.dtChannel,'\si{\micro\second}'
'\Delta t','forward-Euler integration step',p.sim.dtEuler,'\si{\micro\second}'

'Â§','Controller',[],''
'\lambda_r','regulator bandwidth',K.lambdar/2/pi,'\si{\hertz}'
'\lambda_e','observer bandwidth',K.lambdae/2/pi,'\si{\hertz}'
'R_{\mathrm{LQR}}','LQR control weight (bisected to meet $\lambda_r$)',K.lqrR,'--'
'V_{\mathrm{TIA}}','fixed bias on the TIA pair',K.Vtia,'\si{\volt}'
'V_{\max}','drive voltage clamp',80,'\si{\volt}'

'Â§','Derived operating point',[],''
'W','plate weight, $mg$',W,'\si{\micro\newton}'
'V_{\mathrm{hov}}','hover voltage on each quadrant pair',u0(1),'\si{\volt}'
'k_{zz}','heave stiffness at fixed voltage, $-2W/g_0$',Kst(3,3),'\si{\newton\per\metre}'
'k_{\theta_x}','roll stiffness at fixed voltage',Kst(4,4),'\si{\newton\metre\per\radian}'
'k_{\theta_y}','pitch stiffness at fixed voltage',Kst(5,5),'\si{\newton\metre\per\radian}'
'k_{xx}','lateral stiffness (passively trapped)',Kst(1,1),'\si{\newton\per\metre}'
'k_{\theta_z}','yaw stiffness (passively trapped)',Kst(6,6),'\si{\newton\metre\per\radian}'
'p_{z}','heave unstable rate, $\sqrt{2g/g_0}$',pr(1),'\si{\radian\per\second}'
'p_{\theta_x}','roll unstable rate',pr(2),'\si{\radian\per\second}'
'p_{\theta_y}','pitch unstable rate',pr(3),'\si{\radian\per\second}'
'\cosh(p_zT)','per-tick velocity growth over one ISR, heave',cosh(pr(1)*K.Tisr),'--'
};

% scale factors so values print in the unit named in the row
scaleOf = containers.Map( ...
 {'\si{\micro\metre}','\si{\milli\gram}','\si{\femto\farad}','\si{\micro\second}', ...
  '\si{\micro\newton}'}, {1e6,1e6,1e15,1e6,1e6} );

out = fullfile(root,'presentation','trap_parameters_table.tex');
fid = fopen(out,'w');
fprintf(fid,'%% Generated by presentation/makeParameterTable.m -- do not edit by hand.\n');
fprintf(fid,'%% Requires \\usepackage{booktabs} and \\usepackage{siunitx}.\n');
fprintf(fid,'\\begin{table}[htbp]\n  \\centering\n');
fprintf(fid,'  \\caption{Simulation parameters. Electrode footprint geometry\n');
fprintf(fid,'           ($X_c$, $Y_c$, $\\Delta X$, $\\Delta Y$) is tabulated separately.\n');
fprintf(fid,'           Stiffnesses are quoted at fixed voltage, so a negative value is\n');
fprintf(fid,'           anti-restoring; $p=\\sqrt{-k/m}$ is the corresponding unstable rate.}\n');
fprintf(fid,'  \\label{tab:trap_parameters}\n');
fprintf(fid,'  \\begin{tabular}{llSl}\n    \\toprule\n');
fprintf(fid,'    Symbol & Description & {Value} & Unit \\\\\n');
for i = 1:size(R,1)
    if strcmp(R{i,1},'Â§')
        fprintf(fid,'    \\midrule\n    \\multicolumn{4}{l}{\\textit{%s}} \\\\\n', R{i,2});
        continue
    end
    v = R{i,3};  u = R{i,4};
    if isKey(scaleOf,u), v = v*scaleOf(u); end
    if v == 0
        s = '0';
    elseif abs(v) >= 1e-3 && abs(v) < 1e5
        s = sprintf('%.4g', v);
    else
        s = sprintf('%.4e', v);
    end
    fprintf(fid,'    $%s$ & %s & %s & %s \\\\\n', R{i,1}, R{i,2}, s, u);
end
fprintf(fid,'    \\bottomrule\n  \\end{tabular}\n\\end{table}\n');
fclose(fid);

fprintf('Wrote %s  (%d rows)\n', out, sum(~strcmp(R(:,1),'Â§')));
type(out)
