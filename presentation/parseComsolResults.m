function data = parseComsolResults(txtFile)
%PARSECOMSOLRESULTS Parses a COMSOL Global Evaluation/Probe table export
%(4 metadata lines + 1 header line, then rows of
%[rz,rx,ry,thx,thy,thz,Fx,Fy,Fz,capacitance_pF,Taux,Tauy,Tauz,Vfloat8,
%Vfloat9]) into a table with consistent naming and input column order
%(rx,ry,rz,thx,thy,thz) -- the raw export orders inputs as rz,rx,ry.

raw = readmatrix(txtFile, 'NumHeaderLines', 5);
if size(raw, 2) ~= 15
    error('Expected 15 columns in %s, got %d -- check the export format hasn''t changed.', ...
        txtFile, size(raw, 2));
end

rz  = raw(:,1);
rx  = raw(:,2);
ry  = raw(:,3);
thx = raw(:,4);
thy = raw(:,5);
thz = raw(:,6);
Fx             = raw(:,7);
Fy             = raw(:,8);
Fz             = raw(:,9);
capacitance_pF = raw(:,10);   % es.Q0_5/es.V0_5 (pF)
Taux    = raw(:,11);
Tauy    = raw(:,12);
Tauz    = raw(:,13);
Vfloat8 = raw(:,14);
Vfloat9 = raw(:,15);

data = table(rx, ry, rz, thx, thy, thz, Fx, Fy, Fz, capacitance_pF, Taux, Tauy, Tauz, Vfloat8, Vfloat9);

end
