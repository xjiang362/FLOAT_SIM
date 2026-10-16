function [Fx,Fy,Fz,Taux,Tauy,Tauz] = trapElectrode3Model(rx,ry,rz,thx,thy,thz,U)
%TRAPELECTRODE3MODEL Force/torque from the quadrant-3 electrode
% (quadrant 1 rotated 180deg about z). Same 7-in/6-out signature as
% trapElectrode1Model.
%
% A 180deg rotation about z is a PROPER transform (det=+1), so both
% polar vectors (rx,ry,rz / Fx,Fy,Fz) and axial vectors (thx,thy,thz /
% Taux,Tauy,Tauz) transform the same simple way -- no extra pseudovector
% sign flip, unlike the mirror-based quadrant 2/4 transforms.

[Fx1,Fy1,Fz1,Taux1,Tauy1,Tauz1] = trapElectrode1Model(-rx, -ry, rz, -thx, -thy, thz, U);

Fx = -Fx1; Fy = -Fy1; Fz = Fz1;
Taux = -Taux1; Tauy = -Tauy1; Tauz = Tauz1;

end
