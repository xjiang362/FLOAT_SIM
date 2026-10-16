function [Fx,Fy,Fz,Taux,Tauy,Tauz] = trapElectrode2Model(rx,ry,rz,thx,thy,thz,U)
%TRAPELECTRODE2MODEL Force/torque from the quadrant-2 electrode
% (mirror image of quadrant 1 across x=0). Same 7-in/6-out signature as
% trapElectrode1Model.
%
% rx,ry,rz are polar vectors (flip sign under mirror along with the
% mirrored axis); thx,thy,thz and Taux,Tauy,Tauz are axial vectors
% (rotations/torques), which pick up an EXTRA sign flip under a mirror
% (improper transform) beyond the mirrored-axis component -- confirmed
% via the standard pseudovector rule v' = det(R)*R*v.

[Fx1,Fy1,Fz1,Taux1,Tauy1,Tauz1] = trapElectrode1Model(-rx, ry, rz, thx, -thy, -thz, U);

Fx = -Fx1; Fy = Fy1; Fz = Fz1;
Taux = Taux1; Tauy = -Tauy1; Tauz = -Tauz1;

end
