function [Fx,Fy,Fz,Taux,Tauy,Tauz] = trapElectrode4Model(rx,ry,rz,thx,thy,thz,U)
%TRAPELECTRODE4MODEL Force/torque from the quadrant-4 electrode
% (mirror image of quadrant 1 across y=0). Same 7-in/6-out signature as
% trapElectrode1Model. See trapElectrode2Model for the mirror
% (pseudovector) sign-flip rule.

[Fx1,Fy1,Fz1,Taux1,Tauy1,Tauz1] = trapElectrode1Model(rx, -ry, rz, -thx, thy, -thz, U);

Fx = Fx1; Fy = -Fy1; Fz = Fz1;
Taux = -Taux1; Tauy = Tauy1; Tauz = -Tauz1;

end
