function [Fx,Fy,Fz,Taux,Tauy,Tauz] = trapElectrode2Model_v2(rx,ry,rz,thx,thy,thz,V)
%TRAPELECTRODE2MODEL_V2 Quadrant-2 electrode (mirror of quadrant 1
% across x=0), applied to the v2 (direct closed-form) quadrant-1 model.
% Same mirror/pseudovector transform as trapElectrode2Model.m.

[Fx1,Fy1,Fz1,Taux1,Tauy1,Tauz1] = trapElectrode1Model_v2(-rx, ry, rz, thx, -thy, -thz, V);

Fx = -Fx1; Fy = Fy1; Fz = Fz1;
Taux = Taux1; Tauy = -Tauy1; Tauz = -Tauz1;

end
