function [Fx,Fy,Fz,Taux,Tauy,Tauz] = trapElectrode3Model_v2(rx,ry,rz,thx,thy,thz,V)
%TRAPELECTRODE3MODEL_V2 Quadrant-3 electrode (quadrant 1 rotated 180deg
% about z), applied to the v2 (direct closed-form) quadrant-1 model.
% Same proper-rotation transform as trapElectrode3Model.m.

[Fx1,Fy1,Fz1,Taux1,Tauy1,Tauz1] = trapElectrode1Model_v2(-rx, -ry, rz, -thx, -thy, thz, V);

Fx = -Fx1; Fy = -Fy1; Fz = Fz1;
Taux = -Taux1; Tauy = -Tauy1; Tauz = Tauz1;

end
