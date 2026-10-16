function [Fx,Fy,Fz,Taux,Tauy,Tauz] = trapElectrode1Model_v2(rx,ry,rz,thx,thy,thz,V)
%TRAPELECTRODE1MODEL_V2 Consolidated quadrant-1 electrode force/torque
% model, replacing the differentiated-fringing-formula approach
% (trapElectrode1Model.m / _crossTerms / _coupledWL) with the family of
% direct closed-form models fit from the dedicated single-axis energy
% sweeps:
%   Fz    = k*V^2/G_eff^2
%   Taux  = +Leqx*Fz,  Leqx = a0 + a1*rz + a2*thx + a3*thy
%   Tauy  = -Leqy*Fz,  Leqy = a0 + a1*rz + a2*thy + a3*thx  (shared a0-a3)
%   Fx    = (g0z/G_eff)*(C0 + h_own(rx) + h_thz(thz))
%   Fy    = (g0z/G_eff)*(C0 + h_own(ry) - h_thz(thz))
%   Tauz  = (g0z/G_eff)*(g_rx(rx) + g_rx(ry) + g_thz(thz))
% G_eff = g0z - a5*rz - a1_fz*thx - a3_fz*thy (fitted from the Fz sweep)
%
% Inputs (vectorized, same size or scalar-expandable):
%   rx,ry,rz     [m]   plate translation
%   thx,thy,thz  [deg] plate rotation (roll/pitch/yaw)
%   V            [V]   electrode voltage
% Outputs: Fx,Fy,Fz [N], Taux,Tauy,Tauz [N*m]
%
% NOTE: Fx, Fy, Tauz are fit from data with substantial noise (see the
% project report) -- their self-fit R^2 (0.97-0.99, after adding the
% shared baseline C0 for Fx/Fy) reflects fit quality against noisy
% single-axis sweep data, not a full validated cross-check against MST
% across all 6 DOF the way Fz/Taux/Tauy are.

persistent loaded pfitFz tauxTauyCoeffs tauzS FxFyS

if isempty(loaded)
    scriptDir = fileparts(mfilename('fullpath'));
    S1 = load(fullfile(scriptDir, 'simpleFzModel.mat'), 'pfit');
    pfitFz = S1.pfit;   % [k, a1(thx), a3(thy), a5(rz scale)]
    S2 = load(fullfile(scriptDir, 'tauxTauyMomentArm.mat'), 'a0','a1','a2','a3');
    tauxTauyCoeffs = [S2.a0, S2.a1, S2.a2, S2.a3];
    tauzS = load(fullfile(scriptDir, 'TauzModel_kV2GeffSum.mat'), 'C0','g_rx','g_thz','g0_z');
    FxFyS = load(fullfile(scriptDir, 'FxFyModel_kV2GeffSum.mat'), 'C0','h_own','h_thz','g0_z');
    loaded = true;
end

g0_z = pfitFz_g0z_const();   % 15e-6, kept as a local constant for clarity
k_fz = pfitFz(1); a1_fz = pfitFz(2); a3_fz = pfitFz(3); a5_fz = pfitFz(4);

thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;
[rx,ry,rz,thx_rad,thy_rad,thz_rad,V] = deal_expand(rx,ry,rz,thx_rad,thy_rad,thz_rad,V);

G_eff = max(g0_z - a5_fz*rz - a1_fz*thx_rad - a3_fz*thy_rad, 1e-9);

%% Fz
Fz = k_fz*V.^2 ./ G_eff.^2;

%% Taux, Tauy (moment arm, shared coefficients)
a0 = tauxTauyCoeffs(1); a1 = tauxTauyCoeffs(2); a2 = tauxTauyCoeffs(3); a3 = tauxTauyCoeffs(4);
Leqx = a0 + a1*rz + a2*thx_rad + a3*thy_rad;
Leqy = a0 + a1*rz + a2*thy_rad + a3*thx_rad;
Taux = Leqx .* Fz;
Tauy = -Leqy .* Fz;

%% Fx, Fy (shared baseline C0 + own/thz shape terms)
Fx = (FxFyS.g0_z ./ G_eff) .* (FxFyS.C0 + FxFyS.h_own(rx) + FxFyS.h_thz(thz_rad));
Fy = (FxFyS.g0_z ./ G_eff) .* (FxFyS.C0 + FxFyS.h_own(ry) - FxFyS.h_thz(thz_rad));

%% Tauz (sum of rx,ry,thz shape terms)
Tauz = (tauzS.g0_z ./ G_eff) .* (tauzS.C0 + tauzS.g_rx(rx) + tauzS.g_rx(ry) + tauzS.g_thz(thz_rad));

end

function g0_z = pfitFz_g0z_const()
g0_z = 15e-6;
end

function varargout = deal_expand(varargin)
%DEAL_EXPAND Broadcasts scalar/vector inputs to a common column-vector size.
n = 1;
for i = 1:numel(varargin)
    n = max(n, numel(varargin{i}));
end
varargout = cell(1, numel(varargin));
for i = 1:numel(varargin)
    v = varargin{i}(:);
    if numel(v) == 1
        v = repmat(v, n, 1);
    end
    varargout{i} = v;
end
end
