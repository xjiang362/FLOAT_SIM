function [Fx,Fy,Fz,Taux,Tauy,Tauz] = trapElectrode1Model(rx,ry,rz,thx,thy,thz,U)
%TRAPELECTRODE1MODEL Force/torque on the plate from the quadrant-1
% electrode, via the fitted fringing-capacitance energy model.
%
%   [Fx,Fy,Fz,Taux,Tauy,Tauz] = trapElectrode1Model(rx,ry,rz,thx,thy,thz,U)
%
%   Inputs (vectorized, same size or scalar-expandable):
%     rx,ry,rz     [m]   plate translation
%     thx,thy,thz  [deg] plate rotation (roll/pitch/yaw)
%     U            [V]   electrode voltage
%   Outputs: Fx,Fy,Fz [N], Taux,Tauy,Tauz [N*m]
%
%   W = k*U^2*Cf(G_eff,W_eff,L_eff), with the 3D fringing-capacitance
%   formula Cf and effective gap/width/length depending on tilt and
%   lateral position -- fit on a 2000-pt full-6D COMSOL LHS sweep at
%   U=30V (k here is that fit's k divided by 30^2, so forces reproduce
%   the original fit exactly at U=30V and scale as U^2 elsewhere).
%   Fx,Fy,Fz,Taux,Tauy,Tauz = +dW/dq (COMSOL's probe sign convention,
%   confirmed empirically against MST for Fz/Taux/Tauy).

g0_z = 15e-6; W0 = 6000e-6; L0 = 6000e-6;
eps0 = 8.854e-12;
floorVal = 1e-8;

% Fitted coefficients (fitFringingCapacitanceModel.m, rz-scale version)
k_fit = 58.828; U_fit = 30;
k = k_fit / U_fit^2;
a1 = 0.00128445; a2 = 0.226905; a3 = -0.00128545; a4 = 0.241927;
b1 = -0.278475;  b2 = -1.28978; b3 = 10177.4;
c1 = 0.0311339;  c2 = -0.65012; c3 = 16472.3;
a5 = 0.971151;

Cf = @(G,Wd,Ld) (eps0*Wd.*Ld./G) .* ...
     (1 + G./(pi*Wd) + G./(pi*Wd).*log(2*pi*Wd./G)) .* ...
     (1 + G./(pi*Ld) + G./(pi*Ld).*log(2*pi*Ld./G));

thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;
[rx,ry,rz,thx_rad,thy_rad,thz_rad,U] = deal_expand(rx,ry,rz,thx_rad,thy_rad,thz_rad,U);
X = [rx, ry, rz, thx_rad, thy_rad, thz_rad, U];

Wfun = @(X) X(:,7).^2 .* k .* Cf( ...
    max(g0_z - a5*X(:,3) - a1*X(:,4) - a2*X(:,4).^2 - a3*X(:,5) - a4*X(:,5).^2, floorVal), ...
    max(W0 - b1*X(:,6).^2 - b2*X(:,1) - b3*X(:,1).^2, floorVal), ...
    max(L0 - c1*X(:,6).^2 - c2*X(:,2) - c3*X(:,2).^2, floorVal) );

% Fixed absolute step sizes (NOT derived from the input batch's own
% range) -- this function must work correctly for a single isolated
% point (e.g. inside an ODE integrator), where a batch-relative step
% would collapse to ~eps and silently zero out the derivative.
% Sized as ~1e-6 of each variable's typical physical range (translation
% +-45um, rz +-14um, tilt up to ~0.13deg, yaw up to ~0.8deg).
n = size(X,1);
h = [1e-11, 1e-11, 1e-11, 2e-9, 2e-9, 1.4e-8];

dWdq = zeros(n,6);
for j = 1:6
    Xp = X; Xp(:,j) = Xp(:,j) + h(j);
    Xm = X; Xm(:,j) = Xm(:,j) - h(j);
    dWdq(:,j) = (Wfun(Xp) - Wfun(Xm)) / (2*h(j));
end

Fx = dWdq(:,1); Fy = dWdq(:,2); Fz = dWdq(:,3);
Taux = dWdq(:,4); Tauy = dWdq(:,5); Tauz = dWdq(:,6);

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
