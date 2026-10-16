function [Fx,Fy,Fz,Taux,Tauy,Tauz] = trapElectrode1Model_crossTerms(rx,ry,rz,thx,thy,thz,U)
%TRAPELECTRODE1MODEL_CROSSTERMS Same as trapElectrode1Model.m, but with
% the same-axis quadratic correction terms (thx^2, thy^2, thz^2, rx^2,
% ry^2) replaced by cross terms (rz*thx, rz*thy, rx*thz, ry*thz):
%   G_eff = g0_z - a5*rz - a1*thx - a2*(rz*thx) - a3*thy - a4*(rz*thy)
%   W_eff = W0 - b2*rx - b1*thz - b3*(rx*thz)
%   L_eff = L0 - c2*ry - c1*thz - c3*(ry*thz)
% Physical motivation: tilting the plate shifts the effective gap/width/
% length by an amount that scales with how far the plate already is in
% the coupled direction, which a cross term captures more directly than
% a same-axis quadratic term.
%
%   [Fx,Fy,Fz,Taux,Tauy,Tauz] = trapElectrode1Model_crossTerms(rx,ry,rz,thx,thy,thz,U)
%   Inputs/outputs: identical convention to trapElectrode1Model.m.

g0_z = 15e-6; W0 = 6000e-6; L0 = 6000e-6;
eps0 = 8.854e-12;
floorVal = 1e-8;

% Fitted coefficients (fitFringingCapacitanceModel.m, cross-term version)
k_fit = 59.8431; U_fit = 30;
k = k_fit / U_fit^2;
a1 = 0.00142034; a2 = -4.59845; a3 = -0.00148876; a4 = 13.7869;
b1 = -0.000409518; b2 = -1.43177; b3 = 170.455;
c1 = -0.00308286;  c2 = -0.618676; c3 = 60.1574;
a5 = 0.968282;

Cf = @(G,Wd,Ld) (eps0*Wd.*Ld./G) .* ...
     (1 + G./(pi*Wd) + G./(pi*Wd).*log(2*pi*Wd./G)) .* ...
     (1 + G./(pi*Ld) + G./(pi*Ld).*log(2*pi*Ld./G));

thx_rad = thx*pi/180; thy_rad = thy*pi/180; thz_rad = thz*pi/180;
[rx,ry,rz,thx_rad,thy_rad,thz_rad,U] = deal_expand(rx,ry,rz,thx_rad,thy_rad,thz_rad,U);
X = [rx, ry, rz, thx_rad, thy_rad, thz_rad, U];

Wfun = @(X) X(:,7).^2 .* k .* Cf( ...
    max(g0_z - a5*X(:,3) - a1*X(:,4) - a2*(X(:,3).*X(:,4)) - a3*X(:,5) - a4*(X(:,3).*X(:,5)), floorVal), ...
    max(W0 - b1*X(:,6) - b2*X(:,1) - b3*(X(:,1).*X(:,6)), floorVal), ...
    max(L0 - c1*X(:,6) - c2*X(:,2) - c3*(X(:,2).*X(:,6)), floorVal) );

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
