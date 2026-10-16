function [C, A, Iyy, Ixx] = shapeCap(S, G, thx, thy, eps0)
%SHAPECAP Capacitance (and area / second moments) of a signed trapezoid list.
%
%   S    N-by-7 rows [sign, ya, yb, pL, qL, pR, qR] as built by trapGeometry
%   G    base gap g0 - a5*rz  [m]   (vectorized)
%   thx,thy  tilts [rad]            (vectorized)
%   eps0 permittivity; omit or pass 1 to get the bare geometric integral
%
% C   = eps0 * SUM_i sign_i * INT INT dA/(G - y*thx + x*thy)
% A   = plate-projected area
% Iyy = INT y^2 dA,  Ixx = INT x^2 dA
%
% Iyy and Ixx are returned because they give the SMALL-TILT effective
% half-width of the force model independently of any fit: expanding
% Fz = INT eps V^2/(2 g^2) dA to second order gives
%       Fz ~ (k V^2/G^2) * (1 + b^2 th^2/G^2)   with   b^2 = 3*Iyy/A,
% so a shape's own geometry predicts the b that the force fit recovers.

if nargin < 5, eps0 = 1; end

C = zeros(size(G(:))); A = 0; Iyy = 0; Ixx = 0;

for i = 1:size(S,1)
    sg = S(i,1); ya = S(i,2); yb = S(i,3);
    pL = S(i,4); qL = S(i,5); pR = S(i,6); qR = S(i,7);

    C = C + sg * capTrapz(ya, yb, pL, qL, pR, qR, G, thx, thy);

    % Area and moments of the trapezoid, in closed form.
    % width w(y) = (pR-pL) + (qR-qL)*y  =  dp + dq*y
    dp = pR - pL; dq = qR - qL;
    m0 = dp*(yb-ya)       + dq*(yb^2-ya^2)/2;          % INT w dy
    m2 = dp*(yb^3-ya^3)/3 + dq*(yb^4-ya^4)/4;          % INT y^2 w dy
    A   = A   + sg*m0;
    Iyy = Iyy + sg*m2;
    % INT x^2 dx from xL to xR, then over y:  (xR^3 - xL^3)/3
    fx  = @(y) ((pR+qR*y).^3 - (pL+qL*y).^3)/3;
    Ixx = Ixx + sg*integral(fx, ya, yb, 'AbsTol',1e-30,'RelTol',1e-12);
end

C = eps0 * C;
end
