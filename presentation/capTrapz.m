function C = capTrapz(ya, yb, pL, qL, pR, qR, G, thx, thy)
%CAPTRAPZ Exact closed-form  INT INT dA / g(x,y)  over a y-trapezoid.
%
% The single primitive behind every electrode in the trap device. The
% region is swept in y from ya to yb, bounded by two straight edges that
% are linear in y:
%       xL(y) = pL + qL*y        xR(y) = pR + qR*y
% which covers rectangles (qL=qR=0), right triangles (one q = +-1 and
% the edges meeting at a vertex), and trapezoids alike. Hollow shapes
% (the central electrodes) are signed sums of these.
%
% Local gap under a plate displaced by rz and tilted by (thx,thy):
%       g(x,y) = G - y*thx + x*thy,      G = g0 - a5*rz
% (sign convention established by the quadrant-electrode force fit: +thx
%  closes the gap at +y, +thy opens it at +x).
%
% Capacitance scales as 1/g, NOT 1/g^2 like the electrostatic pressure,
% so these integrals are logarithmic rather than rational. Writing
%       h(u) = u*(log(u) - 1),     h'(u) = log(u)
% the x-integral gives (1/thy)*log(u) and the subsequent y-integral gives
% (1/(B*thy))*h(A + B*y) along each edge, with
%       A_E = G + p_E*thy,   B_E = q_E*thy - thx.
%
% Returns the GEOMETRIC integral only (units m^2/m = m); multiply by the
% permittivity outside. Vectorized over G, thx, thy.
%
% Numerically delicate limits (thy -> 0, B -> 0) are handled by explicit
% series branches rather than by letting the formula cancel catastro-
% phically; verifyCapPrimitive.m checks every branch against integral2.

[G, thx, thy] = deal_same(G, thx, thy);
C = zeros(size(G));

dy   = yb - ya;
xmax = max(abs([pL + qL*ya, pL + qL*yb, pR + qR*ya, pR + qR*yb]));

for i = 1:numel(G)
    g = G(i); tx = thx(i); ty = thy(i);

    if abs(ty) * xmax > 1e-9 * g
        % ---- generic branch: x-integral is logarithmic ----
        AL = g + pL*ty;  BL = qL*ty - tx;
        AR = g + pR*ty;  BR = qR*ty - tx;
        C(i) = ( edgeInt(AR, BR, ya, yb) - edgeInt(AL, BL, ya, yb) ) / ty;
    else
        % ---- thy ~ 0: gap independent of x, width integrates directly ----
        dp = pR - pL;  dq = qR - qL;
        if abs(tx) * max(abs([ya yb])) > 1e-9 * g
            % INT (dp + dq*y)/(g - y*tx) dy, closed form.
            % Same contact guard as edgeInt: clamp the gap away from zero
            % so a plate that has been stepped into the electrode returns
            % a large finite number rather than NaN.
            w  = @(y) max(g - y*tx, UMIN);
            C(i) = -(dp/tx + dq*g/tx^2) * (log(w(yb)) - log(w(ya))) ...
                   + (dq/tx^2) * (w(yb) - w(ya));
        else
            % both tilts negligible: plain area / gap
            C(i) = (dp*dy + dq*(yb^2 - ya^2)/2) / g;
        end
    end
end
end

% ---------------------------------------------------------------------
function val = edgeInt(A, B, ya, yb)
% INT_ya^yb log(A + B*y) dy = (1/B)*[h(A+B*yb) - h(A+B*ya)],  h(u)=u(log u -1)
% As B -> 0 the bracket and 1/B both blow up while the product stays
% finite, so switch to the expansion
%     (yb-ya)*log(A) + B*(yb^2-ya^2)/(2A) - B^2*(yb^3-ya^3)/(6A^2)
% well before the cancellation bites.
if abs(B) * max(abs([ya yb])) < 1e-8 * abs(A)
    val = (yb - ya)*log(max(A,UMIN)) ...
        + B *(yb^2 - ya^2)/(2*A) ...
        - B^2*(yb^3 - ya^3)/(6*A^2);
else
    % h(u) = u(ln u - 1) tends to 0 as u -> 0, so the primitive is
    % analytically well behaved at contact. In floating point, though,
    % u = 0 gives 0*(-Inf) = NaN, and an exact zero gap is reachable:
    % the EOM integrator can step the plate into the electrode, where a
    % single NaN silently poisons the whole state vector. Clamping the
    % gap to UMIN keeps the result large-but-finite, so the integrator
    % sees a stiff wall and backs off instead of failing silently.
    % Capacitance genuinely diverges at a contact line, so a large value
    % here is physics -- only the NaN is spurious.
    h = @(u) u.*(log(u) - 1);
    val = (h(max(A + B*yb, UMIN)) - h(max(A + B*ya, UMIN))) / B;
end
end

% ---------------------------------------------------------------------
function u = UMIN()
% Smallest gap the formulas are evaluated at [m]. 1 pm is far below any
% physically meaningful gap (the nominal is 15.45 um) yet keeps every
% logarithm finite.
u = 1e-12;
end

% ---------------------------------------------------------------------
function varargout = deal_same(varargin)
n = max(cellfun(@numel, varargin));
varargout = cell(1, nargin);
for k = 1:nargin
    v = varargin{k};
    if numel(v) == 1, v = repmat(v, 1, n); end
    varargout{k} = v(:);
end
end
