function [C, dCdq, info] = trapCapacitance(plant, pose)
%TRAPCAPACITANCE Electrode-to-resonator capacitances and their Jacobian.
%
%   [C, dCdq, info] = trapCapacitance(plant, pose)
%
% Capacitance stage of the levitation pipeline. Returns both the
% capacitances and their pose derivatives, because the force stage needs
% dC/dq (the energy method gives F_q = 0.5*sum_i V_i^2 dC_i/dq) and
% numerically differencing C afterwards would reintroduce exactly the
% noise this model exists to avoid.
%
% INPUTS
%   plant   struct from trapDeviceGeometry()
%   pose    6-by-1 [rx; ry; rz; thx; thy; thz], SI: metres and RADIANS
%
% OUTPUTS
%   C       [F]      20-by-1 terminal-to-resonator capacitance, fringing
%                    corrected
%   dCdq    [F/m, F/rad] 20-by-6, column j is dC/dpose(j)
%   info    struct: .gap [m] 20-by-1 local gap at each centroid
%                   .area [m^2] 20-by-1 overlapped area
%                   .shade [-] 20-by-1 fraction of footprint still covered
%                   .minGap [m] smallest gap over all electrodes
%
% MODEL -- pure parallel plate. No fitted constants and NO fringing term.
%   C_i = eps0 * A_i(pose) / g_i(pose)
%
%   g_i = g0 - rz - thx*Yc_i + thy*Xc_i            local gap at centroid
%   A_i = A_i0 - Lx_i*ox_i - Ly_i*oy_i             overlapped area
%
% FRINGING IS DELIBERATELY NOT MODELLED. The field does bulge past each
% electrode edge, so this form UNDER-PREDICTS the capacitance by a known,
% systematic amount: median -4.1%, RMS 4.3%, -3.4% at equilibrium, across
% all 20 terminals and 60 poses. The deficit is largest for electrodes
% with a high perimeter-to-area ratio (the TIA pair, -6%) and smallest for
% the big actuation squares (-1.5%), because fringing is an edge effect.
%
% This is a modelling choice, not an oversight: the bias is a smooth
% multiplicative offset that a downstream gain calibration can absorb,
% whereas a fringing correction would put a non-geometric constant into
% what is otherwise a pure geometry model. If the absolute level matters,
% the conformal-mapping edge extension delta = (2/pi)*g applied as
% Aeff = A + P*delta removes nearly all of it (RMS 1.69%) and is still
% parameter-free -- see the project history before reinstating it.
%
% The area term is what produces lateral force. An electrode whose
% footprint is flush with the proof-mass boundary loses area as the plate
% slides away; dA/d(edge position) is just the length of the electrode
% edge lying on that boundary, which trapDeviceGeometry precomputes as
% truncX.len / truncY.len. Electrodes that do not reach the boundary (the
% Cp and TIA pairs) have zero edge length and therefore no lateral force.
%
% YAW enters through the same overlap. Rotating the plate by thz displaces
% its edge by -Yc*thz at the x boundary and +Xc*thz at the y boundary, so
%   rx_eff_i = rx - Yc_i*thz,    ry_eff_i = ry + Xc_i*thz
% and the chain rule then gives, with no further assumption,
%   dC/dthz = -Yc*dC/drx + Xc*dC/dry
% i.e. tauz is the moment of the lateral force about the electrode
% centroid. That the lever-arm form falls out of a pure overlap
% derivative is a consistency check on the construction.
%
% SMOOTHING. A hard overlap gives a STEP in dA/drx, so Fx would be
% discontinuous at the transition and the integrator would chatter. The
% cut is smoothed over plant.fringe.w (~25 um, of order the gap, which is
% the length scale fringing smears a conductor edge over), using
%   smax0(z,w) = z * (1 + tanh(z/w))/2
% which has every limit the problem needs:
%   smax0(0)    = 0        the flush pose loses NO area
%   d/dz at 0   = 1/2      the cut advances at half rate at the transition
%   z -> -inf   -> 0       fully covered, no correction
%   z -> +inf   -> z       fully exposed, hard geometric overlap
% The more familiar (z + w*log(2cosh(z/w)))/2 has the right asymptotes but
% equals w*ln2/2 at z = 0, i.e. it removes ~8.7 um of edge even when the
% plate is exactly flush -- a 0.6% area error at the operating point, and
% a bias the fitted k_i gains would then double-count.
%
% VALIDATION STATUS
%   Vertical and tilt channels are validated: the C_i model reproduces the
%   measured capacitance for all 20 terminals at R^2 >= 0.993 (see
%   verifyTrapDeviceGeometry.m).
%   The LATERAL channel is NOT validated against simulation and must not
%   be: all Fx/Fy/tauz simulation data in this project is corrupted by
%   meshing. The lateral terms rest on the geometry derivation above.
%
% See also TRAPDEVICEGEOMETRY, TRAPSENSING.

eps0 = plant.phys.eps0;

pose = pose(:);
if numel(pose) ~= 6
    error('trapCapacitance:pose','pose must be 6-by-1 [rx;ry;rz;thx;thy;thz].');
end
rx = pose(1); ry = pose(2); rz = pose(3);
thx = pose(4); thy = pose(5); thz = pose(6);

idx = plant.idx.electrodes;
n   = numel(idx);
g0  = plant.gap.z0;
dOg = plant.fringe.deltaOverGap;        % 2/pi, analytic

C    = zeros(n,1);
dCdq = zeros(n,6);
gap  = zeros(n,1);
area = zeros(n,1);
shade= zeros(n,1);

for jj = 1:n
    e  = plant.term(idx(jj));
    Xc = e.Xc;  Yc = e.Yc;  A0 = e.A;

    % ---- local gap at the centroid ----
    g = g0 - rz - thx*Yc + thy*Xc;
    if g <= 0
        error('trapCapacitance:contact', ...
            'Electrode %d (%s) has non-positive gap %.3g m; pose is in contact.', ...
            idx(jj), e.name, g);
    end

    % ---- overlapped area ----
    % The electrode edge sits exactly on the plate rim, so the overhang is
    % simply -side*r_eff: positive when the plate has slid away from it.
    rxe = rx - Yc*thz;
    rye = ry + Xc*thz;

    % Smoothing width for the overlap cut. This is NUMERICAL
    % REGULARISATION, not a fringing model: a hard geometric overlap makes
    % dA/drx a step, so Fx would be discontinuous and the integrator would
    % chatter. The gap is the only natural length scale at the edge, and
    % w -> 0 recovers the exact hard overlap. Because w tracks g, its own
    % derivative contributes to dA/drz, dA/dthx and dA/dthy -- small, but
    % omitting it leaves the Jacobian inconsistent with the returned C.
    w = dOg*g;
    [ox, hx, oxw] = smax0(-e.truncX.side*rxe, w);
    [oy, hy, oyw] = smax0(-e.truncY.side*rye, w);
    Lx = e.truncX.len;  Ly = e.truncY.len;
    A  = A0 - Lx*ox - Ly*oy;
    if A <= 0
        error('trapCapacitance:noOverlap', ...
            'Electrode %d (%s) has no remaining overlap; pose is out of range.', ...
            idx(jj), e.name);
    end

    % dA/drx and dA/dry: the edge length times how fast the cut advances.
    % side = 0 (no boundary edge) makes both vanish, as it must.
    dAdrx = Lx * e.truncX.side * hx;
    dAdry = Ly * e.truncY.side * hy;

    % ---- capacitance and derivatives ----
    Ci   = eps0*A/g;
    dCdA = eps0/g;
    dCdg = -eps0*A/g^2;

    dCdrx = dCdA*dAdrx;
    dCdry = dCdA*dAdry;

    % Area change from the smear width tracking the gap, dA/dw * dw/dq.
    dAdw  = -(Lx*oxw + Ly*oyw);
    dCdw  = dCdA*dAdw;

    C(jj) = Ci;
    dCdq(jj,1) = dCdrx;                              % rx
    dCdq(jj,2) = dCdry;                              % ry
    dCdq(jj,3) = -dCdg      - dCdw*dOg;              % rz   (dg/drz  = -1)
    dCdq(jj,4) = -dCdg*Yc   - dCdw*dOg*Yc;           % thx  (dg/dthx = -Yc)
    dCdq(jj,5) =  dCdg*Xc   + dCdw*dOg*Xc;           % thy  (dg/dthy = +Xc)
    dCdq(jj,6) = -Yc*dCdrx + Xc*dCdry;               % thz  (chain rule)

    gap(jj) = g;  area(jj) = A;  shade(jj) = A/A0;
end

if nargout > 2
    info.gap    = gap;
    info.area   = area;
    info.shade  = shade;
    info.minGap = min(gap);
end
end

% =====================================================================
function [v, dv, dvdw] = smax0(z, w)
%SMAX0 Smoothed max(z,0) with width w, plus both partial derivatives.
%   v     = z * s,            s = (1 + tanh(z/w))/2
%   dv/dz = s + z*(1-t^2)/(2w)
%   dv/dw = -z^2*(1-t^2)/(2*w^2)
% Chosen over the log-cosh form because v(0) = 0 exactly: at the flush
% pose the electrode must lose no area. tanh saturates rather than
% overflowing, so no large-argument guard is needed. dv/dw is needed
% because the smear width tracks the gap and so varies with pose.
u    = z/w;
t    = tanh(u);
s    = (1 + t)/2;
sech2= 1 - t.^2;
v    = z.*s;
dv   = s + z.*sech2/(2*w);
dvdw = -z.^2.*sech2/(2*w^2);
end

