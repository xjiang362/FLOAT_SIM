function geo = trapDeviceGeometry(varargin)
%TRAPDEVICEGEOMETRY Single source of truth for the trap device geometry.
%
%   geo = trapDeviceGeometry()
%   geo = trapDeviceGeometry('z0', 30e-6)     % widen the vertical gap
%   geo = trapDeviceGeometry('xy0', 60e-6)    % widen the lateral gap
%
% CHANGING THE GAP. Pass 'z0' rather than editing the returned struct. The
% equilibrium gap propagates into the channel height, the travel limits,
% all three squeeze-film damping coefficients, the proximity enhancement
% factor and the frame bottom-plate gap, and patching geo.gap.z0 afterwards
% silently leaves those stale. Everything below derives from the parameter,
% so one call keeps the device self-consistent.
%
% Note the unstable rate follows the gap as p = sqrt(2g/z0): 182 Hz at the
% nominal 15 um, 129 Hz at 30 um. Widening the gap therefore buys delay
% margin at the cost of drive voltage (~z0) and position sensitivity
% (~1/z0^2). What it does NOT update is geo.ceff (C12, Cg, Cf), which were
% extracted from COMSOL at z0 = 15 um; the model stays internally
% consistent at any gap but those constants are uncalibrated elsewhere.
%
% Returns one struct describing the levitated mass, the channel it floats
% in, and all 22 capacitance-matrix terminals with the footprint data the
% capacitance and force models need. Everything is SI (metres, kilograms);
% the geometry was supplied in micrometres and is converted here once.
%
% WHY A STRUCT RATHER THAN SCATTERED CONSTANTS: the terminal numbering,
% the quadrant ordering and the Cp orientation were all established
% empirically against the COMSOL capacitance matrix, and getting any of
% them wrong is silent -- the model still runs, it just fits badly. Fixing
% them in one place makes them checkable (see verifyTrapDeviceGeometry.m).
%
% ---------------------------------------------------------------------
% FIELDS
%
%   geo.mass        levitated proof mass
%     .L,.W,.t      [m]      length (x), width (y), thickness (z)
%     .halfL,.halfW,.halfT   [m] half-dimensions
%     .rho          [kg/m^3] density (silicon)
%     .m            [kg]     mass
%     .I,.Iinv      [kg m^2] inertia tensor about the centroid, and inverse
%
%   geo.gap         the channel the mass floats in
%     .z0           [m] vertical gap at equilibrium (mass centred)
%     .zChannel     [m] total vertical channel height = 2*z0
%     .zTravel      [m] max |rz| before contact = z0
%     .xy0          [m] lateral gap at equilibrium
%     .xyChannel    [m] total lateral channel width = 2*xy0
%     .xyTravel     [m] max |rx|,|ry| before contact = xy0
%
%   geo.nTerm       [-] 22
%   geo.idx         terminal index groups: .sensing .actuation .central
%                   .electrodes (1..20) .resonator (21) .frame (22)
%
%   geo.term(i)     one entry per terminal i = 1..22
%     .name         char    e.g. 'C_sns1+', 'Cp-', 'resonator', 'frame'
%     .type         char    'sensing' | 'actuation' | 'central' |
%                           'resonator' | 'frame'
%     .shape        char    'triangle' | 'rectangle' | 'hollow' | 'body'
%     .quadrant     [-]     1..4 for terminals 1..16, NaN otherwise
%     .poly         cell    signed polygon list {s_k, V_k}: V_k is an
%                           N-by-2 vertex array [m], s_k = +1 adds and
%                           s_k = -1 subtracts that region. Hollow shapes
%                           are an outer polygon minus an inner one.
%     .A            [m^2]   NET footprint area (outer minus subtracted)
%     .P            [m]     NET perimeter, including hole edges
%     .Xc,.Yc       [m]     NET area centroid, plate-centred coordinates
%     .dX,.dY       [m]     dY is the full y-extent; dX is the width at
%                           mid-height, chosen so dX*dY = A exactly for
%                           every shape, triangles included. NOTE dX is
%                           NOT the width at Yc: a triangle is 2/3 of its
%                           base there, which would overstate A by 33%.
%
%   geo.frame       decomposition of the resonator<->frame coupling into
%                   five parallel plate capacitors
%     .side(1..4)   .name .normal .faceAt [m] .len [m] .height [m] .A [m^2]
%                   .gap [m]  (lateral)
%     .bottom       .name .A [m^2] .gap [m]  (vertical)
%
%   geo.fringe      analytic fringing, no fitted constants
%     .deltaOverGap [-]  2/pi; each edge extends outward by delta = this
%                        times the local gap, so A_eff = A + P*delta
%
%   geo.ceff        constants of the compact sensed-capacitance model,
%                   shared by every sensing electrode by symmetry
%     .C12,.C1g,.C1f  [F]
%
% ---------------------------------------------------------------------
% CONVENTIONS THAT WERE VERIFIED, NOT ASSUMED
%
%   Tilt sign:    the local gap is g = g0 - rz - thx*Yc + thy*Xc, i.e.
%                 +thx closes the gap at +y and +thy opens it at +x.
%   Quadrants:    terminals advance by a 90-degree ROTATION, so quadrants
%                 2 and 4 take their pair in the opposite order to 1 and 3.
%   Cp pair:      terminal 17 (Cp+) is the -x half, 18 (Cp-) the +x half.
%   TIA pair:     terminal 19 (C_TIA+) is the +y half, 20 the -y half.
%
% See also TRAPGEOMETRY (trapezoid-row form used by capTrapz),
% VERIFYTRAPDEVICEGEOMETRY.

um = 1e-6;

%% ===================== physical constants =====================
geo.phys.eps0 = 8.8541878128e-12;   % [F/m]     vacuum permittivity
geo.phys.g    = 9.81;               % [m/s^2]   gravity, world -z

%% ===================== levitated mass =====================
geo.mass.L   = 5950*um;
geo.mass.W   = 5950*um;
geo.mass.t   = 350*um;
geo.mass.rho = 2330;                     % silicon
geo.mass.halfL = geo.mass.L/2;
geo.mass.halfW = geo.mass.W/2;
geo.mass.halfT = geo.mass.t/2;
geo.mass.m   = geo.mass.rho*geo.mass.L*geo.mass.W*geo.mass.t;
Ixx = geo.mass.m/12*(geo.mass.W^2 + geo.mass.t^2);
Iyy = geo.mass.m/12*(geo.mass.L^2 + geo.mass.t^2);
Izz = geo.mass.m/12*(geo.mass.L^2 + geo.mass.W^2);
geo.mass.I    = diag([Ixx Iyy Izz]);
geo.mass.Iinv = diag(1./[Ixx Iyy Izz]);

%% ===================== channel / gaps =====================
ipG = inputParser;
ipG.addParameter('z0',  15*um, @(v) isscalar(v) && v > 0);
ipG.addParameter('xy0', 50*um, @(v) isscalar(v) && v > 0);
ipG.parse(varargin{:});

geo.gap.z0        = ipG.Results.z0;      % equilibrium, mass centred
geo.gap.zChannel  = 2*geo.gap.z0;        % total height
geo.gap.zTravel   = geo.gap.z0;          % |rz| to contact
geo.gap.xy0       = ipG.Results.xy0;     % equilibrium lateral
geo.gap.xyChannel = 2*geo.gap.xy0;       % total lateral width
geo.gap.xyTravel  = geo.gap.xy0;         % |rx|,|ry| to contact

%% ===================== gas damping (TUNABLE) =====================
% Squeeze-film damping in the free-molecular regime. At 0.1 mTorr the mean
% free path is far larger than the gap, so this is proximity-enhanced
% free-molecular damping rather than viscous flow:
%
%       c = 4*a^3*P / (pi^2 * vbar * g)        translational [N*s/m]
%       q = (a^2/12) * c                       rotational  [N*m*s/rad]
%
% with a the plate edge length, P the ambient pressure, vbar the mean
% molecular speed and g the nominal gap. The a/g factor is the proximity
% enhancement: a plate one gap-length away damps far harder than the same
% plate in free space.
%
% PHYSICAL SCOPE. The squeeze mechanism needs motion that COMPRESSES the
% film, so it describes rz, thx and thy. For rx, ry and thz the plate
% slides parallel to the surface and the gas is sheared, not squeezed --
% a different and much weaker mechanism, smaller by roughly the same
% 4a/(pi^2*g) = 161 that enhances the squeeze terms. The coefficients are
% applied to all six axes as requested, so the three lateral values are
% an UPPER BOUND; cShear below is the order-of-magnitude alternative.
geo.env.pressureTorr = 0.1e-3;              % [Torr]
geo.env.pressure     = 0.1e-3*133.322;      % [Pa]
geo.env.vbar         = 476;                 % [m/s] mean molecular speed, 300 K
geo.env.temperature  = 300;                 % [K]

geo.damp.c = 4*geo.mass.L^3*geo.env.pressure / ...
             (pi^2 * geo.env.vbar * geo.gap.z0);        % [N*s/m]
geo.damp.q = geo.mass.L^2/12 * geo.damp.c;              % [N*m*s/rad]
geo.damp.C = diag([geo.damp.c, geo.damp.c, geo.damp.c, ...
                   geo.damp.q, geo.damp.q, geo.damp.q]);
geo.damp.cShear = geo.damp.c * (pi^2*geo.gap.z0)/(4*geo.mass.L);  % [N*s/m]
geo.damp.enhancement = 4*geo.mass.L/(pi^2*geo.gap.z0);            % [-] ~161

%% ===================== readout timing (TUNABLE) =====================
% Capacitance acquisition schedule. These are the knobs to turn when
% exploring how much sensing delay the loop can survive -- change them
% here and every sampler picks up the new values.
%
% The four quadrant pairs are converted ONE AT A TIME, dtChannel apart,
% then the whole set is shipped over the bus taking tI2C. Channel k of n
% therefore reaches the controller already stale by
%       age(k) = tI2C + dtChannel*(n - k)
% so with the defaults channel 1 is 1050 us old and channel 4 is 600 us
% old, a skew of 450 us between them.
%
% WHY THIS MATTERS: the plant is open-loop unstable at p ~ 1144 rad/s, and
% a plant with an unstable pole p and loop delay tau stops being
% stabilisable with useful margin once p*tau approaches 1. At the default
% 1050 us, p*tau = 1.20 -- already past it. The budget at this pole is
% about 874 us, so one of these two numbers has to come down.
geo.readout.dtChannel = 180e-6;    % [s] one conversion slot = Tisr/4
geo.readout.tI2C      = 500e-6;    % [s] bus transfer after the last conversion
geo.readout.Tisr      = 720e-6;    % [s] controller interrupt period

%% ===================== simulation settings =====================
% Fixed step for the forward-Euler integration of the continuous plant
% between ISR ticks. Must divide Tisr evenly: 720 us / 10 us = 72 steps.
geo.sim.dtEuler = 1e-5;             % [s]

%% ===================== footprint primitives (quadrant 1) =====================
rectV = @(x1,x2,y1,y2) [x1 y1; x2 y1; x2 y2; x1 y2]*um;

P1.sensPlus  = [1560 1525; 2975 1525; 2975 2940]*um;   % x>y half
P1.sensMinus = [1525 1560; 1525 2975; 2940 2975]*um;   % x<y half
s = 1450/2;
P1.actPlus   = rectV(2250-s,2250+s, 750-s, 750+s);     % x>y half
P1.actMinus  = rectV( 750-s, 750+s,2250-s,2250+s);     % x<y half
q1poly = {P1.sensPlus, P1.sensMinus, P1.actPlus, P1.actMinus};
q1type = {'sensing','sensing','actuation','actuation'};
q1shape= {'triangle','triangle','rectangle','rectangle'};
q1suf  = {'sns','sns','act','act'};
q1pm   = {'+','-','+','-'};

% Quadrant maps. Terminals advance by a 90-degree ROTATION, so the sign
% flip alone is not enough: quadrants 2 and 4 also take the pair in the
% opposite order. Verified -- the quadrant-1 ordering fits quadrants 1,3
% at R^2 = 0.996 but quadrants 2,4 at only 0.927, and swapping reverses it.
sx = [1 -1 -1  1];
sy = [1  1 -1 -1];
pairOrder = {[1 2 3 4],[2 1 4 3],[1 2 3 4],[2 1 4 3]};

T = struct('name',{},'type',{},'shape',{},'quadrant',{},'poly',{}, ...
           'A',{},'P',{},'Xc',{},'Yc',{},'dX',{},'dY',{}, ...
           'truncX',{},'truncY',{});
for q = 1:4
    for j = 1:4
        src = pairOrder{q}(j);
        V = q1poly{src};
        V = [sx(q)*V(:,1), sy(q)*V(:,2)];
        e = makeTerm(sprintf('C_%s%d%s', q1suf{j}, q, q1pm{j}), ...
                     q1type{j}, q1shape{j}, q, {1, V});
        T(4*(q-1)+j) = e;
    end
end

%% ===================== central pairs =====================
% Cp (larger inner pair): outer rectangle minus the well housing the TIA.
% Terminal 17 (Cp+) is the -x half and 18 (Cp-) the +x half -- verified
% from the coupling pattern, NOT inferred from the sign of Xc.
CpOut = rectV(  75, 1475, -1475, 1475);
CpIn  = rectV(  75,  737,  -725,  725);
mirXV = @(V) [-V(:,1), V(:,2)];
mirYV = @(V) [ V(:,1),-V(:,2)];
T(17) = makeTerm('Cp+','central','hollow',NaN,{1,mirXV(CpOut); -1,mirXV(CpIn)});
T(18) = makeTerm('Cp-','central','hollow',NaN,{1,CpOut;        -1,CpIn});

% TIA (smaller inner pair): outer rectangle minus a notch.
% Terminal 19 is the +y half, 20 the -y half -- also verified.
TiaOut = rectV(-700, 700,  25, 675);
TiaIn  = rectV(-100, 100,  25, 100);
T(19) = makeTerm('C_TIA+','central','hollow',NaN,{1,TiaOut;        -1,TiaIn});
T(20) = makeTerm('C_TIA-','central','hollow',NaN,{1,mirYV(TiaOut); -1,mirYV(TiaIn)});

%% ===================== resonator and frame =====================
T(21) = makeTerm('resonator','resonator','body',NaN, ...
    {1, rectV(-2975,2975,-2975,2975)});
T(22) = makeTerm('frame','frame','body',NaN, {});

geo.term  = T;
geo.nTerm = 22;
geo.idx.sensing    = find(strcmp({T.type},'sensing'));
geo.idx.actuation  = find(strcmp({T.type},'actuation'));
geo.idx.central    = find(strcmp({T.type},'central'));
geo.idx.electrodes = 1:20;
geo.idx.resonator  = 21;
geo.idx.frame      = 22;

% The two sensing electrodes of each quadrant form one differential
% readout pair; C_eff is measured across the pair. Row q = quadrant q.
geo.idx.sensePair = [1 2; 5 6; 9 10; 13 14];

% Actuation pairs, one control voltage each: the four quadrant pairs plus
% the two central pairs. Each pair is driven differentially (+V and -V),
% so both members see the same V^2 and contribute with the same sign to
% the electrostatic wrench.
geo.idx.actPair = [3 4; 7 8; 11 12; 15 16; 17 18; 19 20];
geo.idx.actPairName = {'quad1','quad2','quad3','quad4','Cp','TIA'};

%% ===================== frame decomposition =====================
% Resonator<->frame coupling as four side faces plus the bottom face.
% The parameter-free prediction from these areas and gaps gives 22.37 pF
% against a measured 22.18 pF at equilibrium (0.9%), and the bottom face
% alone is 93% of the total.
sideNames  = {'Crs1','Crs2','Crs3','Crs4'};
sideNormal = {'+x','+y','-x','-y'};
faceAt     = [ geo.mass.halfL, geo.mass.halfW, -geo.mass.halfL, -geo.mass.halfW];
for kS = 1:4
    geo.frame.side(kS).name   = sideNames{kS};
    geo.frame.side(kS).normal = sideNormal{kS};
    geo.frame.side(kS).faceAt = faceAt(kS);
    if mod(kS,2)==1, geo.frame.side(kS).len = geo.mass.W;
    else,            geo.frame.side(kS).len = geo.mass.L; end
    geo.frame.side(kS).height = geo.mass.t;
    geo.frame.side(kS).A      = geo.frame.side(kS).len * geo.mass.t;
    geo.frame.side(kS).gap    = geo.gap.xy0;
end
geo.frame.bottom.name = 'Crb1';
geo.frame.bottom.A    = geo.mass.L*geo.mass.W;
geo.frame.bottom.gap  = geo.gap.z0;

%% ===================== fitted model constants =====================
% FRINGING -- analytic, no fitted constants. Each electrode edge behaves
% as if extended outward by
%       delta = (2/pi) * g
% the conformal-mapping result for a parallel-plate edge, so the effective
% area is A_eff = A + P*delta. Tested against the 63-pose sweep over all
% 20 terminals with nothing fitted: median +0.05%, RMS 1.83%, and +0.20%
% at equilibrium -- slightly BETTER than the two-constant k = a + b*(P/A)
% fit it replaces (RMS 1.89%). Alternatives scored: no fringing 4.31%,
% g/pi 2.64%, g/2 1.98%, g 3.01%.
geo.fringe.deltaOverGap = 2/pi;

% The same length scale smooths the overlap cut: the edge the plate slides
% past is fringed over exactly the distance the field bulges, so there is
% no separate smear parameter to fit. Being proportional to the gap, it
% also tightens correctly as the plate descends.


% Sensed-capacitance model, shared across all four quadrant pairs by
% symmetry (the eight electrodes' fitted gains agree to 0.4%). These
% constants apply to the form used in trapSensing.m,
%
%   Ceff = C12 + Cf + (Cg+C_iR)(Cg+C_jR)/(Cg+C_iR + Cg+C_jR)
%
% where C_iR is the UN-FRINGED eps0*A/g produced by trapCapacitance.
% Cf is the per-electrode fringing allowance: the capacitance stage is
% pure parallel plate and under-predicts by ~4%, and Cf restores that at
% the sensing output. It sits inside the series alongside Cg, so the
% familiar closed form is Ceff = C12 + Cg/2 + Cf/2 + C_iR/2.
% If fringing is ever reinstated upstream, Cf drops to ~0.3 fF.
geo.ceff.C12 = 50.6203e-15;
geo.ceff.Cg  = 105.1301e-15;
geo.ceff.Cf  = 31.794e-15;

end

% =====================================================================
function e = makeTerm(name,type,shape,quad,poly)
%MAKETERM Build one terminal entry, deriving A, P, centroid and extents
%from the signed polygon list. Empty poly (the frame) leaves them NaN.
e.name = name; e.type = type; e.shape = shape; e.quadrant = quad;
e.poly = poly;
nullTrunc = struct('side',0,'edge',0,'len',0);
if isempty(poly)
    [e.A,e.P,e.Xc,e.Yc,e.dX,e.dY] = deal(NaN);
    e.truncX = nullTrunc; e.truncY = nullTrunc;
    return
end
A = 0; Mx = 0; My = 0; P = 0; ys = [];
for k = 1:size(poly,1)
    sgn = poly{k,1}; V = poly{k,2};
    a  = polyarea(V(:,1),V(:,2));
    cx = mean(V(:,1)); cy = mean(V(:,2));      % exact for triangles
    if size(V,1) == 4                           % rectangle
        cx = (min(V(:,1))+max(V(:,1)))/2;
        cy = (min(V(:,2))+max(V(:,2)))/2;
    end
    A  = A  + sgn*a;
    Mx = Mx + sgn*a*cx;
    My = My + sgn*a*cy;
    P  = P  + sum(sqrt(sum(diff(V([1:end 1],:)).^2,2)));   % hole edges add
    ys = [ys; V(:,2)]; %#ok<AGROW>
end
e.A  = A;  e.P = P;
e.Xc = Mx/A;  e.Yc = My/A;
e.dY = max(ys) - min(ys);
e.dX = A/e.dY;          % width at mid-height, so dX*dY = A exactly

% ---- truncation edges, for the lateral (overlap) terms ----
% An electrode whose footprint reaches the proof-mass boundary loses area
% when the plate slides away, and that is the ONLY source of Fx, Fy and
% tauz in this model. For each axis record which side can truncate, where
% that edge sits, and how long it is -- the edge length is dA/d(edge
% position) to first order, which is what the force derivative needs.
%
% A triangle's x-leg is full height (1415 um for the sensing electrodes),
% so x-truncation removes a full-height strip and dX*dY = A is preserved.
% Its y-extremum is a VERTEX, not an edge, so the y-edge length is ~0 and
% y-truncation removes a vanishing area to first order -- correctly
% captured by len = 0 rather than by a bounding-box approximation.
a = 2975e-6;                       % proof-mass half-width
tol = 1e-9;
V = vertcat(e.poly{:,2});
e.truncX = edgeAt(e.poly, 1, a, tol);
e.truncY = edgeAt(e.poly, 2, a, tol);
end

% =====================================================================
function t = edgeAt(poly, dim, a, tol)
%EDGEAT Which side of axis DIM reaches the plate boundary, and the total
%length of polygon edges lying on it. dim = 1 for x, 2 for y.
t = struct('side',0,'edge',0,'len',0);
other = 3 - dim;
for sgn = [1 -1]
    L = 0;
    for k = 1:size(poly,1)
        if poly{k,1} < 0, continue; end          % holes cannot touch the rim
        Vk = poly{k,2};
        Vc = Vk([1:end 1],:);
        for s = 1:size(Vc,1)-1
            p1 = Vc(s,:); p2 = Vc(s+1,:);
            onRim = abs(p1(dim) - sgn*a) < tol && abs(p2(dim) - sgn*a) < tol;
            if onRim
                L = L + abs(p2(other) - p1(other));
            end
        end
    end
    if L > tol
        t.side = sgn; t.edge = a; t.len = L; return
    end
end
end








