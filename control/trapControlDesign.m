function K = trapControlDesign(plant, varargin)
%TRAPCONTROLDESIGN Discrete-time state-space regulator with integral action.
%
%   K = trapControlDesign(plant)
%   K = trapControlDesign(plant, 'lambdar',2*pi*50, 'lambdae',2*pi*300)
%
% Designs the stabilising loop for the three controllable degrees of
% freedom -- heave and the two tilts -- from four sensor readings and five
% actuation forces.
%
% PLANT MODEL
%   x  = [rz; drz; thx; dthx; thy; dthy]        6 states
%   u  = [fp; f1; f2; f3; f4]                   5 forces  [N]
%   ys = [s1; s2; s3; s4]                       4 sensed heights [m]
%
%   dx/dt = A*x + B*u,      ys = Cs*x
%
% The input is FORCE, not voltage. That matters: at constant force there
% is no electrostatic stiffness, which is why A carries only damping on
% the velocity rows. The negative stiffness of the real device
% (-37.76 N/m in heave) lives in the force-to-voltage inversion instead,
% and only cancels there if that inversion uses the CURRENT pose. Freeze
% it at nominal and the unmodelled unstable pole comes back.
%
% u(1) is the Cp pair, u(2:5) the four quadrants. The TIA pair is not a
% control input: it carries a fixed bias voltage, so its force is a
% constant that the voltage mapping subtracts from the gravity offset.
%
% INTEGRAL ACTION. Three integrators are appended, one per controlled
% output (rz, thx, thy), giving a 9-state augmented plant:
%       xi[k+1] = xi[k] + T*(r - Cc*x[k])
% They remove steady-state error against the constant force offsets this
% device always carries -- gravity, and any mismatch in the bias voltages.
%
% NAME/VALUE
%   'lambdar'  [rad/s] continuous regulator bandwidth (default 2*pi*50)
%   'lambdae'  [rad/s] continuous estimator bandwidth (default 2*pi*300)
%   'Tisr'     [s]     sample period (default plant.readout.Tisr)
%   'ell'      [m]     sensing-point offset used in Cs (default: taken
%                      from the geometry, 2250 um)
%
% OUTPUT K
%   .A .B .Cs .Cc            continuous plant and controlled-output map
%   .Ad .Bd .Csd             discretised at Tisr (zero-order hold)
%   .Aa .Ba .Br              9-state augmented plant, reference input
%   .Kx .Ki                  regulator gains, u = -Kx*xhat - Ki*xi
%   .L                       estimator gain, 6-by-4
%   .polesR .polesE          closed-loop poles placed
%   .ctrb .obsv              rank reports
%   .Tisr .lambdar .lambdae
%
% Pole placement uses a Butterworth pattern at the requested bandwidth
% rather than a repeated real pole. Repeated poles of multiplicity above
% rank(B) cannot be placed at all, and even where they can the gain matrix
% is badly conditioned; Butterworth gives the same bandwidth with distinct
% poles and a well-conditioned solution.
%
% See also TRAPFORCETOVOLTAGE, TRAPDEVICEGEOMETRY.

ip = inputParser;
% lambdar = 200 Hz, not the 50 Hz originally specified. At 50 Hz the loop
% is unstable at every sample rate tested, with or without the estimator;
% 200 Hz converges from the full heave travel. lambdae = 500 Hz sits at
% half Nyquist for Tisr = 500 us, which is the practical ceiling.
ip.addParameter('lambdar', 2*pi*100);
ip.addParameter('lambdae', 2*pi*300);
ip.addParameter('Tisr', plant.readout.Tisr);
ip.addParameter('ell', []);
ip.addParameter('Vtia', 0);
% 'exact'  : the true sampled map [1 sinh(pT)/p; 0 cosh(pT)] (default)
% 'double' : the zero-order hold of a double integrator [1 T; 0 1], which
%            is what the design originally used. Kept so the two can be
%            compared directly -- it is wrong, not merely approximate.
ip.addParameter('designModel', 'exact');
ip.parse(varargin{:});
lambdar = ip.Results.lambdar;
lambdae = ip.Results.lambdae;
T       = ip.Results.Tisr;
Vtia    = ip.Results.Vtia;

m  = plant.mass.m;
Ix = plant.mass.I(1,1);
Iy = plant.mass.I(2,2);
cz = plant.damp.c;                 % squeeze-film translational [N*s/m]
qx = plant.damp.q;  qy = plant.damp.q;
a  = 1500e-6;                      % [m] moment arm of the quadrant pairs

% Sensing offset: taken from the verified geometry, where the four pair
% midpoints sit at (+/-2250, +/-2250) um. NOTE this is the full offset,
% not half of it -- a Cs written with ell/2 needs ell = 4500 um to agree.
if isempty(ip.Results.ell)
    pr = plant.idx.sensePair(1,:);
    ell = abs((plant.term(pr(1)).Yc + plant.term(pr(2)).Yc)/2);
else
    ell = ip.Results.ell;
end

%% ---------------- continuous plant ----------------
A = zeros(6);
A(1,2)=1;  A(2,2)=-cz/m;
A(3,4)=1;  A(4,4)=-qx/Ix;
A(5,6)=1;  A(6,6)=-qy/Iy;

% INPUT IS THE WRENCH, not a force per pair.
%   u = [dFz; dtaux; dtauy]
% The five drive pairs are only three independent directions -- rank(B) is
% 3 at the centred pose -- so commanding five pair-forces asked place() to
% resolve a 2-parameter redundancy arbitrarily, and it did so badly
% (cond(ctrb) ~ 6e8). More importantly, pinning each pair's own z-force
% leaves the TILT stiffness only partly cancelled, because a pair's two
% members sit at different coordinates. Commanding the wrench and letting
% trapWrenchToVoltage allocate it across all five pairs delivers all three
% components exactly, so the inversion cancels the stiffness completely
% (beta = 1) and the redundancy is resolved explicitly by a minimum-norm
% allocation instead of implicitly by the pole placement.
B = zeros(6,3);
B(2,1) = 1/m;
B(4,2) = 1/Ix;
B(6,3) = 1/Iy;

Cs = [1 0  ell 0 -ell 0;
      1 0  ell 0  ell 0;
      1 0 -ell 0  ell 0;
      1 0 -ell 0 -ell 0];

Cc = [1 0 0 0 0 0;        % controlled outputs: rz, thx, thy
      0 0 1 0 0 0;
      0 0 0 0 1 0];

%% ---------------- structural checks ----------------
rC = rank(ctrb(A,B));
rO = rank(obsv(A,Cs));
K.ctrb.rank = rC;  K.ctrb.full = (rC == 6);
K.obsv.rank = rO;  K.obsv.full = (rO == 6);
if rC < 6
    error('trapControlDesign:uncontrollable', ...
        'ctrb rank %d of 6; the force set cannot reach every state.', rC);
end
if rO < 6
    error('trapControlDesign:unobservable', ...
        'obsv rank %d of 6; the four sensors cannot resolve every state.', rO);
end

%% ---------------- discretise ----------------
% The regulator is NOT designed on the zero-order hold of a double
% integrator. That model is wrong, and measurably so.
%
% Within one ISR the voltage is frozen, and with b ~ 1/g^2 the plate obeys
%       qddot = p^2*(q - q_k)*(1-beta) + p^2*q*0 ... (see below) + Q/m
% where p is the local unstable rate and beta is the fraction of the
% stiffness the V^2 inversion actually cancels. Solving that over one tick
% gives the EXACT sampled map
%
%   A_d = [ beta + (1-beta)*cosh(pT)   sinh(pT)/p ]
%         [ (1-beta)*p*sinh(pT)        cosh(pT)   ]
%   B_d = [ (cosh(pT)-1)/p^2 ; sinh(pT)/p ] * (row of B)
%
% beta = 1 gives [1 sinh/p; 0 cosh] -- the position term cancels but the
% VELOCITY is still amplified by cosh(pT) every tick, so even "fully
% cancelled" heave is not a double integrator. beta = 0 gives the
% uncancelled [cosh sinh/p; p*sinh cosh]. Both limits are wrong for this
% device: measurement gives beta = 1.0000 (rz), 0.8000 (thx), 0.7237 (thy).
% Heave cancels exactly because both members of a pair share a gap; tilt
% does not, because a pair's two members sit at different coordinates and
% a tilt redistributes force between them.
%
% As T -> 0 every entry collapses to the double integrator, which is why
% the old model looked right and was not: at T = 500 us, cosh(pT) = 1.168,
% a 17% per-tick velocity gain the design model simply did not contain.
% Designing on this map instead lowers the usable regulator bandwidth from
% about 100 Hz to about 20 Hz.
Kst  = trapOpenLoopStiffness(plant, 'Vtia', Vtia);
inert = [m; Ix; Iy];
Ktot  = [Kst(3,3); Kst(4,4); Kst(5,5)];
pRate = sqrt(max(-Ktot,0)./inert);        % local unstable rate per DOF
beta  = ones(3,1);                        % wrench allocation cancels all 3

Ad = zeros(6); Bd = zeros(6,3);
for i = 1:3
    r = 2*i-1;  pp = pRate(i);
    if pp*T < 1e-9 || strcmpi(ip.Results.designModel,'double')
        Ad(r:r+1,r:r+1) = [1 T; 0 1];
        Bd(r:r+1,:)     = [T^2/2; T] * B(r+1,:);
    else
        ch = cosh(pp*T);  sh = sinh(pp*T);
        Ad(r:r+1,r:r+1) = [1, sh/pp; 0, ch];
        Bd(r:r+1,:)     = [(ch-1)/pp^2; sh/pp] * B(r+1,:);
    end
end
Csd = Cs;                                  % zero-order hold leaves C alone

% augment with one integrator per controlled output
nx = 6; ni = 3;
Aa = [Ad,            zeros(nx,ni);
      -T*Cc,         eye(ni)];
Ba = [Bd; zeros(ni,3)];
Br = [zeros(nx,ni); T*eye(ni)];

rCa = rank(ctrb(Aa,Ba));
K.ctrb.rankAug = rCa;  K.ctrb.fullAug = (rCa == nx+ni);
if rCa < nx+ni
    error('trapControlDesign:augUncontrollable', ...
        ['Augmented ctrb rank %d of %d. The integrators are not reachable, ', ...
         'which happens when the plant has a transmission zero at z = 1 ', ...
         'for the controlled outputs.'], rCa, nx+ni);
end

%% ---------------- pole placement ----------------
sR = uniformDecayPoles(nx+ni, lambdar);
sE = uniformDecayPoles(nx,    lambdae);
zR = exp(sR*T);
zE = exp(sE*T);

% ---- SCALE BEFORE PLACING ----
% The raw state mixes metres (~1e-5) with radians (~1e-4) and positions
% with velocities, so ctrb(Aa,Ba) has condition number ~1e12 and place()
% is operating far outside its reliable range: perturbing Aa at the 1e-12
% level moved Kx by 16%, which made the design non-reproducible between
% runs and the closed-loop outcome a coin flip at short Tisr.
%
% Fix: non-dimensionalise first. Each state is divided by a characteristic
% magnitude -- the gap for heave, the tilt that closes the gap at the plate
% corner for the angles, position/T for the rates, position*T for the
% integrators -- and each input by a characteristic wrench. Gains are
% placed on the scaled system and mapped back. The pole locations are
% unchanged by a similarity transform, so this is exact, not approximate.
sz = plant.gap.z0;                       % [m]   heave scale
st = plant.gap.z0/(plant.mass.L/2);      % [rad] tilt that closes the corner
Wt = plant.mass.m*plant.phys.g;          % [N]   weight
sX = diag([sz, sz/T, st, st/T, st, st/T]);
sI = diag([sz*T, st*T, st*T]);
sU = diag([Wt, Wt*a, Wt*a]);
sXI = blkdiag(sX, sI);

Aa_s = sXI \ Aa * sXI;
Ba_s = sXI \ Ba * sU;
K.condCtrbRaw    = cond(ctrb(Aa, Ba));
K.condCtrbScaled = cond(ctrb(Aa_s, Ba_s));

% ---- SYNTHESIS: LQR, not pole placement ----
% uniformDecayPoles asks for 9 poles sharing one real part, i.e. a nearly
% DEFECTIVE spectrum. With 3 inputs and 9 states the gains achieving those
% poles form an 18-parameter family, and place() picks one arbitrarily:
% cond(eigenvectors of the closed loop) came out ~1e8, and two different
% but equally valid picks gave opposite large-signal outcomes from the same
% start (survives vs contact at 2 ms). The poles were right in both cases;
% what differed was robustness, which pole placement does not control.
%
% LQR returns the unique minimiser, so the design is reproducible and the
% gain is the robustness-optimal member of that family. Bandwidth is still
% specified as lambdar: the control weight is bisected until the slowest
% closed-loop mode sits at exp(-lambdar*T), so the bandwidth spec is met
% while the remaining freedom is resolved optimally rather than arbitrarily.
% On the scaled system every state is O(1), so Q = I is already a Bryson
% weighting and needs no further tuning.
rhoTarget = exp(-lambdar*T);
Qw = eye(nx+ni);
lo = 1e-10; hi = 1e10;
for it = 1:200
    rr = sqrt(lo*hi);
    Kt = dlqr(Aa_s, Ba_s, Qw, rr*eye(3));
    if max(abs(eig(Aa_s - Ba_s*Kt))) < rhoTarget
        lo = rr;            % too fast -> penalise control more
    else
        hi = rr;
    end
    if hi/lo < 1+1e-12, break; end
end
Ka_s = dlqr(Aa_s, Ba_s, Qw, sqrt(lo*hi)*eye(3));
Ka   = sU * Ka_s / sXI;
K.Kx = Ka(:,1:nx);
K.Ki = Ka(:,nx+1:end);
K.lqrR = sqrt(lo*hi);

% ---- estimator: designed on the STIFFNESS-RETAINING plant ----
% The regulator above is correctly designed on the stiffness-free A: once
% the estimate has converged, the voltage inversion really does deliver the
% commanded force, so the loop sees a pure double integrator.
%
% The ESTIMATOR does not get that luxury. The inversion cancels the
% electrostatic stiffness at the ESTIMATED pose, so with b ~ 1/g^2,
%       F_delivered = Ftot*b(x)/b(xhat) ~ Ftot*(1 - 2*e_rz/g)
%       F_delivered - W ~ -(2W/g)*e_rz = -37.76*e_rz
% i.e. the uncancelled stiffness acts on the ESTIMATION ERROR. Writing
% e = xhat - x and carrying this through,
%       e[k+1] = (Ad - L*Cs)*e - G*e = (Ads - L*Cs)*e
% where Ads is A with the fixed-voltage stiffness restored. So the error
% plant is Ads, NOT Ad, and L has to be placed on Ads. Designing it on Ad
% leaves the observer modelling a plant it is not watching: at
% Tisr = 720 us that gives max|eig(Ads - L*Cs)| = 1.479, an outright
% divergent estimator, which is what made the closed loop blow up while
% every stage tested exact in isolation.
%
% Note xhat is still PROPAGATED with Ad/Bd below: in x coordinates the
% plant is stiffness-free and the stiffness enters only through e. Only the
% gain placement changes.
As = A;
As(2,1) = -Kst(3,3)/m;
As(4,3) = -Kst(4,4)/Ix;
As(6,5) = -Kst(5,5)/Iy;
sysds = c2d(ss(As,B,Cs,0), T, 'zoh');
Ads = sysds.A;  Bds = sysds.B;

% Same scaling argument applies to the observer: the measurements are
% heights [m] while the states run from metres to rad/s.
sY     = sz;
Ads_s  = sX \ Ads * sX;
Csd_s  = (Csd * sX) / sY;
L_s    = place(Ads_s', Csd_s', zE)';
K.L    = sX * L_s / sY;
K.condObsvRaw    = cond(obsv(Ads, Csd));
K.condObsvScaled = cond(obsv(Ads_s, Csd_s));

%% ---------------- pack ----------------
K.A=A; K.B=B; K.Cs=Cs; K.Cc=Cc;
K.Ad=Ad; K.Bd=Bd; K.Csd=Csd;
K.Aa=Aa; K.Ba=Ba; K.Br=Br;
K.polesR=zR; K.polesE=zE;
K.sR=sR; K.sE=sE;
K.Tisr=T; K.lambdar=lambdar; K.lambdae=lambdae;
K.ell=ell; K.a=a;
K.As=As; K.Ads=Ads; K.Bds=Bds; K.Kst=Kst; K.Vtia=Vtia;
K.beta=beta; K.pRate=pRate;
K.clR = eig(Aa - Ba*Ka);
K.clE = eig(Ads - K.L*Csd);        % the error plant is Ads, not Ad
K.clEfree = eig(Ad - K.L*Csd);     % what the old design reported
end

% =====================================================================
function s = uniformDecayPoles(n, w)
%UNIFORMDECAYPOLES n poles all with real part -w, in conjugate pairs.
%
% Every mode then decays as exp(-w*t), which is what "bandwidth w,
% critically damped" means. A Butterworth pattern would spread the real
% parts from -w down to -0.17*w, leaving the slowest mode six times slower
% than the stated bandwidth -- not what the spec asks for.
%
% A small distinct imaginary part is given to each pair because place()
% cannot resolve a repeated pole of multiplicity above rank(B), and nearly
% repeated poles make the gain matrix badly conditioned. The spread is
% kept at a few percent of w so the damping stays essentially critical.
nPair = floor(n/2);
s = [];
for k = 1:nPair
    im = w * 0.08 * k;                      % few-percent spread, distinct
    s = [s; -w + 1i*im; -w - 1i*im];        %#ok<AGROW>
end
if mod(n,2) == 1
    s = [s; -w];
end
s = s(:);
end
