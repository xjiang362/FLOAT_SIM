function M = trapFullStateModel(plant, varargin)
%TRAPFULLSTATEMODEL Closed-form DT model with per-channel fractional delays,
%                   plus LQI regulator and LQE (dual-LQI) observer gains.
%
%   M = trapFullStateModel(plant)
%   M = trapFullStateModel(plant, 'rhoU',1e-2, 'rhoI',1e2, 'qProc',1, 'rMeas',1)
%
% CONTINUOUS PLANT. Three decoupled second-order modes, input is a
% generalised ACCELERATION u = [a_z; alpha_x; alpha_y]:
%
%   A_c = blkdiag( [0 1; kappa_z -2 sigma_z], [0 1; kappa_x -2 sigma_x],
%                  [0 1; kappa_y -2 sigma_y] )
%   kappa_j = k_j / m_j ,  sigma_j = c_j / (2 m_j)
%
% k_j is the FIXED-VOLTAGE electrostatic stiffness with the sign flipped so
% that k_j > 0 is anti-restoring (A_c(2,1) = +kappa). Taken from
% trapOpenLoopStiffness, so k_z = 2W/g0 = 37.76 N/m exactly.
%
% CLOSED FORM. With omega_j = sqrt(sigma_j^2 + kappa_j),
%
%   E_j(d) = e^{sigma d} [ cosh(wd) - (sigma/w) sinh(wd) ,  -sinh(wd)/w
%                          -(kappa/w) sinh(wd)           ,   cosh(wd) + (sigma/w) sinh(wd) ]
%          = expm(-A_c d)                                   (verified to 2e-16)
%
%   phi_j(d) = E_j(d)(1,1),   psi_j(d) = E_j(d)(1,2),
%   g_j(d)   = [ (1-phi_j(d))/kappa_j ; -psi_j(d) ] = int_0^d expm(-A_c s) B ds
%
%   A_d = blkdiag(E_j(-T)),   B_d = -blkdiag(g_j(-T))
%
% All of these reproduce expm/c2d to machine precision.
%
% ------------------------------------------------------------------
% CORRECTION TO THE SOURCE DERIVATION. The note forms the delayed-channel
% rows as
%       z[k+1] = Phi x[k] + Gamma u[k],   Phi_i = C_i expm(-A_c delta_i)
% These two are inconsistent. expm(-A_c delta) propagates BACKWARD, so
% Phi_i x[k] is the channel sampled at t_k - delta_i, which is before frame
% k begins and therefore cannot carry u[k]. Checked against a brute-force
% integration the note's row is wrong by 13-38% at the delays used here.
%
% A channel delivered at t_{k+1} with age delta_i was sampled at
% t_k + (T - delta_i), inside frame k, so the consistent rows are
%
%       Phi_i   = C_i expm(A_c (T - delta_i)) = Phi_i^{note} A_d
%       Gamma_i = C_i int_0^{T-delta_i} expm(A_c s) B ds
%               = Phi_i^{note} B_d + Gamma_i^{note}
%
% i.e. the note's matrices are right up to one extra step of propagation.
% Both corrected forms reproduce the brute-force sample to 2e-16. The
% alternative repair -- keeping Phi^{note} and pairing it with u[k-1] --
% is also self-consistent but describes a different, longer schedule.
% ------------------------------------------------------------------
%
% DELAYS. Channel k of 4 reaches the controller stale by
%       age(k) = tI2C + dtChannel*(n - k + 0.5)
% Channels 1..3 have age > T, so they land one frame late and are carried
% in the state z with residual delay delta_i = age_i - T. Channel 4 has
% age < T and is read directly in the current frame with delta_4 = age_4.
%
% AUGMENTED SYSTEM. chi = [x(6); z(3); uPrev(3); xI(3)], 15 states:
%
%   chi[k+1] = Achi chi[k] + Bchi u[k] + Brchi r[k]
%   u[k]     = -[Kx Kz Ku KI] chi[k]
%
% with xI[k+1] = xI[k] + T(r[k] - H y[k]) and H the 3x4 left inverse of the
% sensing map, H = [1 1 1 1; 1 1 -1 -1; -1 1 1 -1]/4 scaled by 1/Lc on the
% tilt rows.
%
% GAINS. Both by discrete Riccati, on a NON-DIMENSIONALISED system: the raw
% state mixes metres with radians over six orders of magnitude and the
% Riccati solve is ill-conditioned without scaling (see trapControlDesign
% for the same issue).
%   regulator  dlqr(Achi, Bchi, Q, R),   Q = blkdiag(Qx, 0, 0, QI)
%   observer   dlqr(Abar', Cbar', Qo, Ro)'  -- the LQR dual, i.e. LQE
%
% NAME/VALUE
%   'rhoU'   control weight scale (default 1e-2)
%   'rhoI'   integrator weight scale (default 1e2)
%   'qProc'  process-noise scale for the observer (default 1)
%   'rMeas'  measurement-noise scale for the observer (default 1e-2)
%   'Tisr','dtChannel','tI2C'  override the readout timing
%
% OUTPUT M: .Ac .Bc .Cc .Ad .Bd .Phi .Gam .Phi4 .Gam4 .Cbar .Dbar .H
%           .Achi .Bchi .Brchi .Kx .Kz .Ku .KI .L .Abar .Bbar
%           .delta .age .sigma .omega .kappa .Lc .T .clR .clE
%
% See also TRAPOPENLOOPSTIFFNESS, TRAPCONTROLDESIGN, SIMULATEFULLSTATELQI.

ip = inputParser;
% Defaults chosen from a sweep at the 60 um gap with the full 1130 us
% acquisition delay: the loop holds for rhoU >= 1e4 with rhoI <= 1e2, and
% 1e5/1e1 sits centrally in that region. Faster weights (rhoU <= 1e2) do
% not survive the delay.
ip.addParameter('rhoU',  1e5);
ip.addParameter('rhoI',  1e1);
ip.addParameter('qProc', 1);
ip.addParameter('rMeas', 1e-2);
ip.addParameter('Tisr',      plant.readout.Tisr);
ip.addParameter('dtChannel', plant.readout.dtChannel);
ip.addParameter('tI2C',      plant.readout.tI2C);
ip.parse(varargin{:});
T  = ip.Results.Tisr;
dc = ip.Results.dtChannel;
ti = ip.Results.tI2C;

m  = plant.mass.m;  Ix = plant.mass.I(1,1);  Iy = plant.mass.I(2,2);
Kst = trapOpenLoopStiffness(plant, 'Vtia', 0);
kk  = [-Kst(3,3); -Kst(4,4); -Kst(5,5)];       % >0 means anti-restoring
mm  = [m; Ix; Iy];
cc  = [plant.damp.c; plant.damp.q; plant.damp.q];

sigma = cc ./ (2*mm);
kappa = kk ./ mm;
omega = sqrt(sigma.^2 + kappa);

pr = plant.idx.sensePair(1,:);
Lc = abs((plant.term(pr(1)).Yc + plant.term(pr(2)).Yc)/2);

%% ---------------- continuous plant ----------------
Ac = zeros(6); Bc = zeros(6,3);
for j = 1:3
    r = 2*j-1;
    Ac(r:r+1, r:r+1) = [0 1; kappa(j), -2*sigma(j)];
    Bc(r+1, j) = 1;
end
Cc = [1 0  Lc 0 -Lc 0;
      1 0  Lc 0  Lc 0;
      1 0 -Lc 0  Lc 0;
      1 0 -Lc 0 -Lc 0];

%% ---------------- closed-form blocks ----------------
Ej = @(j,d) exp(sigma(j)*d) * ...
     [cosh(omega(j)*d)-sigma(j)/omega(j)*sinh(omega(j)*d), -sinh(omega(j)*d)/omega(j);
      -kappa(j)/omega(j)*sinh(omega(j)*d),                  cosh(omega(j)*d)+sigma(j)/omega(j)*sinh(omega(j)*d)];
phij = @(j,d) exp(sigma(j)*d)*(cosh(omega(j)*d)-sigma(j)/omega(j)*sinh(omega(j)*d));
psij = @(j,d) -exp(sigma(j)*d)*sinh(omega(j)*d)/omega(j);
gj   = @(j,d) [(1-phij(j,d))/kappa(j); -psij(j,d)];

Ad = blkdiag(Ej(1,-T), Ej(2,-T), Ej(3,-T));
Bd = zeros(6,3);
for j = 1:3, Bd(2*j-1:2*j, j) = -gj(j,-T); end

%% ---------------- per-channel delays ----------------
n   = 4;
age = ti + dc*(n - (1:n)' + 0.5);
delta      = zeros(n,1);
delta(1:3) = age(1:3) - T;      % land one frame late, residual delay
delta(4)   = age(4);            % read in the current frame
if any(delta < 0) || any(delta(1:3) >= T) || delta(4) >= T
    warning('trapFullStateModel:schedule', ...
        ['Channel delays fall outside the assumed schedule ' ...
         '(0 <= delta < T for all four). delta = %s us.'], ...
        mat2str(delta'*1e6, 5));
end

% Backward-propagated rows, exactly as the source derivation writes them:
%   z[k+1] = C_i x(t_{k+1} - age_i) = C_i x(t_k - delta_i)
% since age_i = T + delta_i for channels 1..3. The interval
% [t_k - delta_i, t_k] lies inside frame k-1, so these rows carry u[k-1],
% NOT u[k] -- that is the one index the derivation has wrong, and it is the
% only thing wrong with Phi/Gamma.
bwd = @(d) blkdiag(Ej(1,d), Ej(2,d), Ej(3,d));
Phi = zeros(3,6); Gam = zeros(3,3);
for i = 1:3
    Gblk = zeros(6,3);
    for j = 1:3, Gblk(2*j-1:2*j, j) = gj(j, delta(i)); end
    Phi(i,:) = Cc(i,:) * bwd(delta(i));
    Gam(i,:) = -Cc(i,:) * Gblk;
end
% Channel 4 has age < T, so its interval also lies in frame k-1 and its row
% enters y[k] through Dbar*u[k-1] with delta_4 = age_4.
G4 = zeros(6,3);
for j = 1:3, G4(2*j-1:2*j, j) = gj(j, delta(4)); end
Phi4 = Cc(4,:) * bwd(delta(4));
Gam4 = -Cc(4,:) * G4;

%% ---------------- measurement and integral maps ----------------
Cbar = [zeros(3,6), eye(3);
        Phi4,       zeros(1,3)];
Dbar = [zeros(3,3); Gam4];
H    = 0.25*[1 1 1 1; 1/Lc 1/Lc -1/Lc -1/Lc; -1/Lc 1/Lc 1/Lc -1/Lc];

% Observer plant. z[k+1] takes u[k-1], x[k+1] takes u[k], so the two inputs
% enter through separate matrices:
%   xbar[k+1] = Abar xbar[k] + Bbar u[k] + BbarPrev u[k-1]
Abar     = [Ad, zeros(6,3); Phi, zeros(3,3)];
Bbar     = [Bd; zeros(3,3)];
BbarPrev = [zeros(6,3); Gam];

%% ---------------- augmented system ----------------
Cbx = [zeros(3,6); Phi4];
Cbz = [eye(3); zeros(1,3)];
nx = 6; nz = 3; nu = 3; ni = 3;
% z[k+1] = Phi x[k] + Gamma u[k-1]: the sample interval lies in frame k-1,
% so Gamma multiplies the PREVIOUS input and belongs in the uPrev column of
% Achi, which is where the source derivation puts it.
Achi = [ Ad,            zeros(nx,nz), zeros(nx,nu), zeros(nx,ni);
         Phi,           zeros(nz,nz), Gam,          zeros(nz,ni);
         zeros(nu,nx),  zeros(nu,nz), zeros(nu,nu), zeros(nu,ni);
        -T*H*Cbx,      -T*H*Cbz,     -T*H*Dbar,     eye(ni)];
Bchi  = [Bd; zeros(nz,nu); eye(nu); zeros(ni,nu)];
Brchi = [zeros(nx,ni); zeros(nz,ni); zeros(nu,ni); T*eye(ni)];

%% ---------------- scaling ----------------
sz = plant.gap.z0;                      % [m]   heave scale
st = plant.gap.z0/(plant.mass.L/2);     % [rad] tilt that closes the corner
Sx = diag([sz, sz/T, st, st/T, st, st/T]);
Sz = sz*eye(nz);
Su = diag([sz/T^2, st/T^2, st/T^2]);    % accelerations
SI = diag([sz*T, st*T, st*T]);
Schi = blkdiag(Sx, Sz, Su, SI);

Achi_s = Schi \ Achi * Schi;
Bchi_s = Schi \ Bchi * Su;

%% ---------------- LQI regulator ----------------
rhoU = ip.Results.rhoU;  rhoI = ip.Results.rhoI;
Q = blkdiag(eye(nx), zeros(nz), zeros(nu), rhoI*eye(ni));
R = rhoU*eye(nu);
Ka_s = dlqr(Achi_s, Bchi_s, Q, R);
Ka   = Su * Ka_s / Schi;
M.Kx = Ka(:,1:nx);
M.Kz = Ka(:,nx+1:nx+nz);
M.Ku = Ka(:,nx+nz+1:nx+nz+nu);
M.KI = Ka(:,nx+nz+nu+1:end);

%% ---------------- LQE observer (the LQR dual) ----------------
Sxz = blkdiag(Sx, Sz);
Sy  = sz*eye(4);
Abar_s = Sxz \ Abar * Sxz;
Cbar_s = Sy \ (Cbar * Sxz);     % yhat = Sy^-1 y, xhat = Sxz^-1 x
Qo = ip.Results.qProc * eye(nx+nz);
Ro = ip.Results.rMeas * eye(4);
L_s  = dlqr(Abar_s', Cbar_s', Qo, Ro)';
M.L  = Sxz * L_s / Sy;

%% ---------------- pack ----------------
M.Ac=Ac; M.Bc=Bc; M.Cc=Cc; M.Ad=Ad; M.Bd=Bd;
M.Phi=Phi; M.Gam=Gam; M.Phi4=Phi4; M.Gam4=Gam4;
M.Cbar=Cbar; M.Dbar=Dbar; M.Cbx=Cbx; M.Cbz=Cbz; M.H=H;
M.Abar=Abar; M.Bbar=Bbar; M.BbarPrev=BbarPrev;
M.Achi=Achi; M.Bchi=Bchi; M.Brchi=Brchi;
M.delta=delta; M.age=age; M.sigma=sigma; M.omega=omega; M.kappa=kappa;
M.k=kk; M.Lc=Lc; M.T=T; M.dtChannel=dc; M.tI2C=ti;
M.rhoU=rhoU; M.rhoI=rhoI; M.qProc=ip.Results.qProc; M.rMeas=ip.Results.rMeas;
M.clR = eig(Achi - Bchi*Ka);
M.clE = eig(Abar - M.L*Cbar);
M.condRegRaw    = cond(ctrb(Achi, Bchi));
M.condRegScaled = cond(ctrb(Achi_s, Bchi_s));
M.condObsRaw    = cond(obsv(Abar, Cbar));
M.condObsScaled = cond(obsv(Abar_s, Cbar_s));
end
