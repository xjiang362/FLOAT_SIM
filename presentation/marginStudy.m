%% marginStudy.m
% Gain and phase margins of the levitation controller.
%
% WHY THIS IS NOT A FORMALITY HERE. State-feedback LQR has guaranteed
% margins at the plant input -- at least 6 dB gain reduction, infinite gain
% increase, 60 deg of phase. Those guarantees apply to the FULL-STATE loop
% only. Doyle's 1978 result is that LQG, i.e. the same LQR gain closed
% through an observer, has NO guaranteed margins at all: the observer can
% erode them arbitrarily. This device is open-loop unstable and the loop
% carries a large acquisition delay, so the erosion is worth measuring
% rather than assuming.
%
% WHAT IS COMPUTED
%   1. Loop broken at the PLANT INPUT, one channel at a time with the other
%      loops closed (classical loop-at-a-time margins via allmargin).
%   2. The same for the state-feedback loop with the observer removed, to
%      isolate what the observer costs.
%   3. Peak sensitivity Ms = ||S||_inf and complementary Mt = ||T||_inf,
%      which bound the margins simultaneously in all loops:
%           GM >= Ms/(Ms-1),    PM >= 2 asin(1/(2 Ms))
%      These are MIMO-valid, unlike the loop-at-a-time numbers.
%
% The loop transfer is verified before use: the poles of feedback(L,I) must
% reproduce the known closed-loop eigenvalues.
%
% See also TRAPCONTROLDESIGN, TRAPFULLSTATEMODEL.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

p = trapDeviceGeometry();
K = trapControlDesign(p);
T = K.Tisr;
nx = 6; ni = 3; nu = 3;

fprintf('Design point: Tisr = %.0f us, lambda_r = %.0f Hz, lambda_e = %.0f Hz\n', ...
    T*1e6, K.lambdar/2/pi, K.lambdae/2/pi);
fprintf('Open-loop plant is UNSTABLE: max|eig(Ad)| = %.4f\n\n', max(abs(eig(K.Ad))));

%% ---------------- plant and controller as DT systems ----------------
G = ss(K.Ad, K.Bd, K.Csd, zeros(4,nu), T);          % u -> y

% observer-based controller, y -> u, with the sign convention u = -Kc(z) y
Ak = [K.Ad - K.Bd*K.Kx - K.L*K.Csd,  -K.Bd*K.Ki;
      -T*K.Cc,                        eye(ni)];
Bk = [K.L; zeros(ni,4)];
Ck = [K.Kx, K.Ki];
Kc = ss(Ak, Bk, Ck, zeros(nu,4), T);

Lin = Kc*G;                                          % loop at the plant input

%% ---------------- verify the loop transfer ----------------
% The loop above is built on G = (Ad, Bd, Csd), the stiffness-FREE design
% plant, so its closed-loop poles are eig(Aa - Ba Ka) together with
% eig(Ad - L Csd). K.clE is the error plant WITH the stiffness retained and
% is the wrong reference here.
clRef = sort(abs([eig(K.Aa - K.Ba*[K.Kx K.Ki]); K.clEfree]),'descend');
clLoop = sort(abs(pole(feedback(Lin, eye(nu)))),'descend');
nCmp = min(numel(clRef), numel(clLoop));
err = max(abs(clRef(1:nCmp) - clLoop(1:nCmp)));
fprintf('loop-transfer check: max|pole mismatch| vs the known closed loop = %.3e\n', err);
if err > 1e-6
    warning('marginStudy:loop','Loop transfer does not reproduce the closed loop; margins below are not trustworthy.');
end
fprintf('closed-loop max|z| = %.6f\n\n', max(clLoop));

%% ---------------- 1. loop-at-a-time margins, observer in the loop ----
fprintf('=== margins at the PLANT INPUT, observer-based (LQG) ===\n');
% For an UNSTABLE plant the gain margin is an interval: too little gain is
% as destabilising as too much, so both edges matter. Report the range.
fprintf('%10s %22s %12s %14s %16s\n','channel','gain margin [dB]','PM [deg]','xover [Hz]','delay margin[us]');
nmU = {'a_z','alpha_x','alpha_y'};
for i = 1:nu
    Li = looptAtATime(Lin, i, nu);
    am = allmargin(Li);
    [glo, ghi] = gmRange(am.GainMargin);
    pm = pickFinite(am.PhaseMargin);
    wc = pickFinite(am.PMFrequency);
    dm = pickFinite(am.DelayMargin);
    fprintf('%10s %22s %12s %14s %16s\n', nmU{i}, ...
        sprintf('[%s, %s]', fmt(20*log10(glo),'%+.2f'), fmt(20*log10(ghi),'%+.2f')), ...
        fmt(pm,'%.1f'), fmt(wc/2/pi,'%.1f'), fmt(dm*T*1e6,'%.0f'));
end

%% ---------------- 2. same loop, observer removed --------------------
fprintf('\n=== margins at the PLANT INPUT, full-state feedback (LQR) ===\n');
Ga = ss(K.Aa, K.Ba, eye(nx+ni), 0, T);               % augmented, full state out

Lsf = [K.Kx K.Ki] * Ga;
fprintf('%10s %14s %14s %16s\n','channel','GM [dB]','PM [deg]','gain xover[Hz]');
for i = 1:nu
    Li = looptAtATime(Lsf, i, nu);
    am = allmargin(Li);
    fprintf('%10s %14s %14s %16s\n', nmU{i}, ...
        fmt(20*log10(pickFinite(am.GainMargin)),'%.2f'), ...
        fmt(pickFinite(am.PhaseMargin),'%.1f'), ...
        fmt(pickFinite(am.PMFrequency)/2/pi,'%.1f'));
end

%% ---------------- 3. sensitivity peaks (MIMO-valid bounds) ----------
fprintf('\n=== sensitivity peaks and the margins they guarantee ===\n');
S = feedback(eye(nu), Lin);          % S = (I+L)^-1
Tc = feedback(Lin, eye(nu));         % T = L(I+L)^-1
Ms = getPeakGain(S, 1e-6);
Mt = getPeakGain(Tc, 1e-6);
fprintf('  ||S||_inf = %.3f (%.2f dB)   ||T||_inf = %.3f (%.2f dB)\n', ...
    Ms, 20*log10(Ms), Mt, 20*log10(Mt));
fprintf('  guaranteed simultaneously in ALL loops:\n');
fprintf('    gain margin  >= %.2f  (%.2f dB)\n', Ms/(Ms-1), 20*log10(Ms/(Ms-1)));
fprintf('    phase margin >= %.1f deg\n', 2*asin(1/(2*Ms))*180/pi);
fprintf('  (a well-conditioned design keeps Ms below about 2)\n');

%% ---------------- 4. delay already in the loop ----------------------
fprintf('\n=== how much delay is already present ===\n');
tOld = p.readout.tI2C + p.readout.dtChannel*3.5;
fprintf('  acquisition delay (oldest channel) %.0f us = %.2f ISR periods\n', tOld*1e6, tOld/T);
fprintf('  one ISR of zero-order hold           %.0f us\n', T*1e6);
fprintf('  NOTE the margins above are for the loop WITHOUT the acquisition\n');
fprintf('  delay; compare the delay margin in section 1 against %.0f us.\n', tOld*1e6);

% =====================================================================
function Li = looptAtATime(L, i, n)
% break channel i, leave the other loops closed
others = setdiff(1:n, i);
if isempty(others)
    Li = L; return
end
Li = feedback(L, eye(n-1), others, others, +1);  % close the others (neg fb)
Li = Li(i,i);
end

function v = pickFinite(x)
if isempty(x), v = Inf; return; end
x = x(isfinite(x) & x > 0);
if isempty(x), v = Inf; else, v = min(x); end
end

function s = fmt(v, f)
if ~isfinite(v), s = 'Inf'; else, s = sprintf(f, v); end
end

function [glo, ghi] = gmRange(gm)
% gain-margin interval around 1 (unity). Unstable plants have a lower edge.
gm = gm(isfinite(gm) & gm > 0);
if isempty(gm), glo = 0; ghi = Inf; return; end
lo = gm(gm < 1); hi = gm(gm > 1);
if isempty(lo), glo = 0; else, glo = max(lo); end
if isempty(hi), ghi = Inf; else, ghi = min(hi); end
end
