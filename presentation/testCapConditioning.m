%% testCapConditioning.m
% The closed form is a difference of logarithms divided by a small
% number, TWICE over:
%
%   C = (eps/thy) * [ B_R^-1 h(A_R+B_R y) - B_L^-1 h(A_L+B_L y) ]
%
%   inner cancellation:  1/B  with B = q*thy - thx   (B -> 0)
%   outer cancellation:  1/thy                        (thy -> 0)
%
% Both are 0/0. capTrapz switches to series expansions below fixed
% relative thresholds, but a threshold is only correct if it sits where
% the SERIES error and the CANCELLATION error cross. Below the switch the
% error is the series truncation (grows with the parameter); above it the
% error is lost significant digits (shrinks with the parameter). So the
% total error must PEAK at the threshold. That is a falsifiable
% prediction about the implementation -- test it by sweeping the
% offending parameter across many decades.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

g0 = 15e-6; um = 1e-6;
ya = 25*um; yb = 1475*um; pL = 1525*um; pR = 2975*um;   % one actuation square
xmax = pR;
thx = 0.1289*pi/180;

ref = @(G,tx,ty) integral2(@(yy,xx) 1./(G - yy*tx + xx*ty), ya, yb, ...
        @(yy) pL+0*yy, @(yy) pR+0*yy, 'AbsTol',1e-18,'RelTol',1e-13);

%% ---- Outer cancellation: sweep thy across 14 decades ----
fprintf('============ OUTER CANCELLATION: 1/thy ============\n');
fprintf('capTrapz switches branch when |thy|*xmax <= 1e-9*G,\n');
fprintf('i.e. at thy = %.3e rad. Error should peak near there.\n\n', 1e-9*g0/xmax);
fprintf('%-12s %-14s %14s %14s %11s\n','thy [rad]','thy*xmax/G','closed form','quadrature','rel err');

tys = 10.^(-16:1:-2);
errs = zeros(size(tys));
for i = 1:numel(tys)
    ty = tys(i);
    Cc = capTrapz(ya,yb,pL,0,pR,0, g0, thx, ty);
    Cn = ref(g0, thx, ty);
    errs(i) = abs(Cc-Cn)/abs(Cn);
    fprintf('%-12.0e %-14.2e %14.8e %14.8e %11.2e\n', ty, ty*xmax/g0, Cc, Cn, errs(i));
end
[~,ip] = max(errs);
fprintf('\nerror peaks at thy = %.0e  (predicted switch near %.0e)\n', tys(ip), 1e-9*g0/xmax);

%% ---- Inner cancellation: drive B -> 0 ----
% B = q*thy - thx. For a rectangle q=0, so B = -thx: drive thx -> 0.
fprintf('\n============ INNER CANCELLATION: 1/B ============\n');
fprintf('capTrapz uses the series when |B|*ymax < 1e-8*|A|,\n');
fprintf('i.e. at thx = %.3e rad.\n\n', 1e-8*g0/yb);
fprintf('%-12s %-14s %14s %14s %11s\n','thx [rad]','|B|*ymax/A','closed form','quadrature','rel err');
ty_fixed = 1e-3;
txs = 10.^(-16:1:-3);
for i = 1:numel(txs)
    tx = txs(i);
    Cc = capTrapz(ya,yb,pL,0,pR,0, g0, tx, ty_fixed);
    Cn = ref(g0, tx, ty_fixed);
    fprintf('%-12.0e %-14.2e %14.8e %14.8e %11.2e\n', tx, tx*yb/g0, Cc, Cn, abs(Cc-Cn)/abs(Cn));
end

%% ---- Is the threshold optimally placed? ----
% Series error ~ d^3 (three terms kept); cancellation error ~ 1e-16/d.
% They cross where d^4 = 1e-16, i.e. d = 1e-4 -- so the 1e-8 threshold on
% B is CONSERVATIVE by four decades. Harmless (it only means the series
% is used where the exact form would also have been fine), but worth
% recording so nobody later "fixes" it in the wrong direction.
d = 10.^(-10:0.5:-2);
fprintf('\n============ THRESHOLD PLACEMENT ============\n');
fprintf('%-10s %12s %12s\n','d=|B|y/A','series err','cancel err');
for dd = [1e-9 1e-8 1e-6 1e-4 1e-3 1e-2]
    fprintf('%-10.0e %12.2e %12.2e\n', dd, dd^3, 1e-16/dd);
end
fprintf('crossover at d ~ 1e-4; current B threshold is 1e-8 (conservative, safe)\n');

%% ---- Near-contact behaviour: does h(u) survive u -> 0? ----
fprintf('\n============ NEAR CONTACT (u -> 0) ============\n');
fprintf('h(u) = u(ln u - 1) -> 0 as u -> 0, so the primitive itself is\n');
fprintf('well behaved. But in floating point 0*(-Inf) = NaN, so an exact\n');
fprintf('zero gap would poison the result. Check the approach:\n\n');
fprintf('%-14s %14s %14s\n','min gap [m]','C','finite?');
for gmin = [1e-6 1e-8 1e-10 1e-12 0]
    % tilt chosen so the gap reaches gmin exactly at y = yb
    tx = (g0 - gmin)/yb;
    Cc = capTrapz(ya,yb,pL,0,pR,0, g0, tx, 0);
    fprintf('%-14.0e %14.6e %14s\n', gmin, Cc, string(isfinite(Cc)));
end
fprintf('\n(capacitance genuinely diverges as a contact line forms --\n');
fprintf(' a growing value here is physics, not a numerical failure.)\n');
