function K = ackermannSISO(A, B, poles)
%ACKERMANNSISO Ackermann's formula: state feedback gain K such that
%eig(A - B*K) equals the given desired poles, for a single-input pair
%(A,B). Also used (via the A',C' duality) to design an observer gain L
%such that eig(A - L*C) hits desired poles: call with (A', C') and
%transpose the result.
%
%   Implemented directly (no Control System Toolbox dependency, matching
%   the rest of this repo) since MATLAB's own place()/acker() refuse
%   exactly-repeated poles for a single-input system, which is exactly
%   the bandwidth-driven repeated-real-pole design used here.

n = size(A, 1);

p = poly(poles);   % monic char. poly coefficients, descending powers, length n+1

phiA = zeros(n);
for i = 0:n
    phiA = phiA + p(i + 1) * A^(n - i);
end

Cm = zeros(n);
Apow = eye(n);
for i = 1:n
    Cm(:, i) = Apow * B;
    Apow = Apow * A;
end

lastRow = zeros(1, n);
lastRow(end) = 1;

K = lastRow * (Cm \ phiA);

end
