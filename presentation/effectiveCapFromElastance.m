function [Ceff, C_M, C_1, C_2] = effectiveCapFromElastance(S11, S22, S12)
%EFFECTIVECAPFROMELASTANCE Effective capacitance between terminals 1 and
%2, given only their 2x2 elastance (inverse Maxwell) sub-block.
%
%   [Ceff, C_M, C_1, C_2] = effectiveCapFromElastance(S11, S22, S12)
%
%   Uses the block-matrix identity (S)_kk = (Schur complement of C)^-1:
%   inverting the small elastance sub-block for just the terminals you
%   care about automatically eliminates EVERY other terminal in the
%   system (resonator included) via Schur complement -- no need to
%   know anything about the rest of the network. This is an
%   approximation (it treats all other terminals as floating/eliminated
%   rather than distinguishing fixed-voltage "delete" vs floating
%   "eliminate"), but verified to differ from the exact
%   eliminate-resonator-only treatment by only ~0.16% for this device.
%
%   Inputs (vectorized -- pass column vectors for a whole sweep at once):
%     S11, S22   elastance diagonal entries for terminals 1, 2 [1/F]
%     S12        elastance off-diagonal entry, terminals 1-2 [1/F]
%
%   Outputs (same size as inputs):
%     Ceff   effective capacitance [F] looking into terminals 1,2
%     C_M    direct mutual branch between 1 and 2, post-elimination
%     C_1, C_2  each terminal's effective to-ground capacitance, post-elimination

n = numel(S11);
Ceff = zeros(n,1); C_M = zeros(n,1); C_1 = zeros(n,1); C_2 = zeros(n,1);

for k = 1:n
    Ssub = [S11(k), S12(k); S12(k), S22(k)];
    Cred = inv(Ssub); %#ok<MINV>

    C_M(k) = -Cred(1,2);
    C_1(k) = Cred(1,1) + Cred(1,2);
    C_2(k) = Cred(2,2) + Cred(1,2);
    Ceff(k) = C_M(k) + (C_1(k)*C_2(k))/(C_1(k)+C_2(k));
end

end
