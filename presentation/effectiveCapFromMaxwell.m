function [Ceff, C_M, C_1, C_2] = effectiveCapFromMaxwell(C11, C22, CRR, C12, C1R, C2R)
%EFFECTIVECAPFROMMAXWELL Effective capacitance between terminals 1 and 2,
%eliminating a single floating terminal R (the resonator), given the 6
%relevant MAXWELL capacitance matrix entries directly (no mutual
%conversion needed -- these are assumed to already be in Maxwell
%convention, i.e. off-diagonal entries negative).
%
%   [Ceff, C_M, C_1, C_2] = effectiveCapFromMaxwell(C11, C22, CRR, C12, C1R, C2R)
%
%   Inputs (all Maxwell-convention, [F]; vectorized -- pass column
%   vectors, one row per pose/sweep point, for a whole sweep at once):
%     C11, C22, CRR   diagonal entries (self-terms) for terminals 1, 2, R
%     C12             off-diagonal, terminals 1-2 (negative)
%     C1R, C2R        off-diagonal, terminal 1-R and 2-R (negative)
%
%   Outputs (same size as inputs):
%     Ceff   effective capacitance [F] looking into terminals 1,2
%            (e.g. f = 1/(2*pi*sqrt(L*Ceff)) for an LC resonance)
%     C_M    direct mutual branch between 1 and 2, post-elimination
%     C_1, C_2  each terminal's effective to-ground capacitance, post-elimination

% Schur complement, eliminating R
C_red11 = C11 - C1R.^2 ./ CRR;
C_red22 = C22 - C2R.^2 ./ CRR;
C_red12 = C12 - (C1R.*C2R) ./ CRR;

% Convert reduced 2x2 Maxwell submatrix back to mutual form
C_M = -C_red12;
C_1 = C_red11 + C_red12;
C_2 = C_red22 + C_red12;

Ceff = C_M + (C_1.*C_2) ./ (C_1+C_2);

end
