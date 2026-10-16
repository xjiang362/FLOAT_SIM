function [Ceff, C_M, C_i, C_j] = effectiveMutualCapacitance(Cm, term_i, term_j, floatingTerminals)
%EFFECTIVEMUTUALCAPACITANCE Effective 2-terminal capacitance between two
%terminals in a multi-terminal capacitive network, accounting for
%floating terminals via Schur-complement (Kron) reduction.
%
%   [Ceff, C_M, C_i, C_j] = effectiveMutualCapacitance(Cm, term_i, term_j, floatingTerminals)
%
%   Inputs (1-INDEXED terminal numbers, matching COMSOL/simulation
%   convention -- terminal k corresponds directly to row/col k of Cm):
%     Cm                 NxN mutual (lumped/partial) capacitance matrix [F]
%     term_i, term_j     the two terminals to find effective capacitance between
%     floatingTerminals  vector of OTHER floating terminal numbers to
%                        eliminate (voltage not fixed by anything --
%                        must NOT include term_i or term_j)
%
%   Any terminal NOT in {term_i, term_j, floatingTerminals} is assumed
%   to be held at a fixed voltage (driven or grounded) and is simply
%   dropped from the reduction -- its coupling doesn't affect the AC/
%   resonance behavior between term_i and term_j (a fixed-voltage
%   terminal's own voltage never changes, so it contributes no AC
%   current, only a constant DC charge offset).
%
%   Outputs:
%     Ceff   effective capacitance [F] as seen looking into term_i/term_j
%            (e.g. for an LC resonance: f = 1/(2*pi*sqrt(L*Ceff)))
%     C_M    direct mutual branch between term_i and term_j (post-reduction)
%     C_i, C_j  each terminal's effective "to ground" capacitance (post-reduction)
%
%   Method: Cm is converted to the Maxwell (admittance-equivalent)
%   convention, since Schur-complement/Kron-reduction network
%   elimination is defined on that convention (Y=j*omega*C_Maxwell).
%   Floating terminals are eliminated via the Schur complement of their
%   self-block; the remaining 2x2 Maxwell submatrix (now an isolated
%   2-terminal system) is converted back to mutual form to get C_M,
%   C_i, C_j, then combined via the direct-path-parallel-with-
%   through-ground-series-path formula.

N = size(Cm,1);

% --- Step 1: Cm (mutual) -> C (Maxwell) ---
% C_ij = -Cm_ij (i~=j);  C_ii = sum_k Cm_ik (row sum, INCLUDING diagonal)
Cm_offdiag = Cm - diag(diag(Cm));
C = diag(Cm * ones(N,1)) - Cm_offdiag;

% --- Step 2: index bookkeeping (1-indexed terminals map directly to rows/cols) ---
idx_i = term_i;
idx_j = term_j;
idx_float = floatingTerminals(:)';

keep = [idx_i, idx_j];
elim = idx_float;
% Everything else (not keep, not elim) is fixed-voltage/grounded --
% simply dropped by never including it in the restricted submatrix below.

sub = [keep, elim];
Csub = C(sub, sub);
nk = numel(keep);

C_kk = Csub(1:nk, 1:nk);
C_ke = Csub(1:nk, nk+1:end);
C_ek = Csub(nk+1:end, 1:nk);
C_ee = Csub(nk+1:end, nk+1:end);

% --- Step 3: Schur complement (Kron reduction), eliminating floating terminals ---
if isempty(elim)
    C_reduced = C_kk;
else
    C_reduced = C_kk - C_ke * (C_ee \ C_ek);
end

% --- Step 4: reduced 2x2 Maxwell submatrix -> mutual form (isolated 2-terminal system) ---
C_M = -C_reduced(1,2);
C_i = C_reduced(1,1) + C_reduced(1,2);
C_j = C_reduced(2,2) + C_reduced(1,2);

% --- Step 5: effective capacitance (direct path parallel with through-ground series path) ---
Ceff = C_M + (C_i*C_j)/(C_i+C_j);

end
