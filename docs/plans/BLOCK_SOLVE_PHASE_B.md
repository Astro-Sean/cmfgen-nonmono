# Phase B: production-grid census and method sizing

24-Nov-2026. Grid facts from the m15 lamesF arms (OUTGEN census):
ND_NM = 202 depth points on the nonmono grid; total moment cells
~ 970K per pass census -> NCF ~ 4800 frequencies. Unknowns per
direction set n = ND_NM x NCF ~ 9.7e5.

## Method cost table

Dense assembled direct solve: n^2 storage ~ 7.5e12 doubles -
absurd, rejected without further analysis. Paper II's ceiling
(64 layers x 1000 wavelengths, memory-bound) confirms: assembly
scales badly and production grids are 5x beyond it in frequency.

Banded dense LU: bandwidth ~ 2 x ND_NM ~ 404; banded storage
n x 2bw x 8B ~ 6 GB, factorisation fill and flop cost scale
n x bw^2 ~ 1.6e14 flops - hours per solve, per refresh, per
iteration. Rejected.

Matrix-free GMRES + block-Jacobi-in-frequency preconditioner:
- operator action storage: TA/TB/TC + coupling coefficient + edge
  index per row ~ 40 MB total
- Krylov basis at restart m=40: 40 x n x 8B ~ 310 MB - fits a node
- preconditioner: per-frequency tridiagonal LU, same coefficients
  as the action; application is exactly the existing DGTSV chain,
  parallel across frequencies
- cost per GMRES iteration ~ one sweep-pass worth of tridiagonal
  work; current refreshes already spend 20-40 passes stall-bound,
  so the solve cost is the SAME ORDER as present spend - the
  difference is convergence instead of stall

## Ring/bandwidth structure

The coupling graph is functional (one edge per row): expanding
rows read k-1, compressing rows read k+MOM_NM_J. Sign flips along
frequency create rings of length ~ 2 x MOM_NM_J; production census
counted 2-6K rings per pass of a few members each. In matrix terms
the frequency-direction structure is block-banded with bandwidth
set by MOM_NM_J x ND_NM, not the full n - reinforcing that sparse
iterative, not dense direct, is the only sane method.

## Verdict

Method choice confirmed: matrix-free restartless/restarted GMRES,
left-preconditioned by the per-frequency tridiagonal inverse,
one direction set at a time. Memory ~ 350 MB peak, cost per solve
comparable to one existing refresh. No assembly ever formed.
Phase B exit criterion met - proceed to phase-D design doc.
