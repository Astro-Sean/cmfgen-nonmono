# Phase C prototype verdict: block solve contracts where lagged diverges

24-Nov-2026. Prototype: unit_tests/python/proto_block_gmres.py -
dense toy system, 24 frequencies x 16 depths, per-frequency
diffusion-like tridiagonal blocks, sign-keyed upwind frequency
coupling with a compression band, diagonal upwind absorption per
the Baron-Hauschildt chi-hat form. Lagged iteration modelled as
block-Jacobi in frequency (exactly what the production sweep does:
each frequency's depth system solved exactly, neighbours read
lagged). GMRES: restartless, left-preconditioned with the
frequency block-diagonal inverse - the matrix-free analogue of
the existing per-frequency DGTSV solves.

## Result 1: full absorption (xi=1, the B&H setup)

coupling  rho_lagged  lagged iters to 1e-10  GMRES+blockJacobi
   4.0     0.688            61                   24
   8.0     0.815           112                   26
  16.0     0.898           214                   29
  32.0     0.946           417                   36

Even in the benign regime the lagged iteration stalls as coupling
grows (rho -> 1), while preconditioned GMRES stays nearly flat.
This is the lamesF stall phenotype reproduced in miniature.

## Result 2: partial absorption (production's descaled rows)

Production rows do not always carry the full upwind weight on the
diagonal: LAG_BND descales lagged content and compressing rows sit
at excess ~ 0. Modelling this as absorption factor f at coupling 8:

f     rho_lagged   GMRES+blockJacobi
1.0    0.815           26
0.8    0.974           48
0.6    1.210           48    <- lagged DIVERGES, GMRES contracts
0.4    1.596           46    <- lagged DIVERGES, GMRES contracts

## Verdict

The coupled wavelength-depth system solved by preconditioned GMRES
contracts in exactly the regime where the lagged iteration provably
diverges, at grid sizes where the iteration count is insensitive to
the divergence strength. Phase C exit criterion met: the route
lives. The preconditioner is the per-frequency tridiagonal inverse -
already in production as the DGTSV chain - so the integration cost
is the operator action and the GMRES shell, not new linear algebra.

Caveats for phase D: dense toy, no rings resolved per depth, no
nonlinear source feedback (the production H-store loop), no
opacity/closure outer iteration. The toy proves the linear algebra
route; it does not prove the outer iteration. That is what the
phase-D design must scope honestly.

Next: phase B grid census (matrix-free memory/iteration estimates
at production NCF x ND) then the phase-D design doc for user
sign-off before any integration.
