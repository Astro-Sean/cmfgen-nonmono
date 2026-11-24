# Coupled wavelength-depth block solve: scoping investigation

Opened 24-Nov-2026, user-directed. Question: can a Baron &
Hauschildt (2004, A&A 427, 987) style coupled wavelength-depth
direct solve replace the lagged frequency iteration in
MOM_J_CMF_NONMONO_V1, and at what cost. This document is the
scope - each phase is a bounded investigation with its own exit
criterion. No implementation before phase D reports.

## Why this route (evidence recap)

Three arms + audit: lagged splitting cannot contract in
compressing advection-dominated regions (no positive physical
source; H-store feedback carries inflated commits into the next
RHS; matrix side proven innocent - a < 1 strictly, NMARG=0).
The frequency-lagged iteration is the defect; the block solve
removes it rather than mitigating it. Proven in principle:
PHOENIX supernova atmospheres with non-monotonic velocity.

## Phase A: literature reread (1-2 days agent time)

Extract from Baron & Hauschildt 2004 and follow-ups (incl.
Hauschildt & Baron PHOENIX CMF papers, Knop et al. 2009 for our
own discretization lineage):
- exact unknown vector layout (per-depth frequency blocks or
  per-frequency depth blocks), bandwidth of the assembled operator
- how they treat the closure (H reconstruction) inside the solve:
  linear system in J only, or J+H coupled
- boundary-condition embedding at both depth and frequency edges
- what stays outer-iterated (opacities, EPS, source functions)
- their reported convergence behaviour in compressing regions
Exit: a written discretization mapping their equations onto our
row assembly (TA/TB/TC/XM, PSI coupling, VB/VC H coefficients),
with equation-number citations verified by rt_literature.

## Phase B: real-grid census (needs arm 13654121 data + static counts)

- NCF x ND per direction set; sign-flip count along frequency per
  depth (the ring density -> true bandwidth of the coupled system)
- memory/flop estimates: banded dense LU vs sparse direct vs
  block-Jacobi-preconditioned GMRES at production grid sizes
- how much of the operator is frequency-local (preconditioner
  candidate: exact depth-tridiagonal inverse per frequency)
Exit: cost table at m15 production grid; go/no-go on direct vs
iterative method choice.

## Phase C: toy 2D solver prototype

Standalone driver (Fortran, unit_tests/f90) on a manufactured
non-monotonic velocity field with analytic opacities:
- assemble the full wavelength-depth operator in sparse form
- verify: (1) reduces to the per-frequency tridiagonal solves in
  the purely-expanding limit; (2) contracts in a compressing
  manufactured case where the emulated lagged iteration diverges
  (reuse the iterated-probe infrastructure of drv_mom_nu iter);
  (3) GMRES iteration counts with depth-block preconditioner
Exit: measured contraction in the regime where the current scheme
provably fails; iteration-count scaling table.

## Phase D: integration design (decision gate, user sign-off)

Written design covering:
- what enters the block operator implicitly (advection/aberration
  frequency coupling - the defect) vs what stays outer-iterated
  (opacities, EPS scatter, H closure, source functions)
- how the refresh/pass/sweep machinery reshapes: the block solve
  replaces the per-frequency DGTSV chain inside the sweep; the
  outer population iteration is untouched
- VADAT surface (one new flag, default F, per AGENTS 9.1 proposal)
- memory plan: assembly once per refresh vs per sweep
Exit: user-reviewed design doc; implementation go/no-go.

## Phase E: implementation (only after D sign-off)

Staged: sparse assembly module -> toy-verified solver import ->
integration behind flag -> unit batteries -> arm on eps1em3 seed.
Each stage a bounded investigation; full testing ladder per
AGENTS.md before any model submission.

## Cost honesty

Realistic estimate: weeks-scale. The frequency-loop rearchitecture
touches sweep, refresh, boundary closures, and store ownership.
The toy prototype (phase C) is the cheap decisive evidence - if
the block solve does not contract the manufactured compressing
case, the route dies there at low cost. Phase C is therefore the
first implementation target after A/B report.


## Implementation status (24-Nov-2026)

Phases A-D are complete and the phase-D design was implemented as
`MOM_BLKGMRES`, default F. E1-E3 are green: the matrix-free production
solver matches dense references to machine precision on manufactured
chains, rings, marginal rings, mixed-sign structures and inactive
columns; the integrated moment-solver and VADAT batteries are green.
E4 is in progress on the eps1em3 seed through CMFGEN job 13655952.
The production verdict remains pending, so this document is closed as
a scoping plan but not as a solver validation.
