# KNOP second-pass audit - open items

Follow-ups from `KNOP_SECOND_PASS_AUDIT.md` (11-Nov-2026).  Ordered by
decision value.  None of these relax convergence criteria or raise
watchdog budgets; they isolate the failure mechanism before any solver
change is considered.  Do not port SOB population-solver knobs
(`BA_DEM_PFLOOR`, `BA_EX_RELFRAC`) until P5-T1 settles the conditioning
question.

## P1 - field non-contraction (highest value)

- **P1-T1** Finite-difference `I_NEXT` sensitivity — **DONE 11-Nov-2026**
  on a manufactured `|Q|->1` compressing model (`p1t1` driver case).
  Confirmed by **two independent methods**: (a) `p1t1_analyze.py`
  re-march validated to ~1e-8 against production `I_M`/`I_P`, and (b) a
  direct production-solver store perturbation (`FREEZE_I_NEXT=T` +
  `NEW_FREQ=F` re-solve, two amplitudes).  Result: outward
  formal-solve gain `g=dI/dI_NEXT -> ~0.97` at `|Q|=0.998` (vs <=7e-2
  all other channels, inward `~1.3e-3`); effective store gain
  `(1-relax)+relax*g ~ 0.985` at `INEXT_RELAX=0.5` — essentially
  non-contractive.  Strong-coupling hypothesis **supported** for the
  outward `I_P` channel.  L=1 (bluest freq, `I_*_PREV=0` boundary) shows
  a non-local outward response (peak depth I=3) and an
  amplitude-independent inward step — flagged as a boundary artefact,
  not the smooth mode.  See `KNOP_P1_T1_FIELD_SENSITIVITY.md`.  Remaining
  optional refinement: repeat on production `FRESH_knop` argmax cells
  (`d=3-5, f~1249`) with a warm store to confirm `g` there.
- **P1-T1b** Global lagged-map eigenvalue.  The P1-T1 gain `~0.985` is a
  local single-cell estimate.  Build the dominant-eigenvalue estimate
  over the full `(compressing cell × frequency × ray)` store — either a
  matrix-free power iteration through the production sweep or a
  coherent multi-cell perturbation — to test whether the *network*
  gain (not just the local cell) approaches or exceeds 1 once
  population feedback and cross-frequency store sharing are included.
  **DONE 11-Nov-2026** (`KNOP_P1_T1B_LAGGED_MAP.md`): dense `30x30`
  central-FD Jacobian of the full lagged-store map (`knray p1t1b`,
  `DRV_P1T1B`).  The transport Jacobian is **nilpotent** — every
  compressing read targets `write_slot + NM_RAY` (`NM_RAY >= 1`), an
  acyclic frequency-ladder shift, so all eigenvalues are `0` and the
  relaxed map has spectral radius `1 - INEXT_RELAX = 0.5` (measured
  contraction `0.500015` in direct iteration).  The full production map
  (sweep + `FG_KNOP_REFRESH_V1`) is exactly seed-independent (`J = 0`);
  the store is slaved to the fixed opacities.  **The network gain does
  not approach 1; the radiation `I_NEXT` map cannot sustain drift.**
  Caveat: populations held fixed, so the outer population-opacity map
  (P5) remains the drift suspect.
- **P1-T2** Add a per-cell coefficient dump (argmax cell's
  `V/R, SIGMA, MU, GAM, CHI, CHI_HAT, DLN, Q, I_NEXT, S`) to the
  `KNOP FIELD` block or a side file, gated by a debug flag, so the
  worst-cell physics is inspectable without a custom build.
- **P1-T3** Confirm the depth-index convention (`NEW_SN_R_GRID` row vs
  computational depth) so argmax depths map unambiguously to `SIGMA`.

## P2 - `KN_NO_GAIN` physical cost

- **P2-T1** Quantify what the floor removes.  A `KN_NO_GAIN` variant
  that caps gain (albedo <= 1) rather than zeroing net absorption, or a
  run with smaller `INEXT_RELAX`, would let the physical cost be
  measured against the `8.6x/iter -> L~1E29` divergence it prevents.
  Record emergent luminosity, continuum flux, energy balance, population
  residuals and max `chi_scat/chi`.
- **P2-T2** Verify the "iterate artifact" rationale for the suppressed
  inversions (high-n Rydberg / intra-n far-IR) against a converged or
  near-converged population state - currently a code-comment assertion.

## P3 - ordered march

- **P3-T1** Add an `ORD_SWEEP`/`HAVE_PACK`/`NI` census line to the
  `KNOP FIELD` block and a matching ITER_LEDGER branch column so the
  executed branch is auditable post hoc.
- **P3-T2** Controlled `KN_RAY_ORDERED=T` vs `F` A/B on an identical
  warm-started checkpoint (same model, atomic data, `I_NEXT` seed,
  tolerances, pass limits, OMP).  Record post-refresh residual,
  per-ray convergence rate, pass count, failed rays, store drift,
  wall-clock.
- **P3-T3** Preserve `OUTGEN` (or tee the `KNOP FIELD` block to a small
  dedicated file) in run dirs - restart prep currently consumes it,
  destroying the only branch-level evidence.

## P4 - transfer equation

- **P4-T1** Quantify the `+3 a_nu I` omission in the compressing regime.
  Compare the implemented equation against the full nu-form on a
  representative compression profile; bound the error it introduces.
- **P4-T2** Fix the comment contradiction: header (~89-101) says the
  term "cancels identically" while the detailed block (~3501-3537) says
  it is "deliberately not added".  Keep the detailed (correct)
  justification; correct the header.  Also reconcile the `+/- GAM`
  `dI/d(ln nu)` sign between ~3474 and ~4023.  Comment-only change.

## P5 - population conditioning

- **P5-T1** Warm-restart `VAL_knop_formal` (or a shortened reproduction
  reaching `maxraw=1E7`) with `POP_AUDIT=T` plus the already-armed
  `POP_RUNAWAY_DEM=1E3`, `TOPN=10`, and **preserve OUTGEN**.  At the
  `maxraw=1E7` pin, correlate `RUNAWAY it=` channel lines (`fl=1` floor
  flag, `dem%`) with per-depth `RUNAWAY d= nflr= colmin_u= rcond=`.
  Confirms or refutes H1 on the production model.
- **P5-T2** Only after P5-T1: evaluate whether the floor-column null
  space needs a conditioning fix (column scaling / truncation / pivot
  fallback in `CMF_BLKBAND`), distinct from the demand-metric
  exemptions.

## P6 - observer

- **P6-T1** Controlled `NUM_ES = 2, 3, 5` comparison on identical
  populations once a converged `I_NEXT` exists; record post-refresh
  residual and integrated luminosity.  Final observer render must use a
  converged `I_NEXT`.
- **P6-T2** Apply OBS fixes C1-C5 from `CMF_FLUX_KNOP_CONVERGENCE_AUDIT`
  before any physical SOB-vs-KNOP spectral comparison; the `~9.5x`
  suppression is a cold-store artefact until then.

## P7 - documentation

- **P7-T1** Amend `KNOP_FORMAL_SOLVER_AUDIT.md` items 19-20 to state
  explicitly that the "end-to-end integration test" is execution-only
  (output produced), not convergence - matching F-23 in the production
  audit.
- **P7-T2** Done in this audit: future-dated `20-Nov-2026` comments at
  `fg_j_cmf_nonmono_v1.f` ~3747, ~5483, ~9814 corrected to `10-Nov-2026`.
