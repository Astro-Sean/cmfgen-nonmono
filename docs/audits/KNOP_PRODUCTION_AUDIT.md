# KNOP solver production-readiness audit - 11-Nov-2026

Scope: `DO_NONMONO_KNOP` end to end - the `FG_J_CMF_NONMONO_V1`
formal solver (`plane/`), the lagged `I_NEXT` store and
`FG_KNOP_REFRESH_V1` machinery, the GS freeze/maturity/watchdog
stack in `SOLVE_FOR_POPS`, the population/BA/resacc coupling, the
accelerator keys (`DO_ES_ALI`, `DO_ES_ALI_BA`, `DO_BA_RES_ACC`,
JFNK), the config/gating surface, and the `CMF_FLUX` observer
path.  This document synthesises the existing audit chain rather
than re-deriving it; primary sources are listed per section.

## Verdict

**NOT PRODUCTION READY.**  The discretisation is audit-verified
against Knop Papers I/II and the stock path is untouched, but no
KNOP run has ever converged a model end to end, the observer
render is a cold-start diagnostic with confirmed defects, and the
newest hardening carries "no model run yet" tags throughout.
Per user directive 11-Nov-2026, `DO_NONMONO_SOB` +
`USE_LAM_ES=T` is the production route; KNOP is retained for
reference and diagnostic comparisons only.

## 1. Convergence evidence - the decisive gap

The standing acceptance test is "a converged full non-monotonic
model run end to end through `run_model_grid.py`."  No run has
met it:

| run | recipe | outcome |
|---|---|---|
| `FRESH_knop` (job 13005611) | cold, `KN_NO_GAIN=T`, `DO_ES_ALI=T`, `DO_ES_ALI_BA=T`, `DO_BA_RES_ACC=T` | **STOP non-convergent it124** - `GS_REL_NOPROG_LIM=10`: "the release/refreeze cycle never matured the lagged field" |
| `FRESH_knop_d1a` (job 13005609) | same, D1a-instrumented build | running, it70 / 50% at audit time - same freeze/release signature |
| `VAL_knop_formal` | warm from SOB `tauto` it1445 | **STOP it1489** - "luminosity diverging, L_out=1.28E+15"; `maxraw` pinned at the 1E7 cap for the last iterations; `KN_NO_GAIN` floored 8447 cells (max chi_scat/chi=5.1) without saving it |
| `VAL_knop_mom` | `USE_LAM_ES=F` moment route | repeated "Possible error converging f" in the continuum formal solve; route closed anyway (`USE_LAM_ES` forced T under KNOP since 08-Nov) |
| `KNOP_opt`, `VAL_knop_jfnk` | optimum recipe | died at setup - `GENOSC_V9` "incorrect number energy levels" (data error, never exercised the solver) |
| `KNOP_fresh` (02-Nov) | `KN_NO_GAIN=F` | diverged it3->it14, `L_out`->1E29 |
| `KNOP_l14_freet` (29-Oct) | warm restart | diverged it1162 - cold `I_NEXT` rebuild + `N_LAM_RWARM` released BA on an immature field |

Two distinct kill mechanisms appear in the live evidence:

- **Field-maturity stall** (`FRESH_knop`): the GS freeze/release
  cycle ran 27 population releases in 124 iterations; demanded
  correction fell 1.6E12 -> ~38 and then oscillated 38-52 without
  progress for ~60 iterations.  `GS_HOLD_POPS` held populations on
  every LAMBDA_GS pass; only the periodic LAMBDA_TFIX releases
  moved them.  The watchdog stopped the run rather than let it
  grind - the guard worked as designed, but the underlying map
  did not mature.
- **Population runaway / luminosity divergence**
  (`VAL_knop_formal`): `maxraw` at the 1E7 trust cap, consistent
  with the floor-pinned trace-population null-space mode
  instrumented in `KNOP_ITEM1_RUNAWAY_2026-11-08.md` (hypothesis
  H1: per-depth SE block loses rank at the population floor, so
  `SOLVEBA_V13` demands ~1E7 corrections forever).  H1 remains
  theorised; the `POP_RUNAWAY_*` census exists but no confirming
  run data is recorded.

Cold starts diverge; every serious attempt has needed a warm
start from a converged SOB checkpoint, and even that has not
produced convergence.

## 2. Formal solver correctness - verified, with documented deviations

`KNOP_AUDIT_2026-10-30.md`, `KNOP_FORMAL_SOLVER_AUDIT.md`,
`KNOP_MATH_AUDIT.md` (vault), `NONMONO_COMPREHENSIVE_AUDIT_2026-10-30.md`:

- Term-by-term verification vs Knop Paper I (nu-form
  `chi_hat = chi + |GAM|/dln nu`, signed-Q upwind source) and
  Paper II (ordered march) - PASS WITH LIMITATIONS.
- Sweep/refresh operator identity verified to 2.6E-8 on a
  manufactured ray; `FGINEXT_GS` checkpoint round-trip bit-exact;
  1/2/4 OMP threads bit-identical after the `FG_JUMP_MAX` race fix.
- Documented deviations from the papers: diagonal ALO only (the
  paper's Lambda* is tridiagonal in wavelength), no `3*a_nu*I`
  aberration term (inherited from stock V13), red-edge absorbing
  BC, nu-form frequency-spacing approximation, `NONMONO_INEXT_CAP`
  overflow guard, `KN_CHI_EXT/KN_ETA_EXT` single-precision capture
  (~1E-7 rel), zero-GAM `MOM_GAM_FAC_MIN` floor bias on the
  (closed) moment route.
- `KN_RAY_ORDERED=T` ordered march - the main paper-conform
  convergence fix - is **theory-derived and runtime-unverified**;
  the A/B vs `KN_RAY_ORDERED=F` on a compressing model is still
  pending (vault Open issues).

## 3. Lagged-field machinery - hardened at unit level, unproven at model level

`FG_FIELD_DRIFT_AUDIT.md` (authoritative; supersedes
`FIELD_DRIFT_AUDIT.md` for status) plus the 10-Nov commit cycle
(`b49b32d3..280c4a44`):

- Fixed and unit-pinned: per-ray contraction estimator (A1),
  divergence-check ordering (A3), shared refresh/sweep DTAU
  operator (B1), commit-census metrics (C1), `HEAL_STORE` census
  (C2), jump-rejected false maturity via `FG_RF_PEND` (C4), JFNK
  eval telemetry bracketing (D1), OMP race, unsolved-ray census.
- The drift plateau itself is **genuine non-contraction of the
  discretised map**, not a measurement artefact.  CORRECTED
  11-Nov-2026 by `KNOP_SECOND_PASS_AUDIT.md` Pillar 1: the marginal
  (`|Q| -> 1`) cells are the **strong-coupling compressing** cells
  (`|GAM|/(DLN*CHI) -> infinity`), not the `|1+sigma*K| -> 0`
  near-zero-`GAM` cells - those have `|Q| -> 0` and contract at the
  relaxation rate.  The earlier phrasing here ("marginally
  contractive ... `|1+sigma*K| -> 0` resonance cells") named the wrong
  limit and is retracted.  Whether the *outer* population-coupled map
  is also expansive (A2) is open; the audit's own words: "if a model
  still drifts at measured rho << 1, the outer map itself is
  expansive and the remedy is structural."
- `INEXT_RELAX` adaptive omega now pools per-ray P90
  contraction - `#status/partial`, "no model run yet" at commit
  time; FRESH_knop is the first real exercise and it stalled.
- `FG_RF_DCHG` reset defect (FIELD_DRIFT D1) fixed; metric
  relabelling (D2) cosmetic-open.

## 4. Population-solver coupling - the runaway mode is unresolved

- `KNOP_STALL_DIAGNOSIS_2026-10-31.md`: the lagged field showed a
  bounded oscillation (`fg d ~ 6e-2..1.5e-1`, hundreds of K bad
  cells) while the population channel pinned at the 1E7 trust cap
  - "a stable period-2/limit-cycle of the coupled map, not
  approaching a fixed point."
- `KNOP_ITEM1_RUNAWAY_2026-11-08.md` + `ITEM2_POPULATION_FLOOR_AUDIT.md`
  + `POPULATION_SOLVER_FIXES_2026-11-09.md`: floor-pinned trace
  populations may give the per-depth SE block an approximate null
  space (H1, unconfirmed).  `STEQ_RAW` raw-demand metrics landed;
  the `BA_QUAR` quarantine and floor-pin fixes are open/partial.
- `FULL_ITERATION_CYCLE_AUDIT.md`: F-01/F-02/F-16 resolved
  04-Nov (convergence exit now requires field maturity;
  JFNK gated to `USE_LAM_ES`; `HAS_LAGGED_FIELD` route flag).
  Still open: F-03 (frozen JFNK tier inert under SOB), F-05
  (resacc-off commit on immature field), F-08/F-10 (freeze budgets
  tuned per-model), F-12 (exit is bookkeeping, physical
  self-consistency inferred not proven), F-23 (**no
  full-iteration integration tests exist**; T-08..T-11 not
  started).
- Accelerators on KNOP (`DO_ES_ALI`, `DO_ES_ALI_BA`,
  `DO_BA_RES_ACC`, JFNK): `DO_ES_ALI` verified working on the
  live run (J-loop `unconverged=0`, `dJ/J=1e-4`, zero ALO
  fallbacks after the 02-Nov accessor fix); `KN_NO_GAIN=T` forced
  keeps the DEN-floor albedo assumption valid.  But the dedicated
  `KNOP_esA/esB` A/B never ran, `DO_ES_ALI_BA` showed no benefit
  on SOB and is unmeasured on KNOP, and none of them address the
  two observed kill mechanisms (field maturity, population
  runaway).  They are safe to leave armed; they are not a
  convergence fix.

## 5. Observer path - confirmed defects, fixes unapplied

`CMF_FLUX_KNOP_CONVERGENCE_AUDIT.md` (11-Nov, newest and
authoritative for OBS):

- **D1** `NUM_ES=1` + KNOP accepted; the refresh gate never fires
  -> cold `I_NEXT`.  Only the Python pipeline guards it.
- **D2** SOB populations never write `FGINEXT_GS` -> any KNOP OBS
  render cold-starts the lagged field.
- **D3** `USE_LAM_ES=F` populations have no matched renderer.
- **D4** runtime: KNOP OBSFLUX median 4.63E4 vs stock 4.41E5
  (~9.5x continuum suppression), consistent with the cold-store
  mechanism (chi_hat inflation + missing `Q*I_NEXT` source =
  artificial sink in compressing cells).
- **S6 BLOCKED**: no converged KNOP population exists, so the
  intended production OBS path has never been exercised.
- Proposed fixes C1-C5 (NUM_ES guard, maturity export, adaptive
  refresh, docs, edge-slot census) - **none applied**.

## 6. Readiness matrix

| subsystem | correctness | convergence | runtime evidence | verdict |
|---|---|---|---|---|
| `FG_J_CMF_NONMONO_V1` discretisation | verified vs papers | n/a (per-solve) | manufactured + live sweeps | PASS WITH LIMITATIONS |
| Lagged `I_NEXT` store + refresh | operator identity 2.6E-8, checkpoint bit-exact | **stall observed live; suspected marginal contraction on `|Q|->1` strong-coupling cells (mechanism corrected by second-pass audit, intrinsic non-contraction unproven - P1-T1)** | FRESH_knop it124 STOP | FAIL (convergence) |
| GS freeze/maturity/watchdog | guards fire correctly | catches stalls but cannot fix them | FRESH_knop, VAL_knop_formal | PASS (as guard) / FAIL (as cure) |
| Population SE/BA coupling | ordering audited | runaway/null-space mode open (H1) | VAL_knop_formal L-divergence | FAIL |
| Accelerators (ES_ALI/BA/resacc/JFNK) | math audited | unverified benefit on KNOP | J-loop clean on FRESH_knop | PARTIAL - safe, not a fix |
| OBS `CMF_FLUX` KNOP render | defects D1-D4 verified | cold-start only | 9.5x suppression measured | FAIL |
| Config/gating surface | auto-sets + hard stops in place | per-model tuned | - | PARTIAL |
| Integration tests | unit pins green (41/41 knop_ray etc.) | **none exist for full iteration** | - | FAIL |

## 7. Blockers to production, ranked

1. **No converged model.**  Two kill mechanisms observed live:
   field-map non-maturity (GS_REL_NOPROG_LIM stop) and
   population-side runaway/luminosity divergence.  Until one
   m15-class model converges end to end, nothing downstream is
   exercisable.
2. **Population runaway mode unconfirmed/unfixed.**  H1 census
   instrumentation exists; the `BA_QUAR`/floor-pin remedies are
   open.  This is the same defect family SOB's
   `BA_DEM_PFLOOR`/`BA_EX_RELFRAC` work addressed on the other
   route - the KNOP path may need the equivalent.
3. **Ordered-march convergence rate unverified.**  The main
   paper-conform fix is theory-derived; the runtime A/B vs
   `KN_RAY_ORDERED=F` has never been run.
4. **OBS render unusable.**  D1-D4 verified; fixes C1-C5
   unapplied; and it cannot be validated until blocker 1 clears.
5. **Outer-map expansiveness (A2) unresolved.**  If the
   population-coupled lagged map is genuinely expansive on
   compression cells, no damping schedule fixes it - the remedy
   is structural (BLOCK_GS design, vault `audits/BLOCK_GS_DESIGN.md`,
   still a design doc).
6. **Integration-test vacuum.**  No full-iteration tests
   (F-23/T-08..T-11); every gate is tuned to one model family's
   failure modes.
7. **Cold-start divergence.**  Every viable attempt warm-starts
   from SOB; a production pipeline cannot depend on that.

## 8. Path to production (if KNOP is ever revisited)

Ordered by leverage:

1. Confirm or kill H1 with the `POP_RUNAWAY_*` census on the
   next KNOP attempt; if confirmed, port the SOB floor/relfrac
   remedies (`BA_DEM_PFLOOR`, `BA_EX_RELFRAC`, quarantine) to the
   KNOP route.
2. Run the pending `KN_RAY_ORDERED` A/B on a compressing model -
   it is the only untested convergence mechanism of size.
3. Get one cold-start convergence on m15 (or declare warm-start
   a hard requirement and document it).
4. Apply OBS fixes C1-C5, then render from the converged
   checkpoint and compare against the SOB `blank` render and the
   Dessart reference.
5. Build the F-23 integration tests before claiming any gate
   behaviour generalises.

Until then the honest status tag for `DO_NONMONO_KNOP` is
**#status/unverified - fails to converge on the m15 benchmark**.

## Sources

- `docs/audits/KNOP_AUDIT_2026-10-30.md`,
  `KNOP_FORMAL_SOLVER_AUDIT.md`, `KNOP_STALL_DIAGNOSIS_2026-10-31.md`,
  `FIELD_DRIFT_AUDIT.md` (superseded), `FG_FIELD_DRIFT_AUDIT.md`
  (authoritative), `KNOP_ITEM1_RUNAWAY_2026-11-08.md`,
  `ITEM2_POPULATION_FLOOR_AUDIT.md`, `POPULATION_SOLVER_FIXES_2026-11-09.md`,
  `NONMONO_COMPREHENSIVE_AUDIT_2026-10-30.md`,
  `NONMONO_NUMERICAL_ARCHITECTURE_AUDIT_2026-11-03.md`,
  `CONVENTIONAL_SOLVER_AUDIT_2026-11-03.md`,
  `FULL_ITERATION_CYCLE_AUDIT.md`,
  `JFNK_PRODUCTION_OPERATOR_INVESTIGATION_2026-11-07.md`,
  `CMF_FLUX_KNOP_CONVERGENCE_AUDIT.md`,
  `KNOP_SECOND_PASS_AUDIT.md` + `KNOP_SECOND_PASS_TODO.md`
  (mechanism corrections: `|Q|->1` strong coupling, not near-zero-GAM)
- Vault: `audits/KNOP_CONVERGENCE_STRATEGY.md`,
  `audits/KNOP_MATH_AUDIT.md`, `audits/BLOCK_GS_DESIGN.md`,
  `audits/MOM_KNOP_STIFFNESS_AUDIT.md`, `solvers/Knop formal solver.md`,
  `solvers/KNOP_JFNK.md`, `history/Open issues.md`,
  `history/outstanding_issues.md`, `history/convergence_issues.md`
- Live run evidence: `cmfgen_validation/m15_2p82d_FRESH_knop`,
  `m15_2p82d_FRESH_knop_d1a`, `m15_2p82d_VAL_knop_formal`,
  `m15_2p82d_VAL_knop_mom`, `m15_2p82d_KNOP_opt*`
