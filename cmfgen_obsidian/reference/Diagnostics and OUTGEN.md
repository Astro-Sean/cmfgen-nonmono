# Diagnostics and OUTGEN

What to grep when watching a run. Section hierarchy: double `***`
lines separate `Current great iteration count` blocks; `====` headers
mark standalone major reports; solver, field, population-update, store
and summary sections nested inside an iteration use `----` subheaders.
Output conventions: since commit
`4360051` diagnostics print **description-as-label** tables:

```text
    max SE residual, current call   :  4.0796E+10
    residual before last full step  :  2.6668E+12
    acceptance test armed this call :  T
```

(not the older `res : value\n  (description)` layout). Bordered
`----` tables wrap ITER_SUM, refresh pass/exit blocks, GS sweep
diagnostics, MOM_TRIAL terms, quarantine, and requeue reports. The
frozen-field sweep block uses nested key/value lines:

```text
KNOP inner-GS sweep
  immature-field iteration : 5
  FG / MOM maxchg          : 1.9362E-03 0.0000E+00
  MOM bulk
    nbad5 / nbad50         : 0 0
    negJ / quar            : 0 0
    ntot / pquar           : 0 0
  MOM sweep
    maxchg                 : 0.0000E+00
    nbad5 / nbad50         : 0 0
    negJ / quar            : 0 0
  LAMBDA correction cap
    churn fraction         : 4.6495E-03
    max correction         : 3.0000E+00
```

## ITER_SUM (per iteration)

```text
ITER_SUM: it=  mode= -> ba= frozen= held=
    pops           demand / applied correction
    se residual    all cells / excl. worst cols
    ba             trust / consec rejects
    lum ratio      inner/outer luminosity
    fg             bad= churn= cap=  lagged field drift
    mom            n5= n50=/ntot neg= quar= pquar= cap=
    jfnk           demand vs frz gate, hits frz=N/M   [OFF|DEFER|ARMED]
    jfnk_cpl       demand vs cpl gate, hits cpl=N/M, coupled=T/F
```

## JFNK blocks

`JFNK DEFER` explains each deferral: warmup deferral, GS field not
ready, required mature-sweep count, field-change metric,
Lambda-designated iteration, demand and contraction status.

`JFNK ADMIT` on admission: current demanded population change, frozen
and coupled gates, consecutive hit counts, selected tier
(`jfnk_frz` / `jfnk_cpl`).

## ITER_TYPE banner (per great iteration, 08-Oct-2026)

Once per great iteration (`MAIN_COUNTER`, reset on NEWMOD) the SE
dispatch prints a bordered record of the selected solve type:

```text
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ITER_TYPE it=   12  Running JFNK-FRZ ITERATION (Newton-Krylov, frozen J)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```

Variants: `LAMBDA ITERATION (fixed-point step)`,
`FULL ITERATION (BA linearization step)`, `JFNK-FRZ` (frozen J),
`JFNK-CPL` (KNOP only; J re-solved per eval). Additional SE calls
within the same great iteration (e.g. the T-update solve) reuse the
first dispatch type. The banner records the dispatch decision; if
JFNK subsequently falls back to LAMBDA inside the solve, the fallback
blocks below name the actual outcome.

Also: JFNK context/setup/`F(0)` validation/iterations/results/commit/
completion blocks, and `BA_RES_ACC` acceptance/rejection events.

Once per JFNK solve a formulation block names the unknown, whether
pseudo-time is on, starting dtau, forcing floor. Per accepted Newton
iteration: one line with old/new dtau and |F|. Per-frequency sweep
text goes to `JFNK_EVAL.log`, not OUTGEN.

## Grep targets

```text
Current great iteration count - *** iteration boundary
JFNK gate: deferred -- warmup, call N of M
JFNK gate: deferred -- GS field not ready, nmature= FG maxchg= MOM n50=
JFNK gate: deferred -- LAMBDA-designated iteration
JFNK gate: last demanded MAXCH= vs gate= hits= -- deferred
JFNK admission: MAXCH= hits frz= cpl= coupled=
ITER_TYPE it=                          - per-iteration solve dispatch banner
JFNK: reusing CMAT_RAD (age N)         - frozen-J operator cache hit
JFNK re-armed                          - stagnation disable cooldown expired
BA_RES_ACC              - acceptance / rejection events
SOLVEBA POPULATION UPDATE - held variables, max corrections, DEC_VEC/INC_VEC, scale
FG_KNOP_REFRESH: pass   - per-pass residual history
ray exits: exact= tol= bulk= stall= cap= div=  - per-ray stop-rule census
stall/cap rays: maxres= bad cells=  - share of unconverged cells left by early exits
MOM_KNOP_REFRESH guard  - why the MOM refresh skipped
MOM field marginal      - trial machinery arming
Calling FG_J_CMF_NONMONO_V1 / MOM_J_CMF_NONMONO_V1 / MOM_JREL_V9
no compressing points - GS storage skipped
SOLVE_FOR_POPS: KNOP GS warm-up / KNOP inner-GS sweep / immature-field iteration
SIGMA_TRUE_IN is all zeros   - missing RECOMPUTE_SIGMA_TRUE warning
```

## `iter_reason` tags (ITER_HISTORY / mode lines)

```text
gs_warmup            requeue-warmup forced Lambda
gs_immature_demote   scheduled BA demoted on unsafe field
gs_freeze_forced     freeze-controller Lambda
mom_trial_release    marginal-release BA probe
jfnk_frz / jfnk_cpl  admitted JFNK tier
jfnk_defer           JFNK eligible but gate not met
hold_pops            GS_HOLD_POPS snapshot/restore iter
lam_stable_exit      merit-based Lambda->full transition
lum_drift_hold       luminosity gate hold
maxch_above_enter    demand above LAM_ENTER -> stay Lambda
rel_stable_hold      N_MATURE < GS_REL_STABLE
cnt_lam_escape       NUM_LAM escape leash forcing a full step
ba_demand_hold       BA_ENTER_RAW raw-demand hold (SOB)
```

## Side files

| File | Content |
|---|---|
| `OUTGEN` | main log |
| `ITER_HISTORY` | one row per iteration: mode, frozen/held, raw/applied correction, se/se2, trust, lum, drift - post-run plotting |
| `ITER_SUM` | tagged block per iteration |
| `GS_STATE` | wrapper-facing freeze state each call (`n_frozen`, `n_trial_acc`, ...) |
| `JFNK_EVAL.log` | per-frequency eval chatter during JFNK sweeps |
| `CONVERGENCE_HISTORY` | MAXCH series (frozen iters report the field drift as MAXCH) |
| `.CMFGEN_KNOP_STALLED` | dropped by the wrapper at the stall bound |

## Healthy vs sick signatures

- Healthy: `gs_warmup` -> `gs_freeze_forced` Lambda steps with
  contracting demand -> live BA -> `JFNK gate` defers with growing
  hits -> `jfnk_frz` -> `jfnk_cpl` near the root.
- Sick: `jfnk_cpl` while `frozen=T`; refresh residual pinned at
  ~1e3-1e5 for many iterations; `pops` demand pinned at 1e7 sentinel;
  lum ratio swinging; `n_trial_acc` climbing without accepted trials.

## Related

- [[Checkpoints and state files]]
- [[Failure modes and guards]] - the signatures' meanings
- [[Field maturity and GS freeze]]
