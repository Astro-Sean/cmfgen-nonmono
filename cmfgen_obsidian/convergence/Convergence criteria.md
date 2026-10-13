# Convergence criteria

What the code tests to declare convergence, what it should be testing,
and the diagnostics that reveal stalls before the exit test ever
fires.

## The exit test

In `cmfgen_sub.f`, after `COMP_J_CONV`:

```text
converged when:  MAXCH < EPS   (EPS = 1%)
  AND            LAST_LAMBDA != MAIN_COUNTER   (not a Lambda step)
  AND            MAIN_COUNTER != LAST_NG       (not the NG iteration)
  -> NUM_ITS_TO_DO = 1 : one final loop writes the model data
```

`MAXCH` is the maximum demanded population correction this iteration.
The two guards exist because both a forced-Lambda step and an Ng
acceleration are kicks - the model cannot be called converged on the
iteration that just moved it.

## Caveats on `MAXCH`

- `DO_LEVEL_CHK` reports the 10th-largest change, so up to ~10 levels
  can still move >10% while `MAXCH` reads <1% ([[The great
  iteration]]). The top-10 vectors now respect `BA_DEM_PFLOOR` /
  `BA_DEM_RELFRAC` (13-Oct-2026): before that, dead levels exempted
  from the demand scan could still fill all 10 slots and pin `MAXCH`
  (X_SOB 40d: HeI singlet dead levels at depth 73, `MAXCH` ~9%).
- `BA_DEM_PFLOOR` / `BA_DEM_RELFRAC` (default 0 = off) exempt
  floor-level population rows from the demanded-correction metric;
  their corrections still commit. Use only where the exempted rows
  are demonstrably opacity-inert (see [[BA matrix and SOLVEBA]]).
- Masked-demand false convergence (fixed 04-Oct-2026): the resacc
  trust mask used to erase failed columns' demand from `MAXCH`, so a
  fully masked step reported demand 0 - jfnk_hyb exited "converged"
  at it8 with 69/70 columns masked, `se_res=1.1e14`, `L_outer=0`.
  `MAXCH` now keeps masked demand (`STEQ_MASKED_MAX`); `MAXCH_APP`
  reports what was applied.
- `MAX_CHNG_LIM=1e100` (these configs) disables the runaway-abort
  guard - nothing stops a diverging run by itself. The stall
  watchdogs below are the substitute.

## Raw vs applied vs demanded

Three different numbers answer three different questions; mixing them
was a real bug source:

| Quantity | Meaning |
|---|---|
| `MAXCH` / `BA_LAST_RAW` | what the solver *demanded* (raw, pre-cap) |
| `MAXCH_CTRL` | what was *applied* after limitation (mode control) |
| `MAXCH_APP` | post-mask applied correction (diagnostics) |

Gates that must see the operator read `BA_LAST_RAW` (JFNK admission,
`BA_ENTER_RAW`); gates that must see the committed step read
`MAXCH_CTRL` (Lambda hysteresis).

## Pipeline tracker false-convergence veto

The external convergence tracker does not count a rejected/masked BA
record as evidence that corrections are below tolerance. It also
requires a nonzero correction above `1e-12`; zeroed corrections from a
fully masked step therefore cannot produce `CONVERGED`. The CMFGEN
iteration exit test above remains separate from this post-run
classification. The tracker tests include the real hybrid
`ITER_HISTORY` false-convergence case and a synthetic converged case.

JFNK commits need an additional status check: an improving best iterate,
rebuilt baseline, or clipped correction is progress but not proof that the
nonlinear target passed. `JFNK_UNCONVERGED_COMMIT` blocks CMFGEN_SUB's
MAXCH-based temperature release, staged source/solver ramps, and model
finalization for those cases. An unclipped `JFNK_INFO=0` commit may
pass the outer gates. The final
population-change test still uses `MAXCH`; the 50-percent administrative
handoff gates use `MAXCH_CTRL`.

## Field-side maturity metrics

The `GS_*` / `MOM_SW_*` sweep-drift metrics (`GS_MAXCHG`,
`GS_NBAD`/`NBAD5`/`NBAD50`, `GS_NCONV`, `GS_CAPPED`, `MOM_SW_NEGJ`,
`quar`) are the convergence observables for the lagged field map -
table in [[Field maturity and GS freeze]]. Two rules of thumb from
the run log:

- FG refresh residual and store churn answer different questions (map
  self-consistency vs movement from the field populations saw). `GS_NBAD`
  counts with the same normalized residual used by the ray solver;
  `GS_NCONV` measures store churn. Do not mix them in gate logic or
  reporting. MOM refresh and `MOM sweep` metrics likewise measure
  different fields.
- Quarantined cells are excluded from `nbad5/nbad50/maxchg` -
  otherwise the residual is pinned by cells that can never update.

## Stall detectors

- `GS_FREEZE_LIM` (default 0=off): in-model hard stop after N
  consecutive frozen iterations.
- `KNOP_STALL_MAX` (60): slurm-wrapper stall bound, reads `GS_STATE`.
- Plateau detector: 6-iteration ring of applied correction, second SE
  channel, drift fraction; <1.5x window flags `CONV_PLATEAU`.
- `BA_NREJ_LIM` (20): consecutive-rejection hard stop.
- `JFNK_MAX_STAGNATE` (3): repeated-failure JFNK self-disable.

## What "converged" still requires

From [[convergence_issues]] Sect. 5 (Type IIn 10 d model, open):

1. BA kick magnitudes decreasing across cycles (currently
   increasing 3.4e3 -> 2.7e4).
2. `nbad5` sweep drift staying below the 10% gate - the field mature
   without re-freezing after BA.
3. Bounded or declining `sweep maxchg` (flickers to ~50 on single
   cells).
4. Population corrections under `EPS` in consecutive iterations, not
   just recovering after each kick.
5. An OBS-stage spectrum (requires `USE_LAM_ES=T` for `cmf_flux`
   regardless of the model-run path).

## Wreck signatures

Reading diagnostics for divergence vs stall: `MAXCH` pinned at the
1e7 sentinel = saturated demand; `res` flat with `MAXCH` pinned = the
resacc cliff; `LUM_RATIO_OI` drifting = inner-boundary flux still
swinging. Full list in [[Failure modes and guards]] and
[[Diagnostics and OUTGEN]].

## Related

- [[Convergence strategy]] - the staged plan
- [[Field maturity and GS freeze]] - metric definitions
- [[Failure modes and guards]] - observed failure signatures
- [[Diagnostics and OUTGEN]] - where each number is printed
- [[Open issues]] - which criteria are currently unmet
