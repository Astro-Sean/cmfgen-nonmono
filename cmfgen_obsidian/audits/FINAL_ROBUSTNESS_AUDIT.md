# Final robustness + numerical-stability audit — cmfgen_jfnk

Date: 03-Oct-2026.  Tree: `cmfgen_jfnk` @ `1ac3086` (devin/jfnk-bringup).
Companion to [[JFNK_SOLVER_AUDIT]] (performance/robustness) and the
order-of-operations fixes committed as `9956744`, `e8abf97`,
`900d515`.  Scope: (a) re-review of those fixes against the code,
(b) audit of the newest solver-behavior features not previously
audited (BA_REJ_FRAC masked rollback, BA_DAMP_CAP, BA_ENTER_RAW,
BA_PROBE_LAM), (c) a robustness pass over the dev-diff surface
(NaN handling, alloc guards, sentinels, VADAT validation).

## A. Order-of-operations fixes — verified

| Check | Result |
|-------|--------|
| `JFNK_SAVE/RESTORE_MOD_STEQ_V1` coverage | PASS — `STEQ_T`, `STEQ_ED`, `STEQ_T_EHB`, per-ion `SE%STEQ` (offset-indexed), `SE%STEQ_ADV`, `SE%T_EHB` (flat `(ID-1)*ND` layout, `MIN(size,ND)` guard) |
| Buffer allocation | PASS — `NION*ND` flat buffers, realloc on size change, `MAX(1,...)` floor, dealloc in both init-cleanup lists, zeroed on alloc |
| Save→restore bracket integrity | PASS — `JFNK_RESIDUAL_V1`: both early `RETURN`s precede the save, no exit between save@~1535 and restore@~1739. `JFNK_BUILD_BA_RAD_V1`: only `STOP`s (alloc failures) precede the save. `JFNK_REBUILD_BASELINE_V1` and the driver frozen-split bracket are adjacent pairs |
| Coupled no-callback guard | PASS — identical sentinel to the budget path (`Y=1e20`, `jfnk_resid_nbad=N`, eval counted), returns before any state mutation, forces F(0) abort → LAMBDA fallback |
| `JFNK_DISABLED` path | PASS — sets the LAMBDA flags (matching the JFNK-failure fallback block) then clears `DO_JFNK_NOW`; falls into the shared ELSE `SOLVEBA_V13` call. `FORCE_LAMBDA` dead var removed |
| `ACCESS_F` cursor guard | PASS — park@5064 / restore@5941 both guarded on `JFNK_EDD_START_REC > 0`; the `EDD_CONT=F`+`ACCURATE=F` combination can no longer position a negative record |
| Fixed-form line length | PASS — no line >132 cols in any edited file |
| Default-path equivalence | PASS — every change sits inside `USE_JFNK`-gated code or is comment-only |

## B. Recent solver-behavior features — verified

### BA_REJ_FRAC masked per-depth rollback (`3432b53`, default 0.25)

- Per-column failure uses the same GROW/HARD bounds as the global
  channels, with the `BA_RES_SAT` exemption applied per column and a
  `0.5*HUGE` clamp floor for nonfinite entries.
- Bulk veto (RES2 = K-th largest, K=ND/10) still rolls back the whole
  state; `NFAIL > NINT(FRAC*ND)` also escalates to whole rollback —
  the mask cannot hide a grid-wide failure.
- Rollback ordering is correct: `POPS(1:NT,I)` restored per masked
  column before the linear solve; `STEQ(:,I)` zeroed after the solve
  so masked depths hold the restored populations; `COL_RES(I)` re-set
  to `BAACC_COL_PREV(I)` so re-arming does not inherit a rejected-
  state residual.  Under `DO_NONMONO_KNOP` the lagged `I_NEXT` store
  is restored whole (`RESTORE_FG_INEXT_LATEST`) — conservative and
  documented.
- Alloc chain is failure-consistent: `POPS_BA_SAVE` → `BAACC_COL_PREV`
  → `BAACC_COL_MIN` → `BA_RBK_MASK` chained on `IOS`; any failure
  leaves `BA_SNAPSHOT_MISSING=T` and the test falls back to the
  legacy scalar decision (`LCOL_OK=F`).  `LCOL_OK` ⇒ `POPS_BA_SAVE`
  allocated, so the rollback cannot silently no-op.
- Lambda follow-ups now *evaluate* the armed step (soft reject +
  trust halving) instead of clearing it as an untested acceptance —
  closes the invisible-bad-step channel.
- Live evidence: `SOB_mask` it39-40 shows
  `ba_partial_rej+ba_damped_full` / `ba=MASKED` /
  `BA_RES_ACC PARTIAL REJ` firing as designed.

### BA_DAMP_CAP (`2cb4215`, default 0 = off)

- Factor cap `MIN(MAX_LIN_COR, BA_DAMP_CAP)` applied only on
  non-LAMBDA calls while `BA_LAST_RAW >= BA_SAT_RAW` — a capped
  applied step cannot bypass the demand test since the gate reads
  raw `MAXCH`, not `MAXCH_CTRL`.
- Bounded leash: `BA_SAT_CNT` increments on consecutive saturated
  non-LAMBDA calls (keyed on entry mode `LAM_WAS`); `ba_sat_demote`
  fires before the damped-full branch, so run length in damped-full
  is bounded by `BA_SAT_LIM` (default 8 when the cap is armed).

### BA_ENTER_RAW + BA_PROBE_LAM (`15cc756`, `b4e2b3c`, `b4b7049`)

- `BA_LAST_RAW` recorded on every non-LAMBDA call (unconditional
  telemetry; gate inert without the VADAT key).
- Merit-exit ordering correct: `lum_drift_hold` → demand hold
  (`ba_demand_hold`, resets `CNT_STABLE_LAM` so the `NUM_LAM` escape
  still works) → `lam_stable_exit`.
- Held-gate probe releases exactly one full step after
  `BA_PROBE_LAM` held calls with the same actions as a merit exit;
  demotion resets `CNT_LAM` so probing cannot spam.
- Hysteresis semantics preserved: probes and exits reset
  `CNT_STABLE_LAM`/`CNT_LAM` symmetrically; `T_HELD_PHASE` /
  `N_ACC_FULL` handling identical across exit kinds; masked
  (`BA_PART_WAS`) steps do not advance the T-release counter.

### VADAT validation (`rd_control_variables.f`)

- `BA_REJ_FRAC` clamped to ≤1.0; `<=0` restores legacy behavior.
- `BA_SAT_LIM`, `BA_ENTER_RAW` (default 1e30 = off), `BA_PROBE_LAM`
  (0 = off), `BA_DAMP_CAP` (0 = off) all optional with documented
  defaults; `BA_DAMP_CAP` warns when armed without a `BA_SAT_LIM`.

## C. Robustness pass

- NaN handling: pervasive self-comparison checks (`X.NE.X`) in the
  residual/JVP/preconditioner paths, `jfnk_resid_nbad` component
  masking, `jfnk_nan_rejects`/`JFNK_N_RESIDUAL_NAN` counters,
  `1e20` sentinel convention consistent across eval-budget, NaN,
  and no-callback paths.
- Newton core: `jfnk_x`/`jfnk_fx`/`jfnk_tol` restored as a triple on
  every rejection path (verified in the earlier audit; unchanged).
- `JFNK_INIT_V1` resets counters/readiness flags and reallocates on
  grid or ion-set change.
- Non-JFNK dev code reviewed for the standard traps: no new
  uninitialized-local reads found; alloc chains all check `IOS`;
  division sites in new code are guarded or use documented
  denominators.

## D. Defects found in this audit

- **Doc-level only**: widespread misdated comments (`05-Oct-2026`,
  `07-Oct-2026`, `4-Oct-2026`, `10-Oct-2026` on Devin-dated blocks)
  — all the work is 02-03-Oct-2026.  Fixed in `514f871`, `1ac3086`;
  upstream `Altered` lines (1989/1995) and Sean Brennan's ported
  V6 stamps left untouched.  No code defect found.

## E. Residual risk / open items (unchanged from the order-of-ops audit)

- **L1**: eval-sweep side state (`RLUMST`, opacities, `OBS_FLUX`,
  `STEQ_T_*`, `TWO_PHOT` rates, cursor-adjacent locals) is not
  restored — verified self-healing today (all consumers read before
  `SOLVE_FOR_POPS`); fragile to any future post-solve reader.
- **L4**: a non-fatal return path inside an eval would leave
  `JFNK_IN_EVAL`/`FREEZE_I_NEXT`/`COMPUTE_BA`-hold set.  All current
  exits are `STOP` or clean; only matters if a soft-error path is
  ever added.
- `IOS_EVAL` unchecked at `cmfgen_sub` ~5092/5941 (file-positioning
  read in the eval sweep).
- **Runtime validation still pending**: the resubmitted JFNK_opt /
  JFNK_hyb runs are the first live exercise of the save/restore +
  self-disable fixes.  Do not claim "working" until a coupled
  activation (or a stagnation-triggered disable) is observed in
  OUTGEN.

## Verdict

Source-level: **PASS WITH LIMITATIONS**.  The dev-diff solver surface
(JFNK, BA resacc/masked rollback, merit-exit gates, damped-full) is
internally consistent, fail-safe on its degradation paths, and
flag-gated so the default path is untouched.  The limits above are
documented fragilities, not live defects.  Numerical-stability
claims still await the running models' behavior at JFNK activation.
