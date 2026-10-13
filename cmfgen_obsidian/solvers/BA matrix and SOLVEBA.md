# BA matrix and SOLVEBA

The population-update machinery: the "BA" matrix is CMFGEN's
approximate Jacobian of the statistical-equilibrium block.

## Structure

- Built by `GENERATE_FULL_MATRIX_V3` inside `CMF_BLKBAND_V3`, driven by
  `SOLVEBA_V13` (`new_main/subs/solveba_v13.f`).
- Block per depth of `NT` unknowns: level SE rows, ionization
  equations, species number-conservation rows (`EQ_SPECIES`), charge
  conservation (`NT-1`), radiative equilibrium (`NT`).
- `SOL_METH`: `DIAG` (block-diagonal) or `TRIDIAG` (block-Thomas over
  `NUM_BNDS`). **In this configuration all off-diagonal depth blocks
  are exactly zero** (VJ=0 under `USE_LAM_ES=T`; dZ coupling is
  VJ-multiplied; `INC_ADV=F`; diffusion off) - TRIDIAG is numerically
  identical to DIAG. Verified bitwise (AGENTS.md TRIDIAG audit).

## The two update modes

| Mode | Entry | Update | Cap |
|---|---|---|---|
| Lambda | `LAMBDA_ITERATION=T` | fixed-point relax via the SE residual | `MAX_LAM_COR` (production `1e10`) |
| Full BA | `LAMBDA_ITERATION=F` | solve `A x = F`, apply `POPS *= (1 - x*SCALE)` | `MAX_LIN_COR` |

Mode control (in `SOLVE_FOR_POPS`): exits Lambda only after
`LAM_N_STABLE` consecutive calls with applied correction
`MAXCH_CTRL < LAM_EXIT`; re-enters on one call at `LAM_ENTER`.
`MAXCH_CTRL` is the *applied* (post-cap) correction. `FRESH_BA_EXIT`
rebuilds the BA on exit; `N_FULL_TREL` phases the temperature release
(`T_HELD_PHASE`).

## Residual acceptance ("resacc")

`DO_BA_RES_ACC=T` - **mandatory on KNOP runs**:

1. Snapshot `POPS_BA_SAVE` (+ per-depth residual history
   `BAACC_COL_PREV`/`BAACC_COL_MIN`).
2. Take the BA step.
3. Re-evaluate the SE residual at the next call.
4. Per-depth verdict (`BA_REJ_FRAC`, default 0.25): each depth column
   is tested against its own `BA_RES_GROW`/`BA_RES_HARD` bounds.

   - **Masked rollback** (<= `BA_REJ_FRAC`*ND columns fail, RES2
     quiet): only the failing columns are restored, their fresh
     corrections are zeroed (`BA_RBK_MASK`), the rest of the step
     stands, and the call continues to a normal solve. No LAMBDA
     demotion (`ba_partial_rej`, `ba=MASKED` in ITER_SUM) - the
     observed failures are local stiff cells (H2:lev1@d70, CV@d13),
     and whole-state rejection threw away mostly-good steps.
   - **Whole-state rollback** (RES2 bulk channel fails, or >
     `BA_REJ_FRAC`*ND columns fail): restore all of `POPS`, discard
     the call's solve, demote to LAMBDA - unchanged legacy path.
   - `BA_REJ_FRAC<=0` disables the mask: any failed channel rolls
     back whole (legacy behaviour).

- Trust radius `BA_TRUST_CUR` halves on rejection (whole or masked)
  (`BA_NREJ_CONSEC`); `BA_NREJ_LIM` (default 20) hard-stops a run that
  alternates Lambda/rejected-BA forever. Two 07-Oct-2026 accounting
  fixes bound what the counter may consume (m15 2.82d stopped at it 75
  after 20 counted events, ~15 of which were residual-improving or
  population-frozen calls):
  - **Held calls are not adjudicated.** While `GS_POPS_HELD` restores
    `POPS` after the call (GS freeze / luminosity hold), the whole
    judgement block is skipped and `BAACC_PENDING` stays armed, so the
    armed step is measured at the first call that can actually move
    the populations. Five of the twenty rejections were iterations
    that could not possibly pass - pops were held, so the armed
    residual re-measured the same frozen state.
  - **Progress exemption.** A masked or soft-rejected step that still
    improved either residual channel (`BA_SE_RES2` or `BA_SE_RES` vs
    the armed pre-step values, with the same saturated-max exemption
    as `REJ_MAX`) resets `BA_NREJ_CONSEC` to zero and leaves the trust
    radius unchanged - it is the partially-committed case the mask
    exists for. A masked step with no channel improvement still
    counts and still halves the radius, so `BA_NREJ_LIM` keeps its
    meaning as a stall detector rather than a cap on total armed
    steps.
- Per-equation quarantine (`DO_BA_QUAR`, default F; `BA_QUAR_MAX`,
  default 25): masked acceptance and the trust mask both operate on
  whole depth columns, but the observed veto on the m15 2.82d SOB
  run was a handful of stiff equations - trace inner-cell rows (CV
  1s_2s_1Se / 1s_2p_1Po at d68-69, the NeV ground-term and HII
  ionization rows at d70-75) demanding >=99.999% annihilation every
  armed step, which pinned `BA_TRUST_CUR` at its floor and zeroed
  ~60/76 columns per solve. With `DO_BA_QUAR=T` a measured-failure
  (equation,depth) pair is parked instead of vetoing the column:
  its row is excluded from `COL_RES`/`BA_SE_RES`/`BA_SE_RES2`, and
  its correction is zeroed on non-LAMBDA solves before the trust
  test. The set is populated only by measured failure - every
  over-trust population variable in a trust-masked column, and the
  worst residual equation of each rollback-masked column. ED
  (J=NT-1) and T (J=NT) are never parked. LAMBDA calls still update
  parked variables, so a parked level keeps contracting under
  fixed-point iteration and its Lambda-call demand still counts
  toward MAXCH - no false converged exit. `BA_QUAR_MAX` bounds the
  blind spot (a parked equation is invisible to acceptance) and
  <=0 disables growth; parked demand is reported per call (count +
  max suppressed |STEQ|) in BA13-TRACE. Inert unless DO_BA_RES_ACC;
  default path unchanged.
- `BA_SNAPSHOT_MISSING` demotes an unchecked step
  (`ba_no_snapshot_demote`, KNOP gate only) and vetoes MOM trials -
  the rollback snapshot is part of the promise. Coverage gap found
  in the cycle audit: the demote sits inside the DO_NONMONO_KNOP
  block and only during the initial RD_LAMBDA window, so under
  DO_NONMONO_SOB (and post-handoff under KNOP) a missing snapshot
  is only flagged after the fact via `BA_RAW_REJECTED` - the
  unchecked step itself already committed.
- Rejection is armed only on BA-designated calls. Ordinary Lambda
  follow-ups do not roll the armed step back (undo-after-Lambda can
  re-issue the same correction and form a FULL/LAMBDA cycle); measured
  degradation halves the trust radius and counts `BA_NREJ_SOFT`, also
  incrementing `BA_NREJ_CONSEC` so `BA_NREJ_LIM` bounds committed bad
  steps. Exception: a degraded `BA_PROBE_LAM` trial is a bounded
  measurement inside RD_LAMBDA, so its Lambda follow-up rolls the step
  back rather than soft-rejecting it. Since 13-Oct-2026 the rollback is
  tiered like any other armed step: a bulk failure (RES2 degraded or
  `> BA_REJ_FRAC`*ND columns) restores the full pre-probe population
  and Knop radiation snapshot and skips the stale Lambda solve, while
  a few-column failure takes the masked path - only the failing
  columns are restored, the kept columns' progress stands, and the
  call proceeds to its Lambda solve with masked corrections zeroed.
  The m15 SOB probes were degrading a single innermost column ~20x
  while improving the other ~69 columns ~45%; the blanket full
  rollback discarded that progress every `BA_PROBE_LAM` calls.
  A clean or masked probe still counts as rejected for the
  `BA_ENTER_RAW`/JFNK gates. This guard is shared by
  `DO_NONMONO_SOB` and `DO_NONMONO_KNOP`/JFNK; the stock path is
  unchanged.
- Masked-demand accounting (04-Oct-2026): `MAXCH` is now the max of
  the post-mask demand and the largest column demand the trust or
  roll-back mask erased (`STEQ_MASKED_MAX`, capped at the 1e7
  sentinel). A column whose correction was refused is not a satisfied
  one - previously a fully masked step reported `MAXCH=0` and tripped
  the `MAXCH<EPS` converged exit in `CMFGEN_SUB` (jfnk_hyb stopped
  "converged" at it8 with 69/70 columns masked, se_res=1.1e14,
  L_outer=0). `MAXCH_APP` stays post-mask: it correctly reports what
  was applied.

## Additional guards layered on

- **Saturated-demand demotion** (`BA_SAT_RAW`/`BA_SAT_LIM`): a full
  step whose raw demanded correction is still at the `1e7` sentinel is
  marching, not contracting - demote and halve the trust radius.
  Demand trend is the only pre-cliff signal: resacc can't see the
  cliff because `SE_RES2` stays flat while demand is pinned.
- **`BA_ENTER_RAW`** (non-mono paths): records raw demand from the
  latest non-LAMBDA step; the `ba_demand_hold` gate blocks LAMBDA->full
  while demand is high. With residual acceptance enabled, pending,
  rejected, masked, or soft-rejected BA trials also keep the gate closed
  until a periodic probe is residual-accepted. The same trial status
  blocks JFNK admission hits. Periodic `NUM_LAM` / `BA_PROBE_LAM`
  probes keep recovery bounded.
- **Luminosity drift gate** (`BA_LUM_GATE`/`BA_LUM_DCHG`): when
  `|LUM_RATIO_OI - prev| > BA_LUM_DCHG*|ratio|`, demotes a scheduled
  BA solve before correction, holds ordinary LAMBDA population updates,
  and blocks LAMBDA->BA exit and MOM_TRIAL. The held iteration advances
  transfer at fixed populations until the ratio settles. It is a
  stability test, not a value test; interacting-SN luminosities need
  not be equal at inner and outer boundaries.
- **`FIX_BA`/`N_FIX_BA`**: frozen-Jacobian reuse below `FIX_BA`%;
  `BA_MATRIX_STORED` records a real write so a mid-call flag flip
  cannot recall a matrix that never existed. `FIX_BA=15,N_FIX_BA=0`
  caused the 175-iteration stale-Jacobian plateau - use
  `FIX_BA=1.0,N_FIX_BA=10` (periodic recompute).
- **`FIX_IMPURITY`**: optional trace-ion freeze (<1e-15) during BA
  updates. Force-cleared on KNOP forced-Lambda iterations
  (`FORCE_LAMBDA` block) to preserve the upstream invariant. The
  pipeline profiles keep `FIX_IMP=F`: freezing the rows suppresses the
  population corrections needed to relax trace stages and leaves large
  CORRECTION_SUM bands at zero.
- **Dead-level demand floors** (`BA_DEM_PFLOOR`/`BA_DEM_RELFRAC`,
  11/12-Oct-2026): rows whose population is below an absolute floor,
  or below `BA_DEM_RELFRAC` times the largest population row at that
  depth, are skipped by the demanded-correction metric and by
  `MAXCH_APP`. Corrections on skipped rows still commit; ED (NT-1)
  and T (NT) rows are exempt. Motivation: SCALE-capped collapse of a
  Boltzmann-dead level is an unsatisfiable demand that pins MAXCH
  forever - the same mode for which the O II 2s2p4 superlevels were
  pruned (TypeIIn_KNOP2) and the X_SOB 40d run stalls at ~9% on HeI
  singlet levels at depth 73. **13-Oct-2026 fix**: the floors now
  also gate the `MAX_INC_VEC`/`MAX_DEC_VEC` accumulation. Under
  `DO_LEVEL_CHK` (all corrections <10%, the normal late-stage regime)
  `MAXCH` is the 10th-largest raw STEQ value, and before the fix those
  vectors were built ungated - dead rows filling all 10 slots kept
  `MAXCH` pinned regardless of the floors. Both floors default to 0;
  default path unchanged.

## SOLVEBA internals

- `CMF_BLKBAND_V3` -> `GENERATE_FULL_MATRIX_V3` (+ `FIXPOP_IN_BA_V3`)
  -> `DGEEQU` equilibration -> `DGETRF` LU -> `DGETRS` solve.
- `IFAIL` from LAPACK is tested in `CMF_BLKBAND_V3`; a singular BA is a
  **hard STOP** (V6's per-depth singularity recovery was deliberately
  not ported - surface the failure).
- `STEQ` sign convention: residual `F`, correction
  `POPS*(1 - STEQ*SCALE)`; saturated corrections get zeroed by resacc
  pre-filtering (`|STEQ|>1e4` rows).

## JFNK connection

Under `USE_JFNK=T` the JFNK solve **replaces the SOLVEBA_V13 call** on
admitted iterations; on any failure it falls back to a Lambda-mode
SOLVEBA safety step. The BA matrix survives as the JFNK
**preconditioner** (`JFNK_FACTOR_PRECOND_V1`: per-depth
equilibrated LU of `BA_full + ridge`, null rows -> identity, Jacobi
fallback). See [[JFNK solver]] and [[JFNK admission]].

## Related

- [[The great iteration]]
- [[Field maturity and GS freeze]]
- [[Failure modes and guards]]
- [[VADAT keys]]
