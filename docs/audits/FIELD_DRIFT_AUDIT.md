# FIELD DRIFT / RADIATION FIELD CONVERGENCE AUDIT

Tree: `cmfgen_jfnk` (cmfgen-nonmono). Scope: the reported
field-drift stagnation on the Knop non-monotonic path
(`DO_NONMONO_KNOP`, `DO_SOB_FG_KNOP`, `DO_NONMONO_SOB` moment path).

## 1. Executive summary

The number reported as "field drift" (the `GS_MAXCHG` /
`MOM_GS_MAXCHG` channel that feeds the maturity gate and the
"largest drift" prints) is the **max-norm of the unrelaxed residual
of the lagged frequency-coupling fixed-point map**, i.e.

```
D_inf = max_cells | G(J_stored) - J_stored | / scale
```

where `G` is one full formal solve (or refresh pass) with the
lagged neighbour field held fixed, and `scale` is a per-cell floor
on the stored/fresh magnitude.

A plateau of this number is, for the most part, **the expected
behaviour of the discretized map, not a measurement error**. The
lagged map is only marginally contractive at a small subset of
cells - strong compression rows where the upwind weight |Q| -> 1,
and the |1+sigma*K| -> 0 "GAM floor" resonance cells near the inner
boundary (documented in the solver comments and in
`docs/audits/KNOP_STALL_DIAGNOSIS_2026-10-31.md`). On those cells
the local map multiplier is ~1, so the worst-cell residual cannot
fall below a floor no matter how many outer iterations run.
`AGENTS.md` records this directly: "the lagged map never contracts
- sweep drift ~1.0 every sweep ... refresh stalls at residual
~1.1" (m15 SOB+moment run).

The audit did not find any of the classic measurement defects the
task asked to rule out: no array is compared with itself, no
previous-state buffer is overwritten before it is used, the
production store is advanced every iteration, and the drift is
evaluated against the *unrelaxed* residual (so damping does not
falsely shrink the reported number).

What the audit did find:

- The headline quantity silently **changes definition across the
  iteration**: `FG_GS_MAXCHG` is the sweep drift before the
  refresh and the refresh residual after it. A user watching one
  number is watching two different metrics stitched together
  (section 5, defect D2).
- `FG_RF_DCHG` (the honest "how much did the committed field move"
  channel) is **not reset at sweep entry**, so on non-refresh
  iterations it reports a stale value (defect D1; fixed).
- The adaptive relaxation factor is updated from the **inner**
  (frozen-source) map's contraction ratio but applied to the
  **outer** (population-coupled) write - a category mismatch that
  can hold omega near 0.5 indefinitely while the outer iterate
  could take larger steps (section 7). This can produce *slow*
  convergence that looks like stagnation; it cannot create a false
  plateau by itself.
- JFNK evals cannot corrupt the production field stores
  (`FREEZE_I_NEXT` + scalar save/restore; MOM store unreachable
  because JFNK admission requires `USE_LAM_ES=T`). Coupled-JFNK
  evals use a *frozen* lagged store: the Newton residual is a
  lagged-J operator, not the fully coupled one (section 8).

Diagnosis vs the task's hypothesis list: primarily **D** (genuine
nonlinear fixed-point non-contraction on a minority of cells) with
a contributing **A** (the headline metric is the wrong instrument
for "is the production field converging" - it conflates inner-map
residual with outer progress) and a plausible secondary **C**
(adaptive-omega category mismatch slowing outer contraction).
B (update sequence), F (stale field state), and the
self-comparison failure modes were ruled out at code level; the
JFNK interaction (E) is guarded by design.

## 2. Field-drift data flow (verified, not assumed)

Per great iteration (`MAIN_COUNTER`), all under
`DO_NONMONO_KNOP` (or the opt-in SOB gates):

```
PREVIOUS STORED FIELD
  I_P/I_M_DEPTH_NEXT(kc,ls,s)   (FG, slot s = freq s+1 field)
  J_NEXT_FRQ(i,k), H_NEXT_FRQ(i,k)  (MOM, slot k = freq k field)
        |
        v
DO ML=1,NCF  frequency sweep            (cmfgen_sub.f, ~3200)
  COMP_J_BLANK -> FG_J_CMF_NONMONO_V1   (comp_j_blank.f 856)
    compressing cells read lagged slot FREQ_IDX_ITER+NM_RAY(I)
    writes I_*_DEPTH_NEXT = SNP + INEXT_RELAX*(I_fresh - SNP)
    drift: |I_fresh - SNP| / max(SNP, 1e-4*peak_fresh)   [fg ~4240]
           -> FG_GS_MAXCHG, censuses (NEW_FREQ only)
  COMP_J_BLANK -> MOM_J_CMF_NONMONO_V1  (comp_j_blank.f 1054;
                                        USE_LAM_ES=F only)
    compressing rows read J_NEXT_FRQ(i,k+MOM_NM_J) (lagged)
    writes J_NEXT_FRQ = J + JNEXT_RELAX*(XM - J)
    compressing-row writes SUPPRESSED (refresh owns them)
    drift: |XM - J| / max(|J|, 1e-4*peak_stored)          [mom ~3452]
           -> MOM_GS_MAXCHG, censuses
        |
        v
post-sweep refreshes                    (cmfgen_sub.f 3281-3293)
  FG_KNOP_REFRESH_V1: per-ray (NI x NCF) re-solve, alternating
    direction; inner march VNEW, under-relaxed FLD step
    (KN_RAY_OMEGA); residual per pass = max unrelaxed
    |VNEW-VOLD|/scale; commit WRITE_STORE (unrelaxed) /
    WRITE_STORE_CAP / reject; sets
    FG_GS_MAXCHG<-KN_RF_RESID, FG_SW_* snapshot, FG_RF_DCHG
  MOM_KNOP_REFRESH_V1: alternating-direction passes over the
    whole (ND_NM x NCF) field; JNEXT_RF_READ Jacobi snapshot
    (live GS on descending when MOM_RF_GS_DESC); per-cell
    MOM_DAMP_NM; AA acceleration; best-pass composite via
    ownership masks; commit or restore; sets MOM_SW_* snapshot,
    MOM_RF_DCHG, MOM_GS_MAXCHG<-BEST_CHG (unconverged)
        |
        v
STEQ assembly (two-photon, charge exchange, ...)
        |
        v
SOLVE_FOR_POPS                            (cmfgen_sub.f 4557)
  GET_FG_GS_MATURITY / GET_MOM_GS_MATURITY  (~1028-1042)
  GS_MAX_ALL = max(residual, sweep drift, retained delta)
  GS_RAW_IMMATURE -> GS_FREEZE (KNOP only) -> forced LAMBDA,
     optional GS_HOLD_POPS population-hold
  release audit GS_JUST_REL; counters N_MATURE / N_FROZEN_TOT
        |
        v
SOLVEBA_V13 (LAMBDA or BA step)  or  JFNK_SOLVE_CORR_V1
        |
        v
next great iteration: the stored field above is the "previous" one
```

The order that actually matters for the metric:

```
stored(x_n) --fresh solve--> G(x_n)
   |             |
   |--drift=|G-x|-->|          (unrelaxed; reported as sweep drift)
   |             v
   |--write x_{n+1} = x_n + omega*(G(x_n)-x_n)
   |
   `--refresh: converge inner map -> commit -> retained delta
```

## 3. The existing drift metrics, defined

### FG / ray side (`fg_j_cmf_nonmono_v1.f`)

- `FG_GS_MAXCHG` (in-sweep, ~4213-4330): max over compressing
  cells `(KC,LS)` and frequencies of
  `|I_fresh - SNP| / max(SNP, 1e-4*VMAX_fresh_ray)`.
  `SNP` is the stored column snapshotted at `NEW_FREQ`. Censused
  only on `NEW_FREQ`; the store write happens in the same loop
  *after* the metric, so the metric is the **unrelaxed**
  fresh-vs-stored residual. Argmax -> `FG_GS_ARGK/ARGR/ARGL`.
- `FG_SW_MAXCHG`, `FG_SW_NCONV` (refresh entry, ~6353): snapshot
  of the in-sweep values before `FG_GS_MAXCHG` is overwritten.
- `KN_RF_RESID` -> `FG_GS_MAXCHG` (post-refresh): max over
  non-diverged rays of the last pass residual
  `max(|VNEW-VOLD|/max(VOLD,1e-4*VMAX))` - the unrelaxed inner-map
  residual at exit.
- `FG_RF_DCHG` (~7017-7038): `max |committed - pre-refresh store|`
  with the same normalization - the **retained-field delta**,
  i.e. how much the production field actually moved.
- Censuses: `FG_SW_CAPPED`/`FG_SW_CAPTOT` (sweep capped writes),
  `FG_RF_CAPPED` (refresh), `FG_GS_NCONV`/`FG_GS_NBAD`/
  `FG_GS_NBAD_BLK` (>50% churn), `FG_GS_NDIV` (diverged rays),
  `FG_INEXT_CAPPED` (non-finite / capped upwind reads),
  `FG_GS_NSTAL`/`FG_GS_NPLAT` (stall/plateau cells).
- Per-iteration reset (`INIT` block, ~2256-2309): clears
  `FG_GS_MAXCHG`, `FG_SW_*`, `FG_GS_N*`, `FG_RF_CAPPED` - **but
  not `FG_RF_DCHG`** (defect D1).

### MOM side (`mom_j_cmf_nonmono_v1.f`)

- `MOM_GS_MAXCHG` (in-sweep, ~3426-3786): max over `(i,freq)` of
  `|XM_NM - J_NEXT| / max(|J_NEXT|, 1e-4*peak_stored_column)`; the
  H field is folded into the same max with a sign-aware rule (H
  crosses zero legitimately). Argmax -> `MOM_GS_ARGI/ARGL`.
- Censuses: `MOM_GS_NTOT=ND_NM*NCF`, `MOM_GS_NBAD5/NBAD50`
  (>5%/>50%), `MOM_GS_NHUGE`, `MOM_GS_NQUAR` (negative-J
  quarantine), `MOM_INEXT_CAPPED`. Compressing rows are excluded
  from the sweep write (refresh owns them) but still contribute to
  the drift census - so on those rows the metric reads
  |fresh - refresh-owned store|, which is nonzero by design until
  outer convergence.
- `MOM_SW_*` (refresh entry, ~5133): sweep snapshot; `MOM_SW_NTOT`
  reset each iteration (~1669) so `GET_MOM_GS_MATURITY`
  (~4427-4525) returns the snapshot when a refresh ran and the
  live `MOM_GS_*` otherwise - the earlier stale-census mode is
  closed by the reset.
- `MOM_RF_DCHG` (~6414-6447): retained vs `JNEXT_SAVE` delta, same
  normalization; non-finite -> `1.0E+30`. `MOM_RF_CONV` =
  converged flag; `MOM_RF_RAN` = refresh ran this iteration.
- Per-iteration reset (~1640-1679): full clear including
  `MOM_RF_DCHG`, `MOM_RF_RAN`, `MOM_SW_NTOT`.

### Maturity gate consumption (`solve_for_pops.f` ~1028-1100)

`GS_MAX_ALL = max(KN_RF_RESID or BEST_CHG, SW_MAXCHG, RF_DCHG)`
compared against `GS_UNSAFE_TOL` (0.5) and `GS_MATURE_TOL` (0.05);
bulk tiers use `NBAD50/NTOT` and churn shares. `GS_RAW_IMMATURE`
drives `GS_FREEZE` (KNOP only), which forces `LAMBDA_ITERATION`
and optionally holds populations (`GS_HOLD_POPS`).

## 4. Answers to the specific verification questions

- Compared-against-itself: **no**. Every metric uses a distinct
  snapshot (`SNP`, `J_NEXT_FRQ` pre-write, `VOLD`, `JNEXT_SAVE`).
- Premature overwrite of the comparison base: **no**. `SNP_*` is
  written on `NEW_FREQ`/slot-shift only; `JNEXT_SAVE` precedes all
  refresh writes; `JNEXT_RF_READ` is the Jacobi snapshot.
- Drift on a field not updated this iteration: **partially by
  design** - MOM compressing rows are not written by the sweep
  (refresh-owned), so their sweep drift measures disagreement
  between the fresh solve and a refresh-owned store; it is a
  residual of the *composite* map, not of a write. Correct for a
  convergence test, but the semantics differ from the FG rows.
- Field advanced every outer iteration: **yes** - sweep relaxed
  writes plus refresh commit/keep; `FGINEXT_GS`/`MOMJNXT_GS`
  checkpoint the stores for restart.
- Relaxed vs unrelaxed distinguished: **yes** in the data
  (`SW_MAXCHG`/`GS_MAXCHG` = unrelaxed residual; `RF_DCHG` =
  committed delta), **no** in the user-facing headline, which
  silently picks whichever the getter last wrote (defect D2).
- Update formula: `x <- x + omega*(G(x)-x)` confirmed at fg ~4281
  (`SNP + INEXT_RELAX*(T1-T2)`) and mom ~3623
  (`J_NEXT + T5*(XM - J_NEXT)`), `T5 = JNEXT_RELAX` or per-cell
  `MOM_DAMP_NM` during refresh.

## 5. Defects

### D1 - `FG_RF_DCHG` stale on non-refresh iterations (FIXED)

`FG_RF_DCHG` is initialized once (`DATA`, fg:1007) and reset only
inside the refresh (fg:6199). The per-iteration `INIT` block
(fg:2256-2309) clears `FG_SW_*`, `FG_GS_*`, `FG_RF_CAPPED` but not
`FG_RF_DCHG`; the MOM analogue `MOM_RF_DCHG` is cleared (mom:1672).
On an iteration where the refresh is skipped (hot-pass budget
exhausted, incomplete store), `GET_FG_GS_MATURITY` reports the
previous refresh's retained delta. Severity: low-medium - it feeds
`GS_MAX_ALL` and the release audit, so a stale *large* value can
hold a freeze one iteration too long (conservative direction), and
a stale *small* value could overstate stability on the RF channel.
Fix: one-line reset added to the `INIT` block.

### D2 - headline drift metric changes definition mid-iteration

`FG_GS_MAXCHG` = sweep drift before the refresh and `KN_RF_RESID`
after it (fg ~6353 overwrites); `MOM_GS_MAXCHG` = best-pass
residual after an unconverged refresh. On refresh-skipped
iterations the sweep value stands. Any single-number monitor sees
three different semantics stitched together and can read a
definition switch as a jump or stall. Severity: medium for
interpretation, none for correctness (the gate deliberately wants
the max of the three - it is *reporting* that is ambiguous).
Fix: not code - the `GS_DRIFT_DIAG` line prints all three
quantities explicitly (section 10). Optional cosmetic fix P1:
relabel the print.

### D3 - normalization asymmetry FG vs MOM

FG floors the denominator at `1e-4` of the **fresh** field peak;
MOM floors at `1e-4` of the **stored** column peak. When the field
is moving fast the two metrics weight low-intensity cells
differently. Both are defensible; they are not interchangeable and
cross-path comparisons of the drift value are approximate.
Severity: low. Not changed - changing either denominator alters
gate telemetry.

### D4 - MOM sweep drift includes refresh-owned rows

Compressing rows are excluded from the sweep store write but
included in `MOM_GS_MAXCHG`. On those rows the metric is
|fresh - refresh-owned store|, a residual of the composite map.
This inflates the sweep max relative to a pure write-consistency
metric - by design (the gate wants the composite residual), but it
guarantees a nonzero floor while the outer iterate disagrees with
the refresh. Severity: informational; documented.

## 6. Update-sequence findings (Phase 2 answer)

The update sequence is implemented as designed:

- The store write is `x <- x + omega*(G(x)-x)` with `G(x)` the
  fresh solve; the drift metric reads `|G(x)-x|`, the **unrelaxed
  residual** - so damping cannot cosmetically shrink the reported
  drift. This is the correct convention (the task's requirement
  that raw and relaxed remain distinguishable is met in the data;
  `RF_DCHG` is the relaxed counterpart).
- The refresh commit is unrelaxed (FG `WRITE_STORE`/`WRITE_STORE_CAP`,
  MOM best-pass composite) - correct: the committed field is the
  converged (or best) iterate, not the damped step.
- `FREEZE_I_NEXT` prevents eval-side writes; rejected JFNK trials
  cannot modify production stores (section 8).
- Checkpoint/restart: `FGINEXT_GS`/`MOMJNXT_GS` dumps carry the
  stores plus the drift history; `RESTORE_FG_INEXT_LATEST`
  validates values; the first post-restart comparison is against
  the restored store, not a zeroed buffer (the zeroed case is
  separately flagged via the cold-start `SOURCE_PRIME=SOURCE`).

No ordering defect found.

## 7. Relaxation findings (Phase 4 answer)

Two distinct relaxation systems exist and must not be conflated:

1. **Sweep write factor** `INEXT_RELAX` / `JNEXT_RELAX`: adapted
   once per iteration. V2 (`FG_ADAPT_OMEGA_V2`, `MOM_ADAPT_OMEGA_V2`,
   default T) estimates rho from the median of the last ~4
   same-map refresh pass ratios and sets `omega -> 1/(1+rho)`,
   smoothed 50/50, clamped [0.05,1]. V1 falls back to the
   cross-iteration residual ratio
   (`MOM_GS_MAXCHG / MOM_RF_RESID_PRV`, FG: `KN_RF_RESID_PRV`).
2. **Refresh inner step** `KN_RAY_OMEGA` (FG march) and per-cell
   `MOM_DAMP_NM` (MOM, seeded from `JNEXT_RELAX`, ratchets down on
   sign flips, recovers 1.15x/call, floored at `MOM_DAMP_MIN`).

Findings:

- **Category mismatch (contributing factor, not a bug):** the V2
  estimator measures the contraction ratio of the *inner*
  frozen-source map but applies the resulting omega to the *outer*
  population-coupled write. When the inner residual is pinned at a
  floor (the D-answer cells), rho -> 1 -> omega -> ~0.5 held
  indefinitely: the outer iterate is damped by an inner-map
  statistic that is itself stagnant. This produces *slow*
  convergence with a constant-looking drift and is a plausible
  mechanism for the reported impression. It cannot fabricate a
  plateau: the residual measured is unrelaxed.
- `MOM_DAMP_NM` can ratchet down across iterations toward
  `MOM_DAMP_MIN` on oscillating cells; recovery is slow by design
  (1.15x/call). Consistent with a damping-limited tail.
- Residual histories are same-map within a refresh (correct) and
  cross-call only in the V1 fallback (documented weakness: it
  conflates source/population drift with map contraction).
- No double-relaxation found: the drift metric is pre-write, the
  commit is unrelaxed, the JFNK path saves/restores the omega
  state (`MOM_GS_MAXCHG`, `MOM_RF_RESID_PRV`, `JNEXT_RELAX`,
  `INEXT_RELAX`, `KN_RF_RESID_PRV` - cmfgen_sub.f 5854-5858).

## 8. JFNK and outer-controller findings (Phases 5-6)

- JFNK admission requires `USE_LAM_ES=T` (`solve_for_pops`
  ~2150-2157, ~2497-2506): under LAM_ES=T the MOM solver does not
  run, so `J_NEXT_FRQ`/`H_NEXT_FRQ` are never touched by evals.
  The MOM-path exclusion is deliberate ("no eval freeze").
- `FREEZE_I_NEXT` is set on JFNK entry (jfnk_solve_corr_v1.f:531)
  and cleared on exit (:2392); the FG store is read-only during
  evals. Scalar adaptive state is saved/restored in
  `JFNK_RAD_SWEEP_V1`. Rejected trials cannot write the store.
- Coupled mode: the eval re-runs the sweep with the lagged store
  **frozen** - the residual `F_eval` is a lagged-J Newton operator
  (`F(0)` legitimately differs from the outer STEQ; documented at
  jfnk_solve_corr_v1.f:644-652). Newton residual convergence is
  therefore convergence of a lagged-coupled system, not the fully
  implicit coupled fixed point; a persistent radiation-field
  inconsistency at the fixed lagged field is *not* visible to the
  Newton residual. This is a designed approximation, flagged so it
  is not mistaken for coupled convergence.
- Outer controller: the freeze/release cycle
  (immature -> frozen -> release -> re-immature) is possible and
  the counters (`N_FROZEN_TOT`, `GS_REL_NOPROG_LIM`,
  `GS_FREEZE_LIM`, `GS_PRECOND_LIM`) exist to terminate it. The
  cycle mechanism the task describes (population correction moves
  the field -> field re-immatures -> refreeze) is real and is the
  intended motivation for the release damping (`GS_REL_N`,
  `GS_LAM_CAP`) and release audit (`GS_JUST_REL`). No counter bug
  found in the audited paths.

## 9. Proposed corrections

P1 (cosmetic): relabel the headline "largest drift" print or split
it into `sweep drift` / `refresh residual` / `retained delta` so a
single channel cannot switch semantics. Not required for
correctness; the diagnostic line already separates them.

P2 (deferred refactor): extract the relative-drift kernel
(`|new-old|/max(|old|, eps*peak)` with census outputs) into a pure
shared subroutine so the ten synthetic cases in the task can be
unit-tested numerically against the production formula. Deferred
to keep this audit's diff minimal.

P3 (research direction, not a change): the adaptive-omega category
mismatch (section 7) suggests estimating rho for the *outer* write
from the *outer* residual series (e.g. `SW_MAXCHG` across
iterations, same-map) rather than from inner pass ratios. Any such
change needs the diagnostic channel added here plus model A/B
runs - explicitly out of scope for a measurement audit.

## 10. Diagnostic added (Phase 7)

New VADAT key `T [GS_DRIFT_DIAG]` (default F). When on,
`SOLVE_FOR_POPS` prints one line per great iteration after the
maturity getters with both fields' raw vs relaxed vs retained
quantities:

```
GS_DRIFT it= <n>
  ray:  sw=<FG_SW_MAXCHG>  res=<KN_RF_RESID>  rf=<FG_RF_DCHG>
        om=<INEXT_RELAX>  prv=<KN_RF_RESID_PRV>  cap=<...>
        arg=(K,ray,freq)  nconv/ntot=<...>  ndiv=<...>
  mom:  sw=<MOM_SW_MAXCHG> res=<BEST/MAXCHG>   rf=<MOM_RF_DCHG>
        om=<JNEXT_RELAX>  prv=<MOM_RF_RESID_PRV>  conv=<MOM_RF_CONV>
        ran=<MOM_RF_RAN>  arg=(depth,freq)   nbad50/ntot=<...>
  state: freeze=<GS_FREEZE> justrel=<GS_JUST_REL>
         mature=<N_MATURE>/<GS_REL_STABLE> maxch=<MAXCH>
```

Channels correspond to the task's requested distinction:
raw sweep drift (`sw`), refresh residual (`res`), retained/
relaxed delta (`rf`), relaxation factor (`om`), residual-history
state (`prv`), argmax cell, censuses, and controller state. Two
small accessors (`GET_FG_RELAX`, `GET_MOM_RELAX`) export the
relaxation factors and argmax indices; no internals are exposed
beyond what the maturity getters already return.

## 11. Tests added, and remaining limitations

Added `unit_tests/python/test_gs_drift_diag.py` (static tier,
matching the `test_knop_lagged_field.py` convention): asserts the
flag is declared default-F and parsed, the accessors exist and
export the right module state, the diagnostic block separates
raw/relaxed/retained and prints the argmax cell and controller
state, and `FG_RF_DCHG` is reset at `INIT`.

The ten synthetic numeric cases require driving either the
embedded loops (not callable standalone) or an extracted pure
kernel (P2). Their analytic oracles:

| # | Case | J_new vs J_old | Expected metric |
|---|------|----------------|-----------------|
| 1 | exact | equal | D=0, converged |
| 2 | geometric | gap halves | D_n+1/D_n = 0.5 |
| 3 | over-damped | 1% steps | residual stays ~gap; D does NOT ->0 while x crawls |
| 4 | oscillatory | A,B,A,B | D alternates, no maturity |
| 5 | near-zero | zeros/1e-30 | finite, floor-governed, no false conv |
| 6 | shock-zone | bulk quiet, few cells bad | max catches the bad cells |
| 7 | relaxation | omega=0.3 | residual |G-x| unchanged by omega; committed move = omega*|G-x| |
| 8 | stale buffer | x_old reused | comparison history must detect non-advance |
| 9 | restart | restored store | first post-restore comparison valid, not artificial 0 |
|10 | non-finite | NaN/Inf | invalid state; `MOM_RF_DCHG=1e30` path; no convergence |

Limitations:

- Model-level validation (the five velocity profiles in the task)
  requires interactive-partition sbatch runs per AGENTS.md - not
  executed here; the open item remains in
  `docs/FIELD_DRIFT_TODO.md`.
- The inner-map residual floor is quantified only through the
  cited code comments and the recorded m15 run (~1.0-1.1 stall);
  a dedicated rho-estimate run with `GS_DRIFT_DIAG=T` would pin it
  per cell.
- The audit covers the KNOP/DO_SOB_FG_KNOP ray path and the
  LAM_ES=F moment path. The stock monotonic path (V13/V11 solvers)
  has no lagged store and is out of scope by construction.

## 12. Source map (exact locations)

| Quantity | File | Lines |
|----------|------|-------|
| sweep dispatch | comp_j_blank.f | 834-862 (FG), 1054-1059 (MOM) |
| refresh calls | cmfgen_sub.f | 3281-3293 |
| FG in-sweep drift+write | fg_j_cmf_nonmono_v1.f | 4168-4330 |
| FG per-iteration reset | fg_j_cmf_nonmono_v1.f | 2256-2309 |
| FG refresh (KN_RF_RESID, FG_SW_*, FG_RF_DCHG) | fg_j_cmf_nonmono_v1.f | ~6199-7070 |
| GET_FG_GS_MATURITY | fg_j_cmf_nonmono_v1.f | 5020-5073 |
| MOM per-iteration reset | mom_j_cmf_nonmono_v1.f | 1640-1679 |
| MOM in-sweep drift+write | mom_j_cmf_nonmono_v1.f | 3405-3828 |
| MOM refresh driver | mom_j_cmf_nonmono_v1.f | 4804-6610 |
| MOM_RF_DCHG | mom_j_cmf_nonmono_v1.f | ~6414-6447 |
| adaptive JNEXT_RELAX | mom_j_cmf_nonmono_v1.f | ~6500-6593 |
| GET_MOM_GS_MATURITY | mom_j_cmf_nonmono_v1.f | 4427-4525 |
| maturity getters in gate | solve_for_pops.f | ~1028-1042 |
| JFNK admission / defer | solve_for_pops.f | ~2150-2157, 2497-2542 |
| SOLVEBA / JFNK dispatch | solve_for_pops.f | 2601-2644 |
| FREEZE_I_NEXT | jfnk_solve_corr_v1.f | 531, 2392 |
| JFNK eval sweep state saves | cmfgen_sub.f | 5766-5860 |
