# KNOP Refresh Convergence Audit (P1-T1e, PART 7)

Question: does `FG_KNOP_REFRESH_V1` / `KNOP_RAY_SOLVE_V1` genuinely
converge to a unique seed-independent fixed point, or is the committed
store's seed-independence a stopping-tolerance / iteration-cap /
commit-cap artefact?

Scope: `plane/fg_j_cmf_nonmono_v1.f` (`FG_KNOP_REFRESH_V1` ~6755-7291,
`KNOP_RAY_SOLVE_V1` ~7303-10518 and its `CONTAINS` workers). No
production equation modified. Evidence classes: **[MATH]** code
structure, **[MANUF]** manufactured model, **[PROD]** production model.

**Headline: the refresh is a separate fixed-point operator whose
internal `FLD` iteration reads OPPOSITE frequency directions in one
array (compressing `L+NFM`, expanding `L-1`), so it can host cycles the
in-sweep store map cannot. Seed-independence of the committed store is
real only where the field solve converges cleanly; guarded commit paths
can retain seed memory. A multi-seed / tolerance / cap recheck on a
production checkpoint is specified but NOT yet run.**

---

## 1. Mechanism (traced)

`FG_KNOP_REFRESH_V1` runs once per great iteration after the frequency
sweep (caller `cmfgen_sub.f:3313-3321`, gated
`(DO_NONMONO_KNOP .OR. DO_SOB_FG_KNOP) .AND. DO_KNOP_REFRESH`), and
self-guards on `KN_STORE_READY` and `FREQ_IDX_ITER==KN_NCF-1` (complete
sweep only).

Per ray (`!$OMP PARALLEL DO` over `LS`, rays independent given fixed
opacities):

1. `SEED_FLD` copies the store column into a local field `FLD`
   (compressing slots only; `FLD(:,1)` stays 0).
2. `KNOP_RAY_SOLVE_V1` runs `SWEEP_FLD`/`SWEEP_FLD_ORD`, iterating the
   *field* `FLD` - not the store - reading the iterated `FLD` via
   `SRC_PRIME`/`KN_ORD_SPV`: compressing at `L+NFM` (redder), expanding
   at `L-1` (bluer).
3. On exit, `WRITE_STORE` (unrelaxed, clean convergence) or
   `WRITE_STORE_CAP` (bounded ~`FG_COMMIT_NC`=5% step) commits; diverged
   or jump-rejected rays write nothing; `HEAL_STORE` zeroes poisoned
   slots.

**Structural difference from the in-sweep store map:** inside one `FLD`
array, compressing rows read *forward* in frequency while expanding rows
read *backward*. Spatially alternating expanding/compressing bands can
therefore create directed cycles in the refresh's internal dependency
graph - a structure the acyclic in-sweep store map does not have.
**[MATH]**

## 2. Convergence controls and exit semantics

All VADAT-wired (defaults): `KN_RAY_TOL`=2e-3, `KN_RAY_MAXPASS`=60,
`KN_RAY_BULK_TOL`=6e-2 (KNOP), `KN_RAY_STALL`=0.98,
`INEXT_RELAX` in [1e-3,1], `FG_COMMIT_JUMP`=50, `FG_COMMIT_NC`=0.05,
`KN_RAY_HOTPASSES`=120 / `KN_HOT_TOTAL_PASSES`=480,
`FG_COMMIT_CONV`=T, `FG_SEED_COMMIT`=T, `FG_BULK_CONTR`=T,
`FG_ORD_SWEEP`=T.

Per-pass residual `RES` = max over cells of
`|VNEW-VOLD| / MAX(VOLD, 1e-4*VMAX_L)`; `NBAD`/`NTOT` count significant
cells. Exit codes IXIT: 0=exact, 1=tol, 2=bulk, 3=stall, 4=cap,
5=diverged, 6=plateau. Hot extension continues up to `KN_RAY_HOTPASSES`
on stall/cap.

Commit-side gates: `JUMP = MAX(FG_JUMP_MAX(...))`; `JUMP >
FG_COMMIT_JUMP`(50) rejects unless `FG_SEED_COMMIT` waives a ray's first
commit; `FG_COMMIT_CONV=T` sends non-clean exits to `WRITE_STORE_CAP`
(if `FG_COMMIT_NC>0`) else rejects.

## 3. Seed-independence mechanism and its limits

- **Clean convergence -> seed-independent.** `WRITE_STORE` overwrites
  every compressing slot unconditionally (`STORE=V1`, floored 1e-30,
  capped `NONMONO_INEXT_CAP`), leaving no seed memory. If every ray
  exits IXIT 0/1/2-clean, the committed store depends only on the fixed
  opacities.
- **Guarded paths -> possible seed memory.** Diverged rays write
  nothing (store keeps old values); rejected commits write nothing;
  `WRITE_STORE_CAP` commits only a bounded step; `FG_COMMIT_CONV=F`
  would commit an unconverged last iterate. In these cells the store
  retains dependence on the incoming seed.

So `J_full=0` (measured exactly on the manufactured model, all
perturbations) demonstrates clean convergence for *those* seeds at
`KN_RAY_TOL=2e-3` - genuine seed independence there, **not** a
proof that every seed/tolerance/cap combination converges. **[MANUF]**

## 4. Distinguishing genuine convergence from artefacts

The zero committed-store Jacobian does **not** by itself prove internal
`FLD` convergence: two different seeds could reach the same committed
store through different (one clean, one capped) paths, or both could be
truncated at the same cap. The decisive test compares the **internal
`FLD` arrays and residual histories**, not only the committed store.

Specified recheck (production checkpoint, fixed opacity):

1. Seed `FLD`/store from zero, baseline, scaled baseline, random
   positive, and boundary-focused states.
2. For each seed run `FG_KNOP_REFRESH_V1` (or `KNOP_RAY_SOLVE_V1`
   directly) and record: full per-pass `RES`/`NBAD`/`NTOT` history,
   iteration count, IXIT stopping reason, cap status, final `FLD`,
   committed store, `WRITE_STORE_CAP`/`HEAL_STORE`/floor/cap counters.
3. Compute the fixed-point residual
   `q = ||FLD - REFRESH_OPERATOR(FLD)|| / max(||FLD||, scale_floor)`.
4. Repeat at tighter `KN_RAY_TOL` and larger `KN_RAY_MAXPASS`/`HOT`
   caps; compare complete `FLD` arrays element-wise.

Failure-mode classification to report: non-convergence (IXIT 4/5/6
dominant), slow convergence (residual decreasing but cap-bound),
multiple fixed points (different `FLD` fixed points), tolerance
artefact (agreement only at loose `KN_RAY_TOL`), cap artefact
(`KN_RAY_MAXPASS`/`HOT` binding), commit-cap (`WRITE_STORE_CAP`
dominant), or numerical floor (`NONMONO_INEXT_CAP`/1e-30 binding).

## 5. Status

- Mechanism and exit semantics: **traced [MATH]**.
- Manufactured `J_full=0`: clean convergence demonstrated for tested
  seeds **[MANUF]**.
- Production multi-seed / tolerance / cap recheck: **SPECIFIED, NOT
  RUN** - needs a production checkpoint and the `FLD`-level residual
  instrumentation (the committed-store dump at unit 291 does not expose
  internal `FLD`).

The correct current statement: the refresh commits a seed-independent
store **where its per-ray field solve converges cleanly**; whether that
holds across all production seeds, tolerances and caps is an open
empirical question, not a proven property.
