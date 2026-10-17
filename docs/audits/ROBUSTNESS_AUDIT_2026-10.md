# Robustness Audit - October 2026

Follow-up implementation of the 36-item full-repo audit (W1-W25 plus
the W26-W36 bookkeeping items folded into the phase groups).  Landed
as three local commits on `development` (no push, per request):

| Commit | Phase | Scope |
|---|---|---|
| `5a772803` | A | Critical/high engine guards, W1-W7 |
| `a97d3bb1` | B | Medium guards and metric hygiene, W8-W20 |
| (Phase C) | C | JFNK literature upgrades, W21-W25 |

Verification for every phase: each touched file compile-checked under
the production ifort flags (`-O3 -xHost -qopenmp -align array64byte
-extend-source`), then a full `./build.sh` completed clean with fresh
`exe/cmfgen_dev.exe`, `exe/cmf_flux.exe`, `exe/obs_frame.exe`.  The
NaN-blind guard conversions in W6 use the fail-closed idiom
`.NOT.(ABS(X) .LE. LIMIT)` so a non-finite X trips the guard instead
of comparing false to every threshold.

## Phase A - critical/high

| W | Status | Notes |
|---|---|---|
| W1 | FIXED | `BAACC_RES2_MIN` update at solveba_v13.f moved inside the `BAACC_PENDING` gate so demoted/LAMBDA/held calls cannot ratchet the SE2 floor. The `BAACC_TRM` trailing-minima feed cited in the audit does not exist in this tree (stale line numbers); the live defect was the ungated update itself. |
| W2 | FIXED (telemetry) | Whole-LU restore contract documented; `GS_POPS_HELD` deferral at :518 already prevents restore-over-hold, so the suspected defect is largely refuted - added `BA_HOLD_RESTORE_REVIVAL` (INTEGER*8, CONTROL_VARIABLE_MOD) bumped and printed if a restore ever runs with a hold armed. |
| W3 | FIXED | Sentinel eval paths set `JFNK_EVAL_BROKEN` (existing module flag); `JFNK_REBUILD_BASELINE_V1` skips RSCALE/PTC adoption on a broken eval and the corr caller exits via `JFNK_CATASTROPHIC`, symmetric with the F0 gate. |
| W4 | FIXED | `NBB_LS`/`NTB_LS` added to the sweep dealloc list; every refresh/sweep alloc failure (main block, FLD, FLD_P, VMAXF) now writes a loud LUER line instead of a silent `RESID=1.0` or fake `IXIT=5`. |
| W5 | FIXED | `get_ibound` validates the grid at load (finite/positive/monotone), bounds both walks, guards a NaN NU; `get_jh_at_prev_time_step` gains IOSTAT on all three REC reads plus a bounded tail walk and NaN-NU guard; `add_dust_eta` ray walk bounded. |
| W6 | FIXED | Three zones converted to `.NOT.(ABS(X) .LE. LIMIT)`: LUM_DIV chain (solve_for_pops.f), `MAXCH > MAX_CHNG_LIM` (cmfgen_sub.f), OBS restore `NU_F` grid tests (fg). |
| W7 | FIXED | `LST_ITERATION`-gated rate-table blocks and `WRITE_RATES` now require `.NOT. JFNK_IN_EVAL`; `WRITE_RATES` remains reachable inside evals at `CONV_DIAG_LVL>=3`. |

## Phase B - medium

| W | Status | Notes |
|---|---|---|
| W8 | FIXED | `NEV_CAP = MIN(JFNK_EVAL_LIM, JFNK_MAX_EVALS)` computed once; every probe JVP checks `>=2` eval headroom and skipped probes are named at `CONV_DIAG_LVL>=2`. |
| W9 | FIXED | `JFNK_ENRICH_ACTIVE` disarmed before the Z0 plain telemetry apply, re-armed for Z1; the two rows now genuinely compare plain vs enriched. |
| W10 | FIXED | `SCR_RITE_V2`/`SCR_RITE_NAM_V2` close the unit on all error returns, check POINT1/POINT2 writes, set `WRITFAIL` consistently; call sites use a dedicated `SCRWF` flag (no more NEWMOD aliasing) and report torn checkpoints. |
| W11 | FIXED | Dead `FILENAME.EQ.FILENAME` predicate removed; EDDFACTOR/rd_cont_j REC reads IOSTAT-guarded; rd_cont_j refuses interpolation on invalid grids; `read_bcd_mat` reports open/read/write errors and rejects unknown types. |
| W12 | FIXED | `MIN(P,0.9999*R_GRID(I))` clamp with `N_CLAMP_POINT` on the ray-geometry SQRTs (characteristics_v2 + _gam, cmf_formal_rel_v4 + _gam); `N_STORE=MAX(1,NINT(...))` floors. |
| W13 | FIXED | (a) `SGN_PRV`/`NOSC_FL` reset at each `SWEEP_FLD` entry - the audit's `SGN_ACC` claim was refuted (it resets every pass) but the host-variable sign state genuinely leaked across the I_M/I_P boundary. (b) New exit class `IXIT=6` for plateau-forced exits (`FG_PLATEAU_DET`), exported as `FG_GS_NSTAL`/`FG_GS_NPLAT`; `solve_for_pops` applies a separate `max(50,NTOT/100)` plateau tier. |
| W14 | FIXED | `KN_RES_PASS` sized `2*(KN_RAY_MAXPASS+KN_RAY_HOTPASSES)`; hot-window passes no longer feed the omega estimator; estimator prints rho/window/applied at `CONV_DIAG_LVL>=2`. |
| W15 | FIXED | (a) `BA_DAMP_CAP` in (0,1] now fails at read time. (b) REVIEWED - `MOM_TRIAL`'s double `(.NOT. GS_UNSAFE)` is the intended parity predicate; no typo found, comment retained. (c) Churn branch can no longer re-cap `T1` above `GS_LAM_CAP` while `LUM_BLOCKED` is armed. |
| W16 | FIXED (partial port) | MOM `KEEP_BEST` commit gains the per-cell retained-store jump guard (>50x cell AND above column scale). MOM already carried the `1E3*MIN_RES` whole-field bound; P2 damping / hot extension not ported (FG-specific sweep machinery, modest MOM compute). No record-layout change. |
| W17 | FIXED | `cm_blkband_v3` divides now floor `POPS` first; `N_DIV_FLOOR` counted and FLAG set on any applied floor. |
| W18 | FIXED | `JNU_CLAMP_DIAG.txt` gated on `CONV_DIAG_LVL>=1` with IOSTAT; `IP_FG_DATA` open IOS checked (`WRITE_IP` disarmed on failure); `KN_HOT_SPENT` cap tested inside the CRITICAL via `L_HOT_OVER`. |
| W19 | FIXED (bounded) | `FG_INEXT_CAPPED` widened to INTEGER*8; the serialized `FGINEXT_GS` field stays 32-bit via saturating `INT(MIN(...,2^31-1),4)` to preserve the record layout. |
| W20 | REFUTED | The claimed FRZ-reflex interpolation-stamp defect assumes a `J1W`/`J2W` stamp mechanism; no such stamps exist anywhere in this tree (grep-verified across `jfnk_cmfgen_int_mod.f` and the corr driver). Commit/baseline-change paths already rebuild the basis unconditionally. Nothing to fix. |

## Phase C - literature upgrades

| W | Status | Notes |
|---|---|---|
| W21 | FIXED (toggle) | Eisenstat-Walker forcing was already fully implemented (Choice-2 `eta_k = gamma*(|F_k|/|F_{k-1}|)^alpha`, `gamma*eta_{k-1}^alpha` safeguard, `[eta_min,eta_max]` clamps, preconditioned-norm reference). Added VADAT `[JFNK_FADAPT]` (default T) switching adaptive vs fixed gtol. Refs: Eisenstat & Walker 1996, Knoll & Keyes 2004. |
| W22 | STUB+DESIGN | `[JFNK_FDCPLX]` (default F) registered; `JFNK_JVP_V1` fails closed with an explicit NOT-IMPLEMENTED STOP if enabled. Port plan in `docs/audits/JFNK_COMPLEX_STEP.md`. Refs: Martins & Ning 2021; Arramy et al. 2024 Sec 2.3. |
| W23 | PARTIAL | `[JFNK_ESFR]` (default F) registered; the source-function-response term (Arramy eqs. 55-61) needs a stored dJ/dn channel from the previous committed Newton update, which does not exist in this tree. Enabling the flag prints a loud WARN and runs the plain anchor basis. Basis recompute already keys off `JFNK_ENRICH_KLOC` state on each assembly. |
| W24 | FIXED | (a) `JFNK_RETREAT_RESET_V1` now clears `JFNK_EVAL_LIM`/`JFNK_SEC_HIT`/`JFNK_N_RESIDUAL_EVALS` and re-stamps `JFNK_T_START` (`JFNK_MAX_EVALS` is a user VADAT re-read on NEWMOD - deliberately not reset). (b) Enrichment workspace split into four guarded allocates with full dealloc on failure. (c) Per-depth LU/LSQ factorization gains the synthetic-RHS backward-error micro-test (`b=M*1`, DGETRS/DPOTRS, reject rel>1E-4) matching the banded path. (d) Reach guard prints a one-line caveat when `JFNK_BANDED_LU_READY` since `JZ_VEC` prices the diagonal band only. |
| W25 | FIXED | Per-activation allocations hoisted to `SAVE` buffers with lazy size-checked init: GMRES basis `v(n,m+1)`, `s_abs`/`idx_sorted`, and the GMRES-dimension-experiment save/restore packs. |

## Known limitations / deferred

- W22 is a flag+stub by design; the complex-step port is scoped in
  `JFNK_COMPLEX_STEP.md`, not implemented this phase.
- W23's stored-dJ/dn channel is the only missing piece; everything
  else (flag, warning, recompute-on-basis-change) is in place.
- W16 MOM parity ported the commit-jump guard only; P2 damping and
  the hot-run extension stay FG-side (they exist to tame ray-sweep
  stall modes that the MOM refresh does not exhibit).
- `FGINEXT_GS`/`MOMJNXT_GS`/`RESACC` record layouts are unchanged.
