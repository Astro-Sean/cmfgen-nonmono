# KNOP three-seed probe experiment

Status: LIVE -- arm C lost to designed non-convergence abort
(20-Nov-2026 03:42); arms A, B, A/A still running. Interim analysis
appended in sections 6-10; the experiment is NOT complete.

## 1. Objective and scope

Determine whether the probe ray's P-direction solve converges to the
same retained solution when its initial field at each J_NEXT refresh is
changed. The experiment exercises the existing `KN_SEED_ARM` machinery
on the designed host (`m15_2p82d_KNOP_p8_probe` recipe).

### 1.1 What this experiment tests

- It varies the initial P-direction field of ONE probe ray
  (`KN_NOCP_PROBE_LS=83`).
- The probe ray never commits its candidate to the retained `I_NEXT`
  store (source evidence in section 2.1): normal writes, capped writes,
  heal writes, and the seed-commit marker are all suppressed.
- It therefore isolates the sensitivity of the probe-ray P solve to
  its initial field, per refresh, under a live non-LTE model.

### 1.2 What this experiment does not test

- It does NOT directly test sensitivity of the entire non-LTE model to
  different population or radiation-field initialisations: all four
  arms boot from the identical iso_M it157 checkpoint, and only the
  probe ray's per-refresh seed differs.
- It does NOT establish global KNOP convergence robustness. The store,
  populations, and spectrum can only differ between arms through
  indirect channels (see section 1.4), not through a seeded store.
- The KNOP solver as a whole must not be classified as robust or
  initialisation-sensitive on the basis of this experiment alone.

### 1.3 The A/A arm

Arm A/A (`KN_SEED_ARM=3`) uses the same retained-store seed as arm A
(`KN_SEED_ARM=0`) - it is a repeatability control for the same seed
and the same probe calculation, not a different physical
initialisation. Any difference between A and A/A in initial fields,
per-pass histories, or final probe fields measures nondeterminism of
the run (e.g. OpenMP scheduling order on 8 threads, timing-dependent
refresh boundaries), not physical sensitivity to the seed. Differences
between A/A and A are reported quantitatively and are not
automatically interpreted as physics.

### 1.4 Propagation channels to watch

Although the probe cannot write the store, its measured quantities are
not fully inert:

- `FG_RF_DCHG` is updated with the probe's candidate delta at
  `fg_j_cmf_nonmono_v1.f:10223` and `:10249` (inside the capped/normal
  commit branches, NOT gated on `L_NOPROBE`). `FG_RF_DCHG` feeds
  `RF_DCHG=MAX(FG_RF_DCHG,FG_RF_PEND)` at `:6666`, which enters the
  `GS_MAX_ALL` maturity reading. A probe candidate that persistently
  differs from the retained generation therefore keeps RF_DCHG
  elevated - conservative, and identical across arms only to the
  extent their probe candidates agree.
- The veto counters `FG_NC_CAP`/`FG_NC_REJ` (`:10220`, `:10230`) count
  the probe's would-be verdict; the ray is then re-marked verdict 7.
- Whether these channels can actually shift the convergence path or
  exit reason between arms is a measurement question for the
  comparison phase.

## 2. What the seed actually modifies (Phase-1 finding)

`KN_SEED_ARM` (VADAT key, default -1 = off) selects the initial field
`FLD` handed to the probe ray's P-direction solve in
`KNOP_RAY_SOLVE_V1` at each refresh (`KN_SEED_SETUP_P`,
`fg_j_cmf_nonmono_v1.f` ~line 14304):

- arm 0 (A, reference): retained `I_NEXT` store candidate, the standard
  path.
- arm 1 (B): sibling-median seed via `KN_SEED_SIB` - cells filled from a
  median over neighbouring rays where mappable (`KN_SEED_NMAP/NFILL`
  count filled vs mappable; `NUNM` unmapped, `NRET` retained).
- arm 2 (C): chained seed - the previous refresh's saved final P iterate
  (`KN_SEED_PC`, archived by `KN_SEED_SAVE_P`). First refresh falls back
  (`KN_SEED_USED=4` records the fallback).
- arm 3 (A/A): reproducibility control - same retained-store seed as
  arm 0, records arm bookkeeping only.

The seed modifies the initial formal-solve ray field on probe ray
`KN_NOCP_PROBE_LS=83` only. Population initialisation, the radiation
field on all other rays, and all model inputs are untouched.

### 2.1 Probe-ray no-commit evidence (source)

`KN_NOCP_PROBE_RAY(LS)` (`fg_j_cmf_nonmono_v1.f:2337-2341`, a pure
function: `KN_NOCP_PROBE .AND. LS .EQ. KN_NOCP_PROBE_LS`) sets the
host-local `L_NOPROBE` at `:9796`. Every store-write pathway is
gated on it:

- Capped commit: `WRITE_STORE_CAP` for M and P at `:10214-10218`,
  both under `IF(.NOT. L_NOPROBE)`.
- Normal commit: `WRITE_STORE` for M and P at `:10237-10239`, both
  under `IF(.NOT. L_NOPROBE)`.
- Seed-commit marker: `FG_RAY_SEEDED(LS)=1` at `:10246-10247` -
  the probe measures its first-commit candidate but never marks
  itself seeded, so the jump-veto waiver state is also unchanged.
- Divergence heal: `HEAL_STORE` for M and P at `:10410-10413`,
  under `IF(.NOT. L_NOPROBE)`.
- Verdict bookkeeping: the probe's reached verdict is copied to
  `NOP_VC` and `FG_RAY_VC(LS)` overwritten with probe code 7 at
  `:10264-10268`.

Runtime verification is built in: `KN-NOP store` lines at
`:8580-8590` print `SUM` of the probe ray's M and P store slices
before and after each solve; identical sums (delta 0) prove no
commit, capped commit or heal wrote those slices.

As noted in section 1.4, `FG_RF_DCHG` / `FG_NC_CAP` / `FG_NC_REJ`
are deliberately NOT gated: the probe's measured delta still
counts toward the refresh-change diagnostic and its would-be
verdict still counts in the veto tallies.

Restart/state files carrying information into each run (all staged
identically from the iso_M it157 donor): `POINT1`, `POINT2`, `SCRTEMP`,
`FGINEXT_GS` (I_NEXT store), `MODEL`, `CUR_MODEL_DATA`, `SN_HYDRO_DATA`.

## 3. Controlled initial states (Phase 2)

| Arm | Model dir | KN_SEED_ARM | Seed |
|-----|-----------|-------------|------|
| A   | `m15_2p82d_KNOP_seed_a`  | 0 | retained store (reference) |
| B   | `m15_2p82d_KNOP_seed_b`  | 1 | sibling-median |
| C   | `m15_2p82d_KNOP_seed_c`  | 2 | chained last-candidate |
| A/A | `m15_2p82d_KNOP_seed_aa` | 3 | retained store (control) |

Pairwise `config_resolved.yaml` diffs: the only non-path differences are
`description` and `KN_SEED_ARM`. Velocity law, density, abundances,
atomic models, frequency grid, boundary conditions, transfer settings,
and all convergence tolerances identical. No converged output from one
arm feeds another; each is an independent job chain from the same donor
checkpoint. The seed is deterministic machinery (no RNG); perturbations
are bounded by construction - arm B interpolates existing sibling ray
values and arm C reuses a previously computed iterate, so no negative
populations or unphysical fields are introduced (populations are never
seeded at all).

## 4. Run metadata and checksums (Phase 3 record)

Executable: pinned frozen build
`exe_frozen/844b17b0fb430c6e_7c87795b85115b49/exe/cmfgen_dev.exe`,
md5 `2e2c1522f5ec` (16,661,392 bytes, built 19-Nov-2026 12:01 from
tree at commit `b2663a9f`, dirty x5: uncommitted `DO_J_RES_DIAG`
SIGMA_TRUE allocation guards in `obs/COMP_JCONT_V4.INC` and
`obs/cmf_flux_sub_v5.f`, plus untracked debris `TIMING`, `fort.4`,
`unit_tests/systems/_pv.txt`). The uncommitted changes are diagnostic
write-path guards only; all four arms run the identical binary, so
arm-vs-arm comparisons are internally consistent. Intel ifort
production flags (`-O3 -qopenmp -mkl`, AVX-512).

Initial-state checksums (identical across all four arms, verified
19-Nov-2026):

| File | md5 |
|------|-----|
| POINT1 / POINT2 | `272f74b8211724f59add39cd891e9a42` |
| SCRTEMP | `1b27c3e975bf2cc104ade855f7a4a308` |
| FGINEXT_GS (I_NEXT store, 5.4G) | `893143f3c47405f23d4c7aca57d64080` |
| MODEL | `9f0c4d8bb272b2595f2c5c1d082b764b` |
| CUR_MODEL_DATA | `da389791b426cef7e52e079e64f4bfec` |
| SN_HYDRO_DATA | `ea6c5d5bd2d1c0ae785d70d4c9a56053` |

Job submissions (`run_model_grid.py`, small partition, 8 cpus, 64G,
23:59 limit, requeue-enabled).

### 4.1 First submission (FAILED - telemetry defect)

| Arm | FINAL job | OBS job | Outcome |
|-----|-----------|---------|---------|
| A   | 13334092 | 13334093 | FAILED 6m18s; OBS cancelled |
| B   | 13334165 | 13334166 | FAILED 5m51s; OBS cancelled |
| C   | 13334174 | 13334175 | FAILED 5m51s; OBS cancelled |
| A/A | 13334185 | 13334186 | FAILED 5m49s; OBS cancelled |

All four died identically: `forrtl: severe (61):
format/variable-type mismatch` in `fg_knop_refresh_v` at
`fg_j_cmf_nonmono_v1.f:8603`, on the second `KN-SEED` summary WRITE
at the first refresh. The format
`'(A,2I8,A,I8,A,I3,A,I4,2(1X,ES11.3),1X,F7.4)'` mismatched the
14-item list: `KN_SEED_PCNP` (integer) hit an `A` descriptor, and the
trailing `' schg='`/SCHG pair had no descriptors at all. These jobs
are not experiment runs; they are retained in this record as the
failed attempt.

Salvage: each crashed run had already written its first-refresh
binary record. Parsed values (reader: sequential unformatted records;
14-int header, 4 reals, XREF then X0 as NI x NCF Fortran-order):

| Arm | used | seed action at r0001 | X0 vs XREF |
|-----|------|----------------------|------------|
| A   | 0 | retained seed | identical |
| B   | 1 | nmap=nfill=169128, nunm=3777192 | 169128 elements differ |
| C   | 2 -> 4 | first-refresh fallback (no saved PC) | identical |
| A/A | 3 | control | identical |

So the seed machinery itself worked before the crash: arm B's X0
genuinely diverged from the reference in exactly the fill count, and
arm C took its documented fallback. The KN-NOP store checksums
printed `delta M/P=0` on every arm before the fatal write, giving a
first live confirmation of the no-commit suppression.

### 4.2 Correction and re-submission

Fix (diagnostic WRITE only, no solver arithmetic):
`fg_j_cmf_nonmono_v1.f:8603` format replaced by
`'(A,2I8,A,I8,A,2I4,A,3(1X,ES11.3),A,1X,F7.4)'`, matching all 14
items in order. Commit `e1df4324`.

Isolated runtime check: the exact statement was compiled and run
under ifort with representative values; output
`unmapped/retained=      17       3 dead=    4210 pcref/pcnp= 158  30 dab/rmax/pcrs=   2.841E-01   9.137E-01   6.290E-01 schg=  0.4375`
shows all 14 items in order with correct types and readable field
separation.

New executable: rebuilt at 13:39 19-Nov-2026 from tree commit
`e1df4324` plus uncommitted obs-path work (`obs/COMP_JCONT_V4.INC`
md5 10dfa3176527c118, `obs/cmf_flux_sub_v5.f` 94b86ad611dda6e3,
`obs/jres_validity.f` 407929f99c824e94 - JRESC event-writer guards,
diagnostic path only, owned by a parallel session and deliberately
not committed here). `exe/cmfgen_dev.exe` md5
`3f37f22bd133468fb2854758d19abf12`, 16,661,392 bytes. Build: ifort
production flags, clean except pre-existing `#6075`/`#8291` remarks;
post-build unit suite all pass at time of resubmission.

| Arm | FINAL job | OBS job |
|-----|-----------|---------|
| A   | 13339955 | 13339956 |
| B   | 13339970 | 13339971 |
| C   | 13339980 | 13339981 |
| A/A | 13339988 | 13339989 |

OBS jobs are dependency-held on their FINAL. Per-run `VADAT` is
generated at job start by `prepare_VADAT.py` from
`config_resolved.yaml`; `KN_SEED_ARM` sits under `vadat_overrides`
(verified in all four resolved configs).

### 4.3 First-refresh validation (all arms PASS)

Each arm ran past the previously fatal write. `KN-NOP store` printed
delta M/P=0 on every arm; all three `KN-SEED` lines printed; the
`KNSEED_P_L0083_a<arm>_r0001.bin` record (2,706,216 bytes) exists and
parses under the documented layout:

| Arm | used | dead | pcref/pcnp | dab | rmax | pcrs | schg |
|-----|------|------|-----------|-----|------|------|------|
| A   | 0 | 3 | 1/30 | 0 | 0 | 0.629 | 0 |
| B   | 1 | 3 | 1/18 | 8.258E4 | 2.157E-2 | 2.280 | 1.0 |
| C   | 4 (fallback) | 3 | 1/30 | 0 | 0 | 0.629 | 0 |
| A/A | 3 | 3 | 1/30 | 0 | 0 | 0.629 | 0 |

All four jobs continued running well past the write (a/b/c >11 min,
A/A >7 min; the old crash fired at ~6 min). The experiment is now
LIVE.

Early readout (one refresh only, not a result): arm B's
sibling-median seed is a large perturbation by design
(dab=8.26E4, 169128 cells filled, schg=1.0) and its first probe
P-solve ended at residual 2.280 vs 0.629 on the retained-seed arms.
Whether that gap closes or persists across refreshes is the
experiment question.

## 5. Telemetry captured per run (Phase 4)

- OUTGEN `KN-SEED` summaries: `used`, `sib fill/map`, `pcix`,
  `unmapped/retained`, `dead`, `pcref/pcnp`, `dab/rmax/pcrs`, `schg`.
- Per-refresh binary records `KNSEED_P_L0083_a<arm>_r<refresh>.bin`:
  14-integer header, 4 reals, then `KN_DK_XREF` (pre-seed reference)
  and `KN_DK_X0` (post-seed start) NI x NCF fields.
- Standard iteration histories (population MAX_CHNG, GS_MAXCHG,
  refresh residuals) from OUTGEN.

## 6. Interim scheduler and job status (20-Nov-2026 ~09:15 CET)

| Arm | FINAL job | State | OBS job | State | Latest refresh / iter |
|-----|-----------|-------|---------|-------|-----------------------|
| A   | 13339955 | RUNNING 18:52 elapsed (~5:07 left) | 13339956 | PENDING (afterok) | r399 / it557 |
| B   | 13339970 | RUNNING 18:52 elapsed | 13339971 | PENDING (afterok) | r435 / it593 |
| C   | 13339980 | FAILED exit 6:0 at 03:42:41 after 13:21:22 | 13339981 | CANCELLED (dependency) | r292 / it449 |
| A/A | 13339988 | RUNNING 18:48 elapsed | 13339989 | PENDING (afterok) | r458 / it615 |

All eight jobs ran on vipc2057 (8 requested, 16 allocated, OMP=16).
FINAL jobs are requeue-enabled on wall limit only
(`FINAL_REQUEUE_WHEN=wall`, `MAX=500`): the three survivors are
expected to hit the 23:59 wall and requeue, which is normal
continuation, not failure. Arm C's exit was a CMFGEN `STOP`, not a
wall kill, so the wrapper correctly refused to requeue.

Failure classification for arm C: NUMERICAL, designed abort.
`SOLVE_FOR_POPS` stopped the run at it449 with

    SOLVE_FOR_POPS: release_no_progress: demand= 8.211E+01
      rel_dem= 8.211E+01  raw_drift= 6.364E-02  count=    10
    SOLVE_FOR_POPS:   10 releases with no demand progress reach
      GS_REL_NOPROG_LIM
      the release/refreeze cycle never matured the lagged field;
      stopping as non-convergent

the `GS_REL_NOPROG_LIM=10` (KNOP auto-set) field-maturity guard -
the same stop class as `FRESH_knop` job 13005611 (it124). No
scheduler, filesystem, or diagnostic-write fault is involved.
Artifacts preserved in place:
`cmfgen_validation/m15_2p82d_KNOP_seed_c/cmfgen/2.82days/`
(`OUTGEN` 10.7 MB, `WARNINGS`, `cmf_13339980.{out,err}`,
`FGINEXT_GS`, `CUR_MODEL_DATA`, `SCRTEMP`, all 292 `KNSEED_*` bins,
manifest). The OBS job never started; `obs/` holds only
`batobs.sh`/`CMF_FLUX_PARAM_INIT`. A partially-written `OBSFLUX`
(1.17 MB) exists but is a live progress file, not a final spectrum.

One refresh per great iteration fires on every arm
(`refresh n = iteration 157+n` while running), so `refresh` is a
clean comparison axis across arms.

## 7. Telemetry findings (refreshes parsed from OUTGEN)

Extraction: `scripts/parse_knop_seed.py` (new; parse patterns
adapted from `parse_knop_iso.py`). Table:
`docs/audits/knop_seed_telemetry.csv` (1575 rows).

### 7.1 Probe-ray P-solve outcome per arm

| Arm | used | resP (final probe P resid) | npP | probe exit (ixit) |
|-----|------|----------------------------|-----|-------------------|
| A   | 0 | pinned 0.69614 essentially every refresh | 26-28 | 3 (stall), always |
| B   | 1 | pinned 2.2799 every refresh | 18 always | 3 (stall), always |
| C   | 4->2 | oscillates 0.2615-2.8738, median 0.614 | 10-61, mean 19.5 | mostly 3; 39 plateau exits (ixit=6) |
| A/A | 3 | pinned 0.69614 | 26-28 | 3 (stall), always |

The probe ray is the per-refresh residual argmax on every arm at
essentially every refresh (`q argmax ray=83`), so `rf_resid`,
`fg_rf`, and `GS_MAXCHG` are the probe's own residual. All other
rays' per-ray residual distributions are identical across arms to
printed precision (`q dist p50=5.716-5.732E-4, p95=6.828E-3` at
matched refreshes).

Pass-1 residuals at refresh ~100 show what the seed is worth at
solve start: A retained seed 7.79E+01, B sibling median 5.24E+03
(~67x worse), C chained candidate 2.01E+00 (~40x better). Within
2-3 passes all three collapse to the same O(1) residual regime:
the perturbation mostly decays during the solve, but the stall
exit lands on a different endpoint per arm (A 0.696, B 2.28, C
0.26-2.87 oscillating). Across refreshes there is no systematic
amplification: B's dab is pinned at 8.256-8.258E4 and resP at
2.2799 for all 435 records; C's dab holds ~0.177-0.18 and resP
oscillates without drift.

Arm B bookkeeping: `sib_fill = sib_map = 169128` (all live cells
filled), `nunm = 3,777,192`, `nret = 0`, `ndead = 3`, `schg = 1.0`
every refresh - the sibling-median map is complete and constant.
Arm C chains correctly: r1 took the documented fallback (`used=4`),
r2+ all `used=2` with `pcref/pcnp/pcrs` describing the candidate
saved at the same refresh (i.e. what the next refresh chains).
The chain has limit-cycle epochs: resP=2.8738, npP=10, ixit=6 for
16 consecutive refreshes (r127-142), 15 (r165-179), 8 (r264-271) -
the chained seed recycles a bad candidate until the solve escapes.
`pcnp` range 10-61.

M direction is untouched: npM=2, resM ~3E-9 on all arms.

### 7.2 Arm C abort mechanism (measured, not hypothesized)

`GS_MAX_ALL = MAX(GS_MAXCHG, GS_SW_MAXCHG, GS_RF_DCHG)`; the gate
reads immature when `GS_MAX_ALL > GS_UNSAFE_TOL = 0.5` (plus the
bad-cell share tier). `GS_MAXCHG` is the refresh residual max =
the probe's residual here. Release events require the gate to
sample a mature reading.

- A/A/A: probe resP pinned at 0.69614 > 0.5 on all 396+/455
  refreshes -> gate never matures -> permanently
  `gs_freeze_forced+hold_pops` -> zero releases -> counter never
  arms. B identical situation at resP=2.2799.
- C: probe resP < 0.5 in 53 of 292 refreshes (chained seed
  occasionally lands a better-than-0.5 candidate). Every
  `release_no_progress` event sits in a sub-0.5 resP window
  (verified at refs 48, 154, 184, 204, 209, 241, 246, 253, 260,
  292). 232 frozen iterations vs 59 release/frozen-alternating
  iterations; 10 judged no-progress releases -> STOP at it449.

So the abort is arm-C-specific but NOT a solver defect or a seed
poisoning: the model on all four arms sits in the same
freeze-hold, non-contracting state (population demand ~82, Model
Convergence pinned at 50% on every logged iteration, `applied=0`).
C died because its chained seed made the probe residual oscillate
across the 0.5 threshold, which let the release machinery prove
10 times that the field was not mature. The other arms carry the
same unmatured field silently under permanent freeze.

Underlying non-maturity is the known cold-KNOP regime
(`FRESH_knop` it124, same stop); what the seed changed is only
the gate-sampling trajectory. Attribution caveat: field-level
nondeterminism exists (section 8), but the release trigger is
correlated with the chained-seed oscillation, not with noise -
retained-seed resP is pinned to 5 digits over 396+ refreshes.

### 7.3 Global model state comparison

- Population demand identical across arms to printed precision:
  it158-160 demand 5.998E+01/6.062E+01/6.129E+01 on every arm;
  late-run plateau 8.211E+01 (C) vs 8.212E+01 (A). `applied=0`
  (held) throughout; `conv=50%` pinned.
- Commit verdicts identical: `v_committed`, `v_capped`,
  `v_unsolved`, `solves_po` all equal at matched refreshes.
- Diagnostics that INCLUDE the probe differ, as designed:
  `unconv_cells` (probe cells counted), `retained_max`
  (`FG_RF_DCHG` gets the probe's candidate delta, :10223/:10249
  not `L_NOPROBE`-gated), `omega` (`rho_max` in the adaptive
  INEXT_RELAX estimate is the probe residual - A 0.9895,
  B 1.09-1.19, C 0.99-1.19 -> applied omega differs ~6E-4 mean,
  ~1.3E-3 max).

## 8. A/A repeatability control (first decision gate)

- Telemetry (all 396 common refreshes): every extracted field
  identical to printed precision - resP, resM, npP, npM, rf_resid,
  q_arg_resid, worst_resid, retained_max, omega, unconv_cells,
  exit census, verdicts. Zero differences including integer pass
  counts.
- Binary field level (XREF = retained seed, NI=3 x NCF=56377):
  NOT bitwise identical. On ~88.5K physical cells (|x| <= 1E15,
  excluding ~80.5K uninitialized-scale cells at identical
  positions and ~70 NaN cells at identical positions): rel diff
  p50 = 3.6E-8 (r1) growing to ~6E-5 (r50), then BOUNDED at
  ~5.5E-5 median / p99 ~2E-3 / max ~7E-3 through r396. Same level
  at j=1,2,3. Within each arm X0 == XREF bitwise (verified).
- Interpretation: OpenMP-order floating-point nondeterminism in
  the committed store, present from refresh 1 (identical donor
  checkpoint, divergence occurs inside it158 sweeps), saturating
  at ~5E-5 median. This sets the experiment's noise floor: seed
  effects (O(1) relative, schg=1.0) are ~4 orders of magnitude
  above it.
- Verdict: A/A establishes adequate repeatability to interpret
  B and C. Macroscopic solve behavior (residuals, pass counts,
  verdicts) is deterministic at diagnostic precision; raw field
  values carry bounded ~1E-4-relative noise.

## 9. Binary record integrity audit

- All 1,575 `KNSEED_P_L0083_a<arm>_r<NNNN>.bin` records verified:
  size 2,706,216 B each, sequential-unformatted layout
  (14-int header, 4 real*8, XREF then X0 as NI=3 x NCF=56377
  Fortran-order real*8), record markers consistent, header
  arm/used/refresh fields matching filename and expected codes.
  No truncated or misaligned records.
- Content checks: `ndead=3` in every record (NI*NCF-ndead=169128
  live = arm B's nfill). NaN cells ~66-95 per record at positions
  consistent between arms. ~80.5K cells hold uninitialized-scale
  values (|x| > 1E15) - the probe FLD workspace outside the
  stored/transported cells; positions identical between A and A/A
  (same pass counts -> same workspace reuse) but differ for B/C
  (different pass trajectories). They are excluded from field
  statistics as INVALID, not silently treated as zeros.
- `pcrs` records the residual of the candidate saved at the same
  refresh (`KN_SEED_SAVE_P`), i.e. what the NEXT refresh chains -
  not the injected seed's quality.

## 10. Seed response and propagation (interim)

Measured facts:

- B perturbation: ~46.3K physical cells changed at rel p50 ~1.0;
  ~42K cells toggled between physical and uninitialized-scale
  values in both directions; ~85 cells seeded to NaN.
  Perturbation decays ~3 orders during the first 2-3 passes
  (5242 -> O(1)) but the solve ends on a fixed worse stall point
  (2.2799) with FEWER passes (18 vs A's 26-28). No growth across
  435 refreshes; not converging better on any observed refresh.
- C perturbation: chained candidate is the best start by ~40x
  (pass-1 residual ~2 vs 78), oscillates mid-solve, lands at a
  variable endpoint (0.26-2.87). It is occasionally BETTER than
  the retained-seed floor (min 0.2615) - which is exactly what
  opened the release gate and killed the run. Chained refreshes
  are dependent trials (limit-cycle epochs r127-142, r165-179,
  r264-271); they must not be counted as independent samples.
- Whole-model propagation: no committed-field or population
  propagation detected (probe store delta M/P = 0 on every
  refresh of every arm; verdicts and demand identical). Two
  diagnostic channels do respond, as documented in 1.4:
  probe-inclusive diagnostics (rf_resid, unconv census,
  retained_max, omega) differ, and the adaptive-omega channel
  measurably perturbs the retained field update globally
  (omega differs ~6E-4 mean; retained-seed XREF for B/C vs A
  differs ~1.3E-3 median - ~20x the A/A noise floor, ~1000x
  below the seed perturbation itself). Populations remain
  identical to printed precision.
- First observable point of global divergence: the probe
  residual itself (deterministic per seed), then omega/retained
  drift at ~1E-3 relative. No propagation into populations,
  commit verdicts, or other rays' residual distributions.
- Spectrum: NOT ready - OBS jobs pending on FINAL completion
  (A/B/AA) or cancelled (C). No comparison performed.

Interpretation vs measurement (kept separate): the measured
records show the probe P-solve converges to a seed-dependent
stall endpoint. Whether that endpoint spread means anything for
a COMMITTING ray (the production question) cannot be answered by
this experiment - the probe never commits by design, so no arm
tests what a seeded store would do downstream. Arm C's abort
demonstrates a real indirect sensitivity: seed-dependent probe
residuals modulate the freeze/release gate trajectory
(`GS_MAX_ALL` samples the probe's max residual). That is a
maturity-guard design observation, not evidence the seed
strategies change the physical solution. No arm shows evidence
of a better CONVERGED candidate: all probe solves end in
stall/plateau exits above KN_RAY_TOL.

## 6-10 original placeholder text (superseded by the above)

Results are reported separately as:

A. Probe P-solve convergence sensitivity - per-refresh seed used
   (`KN_SEED_USED`), per-pass NOP histories, `dab/rmax/pcrs/schg`,
   final probe fields and `KNSEED_P_L0083_*.bin` XREF/X0 records,
   per arm.
B. A/A repeatability control - quantitative A-vs-A/A differences in
   initial fields, per-pass histories and final probe fields,
   attributed to run nondeterminism unless evidence says otherwise.
C. Rest-of-model differences - population, J/store and spectral
   comparisons between arms, on identical grids.
D. Propagation plausibility - whether any observed arm difference can
   plausibly have propagated through a shared store or population
   update given the no-commit design (expected channels are limited
   to the diagnostic accumulators in section 1.4 plus convergence
   timing).

Classification: the experiment can support "probe P-solve robust /
sensitive / non-convergent" per arm. It does not license a
whole-solver robustness claim; a wider claim would need genuinely
reinitialised population/field states, which this design does not
produce.
