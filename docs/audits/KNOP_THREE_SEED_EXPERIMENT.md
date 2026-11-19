# KNOP three-seed probe experiment

Status: RUNNING (jobs queued/submitted 19-Nov-2026; results sections pending)

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

## 6-10. Results, comparison, interpretation

Pending - jobs in queue as of 19-Nov-2026.

When telemetry arrives, results are reported separately as:

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
