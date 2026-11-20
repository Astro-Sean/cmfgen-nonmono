# Phase 2.21: live KNOP validation, observer spectrum and NI-guard deployment

Status as of this writing: the live-executable rebuild, the invalid-ray
check, and the controlled pre/post experiment construction are complete.
All live-model and observer jobs are submitted and in flight; their
quantitative comparison is outstanding and is marked as such below.
Nothing in this audit is inferred from pending jobs.

Baseline for this audit: tip `595f9ae7` (commits `3295014c` NI guard,
`e66492ec` Phase 2.20 audit).

## PART A - live executable rebuild

`make cmfgen KEEP_EXE=1` was run in the production tree. The make rule
publishes the executable via `.tmp` + `mv`, so no running job could
observe a truncated binary; no running job was modified.

- `exe/cmfgen_dev.exe`: `58a699d0f37864e6741c5f089bbb12d8d6ee838b628ab993b6e636c0ba61732d`
- `lib/libsubs.a`: `d4e8738db9ecfff8060c335f4b7cc0ffadb240670c1bb55bd28555b85f81fc69`, 630 members
- `lib/libmod.a`: `7c30cc524cf5ee033c28c81bf253ae3bc22e95b65a18cd75b1dec9a8f58b581e`, 28 members
- `fg_j_cmf_nonmono_v1.o` in the archive: `b2d5ded4307ae81c47c83125813a6ae34afc10b4f8a8970db52aca19031a4dcd`

Provenance chain (`provenance_timestamps.txt`): source edit 14:23:22,
guarded object re-archived 14:25, archive rewritten 15:00:37, executable
relinked 15:00:39. The NI bound `NI .GT. SIZE(R_RAY,1)` is in the
compiled object; `KN_CS_REA(LS)=8` is set on that early-return path
(`fg_j_cmf_nonmono_v1.f` ~10240). `ldd` confirms MKL
(`libmkl_intel_lp64/intel_thread/core`) and `libiomp5` production
libraries (`hashes_exe_libs.txt`).

Invalid-ray test (`inval`, driver `unit_tests/bin/knray_p220`, linked
from the same archive member set): `cert census: inv=1`,
`rf_all_conv= F`, `RF_ALLCONV= F`, `RESULT_OK`, no segfault, no
out-of-bounds geometry access (`inval_full_output.txt`). The census
does not print the reason code per ray; reason 8 is assigned in source
on the guarded path, and `csut` INVR rows 7/8/9 were green at this HEAD
in Phase 2.20. KN_CS_INV + reason 8 + FG_RF_ALL_CONV=F: PASS.

`exe/cmf_flux.exe` was rebuilt at 12:28 (`abe47258...`, contains
`rf_not_allconv`) and is what the pending dependent OBS jobs will run.
It predates the NI guard (guard is in the KNOP solver, not the observer
route); the OBS pinned trees record the exact hash.

## PART B - A1/A2 live-model runs

Metadata recorded at snapshot (`slurm_snapshot.txt`):

- A1 `13414391`, A2 `13414398`: RUNNING, 4h13m elapsed at snapshot,
  16 CPUs (scheduler-expanded from 8 requested), 64G, `vipc2181`,
  started 2026-11-20 10:52:11, 23:59h allocation with requeue.
- Executable: resolved configs record `cmfgen_dev.exe`
  `93364cace2fabdd94577f8ae983abbddd9affdb2a453e1c688a2a4f9bdc9d7de`,
  git head `217049ee`.

Important correction to earlier planning assumptions: the A1/A2 launch
executable `93364cac` is the PRE-FLAG build. `git show
217049ee:plane/fg_j_cmf_nonmono_v1.f` contains zero occurrences of
`FG_RF_ALL_CONV`; `strings` on the binary finds only pre-existing
`GS_STATE` output and none of `fg_rf_all_conv`, `rf_all_conv`, or
`rf_not_allconv` (`exe_string_checks.txt`). The A1 config description
also scopes "commits 8b421089..217049ee". A1/A2 are therefore pre-flag
replicate arms by design - identical to each other, not pre-vs-post.

No `cert census`/`rf_all_conv` lines exist in the pre-flag output (the
census did not exist in that build).

Configuration correction recorded 20-Nov: every arm carries
`KN_MAT_FREEZE=T` in its live VADAT (line 341; inherited from the
`m15_2p82d_KNOP_iso_M` donor config - the "failure-isolation arm M:
frozen matter" description). At `217049ee:solve_for_pops.f` ~2216 the
hold condition `... .OR. KN_MAT_FREEZE` snapshots and rolls back the
population update on every call regardless of field state. The observed
`held=T` on all ~120 population calls, the constant `applied correction
0.000E+00` under an attempted `8.1E+01` maximum change, and the exact
`Model Convergence: 50 %` plateau are therefore BY DESIGN, not a stall:
populations cannot move under this arm, and the "ten releases" counter
cannot accrue (`n_raw_releases=1`, `n_rel_noprog=1` in GS_STATE).

What the arm actually measures is the lagged-field trajectory at fixed
populations. Over 122 refreshes in A2 (`it`~286 at extraction):

- refresh residual `residual max`: median 6.4-6.5E-2, range
  [3.4E-2, 7.5E-2], no downward trend across refresh blocks 1-122;
  1 of 122 refreshes below the 5% maturity bound;
- `retained-field delta` ~4.0-4.5E-2 (GS_RF_DCHG - under the 5% bound);
- sweep drift `sw=` ~5.0-6.3E-2;
- ray exits steady at tol~65 / bulk~16 / stall=3 (probe included;
  NX_CNT counts the probe, so production stall = 2), unsolved=1 (the
  NI=1 ray, residual 0);
- commit verdicts: 121 of 122 refreshes `committed=81 capped=2
  unsolved=1`, one `committed=80 capped=3`;
- the KN-CEN top-offender set locked in at refresh 2 and never changed:
  production rays 22, 45, 50 on the inward (ch1) arm, probe ray 83 on
  ch2; flagged cells show the pinned-floor signature (`v1~1E-30` raw vs
  `vms~1E-28` stored, ratio ~1E7).

So the lagged store is at a near-fixed point in drift while a fixed
small set of rays fails to converge its own solve each refresh - the
committed field contains capped/provisional content every cycle. The
two capped-commit rays are production rays within the {22,45,50} set;
OUTGEN does not name them per ray (per-ray verdicts live in
`FG_RAY_VC`/`NXT_LS`/`RES_LS` in memory only).

Diagnostic arm submitted to name them: `m15_2p82d_KNOP_diagcap`
`13429272` - clone of the B1_PREFLAG deck (same `93364cac` pre-flag
tree, same iso_M donor, same VADAT including `KN_MAT_FREEZE=T` and the
ray-83 probe) plus `KN_PSWEEP_DIAG=T`, `KN_PSWEEP_LS=-1` (per-ray,
per-pass full-map defect `R_MAP` capture for every ray,
diagnostic-only). Its `KN_PSWEEP_DIAG.dat` records which production
rays plateau at ~6E-2 vs contract to tolerance.

Counter note: `GS_STATE` shows `n_mature=123` consecutive under
`KN_MAT_FREEZE` while the printed sweep/refresh metrics exceed 5%. The
release counters under the unconditional MAT_FREEZE hold do not read
as the normal gate transcript; treat them as diagnostic bookkeeping
only.

Completion analysis (great iterations, releases, freeze counts,
GS_MAX_ALL, GS_RF_DCHG, MAXCH, no-progress counters, termination)
requires the finished OUTGENs and is outstanding.

## PART C - valid before/after comparison

A genuine pre-flag baseline exists: `exe_frozen/93364cace2fabdd9_484f34b19ee9e205/`,
frozen at git `217049ee`, exe `93364cac` (no `FG_RF_ALL_CONV` in source
or binary). Two controlled arms were constructed and submitted
(`ab_design.txt`):

- B1_PREFLAG `13427111` (OBS `13427112`): pinned cmfdist
  `exe_frozen/93364cac...`, exe `93364cac` - a third pre-flag
  replicate, redundant with A1/A2 but harmless.
- B2_POSTFLAG `13427115` (OBS `13427116`): pinned cmfdist
  `exe_frozen/58a699d0f37864e6_abe4725807474143`, exe `58a699d0`.

The resolved configs are identical apart from the description and the
pinned-executable block. `VADAT_INIT` is byte-identical
(`f238af97...`) across A1/B1/B2; `SN_HYDRO_DATA` is byte-identical
(`347eda19...`, the iso_M donor copy); the warm-restart checkpoint is
the same `m15_2p82d_KNOP_iso_M` it157 donor (POINT1/POINT2, SCRTEMP,
FGINEXT_GS, MODEL, CUR_MODEL_DATA). Both arms: 8 CPUs requested, 64G,
same partition/qos, same thread expansion logic.

The numerical delta between the pre-flag exe and the post-flag exe is
the full certification bundle, `217049ee..595f9ae7` on the solver and
observer sources (`cert_bundle_commits.txt`): Phase 2.14 bulk-exit
semantics, six-class ray certification, `FG_RF_ALL_CONV` and the
`LAM_FIELD_NOT_READY` fold, the INIT stale-flag reset, the extended
`GET_FG_GS_MATURITY` signature, INVR reason codes, and the NI guard.
The compiled `KN_RAY_BULK_TOL` default is `6.0E-2` under KNOP in both
eras, so the comparison isolates certification behaviour, not a
tolerance change. It is correct to describe the experiment as
"pre-certification-era" vs "current certified", and not as flag-only.

Limitation recorded 20-Nov: both arms inherit `KN_MAT_FREEZE=T` (see
Part B), so B1/B2 cannot answer the ten-release population-release
question - under frozen matter the release counter cannot accrue at
all. What they DO provide is a valid pre/post comparison of the
lagged-field trajectory and the ray-exit/commit-verdict census at
fixed populations: whether the certification bundle (all-ray flag,
census classes, INIT reset) changes which refreshes are treated as
valid, and whether the post-flag build marks the capped-commit
refreshes `rf_all_conv=F` that the pre-flag build silently accepts.
That is the scientifically relevant subset of the original design.

For a ten-release population comparison the arms need
`KN_MAT_FREEZE=F`; that is a design decision for the user, not a
silent reconfiguration (the option is recorded in the README artifacts
and B1/B2 remain submittable in a free-matter variant under new
names).

## PART D - observer execution

Status of the pending OBS jobs (all request 200G, well above the 30G
interactive cap that killed job `13422341` at 97 GB RSS in 22 s):

- `13414392` (A1 OBS), `13414399` (A2 OBS): Dependency on the running
  FINALs.
- `13427112` (B1 OBS), `13427116` (B2 OBS): Dependency on the queued
  FINALs.

Standalone run: the user's complete `obs_p219/run_diag_on` staging
(symlinked to the `m15_2p82d_KNOP_seed_b` model: RVTJ, FGINEXT_GS,
POINT1, SCRTEMP, MODEL, VADAT) was cloned to `obs_p221` to avoid
touching the user's active directory. seed_b was produced by exe
`8a5737ba` (pre-flag era) and ended "Iteration budget exhausted before
convergence" - it is a real but non-converged KNOP model, sufficient
for the observer-route and diagnostic-independence checks; it does not
substitute for a converged model in Part F.

Jobs `13427285` (diag on) and `13427286` (diag off) run
`exe/cmf_flux.exe` `abe47258` via `obs_template/batobs.sh` at
general/200G/8 CPU (`obs_p221_setup.txt`, `obs_p221_job.sh`).
`abe47258` contains `rf_not_allconv` and the extended
`GET_FG_GS_MATURITY` signature; the corrected interface is exercised in
the actual OBS route via `RESTORE_FG_INEXT_OBS` in `cmf_flux_sub_v5.f`.

Runtime, exit status, memory, and OBSFLUX checksum will be recorded on
completion. Link/startup alone is not counted as validation.

## PART E - observer diagnostic independence

The two p221 arms differ only in `PRINT_KNOP_DEBUG` (T vs F) in
`CMF_FLUX_PARAM_INIT`; `DO_KNOP_SPEC_GATE=F` in both, and the key gates
WRITE statements only (`fg_j_cmf_nonmono_v1.f` 6101-6232,
`SET_FG_KNOP_DEBUG`). Byte-for-byte comparison of OBSFLUX, obs_fin_10,
OUT_FLUX, and the frequency grid is scripted for completion; any
difference will be traced to the exact pathway rather than assumed
inert.

## PART F - population/radiation consistency

Requires a completed post-flag model (A2-class or B2). Outstanding.
When available: latest refresh number, `FG_RF_ALL_CONV`,
`GS_MAX_ALL`, `GS_RF_DCHG`, population `MAXCH`, release count,
termination reason. The maturity criterion is a 5 percent bound
(`GS_MATURE_TOL` on `GS_RF_DCHG`) plus the refresh verdict - not exact
equality between successive fields - and will be reported as such.

## PART G - mixed-sign validation

Still outstanding by construction: every available deck is the
all-compressing m15 hydro (v monotonically increasing ~48.9 ->
7073.6 km/s, SIGMA >= -1.037 with no genuine mixed-sign inversion).
Nothing in this phase is used as expanding/compressing coupling
evidence. Requirements unchanged: real non-monotonic V(r), active
seeded cells on both GAM signs, verified sign transitions in solver
inputs, successful refresh, population convergence, observer spectrum.

## PART H - deliverables

- This document.
- `docs/audits/knop_phase221/`:
  `hashes_exe_libs.txt`, `provenance_timestamps.txt`,
  `exe_string_checks.txt`, `inval_full_output.txt`,
  `inval_live_exe_archive.txt`, `slurm_snapshot.txt`,
  `ab_design.txt`, `cert_bundle_commits.txt`,
  `dessart_m15_2p82d_KNOP_B1_PREFLAG.yaml`,
  `dessart_m15_2p82d_KNOP_B2_POSTFLAG.yaml`,
  `obs_p221_setup.txt`, `obs_p221_job.sh`, `README_pending.md`.

## FINAL REPORT (interim - live items pending)

1. Live executable rebuilt with NI guard: YES - `58a699d0`, guarded
   object `b2d5ded4`, provenance chain recorded.
2. Invalid-ray test on that build chain: PASS - `inv=1`,
   `rf_all_conv=F`, reason 8 on the guarded path, no invalid access.
3. A1/A2 complete: PENDING - running 5h+/24h. NOTE: both arms carry
   `KN_MAT_FREEZE=T` (frozen-matter diagnostic); held=T every call is
   by design, populations cannot converge under this arm.
4. Ten-release failure recurred: CANNOT OCCUR under MAT_FREEZE arms -
   the release counter cannot accrue while populations are held. The
   question remains open for a free-matter arm.
5. Valid pre/post comparison: CONSTRUCTED but RE-SCOPED - B1
   `93364cac` vs B2 `58a699d0`, identical inputs, both frozen-matter:
   valid for field-trajectory/certification comparison, NOT for the
   ten-release population question (would need KN_MAT_FREEZE=F arms).
6. cmf_flux completes: PENDING - p221 arms 13427285/6 queued at 200G
   after the diagnosed 30G OOM; dependent OBS jobs queued behind FINALs.
7. Diagnostic on/off OBSFLUX byte-identical: PENDING - arms submitted,
   comparison scripted.
8. Final population evaluated against latest verified field: PENDING -
   needs a completed post-flag model; ordering already proven at source
   (Phase 2.19/2.20).
9. 5% maturity criterion: DOCUMENTED - `GS_MATURE_TOL` bounds
   `GS_RF_DCHG`; it is a bounded-drift criterion, not field equality.
    Live confirmation pending.
10. Mixed-sign velocity fields: UNVALIDATED - no genuine inversion
    model exists; requirements documented in Part G.
