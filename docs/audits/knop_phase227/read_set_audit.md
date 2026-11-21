# Phase 2.27C input read-set audit

Status: SOURCE-REVIEWED with a documented residual gap. No runtime
file-access trace exists; project rules prohibit executing the model
diagnostics on the login node and no capture submission is authorized.
Classification is a static closure, not an observed read set.

Audit inputs: 233 `ln -sf` link names in `batch_v2cap.sh`, all `OPEN`
statements under new_main/, plane/, subs/, newsubs/, unix/, com/, tools/,
`GET_ENVIRONMENT_VARIABLE` call sites, and the captures_v4/v5 manifests.

## Result: manifest coverage

- All 233 launcher-created link names are present in the 251-file
  manifest (the one apparent miss was a grep artifact in a comment).
  These resolve to atomic/scientific data under `$ATOMIC` plus
  `ca2_col_data_filtered` under the project tree. They are the names
  the executable opens in its working directory; the code never opens
  the `$ATOMIC` paths directly in this model.
- The remaining 19 manifest files are the core model inputs:
  VADAT, VADAT_INIT, IN_ITS, RVTJ, MODEL, MODEL_SPEC, SN_HYDRO_DATA,
  T_IN, GAMMAS_IN, RDINR, NUC_DECAY_DATA, SPECIES_MASSES,
  EDDFACTOR_INFO, FIN_CAL_GRID, CUR_MODEL_DATA, CFDAT_IN,
  ADJUST_R_DEFAULTS, arnaud_rothenflug.dat, ca2_col_data_filtered.
- Every literal OPEN name in the source inventory resolves to a
  manifest name, a conditional restart file, or an output file.

## Conditional restart/checkpoint reads (STATUS='OLD' + IOSTAT)

POINT1 (NG accel), BAMATPNT (BA pointer), plus AUTO_CHK_*/GREY_CHK/
SN_GREY_CHK/MU_VALUE_CHK/DDT_WORK_CHK/COEF_CHK families. Absent files
produce a guarded skip, consistent with a fresh start. The capture
assumes fresh-start: all 14 restart names are linked into EMPTY private
scratch. This matches the original launch (links were created dangling),
but whether the original SCRATCH contents now present were run outputs
- not inputs - still needs timestamp confirmation against submission.
If the intended capture is a warm restart, the candidate is WRONG.

## Dynamic and environment-dependent paths

- 36 OPEN sites use constructed names (FILNAME/FILENAME/...); the
  species subset actually opened is fixed by MODEL_SPEC but not
  provable statically. Species-dependent missing-link risk remains.
- `SCRATCH_DIR` drives reads in read_ba_data_v3.f, read_bcd_mat.f and
  scratch writes in cmf_blkband_v3.f/store_ba_data_v3.f/cmfgen_sub.f.
  Manifest pins it to the capture-private scratch: covered.
- `CMF_KN_RF_MB` sets the refresh memory budget in
  fg_j_cmf_nonmono_v1.f. It was ABSENT from the v4 manifest (gap);
  v5 pins its value or its absence via environment/environment_absent.
- `SLURM_JOB_ID` (RUN_ID in solve_for_pops.f) is required by the
  runner; PWD/USER in wr_pwd.f feed a log line only.

## Runtime-generated outputs that alter the read set

- All `./data/...` OPEN sites found are writes (GAMMA_RAY_LINES,
  DIAGN_EDEP, gamma_nu_grid.dat, ETA_ISO_*, GAMMA_J_*, ETA_MUAVG_*).
  The model dir has no `data/` directory; if the gamma-transport path
  runs, the first write fails unless work/data exists. v5 creates it.
- POINT1/SCRTEMP/STEQ_VALS/J_COMP/JEW/EWDATA/CFDAT_OUT/CONT_FREQ/
  LINEHEAT/NETRATE/TOTRATE/TRANS_INFO/ADIABAT_CHK are written into
  SCRATCH during the run - outputs, except where a later iteration
  re-reads them. They stay inside private scratch either way.
- TIMING, OUTGEN, MOM_J_ERRORS, ITER_LEDGER and the KN* diagnostic
  files are work-dir outputs.

## Machine-readable classification

read_set_audit.json lists every manifest name with class, role and
coverage; every conditional, dynamic and env-dependent site; and the
explicit gaps. Gaps that keep the candidate BLOCKED:

1. Observed read set unproven - runtime trace requires allocation.
2. Fresh-vs-warm restart assumption unconfirmed.
3. Species-dependent dynamic-name coverage unproven.
