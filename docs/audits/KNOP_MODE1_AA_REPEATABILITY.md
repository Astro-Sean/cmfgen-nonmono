# KNOP MODE 1 A/A repeatability experiment

20-Nov-2026. A clean A/A pair on the MODE 1 probe-isolation
executable. Two identical jobs from the same donor, same probe
configuration, same immutable binary; the comparison answers how
much of the controller state is deterministic under the isolated
probe.

Status: RUNS SUBMITTED - comparison sections pending.

## 1. Source and executable

- Tree: `cmfgen_jfnk`, branch `main`, HEAD `217049ee`
  (`audit: document empty-pool omega gate and results`).
- MODE 1 implementation commits verified present via
  `git merge-base --is-ancestor`:
  `8b421089 970584dc c85b08df d484c4ae 8eda6c23 282a3ed1
  6f2ec135 217049ee`.
- Working tree: no modified tracked files at submission;
  `git_dirty_files=4` are untracked run artifacts
  (`TIMING`, `fort.4`, `unit_tests/systems/_pv.txt`,
  `unit_tests/systems/p29/`) that do not enter the binary.
- Build: `make DEBUG=nodebug all`, Intel ifort production flags,
  clean (pre-existing #6075 warnings only).

### Executable snapshot (immutable)

Pinned via `scripts/cmf_run_manifest.py pin` (created 20-Nov-2026,
`created: true`):

- Snapshot dir:
  `/viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk/exe_frozen/93364cace2fabdd9_484f34b19ee9e205/`
- `cmfgen_dev.exe` sha256
  `93364cace2fabdd94577f8ae983abbddd9affdb2a453e1c688a2a4f9bdc9d7de`
- `cmf_flux.exe` sha256
  `484f34b19ee9e205f6c9d3a485de707a81cab0c495c251ec1c7e79ba1888e9c2`
- Post-copy verification: pinned-copy sha256 identical to
  `exe/` values for both exes (hashes above are the pinned files).
- Both slurm scripts set `cmfdist=` to this snapshot; the
  generated job re-verifies the hash immediately before exec.

## 2. Run configuration

Cloned from `configs/dessart_m15_2p82d_KNOP_seed_aa.yaml` (the
prior A/A arm). The only yaml differences between the two arms
are `model_name`, `description`, and `output.model_dir`;
`diff` of the resolved configs shows nothing else.

| | A1 | A2 |
|---|---|---|
| Config | `configs/dessart_m15_2p82d_KNOP_A1_MODE1.yaml` | `configs/dessart_m15_2p82d_KNOP_A2_MODE1.yaml` |
| Model dir | `cmfgen_validation/m15_2p82d_KNOP_A1_MODE1` | `cmfgen_validation/m15_2p82d_KNOP_A2_MODE1` |
| FINAL job | 13414391 | 13414398 |
| OBS job | 13414392 | 13414399 |
| Probe | `KN_NOCP_PROBE=T`, `KN_NOCP_PROBE_LS=83` | same |
| Seed arm | `KN_SEED_ARM=3` (retained store) | same |
| Donor | `m15_2p82d_KNOP_iso_M/cmfgen/2.82days` (it157 warm restart) | same |
| CPUs | 8 requested, 64G, small partition, requeue on wall | same |
| OBS | general partition, 200G | same |

Initial-state checksums (identical A1 vs A2, verified post-staging;
match the recorded donor values from the prior experiment):

| File | md5 |
|------|-----|
| POINT1 / POINT2 | `272f74b8211724f59add39cd891e9a42` |
| SCRTEMP | `1b27c3e975bf2cc104ade855f7a4a308` |
| MODEL | `9f0c4d8bb272b2595f2c5c1d082b764b` |
| CUR_MODEL_DATA | `da389791b426cef7e52e079e64f4bfec` |
| SN_HYDRO_DATA | `ea6c5d5bd2d1c0ae785d70d4c9a56053` |
| FGINEXT_GS | `893143f3c47405f23d4c7aca57d64080` (donor; both copies 5,367,543,560 bytes) |

VADAT is generated at job start by `prepare_VADAT.py` from the
identical resolved config, so it is identical by construction.

One setup deviation from the seed_* pipeline: the grid script's
inline `SN_HYDRO_DATA` conversion requires a `.swd` that does not
exist for this model; `SN_HYDRO_DATA` was therefore copied from
the donor into `cmfgen/2.82days/` before submission, and the
script skipped conversion (`use_existing_sn_hydro: true`). Same
file, same md5 as every prior arm.

## 3. Telemetry preserved per run

All in `cmfgen/2.82days/` of each model dir:

- `OUTGEN` - per-refresh `KNOP I_NEXT refresh` block:
  `residual max` (production-only `KN_RF_RESID` under MODE 1),
  `omega est:` (`rho_p90`, `rho_med`, `rho_max`, `n_ray`,
  `applied omega`), `adaptive INEXT_RELAX`,
  `diverged rays`, `clamp events`, `same-census cap/floor`,
  `unconverged cells`, `best-pass bad cells`,
  `commit rejects`, `nc rejects`, `nc capped`,
  `retained-field delta: pending` (`FG_RF_PEND`),
  `ray exits` census, `KN-PROD resid:` banner when it fires.
- `KN-NOP ray probe-only` (probe resM/resP, npM/npP),
  `KN-NOP store` (probe store checksum delta),
  `KN-SEED` (arm/used/fill counts),
  `KN_WATCH` (probe residual and pass count).
- `ITER_LEDGER`, `ITER_HISTORY`, `GS_STATE`,
  `CORRECTION_SUM`, `TIMING`.
- `FGINEXT_GS` final retained `I_NEXT` store.
- OBS job: `obs/` output including final spectrum
  (`OBS_FIN`/`SPEC` products).
- Slurm logs `cmf_<jobid>.out` (environment, thread count,
  exe verify).

`KN_RF_RESID_PRV` is not printed per refresh; its correctness is
observable through the `omega est`/`adaptive INEXT_RELAX` chain
and was unit-verified (`prbz`/`prbl`/`prbr` cases in
`test_knop_mode1.py`).

## 4. Expected nondeterminism sources

- OpenMP scheduling order (16 allocated threads on the prior
  arms despite 8 requested; allocation is scheduler-dependent).
  Bounded binary-field noise, same as the prior A/A arm.
- Wall-clock requeue boundaries: `final_requeue_when: wall` can
  split a run mid-iteration; checkpoint restore is deterministic
  but the iteration at which it splits is not.
- `n_ray` in `omega est` is the production pool size under
  MODE 1 (84 production rays of 85 with ray 83 excluded).

## 5. Comparison plan (Phase 4/5 - to fill after runs)

Per refresh, align A1 vs A2 refresh blocks on refresh index
(refresh n = solver iteration - warm-restart base 157):

- `residual max`, `omega est` fields, `INEXT_RELAX`,
  `diverged rays`, all census/counter numerators, `pending`.
- KN-NOP/KN-SEED probe records: resid, store delta, seed counts.
- Freeze/release decisions from `GS_STATE` history.
- Population iteration counts from `ITER_LEDGER`.
- Final `FGINEXT_GS`: cell-wise max-abs / max-rel / RMS on the
  retained `I_NEXT` store (5.4G binary, same parser as
  `parse_knop_seed.py` with record markers consumed).
- Final spectra (OBS jobs): wavelength-grid agreement, max and
  RMS relative flux difference, continuum-normalized difference,
  diagnostic line profiles.

Report max-abs / max-rel / RMS / N for floats, exact agreement
for discrete counters. No silent dropping of unmatched refreshes
or frequency points - alignment method stated with results.

## 6. Results

PENDING - FINAL jobs 13414391 (A1) and 13414398 (A2) queued
20-Nov-2026; OBS jobs 13414392/13414399 dependency-held.

## 7. Pass/fail assessment

PENDING.
