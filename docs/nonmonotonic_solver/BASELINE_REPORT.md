# BASELINE_REPORT.md — CMFGEN Clean Baseline

## Source Archive

| Field | Value |
|-------|-------|
| URL | `https://sites.pitt.edu/~hillier/cmfgen_files/cur_cmf_18jun25.tar.gz` |
| Download timestamp | 2026-08-21 09:26:54 UTC |
| SHA-256 | `fc81d2dcf2df6538058ee46eab27ab231e01152e02c12666906d583c4dae4a00` |
| Size | 45,060,661 bytes (42.9 MB) |
| Top-level directory | `cur_cmf/` |

## Build Environment

| Field | Value |
|-------|-------|
| OS | Linux viper03 5.14.0-570.128.1.el9_6.x86_64 (RHEL 9) |
| Architecture | x86_64 |
| Compiler (primary) | Intel ifort 2024.0.2 (`/mpcdf/soft/RHEL_9/packages/x86_64/intel/2024.0.2/bin/ifort`) |
| Compiler (fallback) | gfortran (`/u/sbrennan/.local/bin/gfortran`) |
| MPI | Not used (OpenMP only) |
| OpenMP | Yes, 4 threads (controlled by Slurm) |
| Build flags | `-O2 -heap-arrays -assume byterecl -traceback` (from Makefile_definitions) |

## Working Tree Status

| Field | Value |
|-------|-------|
| Path | `/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source` |
| Git remote | `https://djohnhillier@github.com/pitt-phyast/cmfsource.git` |
| Git branch | `main` (28 commits ahead of origin/main) |
| Modified files (vs git HEAD) | 13 files (unstaged) |
| Clean baseline path | `/viper/ptmp1/sbrennan/SN2024pba/cmfgen_baseline/clean_cmfgen/cur_cmf` (immutable) |

## Modified Files (Working Tree vs Clean Baseline)

27 files differ (excluding .git, .mod, exe, RELEASE_NOTES).

### Source files (20)

| File | Diff lines | Category |
|------|-----------|----------|
| `new_main/subs/rv_sn_model_sniin.f` | ~80 | REQUIRED |
| `new_main/cmfgen_sub.f` | 3 | REQUIRED |
| `new_main/mod_subs/rd_control_variables.f` | 3 | REQUIRED |
| `plane/fg_j_cmf_v13.f` | ~40 | REQUIRED |
| `plane/fg_j_cmf_v12.f` | ~8 | REQUIRED |
| `plane/mom_j_cmf_v11.f` | ~40 | REQUIRED |
| `subs/sobew.f` | 4 | REQUIRED |
| `subs/sobew_grad.f` | 4 | REQUIRED |
| `subs/sobew_grad_v2.f` | 4 | REQUIRED |
| `subs/sobew_grad_v3.f` | 4 | REQUIRED |
| `subs/sobjbar_sim.f` | 4 | REQUIRED |
| `disp/subs/sobjbar.f` | 4 | REQUIRED |
| `new_main/subs/do_ng_accel_band_v4.f` | 7 | REQUIRED |
| `new_main/mod_subs/set_abund_clump.f` | ~19 | RELATED |
| `new_main/subs/rd_sn_data.f` | ~11 | RELATED |
| `new_main/subs/generate_full_matrix_v3.f` | ~10 | RELATED |
| `new_main/subs/solveba_v13.f` | ~22 | UNRELATED (debug) |
| `subs/ltepop_wld_v2.f` | ~12 | UNRELATED (debug) |
| `subs/regrid_log_dc_v1.f` | ~15 | UNRELATED (debug) |
| `misc/rewrite_old_model_data.f` | 2 | DANGEROUS |

### Build files (4)

| File | Category |
|------|----------|
| `Makefile_definitions` | UNRELATED (build config: HOST=viper, COMPILER=intel) |
| `new_main/Makefile` | UNRELATED (build config: -lgfortran, plane/lte_hydro deps) |
| `disp/Makefile` | UNRELATED (build config) |
| `obs/Makefile` | UNRELATED (build config) |

## Clean Baseline Preservation

The clean baseline at `/viper/ptmp1/sbrennan/SN2024pba/cmfgen_baseline/clean_cmfgen/cur_cmf` is immutable and will not be modified. All implementation work will be done on a copy or via patch application to the working tree.
