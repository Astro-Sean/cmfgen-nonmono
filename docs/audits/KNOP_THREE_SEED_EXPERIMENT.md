# KNOP three-seed convergence-robustness experiment

Status: RUNNING (jobs queued/submitted 19-Nov-2026; results sections pending)

## 1. Objective and scope

Determine whether the KNOP non-monotonic solver converges to the same
retained solution when the probe-ray seed state at each J_NEXT refresh is
changed. The experiment exercises the existing `KN_SEED_ARM` machinery on
the designed host (`m15_2p82d_KNOP_p8_probe` recipe), not a wholesale
reinitialisation of the model: all four arms boot from the identical
iso_M it157 warm-restart checkpoint and differ only in which field seeds
the probe ray's retained store candidate at each refresh.

Distinguishes:

1. Robust convergence to the same solution.
2. Convergence to different solutions dependent on seeding.
3. False convergence via stopping criteria.
4. Failure to converge / strong sensitivity.
5. Agreement in the emergent spectrum despite internal field differences.

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
`KN_NOCP_PROBE=T` suppresses store writes from the probe ray, so the
probe seed cannot contaminate the retained `I_NEXT` store; arm
differences therefore measure convergence and reproducibility of the
probe solve and of the whole-model trajectory around it.

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
23:59 limit, requeue-enabled):

| Arm | FINAL job | OBS job |
|-----|-----------|---------|
| A   | 13334092 | 13334093 |
| B   | 13334165 | 13334166 |
| C   | 13334174 | 13334175 |
| A/A | 13334185 | 13334186 |

OBS jobs are dependency-held on their FINAL. Per-run `VADAT` is
generated at job start by `prepare_VADAT.py` from
`config_resolved.yaml`; `KN_SEED_ARM` sits under `vadat_overrides`
(verified in all four resolved configs).

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
