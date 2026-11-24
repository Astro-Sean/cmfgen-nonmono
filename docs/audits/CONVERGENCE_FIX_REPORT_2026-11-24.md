# CONVERGENCE_FIX_REPORT — 24-Nov-2026

Work order: convergence liveness, NaN safety, certification on the
CMFGEN non-monotonic fork.

## 1. Provenance

- Base SHA: `3fd9152ba6888491136073fd9cdaedddf01a5674` (main HEAD at
  work start).
- Branch: `fix/convergence-liveness-and-guards` in the git worktree
  `scratch/cmfgen_convguard_wt` (the main checkout `cmfgen_jfnk/`
  and its dirty tree were left untouched).
- Commits (one per work item):
  - `977227ca` W1 NG NaN/Inf-safety
  - `d725a18e` W2 Knop retained-store poison heal
  - `7d944068` W3 ES-ALI fixed-point certification
  - `9b28ebb4` W4 JFNK_PC_QUAL default-on
  - `c0720736` W5 ADJUST_CORRECTIONS validation
  - `39c141bd` W6 energy-balance diagnostic
- Toolchain: ifort 2021.11.1 (module `intel/2024.0`, oneAPI
  2024.0.2), `-O3 -mavx512* -qopenmp -mkl`, `viper_intel_nodebug`
  config. gfortran 11.5.0 + MKL 2026.2 for the standalone JFNK core
  check.
- Post-commit test fix: `938e534a` retargets the
  `test_knop_lagged_field` fully-poisoned-row source pin to
  `KN_HEAL_ROW` (the W2 helper now holding the verbatim `NV=0`
  zero-row branch it pins).
- Executable: `exe/cmfgen_dev.exe` linked 24-Nov 13:03, md5
  `da13b2e30a88382eb1e26668d0c0cfc0`.

## 2. Work item summary

| W | Files | Change | Test | Evidence |
|---|-------|--------|------|----------|
| W1 | `do_ng_accel_band_v4.f`, `control_variable_mod.f`, `rd_control_variables.f` | Pre-scan history + post-scan fit for non-finite values; NaN-safe hold and acceptance predicates; `T2` init every depth; `NG_NF_HIST`/`NG_NF_NEW` census + status line; `NG_NF_ABORT` key (default F). | `test_ng_abstain.py` 23/23 (17 baseline + 6 new pins); direct binary run: NaN PREVRJ -> `NG_DONE=F`, `NG_NF_HIST=1`, `POPS` untouched; 1E300-history overflow -> `NG_NF_NEW=1`. `-O2 -Wall -Wuninitialized` clean on the edited block. | VERIFIED |
| W2 | `plane/fg_j_cmf_nonmono_v1.f`, `drv_knop_ray.f90`, `test_knop_lagged_map.py` | Row-scoped sanity helpers shared with `HEAL_STORE`; entry-side heal before `SEED_FLD`; veto-time retained-row heal + one re-evaluation; retained-branch check in `FG_COMMIT_NC=0`; `FG_POISON_HEAL` counter under the `FG_HEAL_UPD` critical section; refresh summary print. | `test_knop_lagged_map.py` 34/34: sdab flips pinned-defect -> pass (VC=0, converged certified solve commits, poison=16 healed); sdac vetoes a genuine 1E6 jump (VC=1, poison=0); sdad slotwise partial heal (poison=6=2x3, siblings kept). OMP 1 vs 4 byte-identical. `test_knop_ray` 49/49. | VERIFIED |
| W3 | `comp_j_blank.f`, `control_variable_mod.f`, `rd_control_variables.f`, `lam_handoff_mod.f90`, `solve_for_pops.f`, `cmfgen_sub.f`, `drv_conv_exit.f90`, `test_conv_exit.py`, `drv_vadat.f90`, `test_vadat_keys.py` | `ES_ALI_PASS` per-frequency budget separate from `J_IT_COUNTER` (<= `ES_ALI_MAX`, clamp raised 20 -> 50); per-call `ES_ALI_CALL_NUNC/NFALL/DJMAX` exported; `ES_FIELD_NOT_READY` derived in `solve_for_pops`; `ES_ALI_GATE` OFF\|WARN\|HARD (default WARN); HARD adds the conjunct to `CONV_EXIT_EVAL` and both `CONV_EXIT_MET` sites; ITER_SUM ES_ALI line prints census + gate. | `test_conv_exit.py` 53/53 (new ES pins: HARD blocks, WARN/OFF do not, nine physical conjuncts unchanged); `test_vadat_keys.py` 95/95 (parse, case-normalize + `[ES_ALI_GATE]` MODEL_SCR echo, invalid->WARN, clamp 50). | VERIFIED (A/B BLOCKED, sec. 5) |
| W4 | `rd_control_variables.f`, `drv_vadat.f90`, `test_vadat_keys.py` | Absent-key auto-set: `DO_JFNK=T` + absent `JFNK_PC_QUAL` -> T + `JFNK auto-set` banner; solver-module mirror re-synced through `jfnk_pc_set_options_v1`. Explicit F = warned opt-out; `DO_JFNK=F` leaves F. | `test_vadat_keys.py` 95/95 (baseline auto-set, explicit F warning, explicit T kept, DO_JFNK=F unchanged). gfortran `test_pc_quality`: Q2c blindspot rejected (itry=-1), Q2d accepted itry=0 — all pass. | VERIFIED |
| W5 | `fiddle_pop_corrections_v2.f`, `drv_fpc.f90`, `test_fpc.py` | Require `0 < RELAX <= 2` (`.NOT.(>0)` also rejects NaN); require `1 <= L_ST <= L_END <= ND`; invalid inputs warn loudly on `ERROR_LU()` with requested+applied values and reset to safe defaults; T_LIM/POP_LIM checks and Lambda override unchanged. | `test_fpc.py` 32/32: relax0 warns+applies 1.0; relaxhalf/relaxrange pin RELAX semantics; badorder (L_ST>L_END) and lbeyond (L_END>ND) warn + full-range defaults; lamadj pins the Lambda override; prior 4 cases unchanged. | VERIFIED |
| W6 | `cmfgen_sub.f`, `control_variable_mod.f`, `solve_for_pops.f`, `test_iter_output.py` | `STEQ_T_EHB` -> per-shell luminosity (same `r^2*T1` + `LUM_FROM_ETA_V2` recipe as `DEP_RAD_EQ`); worst-shell share of `\|RLUMST(ND)\|` + depth exported via `EHB_RES_*`; printed at `CONV_DIAG_LVL>=1` in ITER_SUM and the final budgeted-iteration report; `N/A` when unassembled. No gate. | `test_iter_output.py` all-pass (+3 pins: token present at >=1, N/A branch in both units, absent from `CONV_EXIT`). | VERIFIED (diagnostic-only, no runtime arm) |

## 3. Behaviour changes and migration notes

- **NG acceleration**: a non-finite band history or fit result now
  refuses the update (`NG_DONE=F`, `POPS` untouched, census + one
  log line) instead of silently committing NaN. Old behaviour on
  invalid data is not reproducible (that is the point); for a fail
  -fast variant set `NG_NF_ABORT=T`. Valid-data path is bitwise
  identical.
- **Knop retained store**: a non-finite/non-positive retained row is
  healed (zero offending slots, `HEAL_STORE` convention) at the
  seed read, at the veto with one re-evaluation, and in the
  retained `NC=0` branch. Genuine jumps on a sane store are vetoed
  exactly as before (`FG_POISON_HEAL` census in the refresh
  summary). No opt-out needed — sane stores are untouched.
- **ES_ALI_GATE**: default `WARN` prints
  `ES_ALI gate WARN: <n> freqs unconverged at cap, max dJ/J <x>
  (not blocking exit)` in `Model Convergence`. `ES_ALI_GATE=OFF`
  reproduces the pre-change behaviour bitwise. `HARD` adds
  `es-ali-unconverged` as a tenth `CONV_EXIT_EVAL` conjunct (pass
  count stays the nine physical gates; the label only fails
  `ALL_MET`). `ES_ALI_MAX` may now reach 50 (default still 20); the
  ES loop has its own `ES_ALI_PASS` budget and is no longer cut at
  the FEDD `J_IT_COUNTER` 20.
- **JFNK_PC_QUAL**: absent key under `DO_JFNK=T` now auto-sets T
  (banner) instead of defaulting F; the probe rejects
  invertible-but-worthless factorizations
  (`||I - M^-1 A|| > JFNK_PC_QUAL_TOL=0.5`) at the cost of one
  matvec+norm per factorization. Explicit `JFNK_PC_QUAL=F`
  restores the old default with a one-line warning.
- **ADJUST_CORRECTIONS**: `RELAX=0` and `L_ST > L_END` are now
  loudly repaired (applied values 1.0 / full grid) instead of
  silently no-opping; messages moved from unit 6 to `ERROR_LU()`.
  Valid files are byte-identical; the Lambda-iteration override is
  untouched.
- **Energy balance**: new `CONV_DIAG_LVL>=1` diagnostic line only;
  nothing to restore — no key, no gate.

## 4. Tests run

Commands (all in `scratch/cmfgen_convguard_wt`, ifort 2021.11.1 +
mkl loaded):

```
python3 unit_tests/python/test_ng_abstain.py      -> 23 passed, 0 failed
python3 unit_tests/python/test_knop_lagged_map.py -> 34 passed, 0 failed
python3 unit_tests/python/test_knop_lagged_field.py -> all checks passed (after 938e534a pin retarget)
python3 unit_tests/python/test_knop_ray.py        -> 49 passed, 0 failed
python3 unit_tests/python/test_vadat_keys.py      -> 95 passed, 0 failed
python3 unit_tests/python/test_conv_exit.py       -> 53 passed, 0 failed
python3 unit_tests/python/test_fpc.py             -> 32 passed, 0 failed
python3 unit_tests/python/test_iter_output.py     -> all layout checks passed
python3 unit_tests/python/test_ave_dc.py          -> 9 passed, 0 failed
python3 unit_tests/python/test_alteration_coverage.py -> 3 passed, 0 failed
```

Named tails:

- knop_lagged_map: `sdab refresh 2 VC=0 poison=16`; `sdac VC=1
  jump=1e6 vetoed`; `sdad poison=6`; `34 passed, 0 failed`.
- OMP 1-vs-4 A/B on the lagged-map driver: outputs byte-identical
  (`diff` clean), `FG_POISON_HEAL` totals match.
- W1 direct binary evidence: `PREVRJ` + one NaN -> `NG_DONE=F`,
  `NG_NF_HIST=1`, POPS untouched; overflowed fit ->
  `NG_NF_NEW=1`.
Full `./build.sh` (24-Nov, ~1 h): compile clean, all four exes
linked (`cmfgen_dev`, `cmf_flux`, `obs_frame`, `plt_steq`);
`tests/f90/build_and_run.sh` ran inside it — 21 solver + 9
interface checks all PASS; `run_all_tests.py` ran to completion:
**99/107 test files green**. Of the 8 reds, 6 are the documented
pre-existing gaps (unchanged from the 22-Nov ledger entry):

- `test_auto_ion` — fixture model dir `m15_2p82d_VAL_stock` no
  longer exists on disk (not a code defect).
- `test_es_conv` — `drv_knop_guard` link misses
  `fg_j_cmf_nonmono_v1.o` (`set_fg_q_resid_`), the known obs/
  archive gap; predates this work, identical on main.
- `test_par_opac` — `opac/mod_par_opac.f` has no production object
  in any archive (same gap class).
- `test_phase227_io` — runner CLI gap (`--prepare/--binary/--run`
  not threaded), as documented.
- `test_rv_snhyd` — the four PINNED defects (non-finite kappa,
  non-finite/negative density, missing sigma); unchanged.
- `test_sobew_angular_floor` — fixture `m15_2p82d_FRESH_sob_dieoff`
  RVTJ missing (run dir cleaned); 8/8 checks passed before the
  missing file aborted the read.

The other two were follow-ups resolved this session:

- `test_knop_lagged_field` — W2 moved the row heal into
  `KN_HEAL_ROW`; the source pin scanned `HEAL_STORE`'s delegating
  body. Retargeted in `938e534a`; test now all-PASS.
- `test_ave_dc` — `He2_F_OSCDAT` fixture lives outside the repo
  at `TREE.parent/handoff_smoke/`; resolved with a
  `scratch/handoff_smoke -> ../handoff_smoke` symlink (same
  pattern as `scratch/cmfgen_reference`). 9/9 green.

`build.sh` exits non-zero on the remaining gaps by design
(`set -euo pipefail`; `run_all_tests.py` returns the fail count) —
the same state as main.

gfortran core check (11.5.0 + MKL 2026.2, `-O0 -g -fcheck=all`):

```
test_newton_core  ALL TESTS PASSED
test_guards       ALL GUARD TESTS PASSED
test_ptc_sign     TEST_PTC_SIGN DONE
test_nm_window    ALL NM-WINDOW TESTS PASSED
test_pc_quality   ALL PC-QUALITY TESTS PASSED (Q2c reject, Q2d accept)
test_schur        ALL SCHUR TESTS PASSED
test_gmres_edge   all checks passed
```

`test_precond_blocks` and `test_trace_state`: BLOCKED under
gfortran (`__INTEL_COMPILER` preprocessing and
`set_kind_module.mod` dependency, per the work order's own
expectation).

## 5. A/B status

- **W2** (sdab/sdac/sdad): complete — all three pins pass on
  production objects; OMP 1/4 identical.
- **W3** (`ES_ALI_GATE` OFF vs WARN, budget 20 and 50): BLOCKED —
  per the work order no production Slurm jobs were submitted.
  Prepared configs cloned from `dessart_m15_2p82d_sobacc_x2.yaml`
  (same warm-restart checkpoint, `DO_ES_ALI=T`, only
  `ES_ALI_GATE`/`ES_ALI_MAX` differ; `install_dir` repointed to the
  worktree so the branch build is exercised):
  - `configs/dessart_m15_2p82d_esgate_off20.yaml`
  - `configs/dessart_m15_2p82d_esgate_warn20.yaml`
  - `configs/dessart_m15_2p82d_esgate_off50.yaml`
  - `configs/dessart_m15_2p82d_esgate_warn50.yaml`
  Metrics to record when run: iterations to exit, max final dJ/J,
  `ES_ALI_NUNC`, wall time, final spectrum.

## 6. Known gaps

- W3 A/B (above) is the only outstanding production-scale evidence;
  until it runs, `ES_ALI_GATE` default stays `WARN` (non-blocking)
  and `ES_ALI_MAX` default stays 20.
- Coupled/banded JFNK, moment-path non-contraction, and the
  cold-start initial-guess work remain open as before (untouched by
  this order).
- The energy-balance line is diagnostic-only by design; a future
  hard energy gate still needs a shocked-CSM-suitable scale.
- `sdp` (diverged-ray lagged-map case): with the entry-side heal the
  formerly divergent poison case now heals-then-commits; the test's
  finite+class assertions still hold. K1 (production-scale stall)
  remains open — the manufactured grid cannot reach rho~0.97.

## 7. Status

`PASS WITH LIMITATIONS — usable only within the documented gates.`

(All six work items verified on production objects; the only
BLOCKED item is the W3 production A/B, which the work order itself
forbids submitting without cluster access.)
