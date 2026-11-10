# AUDIT LOG - active audit programme

Tracks the live audits on `main`: when the code was audited, which
revision, which commits landed the corrections, and what is still
owed. Dated entries are git commit dates; audit-doc dates use the
repo DD-Mon-YYYY convention and match the commit history.

Update this file in the same edit that lands audit-driven code.

## Validation status vocabulary

Each audit reports the highest level reached; a higher level is
never inferred from a lower one.

| Level | Meaning | Statuses |
|---|---|---|
| SOURCE AUDIT | Source inspected and findings recorded | NOT STARTED / IN PROGRESS / PASS / FAIL / BLOCKED |
| IMPLEMENTATION | Code changes landed | same |
| UNIT TESTS | Named tests and results recorded | same |
| BUILD | Build verification recorded | same |
| MODEL VALIDATION | Representative physical model runs completed | same |
| SCIENTIFIC VALIDATION | Populations, conservation and emergent spectrum compared against a baseline | same |

## Active audits

| Audit | Doc | Audited rev | Code commits |
|---|---|---|---|
| FG radiation-field drift | [FG_FIELD_DRIFT_AUDIT.md](FG_FIELD_DRIFT_AUDIT.md) | `7220b616` | `b49b32d3`..`e6e8f47a` (10-Nov-2026): A1 per-ray same-direction estimator + pooled P90, A3 res-pass ordering below divergence check, same-census `FG_RF_CAPPED`, `FG_JUMP_MAX` OMP/finite fixes, `KN_OMEGA_POOL_V1`, `FG_RF_PEND`, ray-class split, `HEAL_STORE` census, eval-state bracket, shared `KN_DTAU_INT_V1` |
| FG fail-safe follow-up | same doc, appended cycle | `e6e8f47a` | `8fe9f4dc`..`280c4a44` (10-Nov-2026): `KN_JMP_SIBREF_V1` ISTAT, `FG_RF_FAILED` fail-closed maturity, `FG_JMP_FAIL_INJ` hook, eval-bracket push moved past last early exit |
| SOB great-iteration acceleration | [ACCELERATION_AUDIT.md](ACCELERATION_AUDIT.md) + [ACCELERATION_TODO.md](ACCELERATION_TODO.md) | `280c4a44` | D1 `6b1fa6d4`..`9901e33b`, D1a `1a21f41b`..`c77114b5` (10-Nov-2026): ITER_LEDGER accel columns (contract v2+v3), ITER_SUM `accel` row, `GET_MOM_AA_STATE` + `MOM_AA_RST`, stale-counter gate on `MOM_RF_RAN`, outer AA decoupled from `MOM_RF_AA` |
| Population-solver integrity fixes | [POPULATION_SOLVER_FIXES_2026-11-09.md](POPULATION_SOLVER_FIXES_2026-11-09.md) | `5c59a53f` | FINDING 2 `STEQ_RAW` raw-demand metrics landed (see doc); FINDING 1 quarantine and FINDING 3 floor-pin analysis only |

## Per-audit validation status (10-Nov-2026)

| Audit | SOURCE AUDIT | IMPLEMENTATION | UNIT TESTS | BUILD | MODEL VALIDATION | SCIENTIFIC VALIDATION |
|---|---|---|---|---|---|---|
| FG drift (both cycles) | PASS | PASS | PASS - test_knop_ray 41/41, pin suites green | PASS - member objects compiled and archived each phase | NOT STARTED - run_model_grid job owed; 1 vs N-thread and emergent-spectrum checks still outstanding | NOT STARTED |
| SOB acceleration | PASS - Phases 1-6 audited | PASS - D1+D1a diagnostics landed | PASS - test_iter_ledger 31/31, test_mom_aa_state 62/62, test_knop_ray 41/41 | PASS - incremental member compiles (`solve_for_pops.o`, `iter_ledger_mod.o`, `mom_j_cmf_nonmono_v1.o`, `rd_control_variables.o`); full `./build.sh` not re-run this cycle | NOT STARTED - benchmark arms E1-E7 not submitted | NOT STARTED - Phase-8 validation not started |
| Population fixes | PASS | IN PROGRESS - FINDING 2 core landed; FINDING 1 analysis only; FINDING 3 partially done | NOT STARTED - cap-vs-raw and floor-semantics tests still owed | IN PROGRESS - `ifort -fsyntax-only` clean on `solveba_v13.f`; full `./build.sh` + `run_all_tests.py` pending per the audit doc | NOT STARTED - fast-model A/B/C pending | NOT STARTED |

### Population-solver outstanding items (from POPULATION_SOLVER_FIXES_2026-11-09.md)

- BA_QUAR permanent quarantine: bounded quarantine age,
  periodic raw-residual reassessment, justified release criteria,
  diagnostics for newly-quarantined and released equations,
  protection against immediate re-quarantine.
- POP_TRACE_MAXFAC: raw-demand vs applied-correction diagnostic
  pair; unit tests proving clipping cannot manufacture
  convergence.
- POP_MASK_FLOOR: floor-semantics tests, conservation tests,
  disabled-path regression, decision on a matrix-level pin.
- Full `./build.sh`, aggregate `run_all_tests.py`, and fast-model
  validation - all pending.

## Rules for new entries

- One row per audit cycle (a new audit doc, or a follow-up cycle
  appended to an existing doc).
- `Audited rev` is the revision the audit inspected, not the rev
  that landed the fixes.
- `Code commits` lists the contiguous hash range the cycle landed;
  leave blank for audit-only findings.
- Update the per-audit validation table in the same edit that lands
  audit-driven code; a level moves only when its named evidence
  exists (e.g. MODEL VALIDATION needs the run recorded, not just
  submitted).
- Filenames in this folder: `NAME_AUDIT[_YYYY-MM-DD].md` for
  point-in-time audits, undated for living docs like this log and
  the two active audit files above.
