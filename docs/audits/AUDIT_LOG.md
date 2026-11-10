# AUDIT LOG - active audit programme

Tracks the live audits on `main`: when the code was audited, which
revision, which commits landed the corrections, and what is still
owed. Dated entries are git commit dates; audit-doc dates use the
repo DD-Mon-YYYY convention and match the commit history.

Update this file in the same edit that lands audit-driven code.

## Active audits

| Audit | Doc | Audited rev | Code commits | Status |
|---|---|---|---|---|
| FG radiation-field drift | [FG_FIELD_DRIFT_AUDIT.md](FG_FIELD_DRIFT_AUDIT.md) | `7220b616` | `b49b32d3`..`e6e8f47a` (10-Nov-2026): A3 res-pass ordering, per-ray contraction estimator, same-census `FG_RF_CAPPED`, `FG_JUMP_MAX` OMP/finite fixes, `KN_OMEGA_POOL_V1`, `FG_RF_PEND`, ray-class split, `HEAL_STORE` census, eval-state bracket | Unit + pin suites green (test_knop_ray 41/41). Model-level verification still owed - run_model_grid job on Slurm, 1 vs N-thread comparison, field residuals and emergent spectrum |
| Fail-safe allocation follow-up | same doc, appended cycle section | `e6e8f47a` | `8fe9f4dc`..`280c4a44` (10-Nov-2026): `KN_JMP_SIBREF_V1` ISTAT, `FG_RF_FAILED` fail-closed maturity, `FG_JMP_FAIL_INJ` hook, shared `KN_DTAU_INT_V1` operator, eval-bracket push moved past last early exit | Same as above - bookkeeping proven by tests, drift behaviour needs the production run |
| SOB great-iteration acceleration | [ACCELERATION_AUDIT.md](ACCELERATION_AUDIT.md) + [ACCELERATION_TODO.md](ACCELERATION_TODO.md) | `280c4a44` | D1 `6b1fa6d4`..`9901e33b`, D2 `1a21f41b`..`c77114b5` (10-Nov-2026): ITER_LEDGER accel columns (contract v2+v3), ITER_SUM `accel` row, `GET_MOM_AA_STATE` + `MOM_AA_RST`, stale-counter gate on `MOM_RF_RAN`, outer AA decoupled from `MOM_RF_AA` | Phases 1-6 audited, diagnostics landed and green (test_iter_ledger 31/31, test_mom_aa_state 62/62). Benchmark arms E1-E7 and Phase-8 validation are Slurm work, not started |

## Rules for new entries

- One row per audit cycle (a new audit doc, or a follow-up cycle
  appended to an existing doc).
- `Audited rev` is the revision the audit inspected, not the rev
  that landed the fixes.
- `Code commits` lists the contiguous hash range the cycle landed;
  leave blank for audit-only findings.
- Status states what is proven (which test suite) and what is still
  owed (model-level run, benchmark arm, validation item).
- Filenames in this folder: `NAME_AUDIT[_YYYY-MM-DD].md` for
  point-in-time audits, undated for living docs like this log and
  the two active audit files above.
