# Population solver scientific-integrity fixes (2026-11-09)

Status: FINDING 2 core implemented and syntax-clean; FINDING 1 and
FINDING 3 analysed with concrete next steps. Build/test/fast-model
validation pending at time of writing. Baseline audited:
`5c59a53f` (post Item-2 floor-mask commits).

> **Tracking note (10-Nov-2026):** this audit is indexed in
> [AUDIT_LOG.md](AUDIT_LOG.md) with per-level validation status.
> Outstanding work (none verified complete): BA_QUAR bounded
> quarantine age, periodic raw-residual reassessment, justified
> release criteria, quarantine diagnostics, re-quarantine
> protection; POP_TRACE_MAXFAC raw-vs-applied diagnostic pair and
> cap-vs-raw unit tests; POP_MASK_FLOOR semantics/conservation/
> disabled-path tests and the matrix-pin decision; full
> `./build.sh`, `run_all_tests.py`, and fast-model validation.

Governing requirement: the population solver must converge the
original statistical-equilibrium system, not a modified system that
looks converged only because difficult equations were suppressed.

## FINDING 2 - POP_TRACE_MAXFAC hides raw demand (IMPLEMENTED, core)

Problem: `POP_TRACE_MAXCAP` clipped `STEQ` in place inside the
`DO_BA_RES_ACC .AND. .NOT. LAMBDA_IT` block, so `MAXCH`,
`BA_DEM_EXM`, `BA_DEM_EXR`, `BA_EX_MAXCH` and the level-check vectors
all saw the capped correction, not the demanded one.

Change (`new_main/subs/solveba_v13.f`):

- `STEQ_RAW` is now populated unconditionally on every call right
  after the band solve and before every safeguard (floor mask,
  quarantine, cap, trust). It holds the raw demanded fractional
  correction. Previously it was allocated only under
  `POP_RUNAWAY_DEM>0` for the RUNAWAY census.
- The convergence-facing demand scans now read `STEQ_RAW` instead of
  `STEQ`: the demanded-correction histogram, `DECREASE/INCREASE`
  (-> `MAXCH`), `EX_DEC/EX_INC` (-> `BA_EX_MAXCH`), the `BA_DEM_EXM`
  exempted-demand scan, the `BA_DEM_EXR` exit-exempted scan, and the
  `MAX_INC_VEC/MAX_DEC_VEC` level-check vectors.
- `STEQ` remains the applied correction: floor mask, quarantine,
  `POP_TRACE_MAXFAC`, NaN/RBK/trust column masks still modify it, and
  it is returned to the caller for the `POPS` update. `MAXCH_APP`
  still reports the applied change.
- Sign convention preserved (`STEQ>0` decrease, `<0` increase);
  non-finite raw corrections still hit the NaN-column zero and now
  also surface through `STEQ_RAW`-based metrics rather than being
  clipped away.

Effect: a clipped, masked or quarantined equation's true demand is
now visible to `MAXCH`/`EXM`/`EXR`/`BA_EX_MAXCH`, so clipping cannot
manufacture convergence. Floor-exempt rows (`POPS<BA_DEM_PFLOOR`)
still route to `BA_DEM_EXM` and `BA_EX_RELFRAC` rows to
`BA_DEM_EXR` - the physical negligibility criterion is unchanged.

Still TODO for FINDING 2: an explicit "largest capped demand +
applied correction for the same equation" diagnostic pair
(`POP_TRACE_MAXFAC`-specific), and unit tests for cap-vs-raw.

## FINDING 1 - BA_QUAR permanent quarantine (ANALYSED, not yet implemented)

Confirmed at `solveba_v13.f`: the LAMBDA quarantine block
(`DO_BA_RES_ACC .AND. LAMBDA_IT`, ~line 1398) zeroes `STEQ(J,I)` for
every `BA_QUAR_MASK` row, and the residual measurement rescans skip
quarantined equations - so a quarantined equation is removed from
both the update and the residual/`MAXCH`/`EXM` accounting. There is
no release path; `BA_QUAR_MASK` is permanent for the run. ED/T rows
are already excluded (`J=1..NT-2` only).

With the FINDING-2 change, a quarantined equation's raw demand now
still reaches `MAXCH`/`EXM`/`EXR` via `STEQ_RAW`, so it can no longer
silently satisfy the exit - but the mask itself is still permanent
and there is no recovery/re-evaluation or bounded-duration release.

Next steps (not yet done):
- Add a per-equation quarantine age + periodic reassessment; release
  when the raw residual/demand meets a justified criterion.
- Bound quarantine duration (`BA_QUAR_MAX` stays a safety cap, not
  convergence evidence).
- Diagnostics: quarantined count, newly-added, released, longest
  duration, largest raw residual and raw demand with species/level/
  depth; distinguish solver-exclusion from physical convergence.
- Ensure release cannot immediately re-park the same equation
  (state transition + reason logging).

## FINDING 3 - POP_MASK_FLOOR doc/impl mismatch (PARTIALLY DONE)

`control_variable_mod.f` now documents the shipped POST-SOLVE mask
(the earlier text described an unimplemented matrix-level pin in
`CMF_DIAG_SOLVE`; the matrix pin crashed the coupled block, job
12868544). `POP_FLOOR_STEP` (added in the prior Item-2 commit
`081d3fdd`) holds the level at its existing value - it does NOT
raise `POPS` to `POP_MASK_FLOOR`, and it does not regularise the
matrix. Masked demand is accounted via `BA_MASK_N/MAX/J/I` and
printed as `POP_MASK:` every call.

With FINDING 2, masked levels' raw demand also reaches the
convergence metrics, so a held equation cannot count as converged.

Still TODO: decide whether a genuine matrix-level pin is safe in the
coupled assembly (`FIXPOP_IN_BA` / `CMF_DIAG_SOLVE` / banded path),
or keep the documented post-solve hold; unit tests for floor
semantics, conservation and the disabled path.

## Files changed this round

- `new_main/subs/solveba_v13.f` - unconditional `STEQ_RAW`; metric
  scans read `STEQ_RAW` (raw demand) while `STEQ` stays the applied
  correction.

## Verification status

- `ifort -fsyntax-only` on `solveba_v13.f`: clean.
- `./build.sh` full build + `run_all_tests.py`: pending.
- Fast-model (`cmfgen_validation_fast`) A/B/C: pending.
- Commit: pending (will be a clearly-identified commit on `main`).
