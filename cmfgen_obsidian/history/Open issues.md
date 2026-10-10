# Open issues

Known limitations, unvalidated machinery, and the honest gaps between
the design and what is demonstrated. Consolidated from
[[outstanding_issues]], [[future_fixes]], and the audits.

## Physics / solver open questions

- **MOM-path non-contractive map.** The `MOM_J_CMF_NONMONO_V1` lagged
  fixed-point map is expansive on strongly compressing cells
  (sigma ~ -37..-55): ~3% of the ND x NCF field sign-flips every pass,
  in-sweep too - a discretization-level property, not a plumbing bug.
  The ordered march does not apply (global tridiagonal, no depth
  march). Containment works (best-pass restore, abort, quarantine) and
  production uses `USE_LAM_ES=T`; if MOM convergence is ever needed
  the paper-conform direction is block/frequency-range iteration
  ([[BLOCK_GS_DESIGN]]).
- **Diagonal ALO only.** The paper's Lambda* is tridiagonal in
  wavelength; dev7 drops the upwind-frequency sensitivity. Same fixed
  point, weaker convergence; cross-frequency Jacobian blocks are out
  of scope for now.
- **Red-edge zero-incident BC**: compressing cells at the grid bottom
  keep `|GAM|/dln` absorption without the compensating source - a
  documented boundary-condition choice.
- **Aberration `3 a_nu I` term absent** - inherited from stock V13's
  nu-form equation; deliberate (proposed as future flag
  `KNOP_USE_FULL_NU_PROP_TERM` with checkpoint-version gating).
- **`DO_NONMONO_FREQ_MATCH`** implemented, default F; recommended ON
  for compressing models (drift ~330 bins/cell on m15).
- **SOB path anti-contractive on TypeIIn** even with RH78 non-local
  coupling - escape-probability line transfer is the diverger for
  that sigma structure; Knop+BLANK is the route for TypeIIn-class
  models.

## JFNK open questions

- **Coupled mode never yet produced an accepted Newton step.** The
  demand gate correctly defers it; endgame validation is BLOCKED on
  the model reaching the gate (demand still >> 1e2 in the last
  checked state).
- **JVP noise floor** (audit A1): the iterative sweep makes F noisy
  (~1e-4 relative); the coupled config uses `JFNK_EPS=1e-2` as the
  noise-aware forward-difference floor (opt/hyb configs use 1e-3).
  This is a mitigation, not a noise detector. Stronger fixes on record:
  deterministic eval sweep (fix inner
  counts, freeze repair/quarantine/branch decisions), or frozen-J
  residuals for Krylov directions with coupled evals only for
  verification.
- **Preconditioner is the wrong operator for coupled mode** (depth-
  block BA misses nonlocal dJ/dn): candidates are banded BA
  (`SOL_METH=TRI` + `JFNK_USE_BANDED`, needs the pipeline's DIAG
  forcing relaxed) or a physics/source-iteration preconditioner
  (Janett 2024) - deeper work.
- **Free-T path unverified** (`JFNK_FIXT=F`): latent state-leak note -
  per-frequency quadrature weights computed at perturbed T are not
  restored (heals at next real sweep; verify before experimenting).
  Related leaks closed 03-Oct: `SE%STEQ_ADV` and `SE%T_EHB` are now
  inside the shared `JFNK_SAVE/RESTORE_MOD_STEQ_V1` bracket, and the
  frozen-mode radiative split no longer leaves a coll-only module
  STEQ in place for the LAMBDA fallback (`9956744`, `e8abf97`).
- No backtracking line search on the committed Newton step - the
  admission gates + PTC + best-iterate fallback substitute; add if
  activations fail rather than defer.

## Performance candidates ([[future_fixes]] R2-R5)

- MOM-path coefficient caching (factorization rebuilt per pass).
- SEED_FLD/WRITE_STORE strided store access (~5 MB strided walks).
- KN_CHI_COEF/ETA_COEF single precision (~1 GB; NOT bit-identical -
  needs convergence A/B).
- Frequency-sweep OMP tail (~29% kmp_wait in one profile; unverified).
- Deliberately not done: tighter refresh exit criteria and
  single-precision field/store (they change the iterated object's
  convergence semantics).

## Upstream bugs found (not dev-tree regressions)

- `USE_J_REL=T` + `TRIDIAG`: zero pivot in `CMF_BLKBAND_V3` at depth 1.
- `'ZERO FLUX'` vs `'ZERO_FLUX'` string mismatch in `comp_j_blank.f`
  (fixed here; exists upstream too).
- Uninitialized `SUCCESS` on the NEWMOD path (fixed `7c6a3c3`).
- `rv_sn_model_sniin.f` SIGMA denominator (fixed `f2f4784`).
- `lte_sub.f` calls missing `WR_ABUND_INFO_V2` (excluded from build).
- CHK_NG `MEAN`/`RATIO` bugs (fixed; dormant under FIX_T).

## Validation still outstanding

- A converged full non-monotonic model run end-to-end through
  `run_model_grid.py` (the standing acceptance test).
- Frozen-J JFNK as routine SOLVEBA replacement (Phase 1) at scale.
- Ordered-march runtime A/B vs `KN_RAY_ORDERED=F` on a compressing
  model (theory-derived; pending runtime confirmation of the
  corner-weighted convergence rate).
- OBS spectrum from a *converged* KNOP model (the clean TypeIIn
  spectrum was from an unconverged it62 state - diagnostic only).
- `DO_NONLOC_SOB` production experiment (Jacobian-clean revision never
  cured the TypeIIn divergence; parked).

## Related

- [[Failure modes and guards]]
- [[JFNK admission]]
- [[Run log highlights]]
