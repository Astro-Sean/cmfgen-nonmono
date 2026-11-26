# KNOP remediation WP-9 -- globalise the first population step after a freeze release

Date: 27-Nov-2026. Status: LANDED pending push (flag defaults F;
the armed arm is the human A/B).

## The defect

The release controller overshoot: the first BA step after a freeze
release demanded ~1E7, was capped at MAX_LIN_COR (400%), committed
that step, and re-poisoned the field the freeze had just converged --
the run reverted FULL->LAMBDA and the release/refreeze cycle
restarted. LAMBDA release steps were already damped at GS_LAM_CAP
(solve_for_pops release window); non-LAMBDA releases relied on the
DO_BA_RES_ACC net alone, which judges a step only after it commits at
whatever radius BA_TRUST_CUR happened to hold from the pre-freeze
ratchet (~1E4 on the stall -- effectively no cap).

## The fix

New VADAT flag `BA_REL_TRUST` (default F): when armed, a saved
`BA_REL_SEED` flag sets on the frozen->live transition (the same
branch that arms GS_JUST_REL), survives LAMBDA calls, clears on
refreeze or GS_PRECOND_ACT, and is consumed by the first non-LAMBDA
call that reaches the SOLVEBA_V13 site -- MOM_TRIAL calls are
excluded like the damp window already does. The consume re-seeds

  BA_TRUST_CUR = MAX(1.0_LDP, MIN(BA_TRUST, GS_COMMIT_RES))

i.e. the trusted per-column correction scale equals the committed
field residual the release was judged on (WP-3 KN_GATE_COMMIT_RES
measurement), floored at the same 1.0 the trust region already
enforces and capped at BA_TRUST. With the seed, the trust test
quarantines runaway columns at a field-consistent radius and the
released step commits bounded corrections instead of the 400% cliff.
The step still arms BAACC_PENDING inside SOLVEBA_V13, so it passes
the existing residual-acceptance test with the resacc backtracking
(reject -> rollback -> halve) unchanged -- that half of the spec was
already structurally present. The strict all-ray-FULL release gate
is untouched. GS_COMMIT_RES < 0 (no committed measurement yet) skips
the seed, leaving the pre-freeze radius.

A CONV_DIAG_LVL>=1 line and the ITER_REASON token 'ba_rel_trust'
record each seeding.

## Verification (all actual outputs)

- test_lam_handoff.py: 66/0 (+6 source pins): flag declaration,
  arm on the transition, clear with N_REL_SINCE on freeze, the
  non-LAMBDA consume site, the seeded formula, and BA_REL_TRUST
  parsed default F.
- test_vadat_keys.py: 138/0 (+2): BA_REL_TRUST absent -> F; armed
  -> T on the module.
- Flag-off battery: test_iter_ledger 33/0, test_iter_conv_gate
  18/0, test_conv_exit 53/0, test_knop_ordering 50/0,
  test_knop_frc 35/0, test_knop_osgn 14/0, test_gfm 21/0,
  test_mat5pen 9/0, test_cmf_blkband_pin PASS -- all identical.
- mod_subs rebuild clean (all dependents recompiled after the
  CVM addition).

## Flag-off identity

BA_REL_SEED only arms when BA_REL_TRUST is set, so with the key
absent or F the release path is byte-identical. The consume is a
single MAX/MIN on module state -- no array work, no exported
channels.

## Owed measurement

- BA_REL_TRUST=T on the m15 stall model: the first BA step after
  each release should commit at <= MAX(1.0, GS_COMMIT_RES) per
  column instead of the 400% cap; watch ITER_LEDGER for
  'ba_rel_trust' followed by resacc accepts that grow the radius
  1.25x per clean step, and for the FULL->LAMBDA revert count to
  fall.
