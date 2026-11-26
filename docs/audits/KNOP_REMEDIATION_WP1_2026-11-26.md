# KNOP_REMEDIATION_WP1 -- defaults, invariants, domain census, 26-Nov-2026

WP-1 of the Knop remediation campaign: flip the WP-0-diagnosed
defaults, make the tolerance ladder a hard invariant set, expose
the maturity gate read-only, and domain-mask the all-ray census.
Config-level change only -- no algorithm was touched. WP-0's
fixture (`test_knop_ordering.py`) is the regression gate.

## Changes

### Default flips (`rd_control_variables.f`)

| key | was | now | route |
|-----|-----|-----|-------|
| `KN_RAY_ORDERED` | T | F | all |
| `KN_NI1_FULL` | F | T | all |
| `KN_RAY_BULK_TOL` | 6.0E-2 (KNOP only) | 0.0 | all |

`KN_RAY_ORDERED=F` selects the alternating `SWEEP_FLD` passes:
the WP-0 fixture measured the ordered march's claimed freshness
broken on inserted ray points (WINX>0 reads the deeper row a full
pass stale; plateau at NXIT=6, RESID 4.1e-2 vs the alternating
solve's clean NXIT=1 in 3 passes). `KN_RAY_BULK_TOL=0` removes
the measured-floor certification path (IXIT=2) everywhere -- the
floor was bounded non-contraction of the lagged map, not
convergence (CHECKLIST nbulk arm); every production ray must now
reach `KN_RAY_TOL` to certify. `KN_NI1_FULL=T` ends the permanent
ninth-conjunct deadlock on single-point rays; combined with the
I4 census the vacuous ray drops out of the verdict denominator
entirely.

`FG_COMMIT_NC` stays 0.05 but is now a promoted variable and its
explicit-set state is staged at read time (see I3 below).

### Promoted variables (`control_variable_mod.f`)

`KN_RAY_TOL_V`, `KN_RAY_BULK_TOL_V`, `FG_COMMIT_NC_V`,
`GS_MATURE_TOL_V` are module variables (VADAT echo convention),
not `rd_control_variables` locals. `GS_MATURE_TOL` itself moved
to a CVM `PARAMETER` (5.0E-2) -- `solve_for_pops` already `USE`s
CVM, so its local copy is deleted and the gate has one source.

### `GS_MATURE_TOL` VADAT exposure (read-only)

The key is read in the GS_ key block of `rd_control_variables.f`.
`CHECK_PARAM_CONSISTENCY` stops the run on any non-default value
-- the ladder invariants are arithmetic on the gate, so letting
the gate float would silently rescale every stop bound.

### Hard invariants (`check_param_consistency.f`)

Active on the FG-census routes (`DO_NONMONO_KNOP .OR.
DO_SOB_FG_KNOP`); each is a STOP in `CHECK_PARAM_CONSISTENCY`,
not a warning:

- I1: `KN_RAY_BULK_TOL < GS_MATURE_TOL` -- a ray may not certify
  on a looser bound than the gate reads mature.
- I2: `KN_RAY_TOL <= 0.1*GS_MATURE_TOL` -- the certified field
  must be at least an order tighter than the gate.
- I3: an explicitly-set `FG_COMMIT_NC` must not equal
  `GS_MATURE_TOL` -- a capped ray then pins `GS_MAX_ALL` exactly
  at the gate, a freeze/release limit cycle. Exemption: the
  absent-key legacy default (0.05 == 0.05) is allowed until WP-4
  replaces the fixed cap, so legacy decks still start. Explicit
  equality hard-stops.
- I4 (census, `fg_j_cmf_nonmono_v1.f`): `KN_NACT`/`KN_NFULL`
  count only production rays with a transfer domain
  (`NI_RAY >= 2`). `FG_RF_ALL_CONV` gains an explicit
  `KN_NINV = 0` conjunct over the full production mask, so a
  non-domain ray left `KN_CS_INV` (flag off, or invalid
  geometry) still fails the verdict closed.

### I3 staging detail

`CLEAN_RD_STORE()` deallocates the VADAT store before
`CHECK_PARAM_CONSISTENCY` runs, so explicit-key presence cannot
be queried in the checker. `RD_CONTROL_VARIABLES` captures it at
read time into `FG_COMMIT_NC_SET` (CVM) immediately after
`SET_FG_COMMIT`.

## Verification

`test_vadat_keys.py`: 114/0 (105 baseline + 9 new WP-1 checks:
defaults echo, `GS_MATURE_TOL` default + read-only abort, I1
pass/fail incl. `DO_SOB_FG_KNOP` route, I2 pass/fail, I3
explicit-equality stop + absent-key pass).

`test_knop_ray.py`: 50/1. The I4 census required re-pinning the
NI=1 tests: `ni1` (flag off) now blocks `RF_ALLCONV` through the
explicit `KN_NINV` conjunct (T29a); `ni1f` became the 3-ray mixed
ensemble -- ray 1 certifies vacuous FULL, drops out of
`NACT`/`NFULL`, and the two converged domain rays carry the
verdict TRUE (T29b); new `ni1s` pins the zero-domain edge --
single NI=1 ray under the flag leaves `NACT=0` and the empty
census fails closed (T29c). The one failure (T18, `diag=False`)
is pre-existing, identical signature on the pre-WP-1 run.

Gate-sensitive suite: `test_iter_conv_gate` 18/0,
`test_conv_exit` 53/0, `test_knop_ordering` 39/0 (WP-0 fixture
unchanged -- driver pokes module vars, bypasses rd path),
`test_fg_gs` 36/0, `test_knop_lagged_map` 25/1 (`sdab` KNOWN
DEFECT identical pre-WP-1), `test_knop_lagged_field` rc=1
(`test_seed_arm_probe` crash identical pre-WP-1).

## Verdict

PASS. Defaults, invariants I1-I4, and the domain-masked census
landed with unit coverage on both arms of every change; the only
suite failures are byte-identical pre-existing ones. STOP per
spec: Arm A (`KN_RAY_ORDERED=F, KN_NI1_FULL=T, KN_RAY_BULK_TOL=0`,
warm restart from it122 on `dessart_m15_2p82d_FRESH_knop_v14`)
is the human's to submit.
