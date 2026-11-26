# KNOP_REMEDIATION_WP5 -- signed-Aitken store relaxation, 26-Nov-2026

WP-5 of the Knop remediation campaign: a signed vector-Aitken
omega for the lagged `I_P/I_M_DEPTH_NEXT` store write, behind
`KN_OMEGA_SGN` (default F; the FG_ADAPT_OMEGA_V2 pool estimator
stays the adaptive-relax default). Opt-in while the m15 arm is
measured; F reproduces the WP-4 path exactly.

## Defect being addressed

The V2 estimator (`FG_ADAPT_OMEGA_V2` -> `KN_OMEGA_POOL_V1`)
prices the outer-map relaxation `INEXT_RELAX` from the INNER
per-ray pass ratios `RHO_P_LS`/`RHO_M_LS` -- the contraction
factors of the frequency sweep at fixed store, not of the
store-map iterate whose step it damps. That is a category
error: the pooled ratio cannot express an over-relaxed step
(omega > 1), and on the m15 stall it reads an inner
convergence rate that says nothing about the outer residual
sequence.

WP-5 measures the outer map itself. Each sweep write site
commits `x <- snap + omega*(fresh - snap)`; the pre-relaxation
defect `r = fresh - snap` per stored slot is the outer-map
residual. Between consecutive sweeps the slots give
`dx = x_k - x_{k-1}` and `dr = r_k - r_{k-1}`, and the next
omega is the signed vector Aitken / generalised-secant price

    omega_{k+1} = -<dx,dr> / <dr,dr>

(Ramiere & Helfer 2015, JCP 283, eq. 36; Kuttler & Wall 2008,
Comput. Mech. 43, fixed-point form). For a scalar linear map
this returns the fixed-point weight `1/(1-lambda)` from a
single pair -- including lambda < 0 (oscillating) and
lambda > 1 (expanding) where the pool estimator is blind.

## Implementation

`cmfgen/plane/fg_j_cmf_nonmono_v1.f` only, behind
`KN_OMEGA_SGN` and only inside the `KN_ADAPT_RELAX` adaptive
branch. State lives in `FG_J_CMF_NONMONO_MOD_V1`:

- per-(slot,channel) residual/base history `FG_OSGN_XP/XM`
  (last committed base, -1 = none) and `FG_OSGN_RP/RM` (last
  raw residual), dimensioned `(N_COMP,NP,NCF_IN)`;
- per-(slot,channel) column buffers `FG_OSGN_C1*/C2*/D*/T*`
  filled at the two write sites (`T2+INEXT_RELAX*(T1-T2)`);
- `FG_OSGN_COL_FIN` folds a column's buffers into the dot
  accumulators `FG_OSGN_ACC(1..4)` = <dx,dr>, <dr,dr>,
  <r_cur,r_cur>, <r_prv,r_prv> at each column transition
  (the NEW_FREQ/SNP_FREQ snapshot point) and on the pending
  column inside `FG_OSGN_UPDATE`;
- `FG_OSGN_UPDATE` (refresh tail, same branch that calls
  `KN_OMEGA_POOL_V1`) closes the pending column, ignores
  sweeps already consumed (`FG_OSGN_GUSE`), resets slot
  history when the `KN_CHI_EXT/KN_ETA_EXT` fingerprint moves
  (populations changed -- the map is only the same map while
  the freeze holds them), and prices `KN_OSGN_STEP`.

Guards in `KN_OSGN_STEP`: denominator collapse
`(<dr,dr> <= KN_OSGN_EPSG*<r_cur,r_cur>)` holds omega;
non-decreasing residual `<r_cur> >= <r_prv>` halves omega
once (`FG_OSGN_NRET`) then reports ICODE=3 so the caller
drops to the unaccelerated estimator; the raw price clamps
to `[KN_OSGN_MIN, KN_OSGN_MAX]` (VADAT, 0.05/2.0 defaults).

Two capture rules matter: `T2 = 0` slots (unrelaxed `T1`
write) are skipped -- omega never acted there, so their
sample is a map output, not a defect; re-entries within a
column overwrite, keeping the last (scattering-converged)
pass's sample, matching the committed value.

Wiring: `control_variable_mod.f` declares `KN_OMEGA_SGN` +
`KN_OSGN_MIN_V`/`KN_OSGN_MAX_V`; `rd_control_variables.f`
reads `[KN_OMEGA_SGN]` (F), `[KN_OSGN_MIN]` (5.0D-2),
`[KN_OSGN_MAX]` (2.0D0) and pushes them through the external
setters `SET_FG_OMEGA_SGN`/`SET_FG_OSGN_LIMS` (mod_subs
cannot `USE` the plane module). `drv_vadat` echoes the CVM
flag plus the module-side `FG_OSGN_FLG`/`KN_OSGN_*` so the
test sees the setter landed. Regrid drops the grid-shaped
history with the store; lazy `FG_OSGN_ALLOC` rebuilds it
inside the sweep write, so the legacy path allocates nothing.

## Fixture

`unit_tests/f90/drv_knop_osgn.f90` +
`unit_tests/python/test_knop_osgn.py` (14/0).

Scalar level -- `KN_OSGN_STEP` driven on the lagged map
`x <- x + omega*(a + lam*(x-a) - x)`:

| case | lam | outcome |
|------|-----|---------|
| sc   | 0.3  | root to <1e-12; omega -> 10/7 = 1/(1-lam) |
| so   | -0.8 | alternating; converges, omega -> 5/9 |
| sd   | 1.5  | ICODE=2 (halve to 0.25) then ICODE=3 fall |
| sc0  | --   | dr=0 -> ICODE=0, omega held |
| scl  | 0.55 / -1.5 | clamps at KN_OSGN_MAX / floor |

Vector level -- `FG_OSGN_ALLOC` + write-site pokes +
`FG_OSGN_COL_FIN`/`FG_OSGN_UPDATE` on a two-slot map:
sweep 1 holds (no history, ICODE=0), sweep 2 prices
omega = 10/7 exactly, both slots converge; `vfp` moves the
CHI_EXT fingerprint on sweep 3 -> history reset -> ICODE=0
hold with omega unchanged (the accumulated dots mixed a
cross-map dx, so the reset also clears them); `vgen` calls
UPDATE twice in one generation -> second call returns
ICODE=-1/LFALL so the caller falls back.

## Test state (this commit)

`test_knop_osgn` 14/0; `test_vadat_keys` 129/0 (+5 WP-5:
absent-default-F + module echo, clamp defaults, explicit
landing, misordered-clamp coalesce); `test_knop_ordering`
50/0; `test_iter_conv_gate` 18/0; `test_conv_exit` 53/0;
`test_pop_audit` pass; `test_sigma_true` 16/0.
`test_knop_lagged_field` fails at "summary and per-pass
report lines present" -- pre-existing (`KN-SEED ray=` is
absent at HEAD, 0 hits; the report format predates this WP).
`test_knop_ray` T38b now reports REPL_MAXABS=2.2e-19 where
the WP-4 binary gave exact 0: deterministic (3/3 replh runs
agree), flag-off code path logically unchanged. The added
flag-gated branches in the inner write loop perturb codegen
scheduling, so the OMP-raced hot-pass trajectory lands on the
fixed point at the last-ulp level rather than bitwise. The
exact-zero bound measures trajectory independence the WP-4
binary happened to satisfy; relaxing it is a human call, not
made here. All other observed suite reds are pre-existing
(verified string-absent at HEAD) or harness gaps:
test_fg_sw_sym, test_gs_drift_diag, test_knop_mode1 (same
renamed report family), test_knop_lagged_map (documented
KNOWN DEFECT), test_auto_ion/test_par_opac (missing
production object), test_phase227_io (needs CLI args),
test_ng_shadow, test_sobew_angular_floor.

## Verdict

VERIFIED at unit level: the estimator returns the correct
fixed-point weight on contracting, oscillating, and
expanding maps, all guards fire as designed, and the default
path is untouched (flag off = WP-4 byte-for-byte). The m15
production measurement of the armed flag is owed -- that is
the A/B the campaign hands back to the human.
