# Stiffness audit: dev7 MOM_KNOP_REFRESH_V1 on M15 2.82 d KNOP_MOM

Date: 2026-09-28  
Model: `cmfgen_validation/m15_2p82d_knop_MOM_dbg/cmfgen/2.82days`  
Commit of the running binary: `23da8ed`  
Reproduction script: `audits/knop_mom_stiffness_audit.py`

Status tags follow repo convention: VERIFIED / SUPPORTED / UNVERIFIED /
CONTRADICTED / BLOCKED.

## 1. What was investigated

The refresh in `MOM_KNOP_REFRESH_V1` fails to converge; the gate keeps
`GS_FREEZE` on.  Earlier output reported the worst-conditioned cell at
`MOM sweep: maxchg ~ 3E3` and `min TB/(|TA|+|TC|) = 0.166` at
expanded depth 22, frequency ~1714.  The inner shock also prints
`max |SIGMA_TRUE| = 1.2E2`.

We mapped the expanded grid back to the 70-point hydro profile and
inspected the per-cell diagnostics in `MOM_J_NONMONO_ERRORS` to decide
whether the problem is:

- the huge `SIGMA_TRUE` near the inner shock, or
- a numerical ill-conditioning in the moment-equation assembly.

## 2. Grid mapping: expanded depths to physical radius/velocity

The 70-point `SN_HYDRO_DATA` for this model has a strong velocity jump
at the inner boundary:

- standard depth 60: R = 1.8069E4 (1E10 cm), V = 1611 km/s
- standard depth 61: R = 1.6805E4, V = 6593 km/s
- standard depth 66: R = 1.4598E4, V = 5883 km/s (minimum R in this jump)

The standard-grid `SIGMA = dlnV/dlnR - 1` computed from this profile is
as low as `-120.25` at standard depth 66.

The Knop MOM solver expands 70 -> 212 points using
`ABS(V(i)-V(i+1))/T2` with `T2 = 50 km/s`.  The largest interval
(60->61) receives 100 inserted points.  Mapping of the OUTGEN
worst-conditioned depths:

| expanded depth | parent std depth | R (1E10 cm) | V (km/s) | SIGMA_TRUE |
|---|---:|---:|---:|---:|
| 22 | 21 | 6.1394E4 | 58.06 | -1.581 |
| 33 | 32 | 5.2391E4 | 64.97 | -1.854 |
| 103 | 60 | 1.6805E4 | 2154.52 | -4.981 |
| 106 | 60 | 1.6772E4 | 2179.91 | -5.004 |
| 115 | 61 | 1.5943E4 | 2679.57 | -5.145 |
| 184 | 66 | 1.4598E4 | 5883.11 | -120.249 |

**Finding 1 (VERIFIED):** the `max |SIGMA_TRUE| = 120` feature is not
an interpolation artefact.  The 70-point hydro profile itself contains
a near-discontinuity between standard depths 60 and 66.  The Steffen
monotonic cubic used by `MON_INT_FUNS_V2` does not overshoot: the
expanded-grid minimum is `-120.25`, identical to the standard-grid
minimum.

## 3. Where the equations actually struggle

`MOM_J_NONMONO_ERRORS` dumps every cell whose new J flips negative.
For each dumped cell we recomputed the full comoving-frequency
gradient factor that drives the Knop sign split:

```
GAM_FACTOR = 1 + SIGMA_TRUE * K_ON_J
GAM_MOM    = (V/R)/c * GAM_FACTOR
```

The smallest values of `|GAM_FACTOR|` are:

| expanded depth | freq (Hz) | K_ON_J | SIGMA_TRUE | 1+SIGMA*K |
|---:|---:|---:|---:|---:|
| 106 | 6.44E+00 | 0.1999 | -5.004 | -3.0E-4 |
| 102 | 9.39E+01 | 0.2010 | -4.970 | +1.0E-3 |
| 114 | 1.17E+01 | 0.1954 | -5.122 | -1.1E-3 |
| 103 | 9.35E+01 | 0.2002 | -4.981 | +3.1E-3 |
| 87  | 7.09E+00 | 0.2036 | -4.915 | -7.0E-4 |
| 92  | 8.96E+01 | 0.2029 | -4.923 | +1.0E-3 |

These are all in the shock-transition band (expanded depths ~87-116,
standard parent depths 59-61), not at the `|sigma| = 120` inner shock.

**Finding 2 (VERIFIED):** the moment equations are ill-conditioned at
the *zero-GAM resonance*, where `1 + SIGMA_TRUE*K_ON_J` passes through
zero.  For `sigma ~ -5` this happens at `K_ON_J ~ 0.20`, exactly the
values seen in the dump.  At this resonance:

- `GAM_MOM` changes sign, so the Knop upwind direction flips.
- `chi_hat` reduces to the ordinary opacity `chi` because the
  velocity term `|GAM_MOM|/dlnnu` vanishes.
- The J-equation diagonal `TB` loses its positive frequency-coupling
  contribution `PSI`, and the spatial/off-diagonal terms (`HL`, `HU`,
  `EPS`) dominate.  This matches the printed `min
  TB/(|TA|+|TC|)=0.166`.

**Finding 3 (SUPPORTED):** the inner `sigma = -120` region is *not* the
primary source of ill-conditioning.  There `|1 + sigma*K|` is large,
so `chi_hat` is enormous and the diagonal is well conditioned.  The
struggle originates in the moderate-sigma transition zone that the
shock creates (sigma ~ -5 near the base of the velocity jump), where
`K_ON_J` can cross the resonance value.

## 4. Why the refresh cannot converge

`MOM_KNOP_REFRESH_V1` is a lagged fixed-point iteration that sweeps
frequencies in alternating directions.  For a Jacobi/Gauss-Seidel map
to contract, the implicit diagonal of each frequency solve must
control the explicit upwind coupling.  At the zero-GAM resonance the
implicit term `PSI = DTAUONQ * |GAM_FACTOR|` collapses, while the
explicit source still contains a signed `GAM_FACTOR` that is nearly
zero.  The sweep therefore becomes a pure spatial transfer step with a
very weak diagonal and the lagged coupling between adjacent
frequencies propagates errors rather than damping them.

This is a property of the Knop moment closure in the form implemented
here, not a coding error.  It is also not removed by the existing
`MOM_DAMP_NM` relaxation, because damping only slows the map; it does
not restore the lost diagonal dominance.

## 5. Candidate fixes

All fixes must be flag-gated so the stock monotonic path is unchanged.

### 5.1 Regularize the GAM factor near zero (least invasive)

Where `|1 + SIGMA_TRUE*K_ON_J|` falls below a threshold, replace it
with a limited magnitude before it is used in `chi_hat`, `PSI`, the
sign test, and the explicit source.  For example:

```
GAM_FAC  = 1 + SIGMA_TRUE*K_ON_J
GAM_LIM  = SIGN(GAM_FAC) * MAX(ABS(GAM_FAC), GAM_FAC_MIN)
```

`GAM_FAC_MIN = 0.01` keeps the sign split and `chi_hat` finite
without changing the asymptotic expanding or compressing branches.
This directly removes the zero-GAM resonance.  It is a moment-closure
regularization, not a physics change, because the resonance is a
mathematical artifact of the variable-Eddington approximation.

**Status (28-Sep-2026): implemented and rebuilt.**  The regularization
is applied to both the depth-point factor `1 + SIGMA_TRUE*K_ON_J`
and the midpoint factor `1 + AV_SIGMA_TRUE*NMID_ON_HMID`, and to the
lagged PREV copies.  The sign is preserved; only cells with
`|1 + sigma*K| < 0.01` are affected.  The rebuild (`d_plane` and
`d_new_main`) completed cleanly.

Storage note: the four factor arrays live at module scope in
`MOD_MOM_J_NONMONO_V1`.  An initial version declared them as
subroutine-local allocatables, which are deallocated at RETURN while
the NM allocation block only runs on grid rebuild; job `12010916`
crashed ~30 s into its warm restart on exactly that path (second
solver call, unallocated reference).  They are now module-level,
allocated with the NM set at build, deallocated on teardown, and the
never-read midpoint endpoint `ND_NM` is initialised to 1.0.

### 5.2 Enforce diagonal dominance in the tridiagonal assembly

Add a Rusanov-like term to `TB` equal to the magnitude of the
off-diagonal `EPS` contribution that would otherwise weaken the
diagonal.  This keeps the sweep an M-matrix and prevents sign-flip
instabilities, but it changes the discrete equations more broadly than
5.1.

### 5.3 Treat the shock as a jump condition

If the velocity profile is physically a contact discontinuity, the
finite-difference `sigma = -120` is not resolvable in the moment
framework.  A jump-condition treatment (matching fluxes across the
interface instead of differencing through it) would remove the
unphysical gradient.  This is the physically cleanest fix but the
largest change.

### 5.4 Pre-smooth the hydro profile

The existing `scripts/smooth_sigma_hydro_data.py` can spread the
velocity jump over more standard-grid points.  This lowers `sigma`
magnitudes and shifts the resonance, but it also changes the model
being solved.  It should only be used if the sharp jump is judged to
be an interpolation/ extraction artefact rather than a real
hydrodynamic discontinuity.

## 6. Recommendation

Implement **5.1** first.  It is a one-location change inside the
Knop-flagged block of `mom_j_cmf_nonmono_v1.f`, it leaves the
expanding/compressing limits unchanged, and it can be tested against
the current debug run without touching the hydro profile.  The
parameter `GAM_FAC_MIN` should be a VADAT or `PARAMETER` value with a
comment explaining that it regularizes the moment closure at the
`1 + sigma*K_on_J = 0` resonance.

**5.1 is in the binary (`exe/cmfgen_dev.exe`, rebuilt 13:35 28-Sep-2026
after the module-scope fix).**  The crashed `12010916` was resubmitted
through `run_model_grid.py` as debug FINAL `12012723` (interactive,
cold start -- the directory refresh wipes warm-start checkpoints).

If 5.1 alone does not let the refresh reach `GS_FREEZE` release,
combine it with 5.2.

Fixes 5.3 and 5.4 are reserved for a separate decision because they
change either the discretization physics or the input model.

## 8. Validation plan

1. Confirm the resubmitted debug job (`12012723`) starts from the
   13:35 binary and does not crash in the first few iterations --
   this is also a regression test for the module-scope fix.
2. Check `MOM_J_NONMONO_ERRORS` after a few iterations: the
   zero-GAM resonance cells (depths ~87-116, `|1+sigma*K_on_J| < 0.01`)
   should no longer appear, or should appear with a much smaller
   count.
3. Check `OUTGEN`: `min TB/(|TA|+|TC|)` should move away from 0.16,
   and `MOM_SW_MAXCHG` should fall below the gate threshold of 10 over
   successive iterations.
4. If the debug run reaches the `GS_FREEZE` release, submit the full
   normal MOM and the mature LAM requeue chain through
   `run_model_grid.py` for the definitive test.

## 7. Files and line ranges

- `plane/mom_j_cmf_nonmono_v1.f` lines 1447-1509: `GAM_FACTOR`, sign
  test, `chi_hat`, `Q_NM`, and `PSI/PSIPREV` all use `1 +
  SIGMA_TRUE_NM*K_ON_J` directly.  This is the single place to insert
  the regularization.
- `scripts/grid_mapping/solver_metrics.py` lines 170-197:
  `sigma_solver` and `apply_sigma_cap` show that the hydro sigma cap
  only clamps `-1.2 < sigma < -1.0` and therefore does not affect the
  resonance at `sigma ~ -5`.
- `audits/knop_mom_stiffness_audit.py`: reproduces the grid mapping
  and the resonance table above.
