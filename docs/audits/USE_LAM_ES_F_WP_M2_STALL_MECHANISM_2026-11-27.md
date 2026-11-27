# USE_LAM_ES=F WP-M1/M2 - compressing-branch stall, reproduced and named

Date: 27-Nov-2026. Evidence for `checklists/use_lam_es_eq_f.md`
WP-M1 and WP-M2. The fixture driver is
`unit_tests/f90/drv_mom_solve.f` (modes `rfst`/`rfj`), the
regression gate `unit_tests/python/test_mom_solve.py`
("rfj:" checks). Line numbers cite the working tree at commit
`da4aec56`.

## WP-M1 - the stall reproduces in the unit driver

Fixture (`rfj` mode): production depth structure on a single-cell
velocity jump (sigt=-60 at model cell 12 of 24): a thin,
scattering compressing band (arg10=10 scales CHI inside the band,
arg13=0.99 sets its albedo) inside an absorption-anchored column
(arg9=0). Tight log-frequency spacing (arg14=1E-4) and high outer
velocity (arg15=1E5) put the band's compressing rows in the
production marginal-map regime. Grid expands 24 -> ~4300 points
through `KN_NM_MAX_INS=256` capped intervals (the production CDS
analogue, mom_j_cmf_nonmono_v1.f:1165,:1482-1519).

Command (unit_tests/bin/mom_solve_ncf32, -DNCF_IN=32 build):

```
rfj nonj 3 1 diff none 1 <tag> 0.0 10 0.0 0.0 0.99 1e-4 1e5 0
```

Measured on HEAD:

- `Compressing direction stalled - dropping it after pass 42;
  the expanding direction continues.` - the production OUTGEN
  signature verbatim (m15 noali arm printed it x1221 at pass 22).
- Refresh exits `converged` at pass 54 behind the expanding
  branch (RFST_CONV=1).
- `nu-census` a-value max = 0.964 (RFST_AMAX; production
  9.947E-01); 7 of 131874 captured rows at a>0.9 (RFST_NA9),
  none at a>0.99, 0 xs<=0 rows.

Deviations from production defaults, both schedule knobs only:

- `MOM_RF_BULK_EXIT` off (driver arg16=0): the fixture's
  compressing drifts sit below the >50% count, so the bulk gate
  (mom_j_cmf_nonmono_v1.f:7007-7020) exits at pass 4 and masks
  the stall. Production's >50% counts evidently kept the gate
  closed past the stall floor.
- `NCF_IN=32` build flag: the jump's matched spans (NM_J~7) must
  stay on-grid at most frequencies; NCF=8 leaves only L=1
  in-grid. Production NCF~56k has every read on-grid.
- AA on, JNEXT_RELAX=0.5, ALT+OWNW schedule - production
  defaults throughout.

## WP-M2 - why the compressing map does not contract

Named mechanism: the frozen-snapshot (Jacobi) upwind read on
compressing rows gives the composite fixed-point map a spectral
radius ~= the census a-value (~0.96 fixture, ~0.995 production).
Inside `KN_M_MAXPASS=60` the compressing half cannot converge -
the same-direction ratio detector fires correctly, and the
expanding half then converges alone. It is intrinsic to the
lagged map, not an instrumentation abort.

Mechanics, all verified at HEAD:

1. Each refresh pass writes the field into a pass-start snapshot
   `JNEXT_RF_READ`/`HNEXT_RF_READ`
   (mom_j_cmf_nonmono_v1.f:6185-6186).
2. Compressing rows (Q<0) read J at the upwind bin
   `k+MOM_NM_J` from that snapshot (:2836-2856) and the H
   midpoint likewise (:2904-2921) - a Jacobi step. Live-GS reads
   apply only on descending passes when `MOM_RF_GS_DESC` or
   `MOM_DOM_CMP` is set (:2838,:2906); production runs neither
   under ALT.
3. Under `MOM_RF_OWNW` (default T), a descending pass writes the
   compressing rows but suppresses expanding ones, and vice
   versa (:4051-4063,:4259-4264). The compressing half therefore
   takes one damped Jacobi step per two-pass cycle; the write is
   `JNEXT_RELAX=0.5`-under-relaxed (:4280).
4. The captured row census a = |NU_CUP|/NU_TB (:3610-3613,
   `NUC_AMAX`) is the upwind-coefficient/diagonal ratio: 0.964
   fixture, 9.947E-01 production. The Jacobi map's spectral
   radius is bounded below by ~amax, so contraction per
   compressing step is ~1-a at best.

Budget arithmetic: a damped Jacobi iterate at rho ~ (1-w)+w*a =
0.5+0.5*0.964 ~ 0.98 needs ln(2E-3)/ln(0.98) ~ 308 passes to
reach `KN_M_TOL=2E-3`; even undamped (rho=a) needs ~150.
`KN_M_MAXPASS=60` makes non-convergence unconditional. At the
production a=0.995 the count exceeds ~1200.

Measured drift trend (the checklist's monotone-vs-oscillating
question): oscillating, not monotone-slow. Same-direction
compressing drifts and ratios in the winning run:

```
pass:   2     4     6     8     10    12    14    16    18
chg:  10.3  67.9  407   37.3  68.5  18.2  131   10.4  29.4
pass:  20    22    24    26    28    30    32    34    36
chg:  10.5  53.8  36.4  2.50  29.9  17.6  8.14  9.10  6.40
pass:  38    40    42
chg:   2.28  3.74  16.7   -> ratios 1.64, 4.45 -> stall at 42
```

The expanding half decays geometrically in the same run
(161 -> 1.5E-3 over 27 odd passes, ~0.7/step) and exits
`converged`. The compressing oscillation is not a per-cell
period-2 cycle - the signed same-parity LAMC check (:6636-6658)
halved `MOM_DAMP_NM` on zero cells all run (`damped 0` every
pass); the direction-level max roams across the marginal tail
instead (7 cells of 131874 at a>0.9). Complex spectrum of the
frequency-shift coupling plus a roaming argmax is enough to
cross the 0.9 bound twice.

Isolation arms on the same geometry (driver arg20):

- `DESC` (all-descending, `MOM_DOM_CMP` -> live-GS compressing
  reads): converged in 39 passes, no stall. The same compressing
  rows contract once the upwind read propagates in-pass - the
  freeze, not the row equations, is the non-contractive
  ingredient.
- `ASC` (all-ascending): `converged` but vacuous on the
  compressing side - under OWNW every pass is expanding-owned,
  so compressing rows are never written.

H channel ("unreconstructed H" from the diagnosis): compressing
midpoints read `HNEXT_RF_READ` (:2904-2921) - one-pass snapshot
lag, identical semantics to J. No separate H lag defect; the
H-lag is the same Jacobi mechanism, not lagged a full refresh.

Secondary instrumentation finding (WP-M5 guard list, not the
production defect): the same-direction stall test (:7083-7090)
uses pure ratios `PASS_CHG/RES_PREV` with no absolute drift
floor. In the rfj/ALB=0.6 variant and the rlx=1.0 arm the
expanding direction stalls on a ~1E-2-1E-1 wobble plateau -
effectively converged but flagged cycling. Production never
showed it (expanding x0 over 1221 events), so the production
expanding map decays through 60 passes without plateauing; the
defect is real but latent.

## Verdict

WP-M1 VERIFIED: the fixture fires "Compressing direction stalled
- dropping it" on HEAD at sub-unity albedo while the expanding
branch converges the refresh; 149/149 checks pass in
test_mom_solve.py.

WP-M2 VERIFIED: non-contraction is intrinsic to the
frozen-snapshot Jacobi compressing map at spectral radius ~
census amax (0.964 fixture / 0.995 production) vs the 60-pass
budget - an unconditional non-convergence, correctly detected.
Not an instrumentation abort, not a separate H-lag defect.
Detector has a latent no-floor weakness recorded for WP-M5.
