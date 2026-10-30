# KNOP arbitrary-velocity solver audit - 30-Oct-2026

Branch: coupled-sob-jfnk @ 3c9182b1.  Scope: plane/fg_j_cmf_nonmono_v1.f,
new_main/subs/recompute_sigma_true.f, new_main/cmfgen_sub.f,
mom_j_cmf_nonmono_v1.f (parallel machinery noted, not reworked).

## Architecture map

Per great iteration (cmfgen_sub.f DO 10000 ML=1,NCF):
- comp_j_blank.f:843 calls SET_FG_NU_NEXT(NU(ML+1)) then
  FG_J_CMF_NONMONO_V1 per frequency (DO_NONMONO_KNOP or SOB_FG_KNOP).
- In-sweep: ray geometry + GAM=V/R*(1+sigma*mu^2) built once per grid;
  per frequency chi_hat=chi+|GAM|/DLN, Q=GAM/(DLN*chi_hat),
  S'=eta/chi_hat+|Q|*I_upwind; upwind = I_PREV (GAM>=0) or the lagged
  I_*_DEPTH_NEXT store (GAM<0, slot s=freq s+1).  Relaxed writes fill
  slots 1..NCF-1; slot NCF is a zero phantom (red-edge free outflow).
- Post-sweep: FG_KNOP_REFRESH_V1 (cmfgen_sub.f:3276) re-solves each
  ray's full (space,freq) field via KNOP_RAY_SOLVE_V1
  (SWEEP_FLD alternating or SWEEP_FLD_ORD ordered march), then commits
  to the same store (WRITE_STORE / WRITE_STORE_CAP).
- SOLVE_FOR_POPS reads GET_FG_GS_MATURITY for the release gate;
  DUMP/RESTORE_FG_INEXT_GS checkpoint at every SCR_RITE_V2.

Discrete operator (verified identical between sweep and refresh):
inward I_M(i+1)=I_M(i)*A0+S'(i)A1+S'(i+1)A2+dS(i)A3+dS(i+1)A4;
outward I_P(i)=I_P(i+1)*A0+S'(i+1)A1+S'(i)A2-dS(i+1)A3-dS(i)A4,
with the same chi_hat, signed->|Q| upwind source, monotonic-parabola
dS clamps, and boundaries (I_M(1)=0; I_P(NI)=I_CORE/I_M(NI)).

## Confirmed bugs (fixed this commit)

1. FG_GS_NSTAL/NPLAT were set from NX_BAD BEFORE the exit-census loop
   filled it (fg ~6256 vs ~6292): the stall/plateau share the release
   gate reads was always 0.  Moved after the census.
2. Refresh NI=2 endpoint slope: sweep sets dS(NI)=S(NI-1) (constant
   slope, e2516057); refresh SLOPES forced DS_L(NI)=0 -> different
   operator on the thinnest insert-free rays.  Aligned to the sweep.
3. DTAU min/max clamp applied unconditionally in refresh BUILD_FREQ;
   the in-sweep applies it only outside INSERT/MON modes -> the
   refresh solved a different operator for INT/INS and MON runs.
   KN_INSERT_SAV now captured; clamp mirrors the sweep rule.
4. SIGMA_TRUE_CORE derivative: two-point secant is 1st-order on
   non-uniform grids (assigned stencil-midpoint slope to R(I)) and
   1st-order one-sided at boundaries.  New scheme: sigma =
   d ln|V|/d lnR - 1 via the non-uniform quadratic in lnR where the
   stencil carries one sign and |V| spans <1 decade -- EXACT for any
   power law V ~ R^b on any grid (homologous SN ejecta gives sigma=0
   to roundoff; the log form equals (R/V)dV/dR-1 for either sign of
   V).  Sign-changing or near-crossing stencils (ln|V| kinked) fall
   back to the linear-space quadratic dV/dR.  Absolute 1E-10
   R-spacing guard replaced by relative RTOL*R.  Trade-off, measured:
   within ~3 cells of a sign crossing the log-space cells carry
   O(10%) magnitude error vs truth (classification sign always right);
   a linear-secant-only scheme had the same order of error there.
5. SIGMA_TRUE_CORE stagnation branch substituted T1=MIN|V(neighbour)|
   for V(i), discarding sign(V(i)): at V<0 near-stagnation cells the
   sigma sign flipped vs the algebraic limit, misclassifying the
   cell's compressing direction.  The signed V(i) is now carried in
   the algebraic form; exact V(i)=0 reports N_UNDEF instead of a
   fabricated classification (the sigma parametrization cannot
   express the finite GAM limit ~mu^2 dV/dR/c).
6. SET_FG_GS_TAG was gated on DO_NONMONO_KNOP only, so under
   DO_SOB_FG_KNOP the FG checkpoint always failed its tag test and
   cold-started I_NEXT every requeue.

## Findings - documented, not changed this commit

- Freq-1 refresh field has no store slot (store slot s=freq s+1 ->
  slot 0 would be needed).  Verified harmless: no in-sweep or refresh
  reader ever needs freq-1 lagged memory (compressing upwind is
  always freq k+NM>=2; expanding upwind is FLD(:,L-1) same-iterate).
  The red-edge phantom slot NCF_IN is the correct free-outflow BC.
- KN_NOSLOT fallback (GAM<0, no compressing depth in bracket): keeps
  chi_hat absorption, drops |Q|I_next source - physically the correct
  "no stored upwind" closure; counted per ray in FG_NOSLOT.
- KN_CHI_EXT/KN_ETA_EXT captured in single precision -> refresh
  coefficients carry ~1E-7 relative error vs the sweep's double
  evaluation (documented memory tradeoff, ~300MB).
- FG_JUMP_MAX reads sibling rays' store columns while other OMP
  threads write them -> benign torn reads (values always finite/
  non-negative); deterministic field, nondeterministic guard inputs.
- NCOMP/NMATCH/NMAX telemetry only fills under HAVE_PACK (pack
  fallback reports 0) - diagnostic gap only.
- AA history (2*(M+1)*NI*NCF), Aitken (+16B/cell) and P2DAMP (+3B)
  memory is outside the ~80B/cell KN_RF_TOT_MB estimate.
- FG_GS_ARG* argmax tuple records the sweep argmax while
  FG_GS_MAXCHG reports the refresh residual post-refresh.
- Boundary-closure derivation (Phase 5): blue edge keeps the V13
  INIT convention (expanding cells have no blue inflow AND no
  outflow term at L=1); compressing cells keep full chi_hat+|Q|
  coupling - consistent between sweep and refresh.
- SWEEP_FLD_ORD's KN_ORD_SPV recomputes the match offset inline
  rather than reusing MATCH_FREQ's NFM_L (currently identical
  formula; drift risk only).

## Residual/convergence notes

The refresh residual RES=max|VNEW-VOLD|/MAX(VOLD,1e-4*col peak) IS
the unrelaxed fixed-point residual |g(x)-x| of the discrete transfer
operator (VNEW is the raw map output before OMEGA_RAY relaxation) -
it is the transfer-equation residual, not just an update magnitude.
MAX norm retained; bulk exit additionally gated on contraction
(FG_BULK_CONTR) and residual floor (KN_RAY_BULK_TOL).

## Validation state

ifort -fsyntax-only clean on both edited production files (mod_dbg
holds a stale GS_FOREIGN_SEED .mod; production mod/ rebuilt by
build.sh).  Model-run validation not performed in this session -
compiling is not validated.

Unit tests added/passing:
- test_sigma_true.py: T3 reworked to true-analytic reference on a
  log-spaced grid; T3b (independent hybrid-stencil mirror, exact),
  T3c (exponential profile), T3d (zero crossing: mirror-exact,
  correct sign both sides, bounded near-zone rel err ~12%).  16/16.
- test_rv_hydro_rdin.py: reference stencil updated to the hybrid
  formula; 9/9.
- test_knop_ray.py + f90/drv_knop_ray.f90: FIRST runtime coverage of
  KNOP_RAY_SOLVE_V1 on a manufactured core ray (NI=8, NCF=6, one
  compressing band).  T1 alternating sweep converges (NXIT=1,
  resid 6e-7, 3 passes); T2 ordered sweep; T3 ord/alt committed
  stores agree to 2.6e-8 (cross-discretization consistency); T4
  FREQ_MATCH fires multi-bin slots under strong drift (NMATCH=12,
  NMAX=2); T5 forced stall/cap exits visibly (NXIT=4), never as
  converged; T6 near-fixed-point seed converges immediately.  6/6.
- test_vadat_keys 70/70, test_knop_lagged_field (static), 
  test_mom_sweep_table 8/8, test_sob_tau_floor 14/14 -- all pass
  after the source changes.
- test_knop_ray T7: FGINEXT_GS checkpoint round-trip is bit-exact;
  fingerprint mismatch, tag mismatch and foreign-seed paths are all
  rejected/flagged correctly (foreign seed restores the field but
  stays immature, I_NEXT_GS_OK=F).  T8: FG_KNOP_REFRESH_V1
  end-to-end over the OMP ray loop converges on the manufactured
  problem (KN_RF_RESID=6e-7); 1 vs 4 OMP threads bit-identical.


## Zero-GAM floor analysis (mom_j_cmf_nonmono_v1.f, mandatory item)

Construction (mom :2011-2013, :2132-2134):
  T1 = 1 + sigma_true*K   (K = NMID_ON_HMID at midpoints, K_ON_J at depths)
  GAM_FAC = SIGN(MAX(|T1|, MOM_GAM_FAC_MIN=1e-2), T1)
  exact-zero exception: IF(T1==0) GAM_FAC=0

Limit analysis: as T1 -> 0 the physical GAM = CON_GAM*T1 -> 0
continuously; chi_hat -> chi, Q -> 0, W -> 0, EPS -> 1/(1+W) -> 1.
The unregularised operator is SMOOTH through the resonance - a
zero-strength coupling has no direction, so the sign flip is
immaterial.  The floor exists to stabilise the SIGN classification
against roundoff/noise at |T1|~0 (otherwise tiny T1 flips the upwind
direction between iterates), at the cost of an O(1e-2) systematic
bias in |GAM_FAC| on |T1|<1e-2 cells: chi_hat gains a spurious
|GAM|/DLN ~ 1% of the nominal coupling, and W/Q carry a 1%
residual coupling the true equations set to zero.

Assessment: the floor is a defensible regularisation ONLY while the
biased |GAM_FAC|*CON_GAM stays small against chi_hat - true for the
1e-2 floor.  Two structural defects worth recording:
(a) discontinuity: GAM_FAC jumps -1e-2 -> +1e-2 across T1=0 while
    physics passes smoothly through zero; the matrix operator is
    non-continuous in T1 at every resonance crossing (a W_NM jump
    of |GAMH|*2e-2).
(b) exact-zero exception reintroduces the very zero the floor
    exists to remove: GAM_FAC(0)=0 between +/-1e-2 neighbours -
    a measure-zero spike that re-creates the chi_hat=chi, W=0,
    EPS=1 stock-decoupled state at exactly the resonance cell.
    Harmless in probability (T1==0 exactly is rare in production,
    though it DOES occur on structured profiles, e.g. sigma=-1 with
    K=1), but it means the "resonance regularised" claim is not
    uniform: at exact resonance the solve silently reverts to the
    unregularised stock equation, and the two sides of the jump
    disagree by O(1e-2*CON_GAM/DLN) in chi_hat.
Recommendation (not implemented - needs controlled tests): replace
sign*max(|T1|,floor) with a smooth bounded-below-magnitude map,
e.g. T1 + SIGN(T1)*floor for |T1|<=floor (a 2*floor offset, still
discontinuous in the same place), or better GAM_FAC =
T1*(1+floor^2/(T1^2+floor^2))^0.5 which is C1, sign-preserving,
|T1|*~1 for |T1|>>floor and +-sqrt(2)*floor at T1=0 - removes both
the jump and the need for the exact-zero exception.  Requires the
mom runtime driver (pending) to verify identical convergence class
before adoption.
