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
   1st-order one-sided at boundaries -> nonuniform quadratic (2nd
   order) interior + boundary.  Absolute 1E-10 R-spacing guard
   replaced by relative RTOL*R.
5. SIGMA_TRUE_CORE stagnation branch substituted T1=MIN|V(neighbour)|
   for V(i), discarding sign(V(i)): at V<0 near-stagnation cells the
   sigma sign flipped vs the algebraic limit, misclassifying the
   cell's compressing direction.  Now evaluates (R/V)dV/dR-1 with the
   signed V(i); exact V(i)=0 reports N_UNDEF instead of a fabricated
   classification.
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
