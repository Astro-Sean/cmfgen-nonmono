# Phase 5.3 - Expanding-Only P-Direction Convergence Audit

Status: AUDIT COMPLETE - source inspection and mathematical reconstruction
only. No solver, relaxation, commit-policy, jump-veto, or population-solver
changes were made. Findings are evidence for independent review before any
implementation is authorised.

Date: 2026-11-XX (post Phase-5.2 commit `5cf564e0`)
Scope: `plane/fg_j_cmf_nonmono_v1.f` at source revision `5cf564e0`
(= tested executable `fee96b6a2a7c02fc176b0bf8ed9d0285` plus the
format-descriptor fix, which is inside a probe-gated report block and
does not alter solver arithmetic).

Subject: ray LS=83, channel P (outward-directed intensity), expanding
rows (`GAM(J,83) >= 0` for all marched rows of interest; `MATCH_FREQ`
is never invoked on this ray). Phase 5.2 measured a persistent
oscillatory residual regime (class D on a class-B plateau) with
`OMEGA_RAY` pinned at its 0.45 floor, residual argmax at ray points
j in {1,2} over a broad frequency range, while the M sweep converged
in two passes.

This document reports what the code actually does, traced from array
reads and writes - not from routine names or comments.

------------------------------------------------------------------------
1. SOURCE CALL GRAPH
------------------------------------------------------------------------

Per refresh, per ray LS (inside the `!$OMP PARALLEL DO` ray loop,
`fg_j_cmf_nonmono_v1.f:4189`):

    KNOP_RAY_SOLVE_V1(LS, ...)
      |-- builds ray grid quantities once per ray:
      |     RAY_PNT/R_RAY/Z (lines ~3031-3115), J_PNT (model depth ->
      |     ray point), GAM, NFM_L, DLNM_L, JPR1/JPR2/IKC1/IKC2/WINX,
      |     DRF_L, KN_DLN_P/KN_DLN_N (module-level, checkpointed).
      |-- per frequency L = 1..NCF: fills DTAU_L, CHI_L, ETA_L,
      |     DCHR_L, Q_L, SR_L, EEL/E0L-E3L, A0L-A4L (coefficient
      |     construction ~lines 11900-12150); snapshots them into the
      |     per-(ray,frequency) PK_* arrays (lines 9575-9584).
      |-- allocates FLD(NI,NCF) and FLD_P(NI,NCF), SPW(NCF,4),
      |     IM_NI(NCF), RES_FREQ(NCF), VMAXF(NCF), P2_* state.
      |
      |-- M (inward) solve:
      |     SEED_FLD(FLD, I_M_DEPTH_NEXT)         ! retained store
      |     SWEEP_FLD_ORD(FLD, INWARD=.TRUE., ...)  ! default path
      |       or SWEEP_FLD(...) if KN_RAY_ORDERED=F
      |     on DIVERGED: deallocate, skip ray (stores untouched)
      |     IM_NI(1:NCF) = FLD(NI,1:NCF)            ! line 9665
      |
      |-- P (outward) solve:
      |     FLD_P = 0; SEED_FLD(FLD_P, I_P_DEPTH_NEXT)  ! line 9693-4
      |     SWEEP_FLD_ORD(FLD_P, INWARD=.FALSE., ...)  ! line 9697
      |       or SWEEP_FLD(...) if unordered
      |     commit path: WRITE_STORE / WRITE_STORE_CAP /
      |       FG_RAY_SEEDED / HEAL_STORE (all probe-guarded in 5.2)
      |
      `-- deallocates all per-ray work arrays (line 10176-10190).

Within SWEEP_FLD_ORD (line 10822), one "pass" = one full sweep over
ray points and frequencies:

    pass start:  VMAXF(L) = MAX(MAXVAL(FLD(:,L)), 1e-25)   per L
    boundary:    row j=NI relaxed toward I_CORE/IM_NI (P) or
                 row j=1 set (M)                            lines 11177-11256
    interior:    KN_ORD_SPW_FILL(FLD, NI)                  line 11257
                 DO J = NI-1, 1, -1                        line 11258
                   KN_ORD_SPW_FILL(FLD, J-1)  (preview)   line 11259
                   freq order: L=1..NCF if GAM(J,LS)>=0
                               L=NCF..1 otherwise          lines 11260-4
                   per (J,L):
                     SPJ = KN_ORD_SPV(FLD,J,L)            line 11267
                     SPW(L, MOD(J-1,4)+1) = SPJ            line 11268
                     DSU = KN_ORD_DS(J+1,L)               line 11278
                     DSD = KN_ORD_DS(J,L)                 line 11279
                     VNEW = FLD(J+1,L)*PK_A0(J,L)
                          + SPW(L,MOD(J,4)+1)*PK_A1(J,L)
                          + SPJ*PK_A2(J,L)
                          - DSU*PK_A3(J,L) - DSD*PK_A4(J,L)
                                                          lines 11280-3
                     floor/cap: VNEW -> [1e-30, CAP]       lines 11284-93
                     residual: T1 = |VNEW-VOLD| /
                                MAX(VOLD, 1e-4*VMAXF(L))   line 11295
                     census:   NBAD_CUR/NTOT_CUR, RES_FREQ,
                                KN_SGN_SUM, P2 tracking      lines 11304-42
                     commit:   FLD(J,L) = VOLD +
                                OMEGA_RAY*(VNEW-VOLD)       line 11356
                 (P2 geometric-mean branch exists but
                  FG_ORD_P2DAMP defaults .FALSE., line 1137)
    pass end:    NB_O/NT_O census, RATIO_OM = RES/RES_PREV,
                 OMEGA_RAY adaptation, exit classification
                 (IXIT)                                     lines 11423-11622

Helper routines called inside the P sweep:

    KN_ORD_SPV(FLD,R,L,SPV)      line 10675 - S'(R,L) evaluation
    KN_ORD_SPW_FILL(FLD,R)       line 10728 - fills SPW ring slot
    KN_ORD_SEC(R,L)              line 10742 - SPW secant slope
    KN_ORD_DS(R,L)               line 10759 - monotonic-parabola DS
    KN_Q_ACCUM(T1,INWARD)        - probe diagnostic accumulator
    NOP_CELL(T1)                 - probe per-cell telemetry (5.2)

MARCH (line 12372) is the legacy unordered comparator: same PK_*
coefficients, frequency loop without directional ordering, used only
when KN_RAY_ORDERED=F. Ray 83 runs the ordered path.

------------------------------------------------------------------------
2. DISCRETE EQUATION ACTUALLY SOLVED (P direction, expanding row)
------------------------------------------------------------------------

Unknown updated in place: FLD(J,L) for J = NI..1 (boundary then
interior), L in the row's frequency order. One pass computes, for
each interior point (J,L), J = NI-1..1:

    S'(R,L) =
        PK_SR(R,L) + PK_Q(R,L) * FLD(R,L-1)        GAM(R,LS) > 0, L > 1
        PK_SR(R,L)                                  GAM > 0, L = 1
        PK_SR(R,L)                                  GAM = 0 (exact)
        PK_SR(R,L) + |PK_Q(R,L)| * INEXT_V          GAM < 0 (not on ray 83)

    VNEW(J,L) =
        FLD(J+1,L) * PK_A0(J,L)                     upwind intensity
      + SPW_up(J+1 -> slot read as MOD(J,4)+1) * PK_A1(J,L)
      + SPJ * PK_A2(J,L)                            fresh S'(J,L)
      - DSU * PK_A3(J,L) - DSD * PK_A4(J,L)         parabola corrections

    FLD(J,L) <- VOLD + OMEGA_RAY * (VNEW - VOLD)

where, for the P march (J decreasing = outward propagation):

- FLD(J+1,L): upwind (inner-side) neighbour intensity. Row J+1 was
  updated earlier in the SAME pass -> FRESH (Gauss-Seidel).
- SPW(L,MOD(J,4)+1): the ring slot that row J+1's own L-loop wrote
  with SPJ = S'(J+1,L) during this pass -> FRESH S' at the upwind
  node. (See section 5 for the ring mechanics; the indexing is a
  deliberate shift register, not an off-by-one.)
- SPJ = KN_ORD_SPV(FLD,J,L): S'(J,L) evaluated on the current row
  with frequency donor FLD(J,L-1) already updated this pass for
  ascending L -> FRESH in frequency (Gauss-Seidel in L).
- DSU/DSD: monotonic-parabola slopes from SPW secants at rows
  J+2..J-1; all contributing slots were written by fresh SPJ values
  earlier in this pass (or by the boundary/preview fills at the
  leading edge) -> consistent within the pass.
- Boundary row j=NI: VNEW = IM_NI(L) for non-core rays
  (LS > KN_NC; ray 83 qualifies), relaxed once per pass:
  FLD(NI,L) <- VOLD + OMEGA_RAY*(IM_NI(L) - VOLD).
- Outer edge j=1: no boundary condition is imposed; the last marched
  point is computed by the same interior formula (outflow edge).

Coefficients (per ray, per frequency, filled once per refresh into
PK_*):

    DTAU(J,L): segment optical depth, Euler-Maclaurin form via
               KN_DTAU_INT_V1, |GAM|/DLN folded into CHI_L for
               compressing cells only; for expanding cells CHI_L is
               the plain opacity. Clamped non-negative (lines
               12076-12098).
    EEL = exp(-DTAU); E0L..E3L: exponential-integral moments
               (lines 12100-12127).
    A0 = EEL
    A1 = E0L - 3*E2L + 2*E3L          (upwind-node source weight)
    A2 = 3*E2L - 2*E3L                (downwind-node source weight)
    A3 = DTAU*(E1L - 2*E2L + E3L)
    A4 = DTAU*(E3L - E2L)
    (large-DTAU branch, T1>=40, uses the 1/T1 asymptotics;
     lines 12129-12144)
    Q_L(J)  = GAM(J,LS) / (DLNM_L(J) * CHI_L(J))   line 12055
              gated by KN_RF_COUPLING(J,L); for GAM >= 0 the gate
              always returns .TRUE. (line 12207-12208).
              DLNM_L is the per-row log-frequency span; for
              expanding cells it equals KN_DLN_P(L) =
              ln(nu_{L-1}/nu_L) - the backward span matching the
              L-1 donor (lines 12215-12218 comment, 3843).
    SR_L(J) = ETA_L(J) / CHI_L(J)                    line 12059

Sanity checks on the coefficient block:

- Q_L -> 0 continuously as GAM -> 0+ provided CHI_L stays finite
  and positive; the GAM == 0 branch returns SR only, which equals
  the GAM -> 0+ limit (Q -> 0). No discontinuity in S'. The
  frequency ORDER, however, switches at GAM = 0 (ascending vs
  descending) - a per-row ordering switch, not a coefficient
  discontinuity. Rows with |GAM| tiny but nonzero still pick a
  direction; the two directions share the same fixed point but
  different transient paths.
- For expanding cells KN_RF_COUPLING is unconditionally .TRUE., so
  the coupling is never masked off for ray 83.
- A-coefficients are pure functions of DTAU and the E-moments;
  A0 in (0,1], A1,A2 >= 0 for the standard exponential-integral
  solution (formal-solution weights). A3,A4 carry the parabola
  correction signs used above.

------------------------------------------------------------------------
3. ARRAY DEPENDENCY MAP
------------------------------------------------------------------------

Read/write sets inside one P pass (FLD = FLD_P):

    FLD(J+1,L)      read  (upwind intensity)          fresh this pass
    FLD(J,L-1)      read  via KN_ORD_SPV (GAM>0)      fresh this pass
                    (L-1 precedes L in ascending order)
    FLD(J,L)        read+write (VOLD then damped)     in place
    FLD(J-1,L')     read  via SPW preview/KN_ORD_DS   lagged (see 5)
    FLD(*,L+NNF2)   read  only for GAM<0 rows           n/a ray 83
    IM_NI(L)        read  at boundary j=NI             fixed snapshot
    PK_A0..A4,PK_SR,PK_Q,PK_DTAU  read-only coefficients (per refresh)
    SPW(NCF,4)      read+write ring, per-ray private allocatable
    GAM,NFM_L,DLNM_L,JPR*,IKC*,WINX,DRF_L  read-only ray geometry
    RES_FREQ,VMAXF,RES,NBAD_CUR,NTOT_CUR,NNEG,NCAP,
    KN_SGN_SUM,P2_*,WRC_*          per-pass accumulators

Lifetime/scope:

- SPW, FLD_P, IM_NI, RES_FREQ, VMAXF, P2_* are LOCAL allocatables of
  KNOP_RAY_SOLVE_V1 (declarations ~lines 9300-9360, allocation
  ~9598, deallocation ~10176-10190). Each OMP ray thread has its own
  invocation -> private instances. No cross-ray sharing, no race on
  these arrays.
- Module-level SAVEs touched per ray are either per-ray indexed
  (KN_Q_LS(LS), KN_PA_*(LS,:), KR_DBG_*(...,LS)) or updated under
  `!$OMP CRITICAL` (KN_HOT_CNT, line 9671-9675). KN_NOP_* probe
  telemetry is written only on the probe ray (L_NOPROBE guard).
- IM_NI is a per-call local array filled once at line 9665 from the
  M solution's FLD(NI,:); never written during the P sweep. No
  aliasing with FLD_P (separate allocation).

------------------------------------------------------------------------
4. BOUNDARY-CONDITION ANALYSIS
------------------------------------------------------------------------
4.1 Ray-point orientation (correction to Phase-5.2 wording)

Phase 5.2 described the residual concentration as "the inner
boundary (j=1-2)". Traced geometry says otherwise:

- Line 4144 comment: "J_PNT(1,LS)=1 is the common case -- R(1)
  sits on the outermost ray". Model depth index 1 is the OUTERMOST
  atmospheric layer; depth index increases inward.
- J_PNT(I,LS) brackets model depth I between ray points K,K+1 by
  radius; for ray 83, J_PNT(2,83)=2 (KN-WATCH geom output), so ray
  point j=2 carries model depth K=2 - the second-outermost layer.
- Ray points are ordered j=1 (outermost, largest R_RAY) to
  j=NI (innermost, smallest R_RAY = core side).
- M (inward-propagating) marches J=2->NI with upwind FLD(J-1):
  information flows outer->inner. Its Dirichlet boundary is at
  j=1 (outer edge, zero incident intensity) - consistent.
- P (outward-propagating) marches J=NI-1->1 with upwind
  FLD(J+1): information flows inner->outer. Its boundary is at
  j=NI (inner edge) supplied by IM_NI/I_CORE - consistent.
- CONCLUSION: the P residual argmax at j in {1,2} is at the OUTER
  boundary of the ray - the escape edge - not the inner/core
  boundary. Depth K~1-2 is the top of the atmosphere along this
  ray. This matters for interpretation: the oscillation peaks
  where the outward intensity leaves the domain, i.e. where the
  accumulated frequency-coupling error along the whole column is
  expressed, or where the outer-edge discretisation is weakest.

4.2 Index-map distinctions (as required)

- RAY_PNT(K,LS): interpolation interval on the EXTended radial
  grid used to build ray point K (built ~lines 3031-3115 from
  R_EXT/Z_EXT and inserted points).
- J_PNT(I,LS): maps MODEL depth I to the ray point whose radius
  brackets R(I); used to place J/H moments on ray nodes.
- COMP_DEPTH_IDX(K): compressed/depth indexing for the moment
  solver; not used inside the ray sweep.
- KN_WATCH_K=2 -> J_PNT(2,83)=2 -> ray point j=2 -> model depth
  K=2 (outer layers). The earlier "inner boundary" phrasing
  conflated small ray index with inner radius; j=2 is outer.

4.3 IM_NI trace (Part 2 of the task)

- Creation: IM_NI(1:NCF)=FLD(NI,1:NCF) at line 9665, executed
  AFTER the M sweep returns and BEFORE FLD_P is allocated. It is
  a snapshot of the M channel's inner-boundary intensity, i.e.
  I_M at the core-side ray endpoint, all NCF frequencies.
- Meaning: for non-core rays the physical inner boundary
  condition is I_P(inner) = I_M(inner) - the two directed
  intensities meet at the thermalizing inner boundary
  (ZERO_FLUX/PNT_SRCE family; see lines 11184-11191 for the
  core-ray variants DIFFUSION/ZERO_FLUX/PNT_SRCE/KN_IC).
- Frequency indexing: IM_NI(L) is indexed by the same L as
  FLD(:,L); no remapping, no interpolation - direct copy of the
  M row at j=NI.
- Sign/scaling: none applied; raw intensity units, same as FLD.
- Updated during P sweep? NO - read-only inside SWEEP_FLD_ORD.
  It is fixed for the whole P solve (all passes).
- Aliasing/stale: none - separate local array, filled once per
  ray per refresh.
- Is M genuinely converged before P consumes it? M exits with
  IXIT_M and RES_M; on the probe ray M converged (IXIT=1,
  res~4e-8, np=2) so IM_NI is the M fixed point to sweep
  tolerance. If M had diverged (IXIT=5) the code skips the ray
  entirely (line 9654-9663). If M stalls/caps without diverging,
  IM_NI is whatever M's last damped pass left - formally
  "converged enough" only if M's exit was clean. For ray 83 this
  holds.
- j=NI vs model depth: j=NI is the innermost ray point; the
  nearest model depth is K=ND (innermost). The boundary value is
  applied at the ray endpoint, not interpolated to a model depth -
  correct for a ray boundary condition.

4.4 Outer edge j=1

P imposes no boundary condition at j=1; VNEW(1,L) is computed by
the interior formula from upwind j=2 plus sources. This is the
correct outflow treatment for outward intensity (no incoming
stream at the outer edge). M's j=1 IS a boundary (set toward 0 /
incident intensity) because inward intensity enters there. The
asymmetry is physical, not a defect.

------------------------------------------------------------------------
5. FREQUENCY-MARCH ANALYSIS (expanding-only operator)
------------------------------------------------------------------------
5.1 Discretisation actually used (GAM > 0)

For a row with GAM(J,LS) > 0 the frequency loop runs L = 1..NCF
ascending. At each (J,L):

    S'(J,L) = SR(J,L) + Q(J,L) * FLD(J,L-1),   L > 1
    S'(J,L) = SR(J,L),                          L = 1
    Q(J,L)  = GAM(J,LS) / (DLNM_L(J) * CHI_L(J))

So the frequency coupling is a ONE-SIDED upwind (in frequency)
term: intensity at the LOWER frequency bin L-1 of the SAME ray
point feeds the source at L. The donor FLD(J,L-1) was updated
earlier in the same row's ascending loop -> fresh Gauss-Seidel
value including its own relaxation.

Edge closures:

- L = 1 (lowest frequency bin): no coupling term at all - the
  source is SR only. This is an outflow/zero-inflow closure in
  frequency space: photons cannot arrive at L=1 from below the
  grid. Whether this matches Knop et al. (2009) Paper II's
  boundary treatment needs checking against the paper (see
  references); it is a modelling choice, not obviously wrong.
- L = NCF: donor L-1 = NCF-1 exists; no special closure needed
  at the top for expanding rows. The coupling is unidirectional
  (low->high L in update order), so there is no downwind
  reflection within a pass.
- GAM < 0 rows (not on ray 83's argmax region but possibly
  elsewhere on the ray): descending L order, donor at
  L+NNF2 with |Q| and optional MATCH_FREQ remapping through
  JPR1/JPR2/IKC1/IKC2/WINX and DRF_L/KN_DLN_N. NNF2 =
  CEILING(DRF_L/KN_DLN_N) clamped to [1, NCF-L]; out-of-range
  donors degrade to SR-only (line 10715-10716). KN_RF_COUPLING
  masks the coupling off for compressing cells whose matched
  frequency would leave the grid (lines 12205-12211) - the
  "implicit outflow must be dropped" rule documented at
  12195-12203.

5.2 Consistency with Knop et al. (2009)

The structure S' = SR + Q*I_upwind with Q proportional to the
signed velocity-gradient term GAM/(DLN*chi) is the comoving-frame
frequency-coupling (advection) term of the Knop formal solution:
in an expanding flow a photon redshifts as it propagates, so the
intensity at frequency L is sourced by the intensity the photon
had at a lower comoving frequency upstream. The upwind choice
(L-1 for expanding, matched-higher for compressing) is the
correct characteristic direction IF the sign convention
(GAM>0 = expansion) matches the paper. This audit did NOT
re-derive Paper II's equations; the specific claim to verify
externally is whether the L=1 closure and the |Q| vs signed-Q
asymmetry (5.3) match the published discretisation.

5.3 Signed vs absolute Q (asymmetry worth flagging)

- Expanding branch (line 10690): SPV = SR + PK_Q*FLD(R,L-1)
  uses SIGNED PK_Q. For expanding cells Q = GAM/(DLN*CHI);
  with GAM>0, DLN>0 this is positive iff CHI_L(J) > 0. If CHI_L
  can go negative anywhere (non-monotonic opacity construction),
  the coupling term changes sign - physical only if the
  formulation intends signed coupling.
- Compressing branch (line 10714): SPV = SR + ABS(PK_Q)*INEXT_V
  forces the coupling positive regardless of Q's sign.
- For ray 83 (GAM>=0) only the signed branch executes. If
  CHI_L>0 throughout (normal), signed Q > 0 and the
  asymmetry is moot; if CHI_L<0 occurs, expanding cells get a
  negative coupling while compressing cells would get positive -
  an inconsistency candidate. Minimal test: dump MINVAL(CHI_L)
  per ray; if negative values exist on expanding rows, this is
  a real defect candidate.

5.4 GAM -> 0 behaviour

- Q -> 0 smoothly (GAM in numerator); S' -> SR continuously.
- The exact-zero branch (line 10678) is an IEEE NaN guard
  (0*NaN = NaN), documented in the comment - correct.
- The frequency ORDER switches discontinuously at GAM = 0
  (line 11260: >=0 ascending, <0 descending). For a ray whose
  GAM changes sign along J, adjacent rows use opposite
  frequency orderings; each row's GS sweep is still
  self-consistent, but the pass couples rows with different
  lag structure. For ray 83 the task states GAM>=0; the
  KN-GEO watch cell showed gam_t = 1.09646E-09 (tiny positive)
  - consistent with an expanding row near the zero crossing.
  Rows with GAM ~ 1e-9 have Q ~ 1e-9/(DLN*CHI) ~ negligible
  coupling but still take the ascending branch - benign.
- Conditioning: Q can become large where CHI_L is small
  (optically thin, high-velocity-gradient cells). Large Q makes
  the L->L+1 coupling strong; combined with the radial GS this
  is the most plausible route to a non-contractive iteration
  (section 7).

------------------------------------------------------------------------
6. ITERATION-MAP ANALYSIS (Part 4)
------------------------------------------------------------------------
6.1 The map

One P pass computes x_{n+1} = H(x_n) where

    H(x) = x + OMEGA_RAY * (G(x) - x)

G being the undamped lexicographic Gauss-Seidel formal-solution
map (boundary row, then J = NI-1..1, L in row order), and
OMEGA_RAY the scalar relaxation. Fixed points of H coincide
with fixed points of G (for OMEGA_RAY != 0). The residual
reported is

    RES = max_{J,L} |G(x)_{J,L} - x_{J,L}| /
                MAX(x_{J,L}, 1e-4*VMAXF(L))

i.e. the UNDAMPED defect of G evaluated at the current (damped)
iterate - measured BEFORE relaxation (line 11294-11301 compute
T1/RES, line 11356 applies the damped commit). This is the
correct fixed-point residual of G; it is NOT the norm of the
actual step x_{n+1}-x_n (which is OMEGA_RAY times smaller).

6.2 Ordering classification

The sweep is a MIXED Gauss-Seidel iteration:

- Radial-upwind (J+1): fresh within the pass.
- Radial-downwind S' (J+1 node source via SPW shift register):
  fresh within the pass.
- Frequency-upwind donor (J,L-1 for expanding): fresh within
  the pass.
- Frequency-downwind: not read for expanding rows.
- Boundary j=NI: relaxed toward fixed IM_NI (not solved).
- Preview fills (KN_ORD_SPW_FILL at J-1): write lagged values
  into the slot that the CURRENT row's SPJ will overwrite
  per-frequency as the L-loop proceeds; the lagged values are
  consumed only by KN_ORD_DS secants at frequencies/rows the
  march has not reached - the "lagged C-corner" of Knop Paper
  II Eqs. 31/32 per the code comment at 10725-10726. At the
  leading edge these are one-pass-stale; elsewhere fresh.

So: block Gauss-Seidel in (J,L) with a one-cell lagged
extrapolation stencil at the marching front - NOT Jacobi, NOT
pure GS. In-place updates throughout (FLD written immediately).

6.3 OMEGA_RAY trace

- Initialised to KN_RAY_OMEGA at the head of each directional
  sweep (line 10921 for ORD; 10271 for unordered).
- Adapted at pass end from RATIO_OM = RES_n / RES_{n-1}
  (same-direction residual ratio, lines 11457-11497):
    RATIO_OM > 1.05 and NPO >= 4  -> OMEGA_RAY *= 0.7 (floor 0.45)
    RATIO_OM < 0.5                -> OMEGA_RAY *= 1.15 (cap 1.0)
    |RATIO_OM - 1| < 0.05 (x3, FG_PLATEAU_DET=T) -> plateau flag
- Observed probe history: OMEGA_RAY walked 1.0 -> 0.7 -> 0.805
  -> 0.926 -> 0.648 -> ... -> pinned at 0.45 from ~pass 12 on.
  The controller responds to the UNDAMPED residual ratio while
  controlling the DAMPED step - mathematically these are
  different sequences; the code comments (lines 10416-10470,
  9895) show awareness that with OMEGA<1 the committed field
  is not the raw map output.
- Residual-history consistency: RES measures |G(x)-x|; for a
  period-2 orbit of H with step amplitude d, |G(x)-x| = d/OMEGA
  at both orbit points (alternating sign). Observed RES ~ 0.6-0.7
  at OMEGA=0.45 implies an orbit amplitude d ~ 0.3-0.32 in the
  normalised defect norm - internally consistent with a
  period-2 (or low-period) limit cycle of the damped map, not
  with slow convergence (no downward trend over the last 20+
  passes) and not with divergence (bounded, recurrent values
  quantised to a few levels across 70 refreshes).

6.4 Is the measured map/norm pair valid for ratio-as-contraction?

NO - and the code does not claim it is. RATIO_OM is used only
to modulate OMEGA_RAY and to classify stall/plateau exits; it
is documented (Phase 5.2 report) as a damping heuristic, not a
convergence factor. Any analysis treating R(n)/R(n-1) as a
spectral estimate would be invalid: the norm is a max-norm
defect of G, the iteration is damped, and the map is
state-dependent (limiters, ordering switches).

------------------------------------------------------------------------
7. FIXED-POINT FAILURE vs NUMERICAL LIMIT CYCLE (Part 5)
------------------------------------------------------------------------
Observed: stationary oscillatory regime, RES ~ 0.57-0.80 with
excursions to ~3.3, argmax j in {1,2}, broad L, OMEGA at floor,
exit IXIT=3 (stall). Candidate mechanisms, ranked by
plausibility for this configuration:

A. NON-CONTRACTIVE UNDAMPED MAP (highest prior)
   The P sweep solves an IMPLICIT problem: S' contains
   Q*FLD(J,L-1), i.e. the unknown intensity appears in its own
   source. G is therefore a genuine fixed-point map, not a
   direct solve. Where |Q| is large (small CHI_L, large
   |GAM|/DLN - optically thin high-gradient cells), the
   frequency-coupling feedback can give J_G eigenvalues with
   |lambda| > 1 or real lambda < -(2/OMEGA - 1) ~ -3.4 at
   OMEGA=0.45, producing a period-2 instability that non-linear
   saturation (floor/cap/positivity) converts into the observed
   limit cycle. This is the classic ALI motivation
   (Rybicki & Hummer 1991): unaccelerated Lambda-type
   iteration diverges in scattering-dominated regimes; here the
   analogue is unaccelerated frequency-coupled GS.
   M converges because its store was already at the fixed
   point and/or its coupling structure differs; P's retained
   store was NOT at the fixed point (Phase 5.2 store anomaly),
   so P iterates a genuinely non-contractive map from a
   non-fixed-point start.

B. LAGGED-CORNER / STALE-STATE COUPLING (medium)
   The preview fills inject one-pass-old S' values at the
   marching front for KN_ORD_DS secants. If the leading-edge
   rows (j small = outer edge, where the argmax lives) depend
   sensitively on those secants, a two-pass lag can sustain a
   period-2 mode: pass n uses stale S' from n-1, pass n+1 uses
   the n values, etc. The ring mechanics verified in section 5
   show fresh values dominate; the lag is confined to
   not-yet-marched rows. Cannot be excluded without J_G.

C. NONLINEAR LIMITER/DISCONTINUITY (medium-low; P2 is OFF)
   FG_ORD_P2DAMP defaults .FALSE. (line 1137) so the
   geometric-mean limiter is inactive unless the run enabled
   it - the probe VADAT did not set KN_ORD_P2DAMP. Remaining
   nonlinearities: VNEW<0/NaN -> 1e-30 floor (line 11284) and
   CAP=NONMONO_INEXT_CAP (default 1e20, line 11289). If the
   linear map routinely produces negative VNEW near j=1-2
   (cancellation among A-terms in steep-gradient outer cells),
   the floor makes G nonsmooth and can support a relaxation
   oscillation. Checkable: NNG_P/NCAP counters per pass (not
   currently emitted per-pass in the probe telemetry - a gap).

D. BOUNDARY INCONSISTENCY AT j=NI (low for argmax location)
   P's inner boundary relaxes toward IM_NI = M's converged
   value. If the true P fixed point requires I_P(NI) !=
   I_M(NI) (the equality assumption is a modelling choice),
   the boundary injects a persistent defect that propagates
   outward and could peak at j=1. But the argmax is at the
   OUTER edge, and the boundary row itself (j=NI) is not the
   argmax - the defect would have to amplify along the march.
   Possible but not the simplest explanation.

E. FREQUENCY-BOUNDARY REFLECTION (low)
   L=1 closure drops the coupling term entirely; L=NCF has no
   special treatment for expanding rows. A mismatch between
   the interior coupling and the edge closure could reflect
   error back into the band. The argmax L values (1230-1660,
   excursions to 27683) are interior, not at the grid edges -
   argues against edge reflection as primary.

F. RESIDUAL MIS-MEASUREMENT (ruled out as sole cause)
   RES is the undamped defect measured pre-relaxation - the
   right quantity for a fixed-point audit. It cannot
   manufacture a plateau; it faithfully reports that G(x)-x
   stays O(0.6). Combined with bounded recurrent values this
   is a genuine non-contractive/limit-cycle signature, not a
   telemetry artefact. (The separate KN_Q_RESID diagnostic
   measures a DIFFERENT map - one undamped evaluation on the
   committed field - hence its larger values ~9-12 vs ~0.7;
   both are consistent with a far-from-fixed-point P field.)

G. OPENMP/RACE (ruled out for SPW/FLD_P; see 3)
   All sweep work arrays are per-call locals; module SAVEs are
   per-ray indexed or critical-guarded. No mechanism found for
   cross-thread contamination of the P sweep.

Local linearisation: delta_x_{n+1} = J_G delta_x_n with
J_G = (1-OMEGA)I + OMEGA*DG (DG the GS Jacobian of G). J_G is
not formed explicitly anywhere; estimating its dominant
eigenvalue requires finite-difference probing (test T6/T7
below) or a matrix-free eigen-iteration on a single-ray
reduced model.

Spectral radius vs non-normal growth: a period-2-dominated
spectrum (lambda near -1 in DG, mapped by damping to
1-2*OMEGA ~ 0.1 ... stable) would NOT explain the cycle; the
cycle requires either |lambda_DG| large enough that the
damped factor |(1-OMEGA)+OMEGA*lambda| >= 1 (real lambda
<~-3.4 or complex pair outside the stability disk), or a
nonlinear saturation of a weakly unstable mode. Non-normal
transient growth (Trefethen & Embree) can delay but not by
itself sustain a stationary cycle - it needs an unstable or
marginally-stable eigen-structure plus saturation.

------------------------------------------------------------------------
8. SOURCE-LEVEL BUG SEARCH RESULTS (Part 6)
------------------------------------------------------------------------
Checked and CLEARED (with evidence):

- SPW ring indexing: FILL writes slot MOD(R-1,4)+1; row J's
  fresh SPJ is deposited at MOD(J-1,4)+1 = slot(J-1), which is
  exactly the slot row J-1 reads as its upwind source
  (MOD((J-1),4)+1 in the VNEW formula = slot(J)). Shift
  register verified consistent; KN_ORD_SEC uses
  MOD(R,4)+1 - MOD(R-1,4)+1 = slot(R+1)-slot(R) secant. No
  off-by-one found.
- VNEW stencil: FLD(J+1) upwind + slot(J)=S'(J+1) via A1 +
  SPJ=S'(J) via A2 - DSU*A3 - DSD*A4. Matches the formal
  solution structure (upwind intensity + two-node source +
  parabola corrections).
- IM_NI: read-only in P, snapshot of converged M, no aliasing.
- Bounds: L-1 donor guarded by L>1; compressing donor guarded
  by L+NNF2<=NCF; J-1 preview guarded by J-1>=1; J+1 reads
  bounded by march range; SPW 4-slot ring covers the
  4-row stencil KN_ORD_DS needs.
- GAM sign usage: >=0 ascending/L-1 signed-Q; <0 descending/
  matched |Q|. Consistent with expanding vs compressing
  characteristics.
- OMEGA_RAY: single scalar, reset per direction, floor 0.45
  enforced by MAX(); no path found where it exceeds [0.45,1].
- Residual state: T1 computed from VNEW vs VOLD BEFORE the
  damped commit; VOLD read before overwrite; RES_FREQ
  accumulated per L correctly; census gate
  MAX(VOLD,|VNEW|)>1e-4*VMAXF consistent with NBAD semantics.
- M/P array separation: FLD vs FLD_P distinct allocations;
  commit gated on both solves.
- OpenMP: per-ray locals private; shared counters critical or
  per-ray indexed.

SUSPECTED DEFECTS / CONCERNS (ranked by confidence that they
contribute to the observed cycle):

D1. [HIGH - mechanism, not yet proven defect]
    Unaccelerated frequency-coupled GS is non-contractive for
    this P field. Evidence: implicit map (Q*I term), large Q
    where CHI small, period-2-like oscillation at OMEGA floor,
    M/P asymmetry explained by different retained-store
    distance-to-fixed-point. NOT a code bug per se - a
    formulation/solver limitation. Test: T1 (manufactured),
    T6 (eigen probe).

D2. [MEDIUM]
    Lagged C-corner preview values at the marching front may
    inject a one-pass delay that participates in the cycle.
    Source: KN_ORD_SPW_FILL(FLD,J-1) line 11259 +
    KN_ORD_DS reads. Consequence: effective two-time-level
    coupling at leading-edge rows (j small = outer edge =
    argmax region). Minimal test: T7 (Jacobi vs GS ordering
    probe) / zero the preview lag on a single-ray harness.

D3. [MEDIUM-LOW]
    Sign asymmetry signed-Q (expanding) vs |Q| (compressing),
    lines 10690 vs 10714. If CHI_L(J) < 0 anywhere on an
    expanding row, the coupling sign flips while compressing
    rows keep positive coupling - inconsistent operator.
    Minimal test: dump MINVAL(CHI_L) on ray 83; manufactured
    negative-CHI case.

D4. [LOW]
    L=1 frequency closure drops the coupling term entirely
    (line 10686 guard). If the physical problem has inflow at
    the lowest bin this is an uncontrolled boundary; for ray
    83 the argmax L is interior so impact is likely small,
    but it is an unverified modelling choice vs Paper II.

D5. [LOW - telemetry gap, not solver defect]
    Per-pass NNEG/NCAP (floor/cap hit counts) are accumulated
    but not emitted in the probe pass records; cannot confirm
    from existing data whether the 1e-30 floor fires during
    the cycle (mechanism C). Minimal test: extend probe
    telemetry or read counters at exit.

D6. [INFO]
    j=1/j=NI orientation: Phase-5.2 docs called j=1-2 the
    "inner boundary"; traced geometry shows j=1-2 is the OUTER
    edge (depth K~1-2, top of atmosphere). Corrected here;
    interpretation of "boundary-concentrated" residual shifts
    from core-side to escape-side. No code change implied.

------------------------------------------------------------------------
9. PROPOSED TESTS AND ACCEPTANCE CRITERIA (Part 7 - design only)
------------------------------------------------------------------------
No implementation authorised. Each test lists purpose and a
concrete acceptance criterion.

T1. Manufactured expanding-only transfer problem.
    Purpose: verify the discrete equation against a known
    solution. Construct a 1-ray, few-depth, few-frequency
    problem with prescribed GAM>0, CHI, ETA, DTAU where the
    continuous/fine-grid solution is known (e.g. constant
    coefficients -> analytic exponential-integral solution,
    or Method-of-Manufactured-Solutions forcing).
    Acceptance: one undamped P pass reproduces the analytic
    VNEW to roundoff on interior cells; the fixed point of G
    equals the analytic solution to discretisation error;
    residual -> 0 under iteration.

T2. Zero-GAM limiting case.
    Purpose: verify S' -> SR and the frequency-order switch
    are benign. Set GAM=0 everywhere on a test ray.
    Acceptance: P sweep converges like a pure-absorption
    formal solve; S'(J,L)=SR(J,L) exactly; no donor reads
    (KN_DBG_NDON stays 0); identical result for ascending vs
    forced-descending order at machine precision.

T3. Smooth monotonic expanding velocity field.
    Purpose: isolate whether non-monotonicity is required for
    the cycle. Run a single-ray (or small-ray) P solve on a
    strictly expanding v(r) with realistic CHI/ETA.
    Acceptance: P converges within the pass cap for a range
    of Q magnitudes; if it also cycles, the defect is in the
    expanding operator itself, not the non-monotonic geometry.

T4. Controlled inner-boundary perturbation.
    Purpose: test D4/boundary coupling. Perturb IM_NI(L) by a
    known delta on one frequency, run P passes, verify the
    perturbation decays (or measure its amplification factor).
    Acceptance: boundary perturbation response is linear,
    bounded, and decays with pass count; any sustained
    amplification localises the unstable mode.

T5. Frequency-grid refinement.
    Purpose: test whether the cycle is a grid-scale
    instability. Repeat the ray-83 P solve (single-ray harness
    or short run) at NCF and 2*NCF (or shifted binning).
    Acceptance: if the limit-cycle amplitude/residual shrinks
    with refinement, the mechanism is a discretisation-scale
    instability (supports D1/D2); if invariant, suspect a
    formulation-level defect.

T6. Computed-vs-reconstructed residual.
    Purpose: verify RES semantics. Independently recompute
    |G(x)-x| for a saved iterate (offline replay of one pass
    from dumped FLD/PK state) and compare to reported RES.
    Acceptance: agreement to roundoff; confirms residual is
    the undamped defect and localises any state mismatch.

T7. Pass-order test (GS vs Jacobi discrimination).
    Purpose: verify the claimed ordering. Instrument (or
    reason from a reduced harness) whether row J's VNEW uses
    same-pass-updated FLD(J+1,L) and FLD(J,L-1) (GS) vs
    previous-pass values (Jacobi).
    Acceptance: demonstrate fresh-donor reads; if any donor
    is found lagged beyond the documented C-corner, that is a
    defect.

T8. Deterministic single-ray test, OpenMP disabled.
    Purpose: exclude all parallelism. Run ray 83's P sweep
    alone (or the whole ray solve) with OMP_NUM_THREADS=1.
    Acceptance: bit-identical or negligibly-different limit
    cycle -> parallelism excluded (expected, given private
    locals); any change would indicate a race.

T9. Retained-store sensitivity.
    Purpose: separate store-distance from sweep-contractivity.
    Seed FLD_P from (a) the retained store, (b) the M field,
    (c) a smoothed/relaxed store; compare cycle amplitude.
    Acceptance: if seeding from a nearer-to-fixed-point field
    shrinks the cycle, the sweep is conditionally convergent
    and the store anomaly is upstream; if the cycle persists
    at fixed amplitude, the map itself is non-contractive.

T10. Production-path regression.
    Purpose: guard the default. Run the standard SOB/KNOP
    configuration matrix cell(s) with all diagnostics off.
    Acceptance: bit-identical fields vs pre-audit baseline;
    no VADAT default changes; probe keys remain off by
    default.

Supporting diagnostics to add ONLY after approval:
- per-pass NNEG/NCAP/NNG_P emission (test D5);
- MINVAL(CHI_L)/MAXVAL(|Q|) per ray dump (tests D3, D1);
- single-ray P-sweep harness for T1-T9 without full model;
- finite-difference J_G probe: perturb FLD(J,L) by eps,
  measure Delta VNEW -> column of DG; power-iterate for
  dominant eigenvalue (T6 companion).

------------------------------------------------------------------------
10. REFERENCES AND CONCEPTS USED
------------------------------------------------------------------------
- Knop et al. (2009), Paper I, A&A 496, 295-298,
  DOI 10.1051/0004-6361:200810290 - comoving-frame formal
  solution for non-monotonic flows; used for the S' = SR + Q*I
  coupling structure and the two-direction (inward/outward)
  sweep decomposition. SPECIFIC CLAIMS TO VERIFY EXTERNALLY:
  the exact frequency-direction convention for expanding vs
  compressing cells, and the edge closures.
- Knop et al. (2009), Paper II, A&A 501, 813-820,
  DOI 10.1051/0004-6361/200911793 - the ordered sweep and the
  lagged C-corner construction (code comment cites Eqs. 31/32
  for the preview fills). Used to classify the SPW preview
  mechanism; the L=1 closure and Q-sign conventions need
  verification against the paper's discrete equations.
- Rybicki & Hummer (1991), A&A 245, 171 - accelerated Lambda
  iteration. Concept used: unaccelerated fixed-point
  iteration on the radiation field diverges in scattering-
  dominated regimes; the analogue here is the unaccelerated
  frequency-coupled GS map G. NOT used to claim the code
  implements ALI - it does not; the reference motivates why
  a non-contractive P map is plausible.
- Higham (2002), Accuracy and Stability of Numerical
  Algorithms, 2nd ed. - concepts used: fixed-point iteration
  x_{n+1}=x_n+Omega(G(x_n)-x_n), convergence governed by
  spectral radius of the damped Jacobian, and the distinction
  between residual ||G(x)-x|| and step norm ||x_{n+1}-x_n||.
- Trefethen & Embree (2005), Spectra and Pseudospectra -
  concept used: non-normal transient growth can produce large
  finite-time amplification even when the spectral radius is
  < 1; relevant to distinguishing a true instability from
  transient growth in J_G (section 7). No quantitative
  pseudospectral claim is made.

Caveat: this audit used the papers for conceptual framing
only. No equation or algorithm is asserted to be "from" a
paper without the relevant section having been checked; the
items flagged "verify externally" are exactly the places
where the code's conventions should be compared against the
primary sources before any fix is designed.

------------------------------------------------------------------------
11. WHAT THIS AUDIT ESTABLISHES vs LEAVES OPEN
------------------------------------------------------------------------
Established by source inspection:
- The exact discrete equation, ordering, and donor freshness
  of the P sweep (sections 2-5).
- The SPW ring is a consistent shift register, not an
  off-by-one defect.
- IM_NI is a clean one-way M->P boundary snapshot.
- Residual is the undamped defect measured pre-relaxation;
  OMEGA history is internally consistent with a period-2-like
  limit cycle of the damped map.
- No array-bound, sign, aliasing, or OpenMP defect found in
  the inspected paths.
- j=1-2 argmax is the OUTER edge, correcting Phase-5.2's
  "inner boundary" phrasing.

Left open (for ChatGPT review / proposed tests):
- Whether J_G's dominant eigen-structure actually produces
  the observed cycle (needs T5-T7 or a reduced harness).
- Whether the L=1 closure and Q-sign conventions match Knop
  Paper II exactly.
- Whether the retained-store anomaly is upstream cause or
  downstream symptom of the non-contractive map (T9).
- Whether any production remedy is warranted at all - none
  is proposed or authorised by this audit.

========================================================================
PHASE 5.3A CORRECTIONS - dated addendum
========================================================================

Independent review (against committed source and Knop et al. 2009,
Paper II, arXiv:0909.4264) found errors in the original audit. This
section supersedes the affected parts of sections 2, 5.1, 5.3, 5.4,
6, 7, 8 and 9 above. The original text is retained for provenance;
where the two disagree, THIS section is authoritative.

Status of claims used below:
  [V] verified code behaviour (traced in source)
  [M] measured numerical behaviour (Phase 5.2 runtime data)
  [H] hypothesised cause (plausible, undemonstrated)
  [D] demonstrated defect (none established to date)

------------------------------------------------------------------------
A1. CORRECTED EFFECTIVE-OPACITY ANALYSIS (supersedes sections 2, 5.1,
    5.3, 5.4)
------------------------------------------------------------------------

The original audit stated that expanding cells use plain CHI_L with
the GAM_MAG/DLN contribution restricted to compressing cells. That is
wrong for L > 1. Traced coefficient construction
(fg_j_cmf_nonmono_v1.f, per-frequency block ~lines 11900-12071):

L = 1 branch (lines 11900-11976):
  - CALL MATCH_FREQ(1) runs, then value guards clamp
    CHI_L >= 1e-30, ETA_L, GAM, GAM_MAG, dGAMdR to finite values.
  - Q_L(I) = 0 for all I (line 11927).
  - The opacity increment and Q_L assignment are gated by
        IF(GAM(I,LS) .LT. 0 .AND. KN_RF_COUPLING(I,L))
    (lines 11928, 11957): ONLY compressing cells (GAM < 0) with an
    in-range matched donor receive
        CHI_L += GAM_MAG/DLN          (line 11940)
        DCHR_L -= dGAMdR/DLN           (line 11941)
        Q_L    = GAM/(DLN*CHI_L)       (line 11963, signed GAM < 0)
  - For expanding cells (GAM >= 0) at L = 1: NO opacity increment,
    Q_L = 0, plain CHI_L, SR_L = ETA_L/CHI_L. In the sweep,
    KN_ORD_SPV's L > 1 guard returns SR only - consistent: no
    donor exists below the grid and the coupling is dropped.
  - Physical/numerical meaning: L = 1 is a one-sided frequency
    boundary. For expanding flow the upwind donor would be L = 0,
    outside the grid, so the coupling term is omitted entirely.
    Whether the |GAM|/DLN opacity increment should still apply at
    L = 1 (the paper's chi_hat depends on p_l^|, which exists
    regardless of donor availability) is an OPEN modelling
    question - see U3 below. Not asserted as a defect.

L > 1 branch (lines 11977-12071):
  - CALL MATCH_FREQ(L), then for EVERY cell with
    KN_RF_COUPLING(I,L) true:
        Q_L(I) = GAM_MAG(I,LS)/DLN            (line 11998, temporary)
        CHI_L(I) = CHI_L(I) + Q_L(I)          (line 11999)
    KN_RF_COUPLING (lines 12205-12211) returns .TRUE. immediately
    for GAM >= 0; for GAM < 0 it additionally requires the matched
    frequency to stay inside the grid (IKC checks). GAM_MAG =
    |GAM| >= 0 by construction (lines 78, 3566).
  - CHI_L is then re-guarded: non-finite or < 1e-30 -> 1e-30
    (lines 12002-12009). DLNM_L is guarded > 0 (lines 11992-11996).
  - DCHR_L gets SIGN(GAM)*dGAMdR/DLN for coupling cells (line 12016).
  - Q_L is then OVERWRITTEN (lines 12048-12058):
        IF(KN_RF_COUPLING(I,L)) Q_L(I) = GAM(I,LS)/(DLN*CHI_L(I))
        ELSE Q_L(I) = 0
    using the UPDATED (incremented, guarded) CHI_L.
  - SR_L = ETA_L/CHI_L (line 12059) with the same updated CHI_L.

CORRECTED effective opacity for L > 1, expanding cells (GAM > 0,
KN_RF_COUPLING = .TRUE.):

    chi_hat(J,L) = CHI_L(J) + GAM_MAG(J,LS)/DLNM_L(J)
                 = chi(J,L) + |GAM(J,LS)| / DLN_L
    chi_hat >= 1e-30   (guard, line 12004)
    Q_L(J)     = GAM(J,LS) / (DLNM_L(J) * chi_hat(J,L))   > 0
    DTAU(J,L), PK_A0..PK_A4, PK_SR, PK_Q all built from chi_hat.

So expanding cells DO receive the |GAM|/DLN generalised-opacity
increment for every L > 1 - identical in form to compressing cells.
The asymmetry is confined to L = 1, where expanding cells get plain
chi and no coupling while compressing cells keep the increment.

GAM sign: GAM(J,LS) > 0 defines the expanding branch; GAM_MAG =
|GAM| > 0; DLNM_L > 0 after guard; chi_hat > 0 after guard. Hence
PK_Q = GAM/(DLN*chi_hat) is STRICTLY POSITIVE for every expanding
coupling cell. The signed-Q vs |Q| asymmetry (signed in the
expanding SPV branch, ABS in the compressing branch) cannot
produce a sign defect on ray 83: both evaluations yield a
positive coupling coefficient. See A3.

------------------------------------------------------------------------
A2. KNOP PAPER II TO CODE CORRESPONDENCE
------------------------------------------------------------------------

Paper: Knop, Hauschildt & Baron (2009), A&A 501, 813-820
( arXiv:0909.4264 ), equations verified against the arXiv HTML
rendering. The code is NOT a literal transcription of Paper II:
it uses a different source-function splitting (S' = SR + Q*I_donor
inside the formal solve) rather than the paper's A/B/C matrix +
Delta-I decomposition. The table below records the strongest
supportable correspondence; entries marked UNRESOLVED are places
where no algebraic identity has been established.

Paper quantity            | Code quantity              | Location            | Correspondence
--------------------------|----------------------------|---------------------|----------------------------------
a_l (signed coupling coef)| GAM(I,LS) (signed)         | fg_j_cmf_nonmono_v1  | STRUCTURAL: both are the signed
                          |                            | ~line 3566           | velocity-gradient coupling term;
                          |                            |                     | a_l>0 <-> GAM>0 <-> upwind donor
                          |                            |                     | at the lower grid index. Sign
                          |                            |                     | convention consistent IF the
                          |                            |                     | paper's wavelength ordering
                          |                            |                     | (lambda ascending) maps to the
                          |                            |                     | code's frequency ordering
                          |                            |                     | (nu ascending) with the
                          |                            |                     | appropriate reversal. UNRESOLVED:
                          |                            |                     | exact coefficient identity
                          |                            |                     | (a_l vs V/R*(1+sigma*mu^2)
                          |                            |                     | factors).
p_l^- (expanding donor wt) | Q*FLD(J,L-1) term in SPV    | KN_ORD_SPV lines     | STRUCTURAL: paper B^- multiplies
                          | (L>1 branch)                | 10684-10691          | I_{l-1} for a_l>=0 with
                          |                            |                     | coefficient ~ -(a_l/chi_hat)
                          |                            |                     | *p_l^- *beta; code uses
                          |                            |                     | Q*FLD(J,L-1), Q =
                          |                            |                     | GAM/(DLN*chi_hat) > 0.
                          |                            |                     | For log-spaced grids
                          |                            |                     | |p_l^-| ~ 1/DLN, so the
                          |                            |                     | magnitudes are compatible.
                          |                            |                     | UNRESOLVED: exact algebraic
                          |                            |                     | identity including the
                          |                            |                     | alpha/beta/gamma weighting
                          |                            |                     | the paper applies.
p_l^| (diagonal coupling)  | GAM_MAG/DLN added to CHI_L  | lines 11998-11999,    | STRUCTURAL: paper chi_hat =
                          | (opacity increment)         | 12055               | chi + xi*a_l*p_l^| with
                          |                            |                     | p_l^| = lambda_l/(lambda_l -
                          |                            |                     | lambda_{l-1}) > 0 for
                          |                            |                     | a_l>=0; code adds
                          |                            |                     | |GAM|/DLN > 0. Both are
                          |                            |                     | positive gradient terms in
                          |                            |                     | the effective opacity.
                          |                            |                     | UNRESOLVED: whether
                          |                            |                     | |GAM|/DLN == a_l*p_l^|
                          |                            |                     | algebraically (depends on
                          |                            |                     | grid-spacing convention).
chi_hat (generalised opac.)| CHI_L after += GAM_MAG/DLN  | lines 11999-12009,   | VERIFIED correspondence at
                          | (chi_hat)                   | used in DTAU/Q/SR    | the level of "effective
                          |                            | (lines 12055,12059,  | opacity = chi + positive
                          |                            | 12076-12088)         | coupling term": code's
                          |                            |                     | chi_hat enters DTAU
                          |                            |                     | (via KN_DTAU_INT_V1),
                          |                            |                     | SR = ETA/chi_hat and
                          |                            |                     | Q = GAM/(DLN*chi_hat),
                          |                            |                     | matching the paper's role
                          |                            |                     | of chi_hat in d tau =
                          |                            |                     | chi_hat ds and in the
                          |                            |                     | source denominators.
PK_Q (off-diagonal wl      | B^-, B^+ coefficients      | KN_ORD_SPV, Q_L      | STRUCTURAL ONLY: the paper's
coupling)                  | (Eqs. 19-21,26-28)          | construction          | wavelength-coupling matrix
                          |                            |                     | elements carry alpha/beta/
                          |                            |                     | gamma spatial weights AND
                          |                            |                     | p_l^+/- factors; the code
                          |                            |                     | collapses this to a single
                          |                            |                     | Q multiplying the donor
                          |                            |                     | intensity inside S'.
                          |                            |                     | UNRESOLVED: whether the
                          |                            |                     | omission of the spatial
                          |                            |                     | quadrature weights is
                          |                            |                     | exact for the piecewise-
                          |                            |                     | linear source assumed, or
                          |                            |                     | an approximation.
B^diag (self-coupling)     | (no explicit equivalent)    | VNEW formula lines   | DIFFERENCE: Paper II Eqs.
(1-B^diag)^-1 prefactor    |                            | 11280-11283          | 31-34 carry an explicit
                          |                            |                     | (1-B^diag)^-1 where
                          |                            |                     | B^diag contains the
                          |                            |                     | S-tilde self-term
                          |                            |                     | (~ -4 a_l I_l/chi_hat for
                          |                            |                     | xi=1). The code's S'
                          |                            |                     | contains NO FLD(J,L)
                          |                            |                     | self-term for expanding
                          |                            |                     | rows, so its analogue of
                          |                            |                     | B^diag is zero and the
                          |                            |                     | prefactor is 1. UNRESOLVED:
                          |                            |                     | whether the code's
                          |                            |                     | splitting (self-coupling
                          |                            |                     | folded entirely into
                          |                            |                     | chi_hat) is algebraically
                          |                            |                     | equivalent to the
                          |                            |                     | paper's xi choice, or a
                          |                            |                     | different (possibly
                          |                            |                     | Paper-I-style) splitting.
PK_A0..PK_A4 (formal        | alpha/beta/gamma-type      | lines 12129-12144,    | STRUCTURAL: both are
solution coefficients)     | exponential-integral       | EEL/E0-E3 moments    | exponential-integral
                          | weights over DTAU          |                     | weights over segment
                          |                            |                     | optical depth; the paper
                          |                            |                     | cites Olson & Kunasz
                          |                            |                     | (1987)/Hauschildt (1992)
                          |                            |                     | for the same family.
                          |                            |                     | UNRESOLVED: term-by-term
                          |                            |                     | mapping (the code's A1/A2
                          |                            |                     | two-node source weights
                          |                            |                     | vs the paper's
                          |                            |                     | alpha/beta/gamma on
                          |                            |                     | i-1,i,i+1 nodes) - the
                          |                            |                     | stencils are organised
                          |                            |                     | differently (parabola
                          |                            |                     | correction vs linear
                          |                            |                     | interpolation).
Lagged C-corner (i+1,l-1)   | SPW preview + KN_ORD_DS     | KN_ORD_SPW_FILL      | VERIFIED correspondence:
(n) in Eqs. 31-32           | slopes at marching front    | line 11259,          | Paper II Eqs. 31-32
                          |                            | KN_ORD_DS 10759      | explicitly keep the
                          |                            |                     | spatial-downwind /
                          |                            |                     | wavelength-upwind corner
                          |                            |                     | at iteration (n); the
                          |                            |                     | code's comment at
                          |                            |                     | 10808-10816 says the
                          |                            |                     | same. This is the one
                          |                            |                     | place lag is BY DESIGN
                          |                            |                     | in both formulations.
L = 1 boundary              | L=1 branch: Q=0, plain chi  | lines 11900-11976,   | DIFFERENCE/OPEN: the
                          | for expanding; increment    | KN_ORD_SPV L>1 guard | paper's p_l coefficients
                          | kept for compressing        | line 10686           | are defined for interior
                          |                            |                     | l; edge treatment is not
                          |                            |                     | specified in the excerpted
                          |                            |                     | equations. Whether the
                          |                            |                     | code's asymmetric L=1
                          |                            |                     | handling matches the
                          |                            |                     | paper's boundary closure
                          |                            |                     | is UNRESOLVED.

Algebra status: no claim is made that the code's discrete
equation is algebraically identical to Paper II Eqs. 2, 11-17,
19-30. The correspondence established is at the level of
structures (signed coupling -> upwind donor; positive gradient
term -> effective opacity; lagged downwind corner; exponential
weights). The term-by-term algebraic identity is UNRESOLVED and
is listed as open question U2.

------------------------------------------------------------------------
A3. SIGNED-Q HYPOTHESIS - REASSESSED (supersedes section 5.3 and
    defect candidate D3)
------------------------------------------------------------------------

Original candidate D3 proposed that CHI_L could go negative on an
expanding row, flipping the sign of the signed-Q coupling while
compressing rows keep |Q|. Inspection of the guards shows this
cannot happen at the point where PK_Q is computed:

- For L > 1: CHI_L is incremented by GAM_MAG/DLN (line 11999),
  then IMMEDIATELY re-guarded: any non-finite or CHI_L < 1e-30 is
  reset to 1e-30 (lines 12002-12009) BEFORE Q_L is recomputed at
  line 12055. DLNM_L is likewise guarded > 0 (lines 11992-11996).
- For L = 1: CHI_L guarded >= 1e-30 at lines 11903-11908 before
  any use.
- Therefore at every Q_L evaluation, chi_hat >= 1e-30 > 0,
  DLN > 0, and for expanding rows GAM > 0 -> PK_Q > 0 strictly.
- GAM_MAG = |GAM| >= 0 by construction (lines 78, 3566), so the
  opacity increment is non-negative; chi_hat >= chi whenever the
  increment applies (chi itself could in principle be small but
  the guard floors the sum).

CONCLUSION: negative CHI_L at the PK_Q computation is impossible
under the current guards. PK_Q is strictly positive on every
expanding coupling cell. The signed-vs-ABS asymmetry between the
two SPV branches is therefore MOOT for ray 83 (and for any model
with chi_hat > 0): both branches produce a positive coupling
coefficient. D3 is REMOVED as an active defect candidate. It is
retained only as a note that IF a future code path computed Q
before the chi_hat guard, the sign asymmetry would matter.

------------------------------------------------------------------------
A4. PERIOD-2 CLAIM - REASSESSED (supersedes sections 6.3, 7 mechanism
    framing)
------------------------------------------------------------------------

Phase 5.2 measured a persistent bounded oscillatory residual
plateau (RES ~ 0.57-0.80, excursions ~3.3, quantised across 70
refreshes). That is the MEASURED result [M]. The original audit
then asserted a "period-2-like limit cycle of the damped map" and
used RES ~ 0.6 at OMEGA = 0.45 to infer an orbit amplitude ~0.3.
That inference is NOT established:

- RES is a max-norm scalar of the undamped defect |G(x)-x|. A
  scalar residual history cannot establish periodicity of the
  full intensity vector x in (NI x NCF) space.
- The orbit-amplitude estimate assumed a strict period-2 orbit
  with alternating sign of G(x)-x; a higher-period orbit, a
  quasi-periodic orbit, or a noisy bounded non-convergent
  trajectory can produce the same scalar history.
- Correct distinction: for a candidate period-p orbit, the
  relevant quantities are the FIELD-LEVEL distances
        d1_n = ||x_n - x_{n-1}||
        d2_n = ||x_n - x_{n-2}||
        (and d_p for the candidate period p)
  Periodicity requires d_p -> 0 while d1 stays bounded away
  from 0. The orbit amplitude about the midpoint is
  (x_n - x_{n+p/2})/2 for even p - NOT the same as the
  pass-to-pass distance d1. The original text conflated these.

REQUIRED TELEMETRY (design only, not implemented):
- Snapshot FLD_P (or a reduced probe set) at pass end for the
  probe ray; compute ||x_n - x_{n-k}|| for k = 1,2,3,4 under a
  FIXED normalisation (e.g. divide by the pass-start VMAXF(L)
  snapshot from the SAME pass index, or by a stored reference
  field - VMAXF(L) itself changes between passes because it is
  recomputed from the current FLD each pass, line ~10930).
- Compare update-vector signs: sign(VNEW - VOLD) per cell across
  consecutive passes; a period-2 mode shows systematic sign
  alternation in the same cells.
- Localise the update vector: report argmax of |x_n - x_{n-1}|
  in (J,L) and compare with the RES argmax (j in {1,2}).
- Acceptance for calling it period-2: d2_n / d1_n -> 0 while
  d1_n stays O(0.1-1) of its running mean, sustained over >= 10
  passes, on the full field (not just the argmax cell).

Until such telemetry exists, the correct statement is:
"persistent bounded oscillatory residual plateau; periodicity
class unresolved; period-2 is a hypothesis, not a measurement."

------------------------------------------------------------------------
A5. M-STORE CLAIMS - SEPARATED (supersedes section 7 mechanism A
    wording)
------------------------------------------------------------------------

The original audit wrote that "M converges because its retained
store was already at the fixed point and/or its coupling
structure differs; P's retained store was NOT at the fixed
point." Only part of that is supported:

FACTS [M]:
- M converges in 2 passes every refresh (IXIT_M = 1, res ~ 4e-8).
- M's retained store shows no identified anomaly comparable to
  the P-channel store discrepancy documented in Phase 5.2.
- IM_NI (M's inner-boundary row) is consumed by P as a fixed
  boundary; its accuracy inherits M's convergence.

NOT ESTABLISHED:
- M's initial distance from its own fixed point is NOT
  independently known. Two passes to tolerance is consistent
  with EITHER a store already near the fixed point OR a
  strongly contractive map OR a loose effective tolerance for
  that channel. The data do not distinguish these.
- No measurement demonstrates that the P store anomaly is the
  CAUSE of the P non-convergence rather than a co-symptom.

Corrected mechanism wording: the M/P asymmetry is an OBSERVED
fact; its explanation is OPEN. Candidate explanations
(store proximity, coupling structure, boundary role) are
hypotheses to be tested (T9 revised below), not conclusions.

------------------------------------------------------------------------
A6. OUTER RADIAL ENDPOINT ANALYSIS (new candidate mechanism -
    distinct from L=1 frequency closure)
------------------------------------------------------------------------

The P residual argmax sits at j in {1,2} - the OUTER/escape edge
(corrected orientation, original section 4.1). The outermost
point j=1 is computed by the interior VNEW formula; its source
reconstruction differs from interior rows through KN_ORD_DS's
R=1 endpoint branch (lines 10762-10772):

    SA = KN_ORD_SEC(1,L) = (S'(2,L) - S'(1,L)) / DTAU(1,L)
    SB = KN_ORD_SEC(2,L) = (S'(3,L) - S'(2,L)) / DTAU(2,L)
    TDT = DTAU(1,L) + DTAU(2,L)
    DD = SA + (SA - SB) * DTAU(1,L)/TDT        (slope extrapolation)
    DS(1,L) = (SIGN(SA)+SIGN(DD)) * MIN(|SA|, 0.5*|DD|)

- SA, SB are secant slopes of S' between rows (1,2) and (2,3);
  all three S' values are fresh within the pass at the time
  J=1 is processed (rows 2,3 solved earlier; SPJ(1,L) computed
  in row 1's own L-loop before KN_ORD_DS(1,L) is called).
- DD extrapolates the slope trend one half-cell beyond row 1;
  the 0.5 factor damps the extrapolated magnitude; the
  (SIGN(SA)+SIGN(DD)) factor ZEROES DS when the extrapolated
  slope reverses sign relative to SA (opposite signs -> 0).
- When signs agree, |DS(1)| = min(2|SA|, |DD|) - up to twice
  the first secant magnitude. The parabola correction enters
  VNEW(1,L) as -DSD*A4(1,L); A4 is O(DTAU) for thin cells, so
  the correction is bounded by ~2*secant*DTAU ~ O(S' jump)
  across the first two cells.
- At J=2: DSD = KN_ORD_DS(2,L) uses the INTERIOR branch
  (secants SEC(1), SEC(2), weighted DD, three-way min with
  0.5|DD|); DSU = KN_ORD_DS(3,L) likewise. So the j=1/j=2
  residual asymmetry is plausible: j=1 carries the endpoint
  extrapolation while j=2 carries the interior clamp.
- Error entering from j=2 can be AMPLIFIED at j=1 only through
  the S'(1),S'(2),S'(3) secants and the A-coefficient weights;
  the sign/min clamps bound but do not eliminate this path.
  Whether the endpoint formula is consistent with Paper II's
  boundary treatment: Paper II does not specify the spatial
  endpoint closure in the excerpted equations (the alpha/beta/
  gamma coefficients at i=1,i=ni come from Olson & Kunasz 1987 /
  Hauschildt 1992); the monotonic-parabola DS is CMFGEN's own
  construction in the spirit of SLOPES. CONSISTENCY UNRESOLVED.

This is registered as a DISTINCT candidate mechanism (M-ENDPOINT
below), not merged with the L=1 frequency-boundary question.

------------------------------------------------------------------------
A7. REVISED TEST ACCEPTANCE CRITERIA (supersedes T1, T5, T6, T9)
------------------------------------------------------------------------

T1 (manufactured expanding-only problem) - REVISED:
  Build a discrete manufactured problem on the SAME grid and
  SAME coefficient machinery (PK_A*, PK_Q, PK_SR, DTAU, SPW/DS
  stencils). Two acceptance levels:
    (a) Discrete manufactured solution: choose FLD* satisfying
        the discrete equations exactly (solve the linear system
        or construct forcing), verify one undamped pass returns
        FLD* to roundoff and RES -> 0 under iteration.
    (b) Continuous analytic solution: if used, acceptance is
        agreement to the DISCRETISATION error (measured by grid
        refinement), NOT roundoff. Roundoff-level agreement with
        a continuous solution is impossible where spatial or
        frequency truncation error exists and would falsely
        indicate a defect.

T5 (frequency-grid refinement) - REVISED:
  Refinement alone cannot identify a mechanism. Required
  accompanying measurements at each resolution:
    - chi_hat(J,L) and Q(J,L) fields (does the coupling
      strength per bin stay constant, grow, or shift?)
    - effective optical depth per frequency bin
    - residual SHAPE: argmax (J,L) distribution, not just the
      scalar max
    - convergence behaviour: passes to exit, OMEGA trajectory,
      plateau level
  Acceptance: a mechanism hypothesis is supported only if the
  predicted scaling (e.g. cycle amplitude ~ f(NCF)) matches;
  otherwise refinement is diagnostic, not discriminating.

T6 (computed vs reconstructed residual) - REVISED:
  A finite-difference Jacobian-vector product must replay the
  COMPLETE ordered sweep from perturbed initial field states:
  perturb FLD_P by eps*direction, run the full SWEEP_FLD_ORD
  (boundary + all rows + all frequencies, fixed PK coefficients
  and fixed IM_NI), measure ||G(x+eps*d)-G(x)||/eps. Perturbing
  one local VNEW calculation measures a single stencil row, not
  the sweep Jacobian - the GS ordering makes the Jacobian
  non-local through the fresh-donor dependencies.
  Acceptance: Jv product consistent with two independent eps
  values (eps and eps/10 agree to truncation level); residual
  reconstruction agrees with reported RES to roundoff.

T9 (retained-store sensitivity) - REVISED:
  Different seeds probe sensitivity to initial condition and
  possible basins of attraction; they CANNOT by themselves
  establish the causal origin of the historical store anomaly.
  Required for the causal claim:
    - an INDEPENDENT reference solution for the P field on ray
      83 (e.g. a fully-converged single-ray solve by a
      different method - dense direct solve of the discrete
      system, or many more undamped iterations with a
      verified-contractive variant), against which both the
      retained store and the sweep endpoint are compared;
    - demonstration that the store anomaly persists when the
      sweep is started from the reference solution (isolating
      store-vs-sweep causality).
  Without the reference solution, T9 output is labelled
  "seed sensitivity", not "cause".

------------------------------------------------------------------------
A8. UPDATED MECHANISM RANKING
------------------------------------------------------------------------

Confidence is assigned to the HYPOTHESIS, not to plausibility.
No mechanism is demonstrated; no code defect is proven.

M1. Non-contractive frequency/radial GS map [H, leading]
    The S' = SR + Q*I_donor term makes G implicit; where Q is
    large the undamped map can have spectral radius > 1 or
    real eigenvalues < -(2/OMEGA-1), producing bounded
    oscillation under damping + nonlinear floors. Consistent
    with: persistent plateau, OMEGA floor, broad bad-cell
    census. NOT demonstrated: no J_G estimate exists.
    Tests: T1, T3, T6, T9.

M2. Lagged C-corner coupling at the marching front [H]
    By-design lag (Paper II Eqs. 31-32 analogue) injects
    one-pass-old S' at leading-edge rows; the argmax sits at
    the outer edge where the front terminates. Could sustain
    a two-pass mode but is also present in convergent M/P
    solves elsewhere. Tests: T7, T6.

M3. Outer radial endpoint reconstruction [H, NEW]
    KN_ORD_DS R=1 extrapolation can produce |DS| up to 2x the
    first secant; interacts with A4 weighting at j=1. A
    plausible amplifier for error arriving from j=2 but the
    clamps bound it. Consistency with Paper II boundary
    closure unresolved. Tests: manufactured endpoint problem,
    T5 residual-shape.

M4. Effective-opacity / wavelength-coupling discretisation
    mismatch [H, downgraded]
    chi_hat = chi + |GAM|/DLN is applied for ALL L>1 coupling
    cells (corrected A1); the remaining question is whether
    |GAM|/DLN == a_l*p_l^| and whether the single-Q source
    splitting (vs paper's A/B/C + Delta-I) is exact. This is
    now a FORMULATION question, not a sign/opacity bug.
    Tests: algebraic comparison (U2), T1.

M5. L=1 frequency boundary closure [H, narrowed]
    Expanding L=1 drops both the donor AND the opacity
    increment; compressing L=1 keeps the increment. Asymmetric
    but arguably consistent (no donor exists). Effect on ray
    83 likely small (argmax L interior). Tests: T2 zero-GAM,
    edge-sensitive manufactured case.

M6. Nonlinear positivity floor / upper cap [H, low]
    VNEW < 0 -> 1e-30 and CAP = 1e20 are nonsmooth; if they
    fire during the cycle they can sustain relaxation
    oscillation. P2_DAMP is OFF by default. UNVERIFIED:
    per-pass NNEG/NCAP not in probe telemetry (gap D5 stands).

M7. Inner boundary condition via IM_NI [H, low]
    I_P(NI) relaxed toward I_M(NI) assumes equality of the two
    streams at the inner boundary. If the true fixed point
    requires I_P != I_M there, a persistent boundary defect
    propagates outward. Argmax is NOT at j=NI, so any effect
    must amplify along the march. Tests: T4 perturbation.

Ordering rationale: M1-M3 are the only mechanisms that act
AT the observed argmax location (outer edge) with the right
qualitative signature (bounded oscillation). M4-M7 are either
upstream formulation questions or weaker local effects.
NONE is assigned HIGH confidence; all require the proposed
tests.

------------------------------------------------------------------------
A9. UNRESOLVED MATHEMATICAL QUESTIONS
------------------------------------------------------------------------

U1. Does J_G (the sweep Jacobian, including fresh-donor GS
    couplings and the lagged C-corner) have spectral radius
    or damped-iteration stability boundary compatible with the
    observed plateau? Requires T6 finite-difference probing or
    a reduced single-ray dense solve.

U2. Algebraic identity between the code's discretisation and
    Paper II Eqs. 2, 11-17, 19-30: specifically whether
    |GAM|/DLN == a_l*p_l^|, whether Q == the paper's
    wavelength-coupling coefficient including alpha/beta/gamma
    spatial weights, and whether omitting the S-tilde self-term
    (B^diag = 0 in the code) is an exact splitting or an
    approximation. Requires term-by-term algebra or a
    manufactured comparison.

U3. Should the L=1 expanding cell retain the |GAM|/DLN
    opacity increment even though the donor is absent? The
    paper's chi_hat depends on p_l^| which exists at the edge;
    the code drops it for expanding L=1 but keeps it for
    compressing L=1. Consistency with Paper II edge treatment
    unresolved (paper excerpt does not specify l=1 closure).

U4. Is the j=1 endpoint slope formula (KN_ORD_DS R=1 branch)
    consistent with the boundary closure the paper's
    alpha/beta/gamma coefficients imply at i=1? The paper
    cites Olson & Kunasz (1987)/Hauschildt (1992) for those
    coefficients; the code's monotonic-parabola DS is a
    different (CMFGEN-native) construction.

U5. Periodicity class of the P orbit: period-2 vs higher
    period vs quasiperiodic vs noisy bounded trajectory.
    Requires the field-level telemetry defined in A4.

U6. Causal direction between the retained-store anomaly and
    the sweep non-convergence: does the store anomaly cause
    the cycle, does the cycle cause the anomaly, or are both
    symptoms of an upstream defect (e.g. coefficient or
    boundary construction)? Requires the independent
    reference solution defined in T9.

U7. Sign convention mapping between the paper's wavelength
    grid (lambda ascending) and the code's frequency grid
    (nu ascending): whether "upwind at l-1" in the paper
    corresponds to "donor at L-1" in the code for expanding
    flow, including the lambda<->nu reversal. Structural
    consistency is established; the index-direction algebra
    is not yet written out.

------------------------------------------------------------------------
A10. WHAT THE CORRECTED AUDIT ESTABLISHES
------------------------------------------------------------------------

Verified code behaviour [V]:
- chi_hat = chi + |GAM|/DLN applied to ALL L>1 coupling cells
  (expanding and compressing); plain chi + no coupling only
  for expanding L=1.
- PK_Q strictly positive on expanding coupling cells; the
  signed/ABS asymmetry cannot flip sign under current guards.
- SPW ring is a consistent shift register; residual is the
  undamped defect pre-relaxation; OMEGA floor 0.45.
- KN_ORD_DS endpoint formulas and their clamps, as written
  in A6.

Measured behaviour [M]:
- Persistent bounded oscillatory residual plateau ~0.6-0.7
  with excursions ~3.3; argmax j in {1,2} (OUTER edge);
  broad L distribution; M converges in 2 passes.

Hypothesised [H]: all seven mechanisms M1-M7; none demonstrated.
Demonstrated defect [D]: NONE.

The audit does NOT establish that any production change is
warranted. The next authorised step is the test programme
(T1-T10 revised) after ChatGPT reviews this corrected
correspondence.

========================================================================
PHASE 5.3B - WAVELENGTH DISCRETISATION AND JACOBIAN AUDIT
========================================================================
Dated addendum to the Phase 5.3A corrections (commit 15ddf6c3).
Audit-only: no Fortran, parameter, or behaviour changes. Claim
labels as before: [V] verified code behaviour, [M] measured
numerical behaviour, [H] hypothesised cause, [D] demonstrated
defect. No defect is demonstrated by anything below.

------------------------------------------------------------------------
B1. EXACT WAVELENGTH-COEFFICIENT COMPARISON (Task 1)
------------------------------------------------------------------------

Paper II (arXiv:0909.4264, verified against the arXiv HTML
rendering), for a_lambda >= 0 on an ascending-wavelength grid
lambda_{l-1} < lambda_l < lambda_{l+1}:

    p_l^-  = -lambda_{l-1}/(lambda_l - lambda_{l-1})   (< 0)
    p_l^|  =  lambda_l    /(lambda_l - lambda_{l-1})    (> 0)
    p_l^+  =  0
    chi_hat_paper = chi + xi * a_l * p_l^|

Code (fg_j_cmf_nonmono_v1.f): expanding cells (GAM > 0) use

    DLN        = ln(nu_{L-1}/nu_L) = ln(lambda_L/lambda_{L-1}) > 0
    chi_hat_cd = chi + |GAM|/DLN                    (lines 11998-11999)
    Q_cd       = GAM/(DLN * chi_hat_cd)             (line 12055)
                 -> donor term Q_cd * FLD(J,L-1)   (KN_ORD_SPV, def. 10675)

Exact transformations (verified algebraically):

    p_l^|   = lambda_l/(lambda_l-lambda_{l-1})
            = 1/(1 - lambda_{l-1}/lambda_l)
            = 1/(1 - exp(-DLN))
    |p_l^-| = lambda_{l-1}/(lambda_l-lambda_{l-1})
            = exp(-DLN)/(1-exp(-DLN))
            = 1/(exp(DLN) - 1)
    (identity check: p_l^| - |p_l^-| = 1 exactly)

Comparison ratios requested:

    R_chi   = [1/DLN] / [1/(1-exp(-DLN))] = (1-exp(-DLN))/DLN
    R_donor = [1/DLN] / [1/(exp(DLN)-1)]   = (exp(DLN)-1)/DLN

Small-DLN expansions:

    R_chi   = 1 - DLN/2 + DLN^2/6 - DLN^3/24 + ...
    R_donor = 1 + DLN/2 + DLN^2/6 + DLN^3/24 + ...

The two ratios are NOT identical and have OPPOSITE first-order
deviations: the code's 1/DLN under-estimates the paper's p_l^|
opacity weight by ~DLN/2 and over-estimates the |p_l^-| donor
weight by ~DLN/2. They do not share a common correction factor.

Interpretation: the paper's p_l coefficients are the exact
coefficients of a FORWARD DIFFERENCE of d(lambda I)/d lambda on
the linear lambda grid: [lambda_l I_l - lambda_{l-1} I_{l-1}]
/(lambda_l - lambda_{l-1}) = p_l^| I_l + p_l^- I_{l-1}. The
code's 1/DLN weights correspond to replacing the linear spacing
lambda_l - lambda_{l-1} by lambda_l * DLN = lambda_l *
ln(lambda_l/lambda_{l-1}), i.e. a first-order LOGARITHMIC-grid
discretisation of the same derivative (equivalently, a
d/dln(lambda) difference). This is a controlled approximation
of Paper II's linear discretisation, exact only as DLN -> 0.

GRID STATISTICS (measured, not invented):
Source: CONT_FREQ of run
cmfgen_validation/m15_2p82d_KNOP_p8_probe/cmfgen/2.82days/
(NCF = 56377, header !NCF,56377; ML 1..56377 monotonic
decreasing frequency 200 -> 3.0e-3 in 1e15 Hz, i.e. wavelength
~1.5 -> ~1e5 Ang ascending with index).

    DLN = ln(nu_{L-1}/nu_L) over all 56376 interior bins:
      MIN    = 1.6592e-05
      MAX    = 1.9368e-03
      MEDIAN = 1.6680e-04
      MEAN   = 1.9702e-04
      92.8% of bins have DLN > 1e-4;  1.5% have DLN > 1e-3.

    R_chi   in [0.999032, 0.999992]   -> code opacity increment
                                          LOW by <= 0.097%
    R_donor in [1.000008, 1.000969]   -> code donor coefficient
                                          HIGH by <= 0.097%
    Worst bins at DLN = 1.94e-3 (the coarsest ~1.5% of the grid,
    around nu ~ 1.1e2 x1e15 Hz and the grid edges).

CONCLUSION for Task 1: classification B-with-caveat.
  - Versus Paper II the code uses a controlled first-order
    logarithmic discretisation whose coefficient error is
    <= ~0.1% on this grid, with opposite-sign errors in the
    opacity and donor coefficients.
  - CAVEAT: the code cites "Knop Paper I Eq 7" (chi_hat =
    chi + |GAM|/dLOG_NU) and "Eq 9" (signed Q) in its own
    diagnostics (OUTGEN lines 702-703). If Paper I defines
    exactly this logarithmic form, then versus Paper I the
    code may be EXACT and the Paper-II difference is a
    formulation choice between the two papers (Paper II
    introduced the xi Crank-Nicholson mixing and the
    fully-implicit variant). Paper I's text was NOT available
    for this audit; verifying Eq 7/9 of Paper I is required
    to promote this from B to A-or-C. See U8.
  - xi question: the code has no [1-xi] remainder terms; its
    unconditional chi_hat increment and signed-Q source are
    consistent with xi = 1 (fully implicit wavelength
    discretisation) in Paper II's notation, modulo the
    discretisation difference above.
  - a_l vs GAM: GAM = 3.33564e-6 * (V_RAY/R_RAY) *
    (1 + SIGMA_RAY*MU^2) (lines 3563-3566) is the signed
    comoving-frame frequency-gradient coefficient. Structural
    correspondence to a_l is established (signed, zero at
    turning points, enters opacity via |.| and source via
    signed Q); the exact algebraic identity including all
    geometrical factors is UNVERIFIED pending Paper I's a_l
    definition (U2/U8). Even if a_l = c*GAM for some constant
    c, that factor rescales both opacity and donor terms
    coherently and cannot by itself produce the observed
    O(0.6) residual.
  - QUANTITATIVE CONTEXT: the coefficient discrepancy is
    <= 0.1% while the measured P residual plateau is ~60-70%.
    The discretisation difference is ~600x smaller than the
    observed defect and CANNOT be its cause. This downgrades
    mechanism M4 accordingly.

------------------------------------------------------------------------
B2. FREQUENCY INDEX DIRECTION - U7 RESOLVED (Task 2)
------------------------------------------------------------------------

Code frequency grid (CONT_FREQ, same run): ML = 1 has the
HIGHEST frequency (200e15 Hz) and ML = NCF the lowest
(3e-3e15 Hz): nu DESCENDS with index L. Therefore

    lambda_L = c/nu_L  ASCENDS with L
    (lambda_1 ~ 1.5 Ang ... lambda_NCF ~ 1e5 Ang)

so the code's frequency index L orders wavelengths in the SAME
direction as Paper II's wavelength index l: code L <-> paper l
with lambda_L = lambda_l (up to grid construction details).
The solver's L is identified with CONT_FREQ's ML: both index
the same NCF-element model frequency array (NCF = 56377 in the
cited run), and KN_DLN_P(L) = ln(nu_{L-1}/nu_L) > 0 confirms
the descending-nu ordering is the one the coefficients use.

Log-interval transformation:

    DLN_code = ln(nu_{L-1}/nu_L)
             = ln( (c/lambda_{L-1}) / (c/lambda_L) )
             = ln(lambda_L/lambda_{L-1})
             = DLN_lambda  > 0   (since lambda_L > lambda_{L-1})

Sign of the wavelength derivative: for expanding flow
(GAM > 0 <-> a_lambda > 0) the upwind direction in wavelength
is toward SMALLER lambda (the photon redshifts to larger
lambda as it propagates; upstream it had smaller lambda).
Paper II's a_l >= 0 branch takes the donor at l-1
(p_l^- != 0, p_l^+ = 0). The code's expanding branch takes the
donor at L-1, which is the same physical bin:
lambda_{L-1} < lambda_L. CONSISTENT.

Therefore: code donor index L-1 <-> paper donor index l-1;
code ascending-L frequency march <-> paper ascending-l
wavelength march for expanding rows; DLN_code = DLN_lambda.
U7 RESOLVED: no index-direction inconsistency found. The
remaining caveat is only the unverified a_l = GAM algebraic
identity (U2/U8), not the direction.

------------------------------------------------------------------------
B3. CORRECTED JACOBIAN EXPERIMENT - T6 REVISED (Task 3)
------------------------------------------------------------------------

Maps distinguished:

    G(x) = undamped one-pass ordered sweep (VNEW field before
           relaxation; mathematically the map whose fixed
           point the sweep seeks)
    H(x) = x + OMEGA_RAY * (G(x) - x)   (one damped pass)

The production sweep applies relaxation AND adaptive OMEGA_RAY
(RATIO_OM feedback, lines 10427-10471) plus optional P2
limiters. A finite-difference probe through the normal path
measures a state-dependent, non-autonomous map - not a clean
Jacobian of either G or H.

EXPERIMENT A - UNDAMPED MAP JACOBIAN J_G:
  - Hold PK_A0..A4, PK_SR, PK_Q, PK_DTAU fixed (per-refresh
    coefficients).
  - Hold IM_NI fixed.
  - Identical boundary data and sweep ordering for all
    perturbations (same KN_RAY_ORDERED path, same SPW preview
    behaviour).
  - Set OMEGA_RAY = 1.0 for the probe pass; disable P2
    (FG_ORD_P2DAMP = F); prevent the adaptive-relaxation
    update from running (or snapshot/restore OMEGA_RAY around
    the probe pass).
  - Evaluate exactly ONE complete ordered P sweep from x and
    from x + eps*d; record full FLD_P output.
  - J_G d ~ [G(x+eps d) - G(x)]/eps.

EXPERIMENT B - DAMPED ITERATION JACOBIAN J_H:
  - Same fixed coefficients, IM_NI, ordering.
  - Hold OMEGA_RAY fixed at a specified value (e.g. 0.45, the
    observed floor, and 1.0); do NOT let RATIO_OM adapt it
    during the probe.
  - Evaluate one complete damped sweep from x and x+eps d.
  - J_H d ~ [H(x+eps d) - H(x)]/eps.

Identity (when differentiable and the same fixed Omega is used):

    J_H = (1 - OMEGA_RAY) * I + OMEGA_RAY * J_G

Verification: compute J_H d via Experiment B and
(1-OMEGA) d + OMEGA (J_G d) via Experiment A on the same
direction d; agreement to finite-difference truncation level
validates both the implementation of the probe and the
linearity assumption locally.

Protocol for both:
  - At least two eps values (e.g. eps and eps/10 scaled to the
    field norm); require consistency between them.
  - Record which cells hit the VNEW<0->1e-30 floor, the
    NONMONO_INEXT_CAP, or the KN_ORD_DS sign/min clamps in
    each evaluation; a derivative measured across an
    active-set change is NOT a smooth Jacobian element and
    must be flagged, not averaged.
  - d should be a smooth probe direction (e.g. a low-order
    polynomial in (J,L), or the previous-pass update vector
    x_n - x_{n-1}) - not a single-cell spike, which measures
    one stencil row only.
  - Power iteration on J_G (or J_H) with the probe sweep as
    the matrix-free operator estimates the dominant
    eigen-direction; the Rayleigh quotient gives |lambda|.
    Non-normal growth is assessed separately by tracking
    ||delta_n|| over many probe iterations (Trefethen &
    Embree transient-growth diagnostic).

This is a DESIGN specification only; nothing is implemented.

------------------------------------------------------------------------
B4. CORRECTED PERIODICITY TELEMETRY - T4/A4 REVISED (Task 4)
------------------------------------------------------------------------

Fixed normalisation (replaces the per-pass VMAXF ambiguity):

  - At experiment start, save a reference per-frequency scale
        s_L = VMAXF(L) computed once from the initial FLD_P
    (or from a stored reference field). Use THE SAME s_L for
    every pass-to-pass comparison. Do not recompute VMAXF per
    pass for the norm (VMAXF(L) = max_J FLD(J,L) changes as
    the field evolves and would make norms incomparable).
  - Norm over a FIXED index set: all marched ray points
    J = 1..NI_83 and all frequencies L = 1..NCF (or the
    documented hot subset, stated explicitly).
  - Report BOTH:
        RMS:   ||x_n - x_{n-k}||_2 / sqrt(NI*NCF)  with each
               cell weighted by 1/s_L^2 (or unweighted, stated)
        MAX:   max_{J,L} |x_n-x_{n-k}| / s_L
  - Keep floor-dominated cells identifiable: flag cells where
    VOLD <= 1e-4*s_L (the residual denominator floor) or where
    the 1e-30 positivity floor / CAP was active; report the
    norm both including and excluding them.

Distances, for k = 1,2,3,4:

    d_k(n) = ||x_n - x_{n-k}||   (fixed s_L, fixed index set)

Periodicity classification: candidate period p requires
d_p(n) -> 0 (below a stated fraction of d_1's running mean)
while d_1(n) remains finite and bounded away from zero.
For a period-2 candidate additionally require, sustained
over >= 10 consecutive passes:
  - alternating update-vector signs: sign(VNEW - VOLD) at the
    argmax cells flips systematically between even/odd passes;
  - stable spatial-frequency support: the (J,L) set carrying
    the largest |x_n - x_{n-1}| is the same across passes
    (compare with the RES argmax at j in {1,2});
  - repeated field values at even/odd pass separation:
    x_{n} ~ x_{n-2} elementwise to within the stated norm,
    not merely the scalar residual.

Scalar RES history alone remains insufficient evidence for
any periodicity class.

------------------------------------------------------------------------
B5. CORRECTED INDEPENDENT REFERENCE SOLUTION - T9 REVISED (Task 5)
------------------------------------------------------------------------

The ordered sweep is NONLINEAR: KN_ORD_DS sign/min clamps,
the VNEW<0 -> 1e-30 positivity floor, the NONMONO_INEXT_CAP
upper cap, and state-dependent active sets. A dense direct
linear solve of G(x)-x=0 is therefore NOT automatically valid.

Reference formulation: solve the nonlinear fixed-point
problem

    F(x) = G(x) - x = 0

on the single-ray (or reduced) P system with:
  - the same PK coefficients (PK_A0..A4, PK_SR, PK_Q,
    PK_DTAU) from the same refresh;
  - the same frequency grid and ray geometry (NI_83, NCF,
    GAM, DLNM_L, SPW stencil);
  - the same IM_NI boundary and j=1 outflow edge;
  - the same donor interpolation and endpoint (KN_ORD_DS)
    treatment;
  - the same floor/cap rules - unless an explicitly separate
    LINEARISED reference (floors/caps removed, clamps frozen
    at their active-set values) is being studied for
    comparison.

 admissible methods, in order of preference:
  1. Nonlinear Newton or trust-region root solve on F(x)=0
     with matrix-free Jacobian-vector products from B3
     Experiment A (the same complete-sweep probes).
  2. Piecewise-linear active-set solve: iterate the linear
     solve, update the active set (floor/cap/clamp states),
     re-solve; converge when the active set is unchanged
     between iterations AND the residual is below tolerance.
  3. Dense direct solve of the linearised operator ONLY
     after demonstrating that (a) no floor/cap/clamp changes
     sign anywhere in the domain of interest and (b) the
     active operator is fixed - i.e. the problem is provably
     linear on the relevant set.
  4. Another independently verified convergent method that
     preserves the same discrete equations (e.g. a
     time-dependent relaxation to steady state with a
     contractive splitting, or a Krylov solve of the
     linearised system preconditioned by the sweep).

Validation of the reference:
  - ||F(x_ref)|| below a stated tolerance under the SAME
    residual norm used by the sweep;
  - cross-check between at least two of the methods above
    (agreement to discretisation/truncation level);
  - spot-check that the reference satisfies the boundary
    conditions and the discrete equation residual cell-by-
    cell, not just in norm.

Only with x_ref in hand can seed-sensitivity experiments
(T9 original) address causality: if the retained store is
far from x_ref AND the sweep from x_ref converges, the store
is upstream; if the sweep from x_ref also cycles, the map
itself is non-contractive. Different seeds alone still do
not prove the historical anomaly's origin.

------------------------------------------------------------------------
B6. U1-U7 STATUS UPDATE
------------------------------------------------------------------------

U1 (J_G spectrum / stability boundary): STILL OPEN. Now has a
    concrete measurement design (B3 Experiments A/B). No
    eigenvalue estimate exists yet.
U2 (algebraic identity code<->Paper II): PARTIALLY RESOLVED.
    The structural correspondence and the index mapping are
    established (B1, B2); the coefficient-level difference is
    quantified (<=0.1% on this grid). The term-by-term
    algebraic identity remains unresolved pending Paper I's
    a_l definition and Eq 7/9 text (now U8). Impact on the
    observed residual is bounded at ~0.1% - cannot explain
    O(0.6) residuals.
U3 (L=1 opacity-increment question): STILL OPEN, unchanged.
U4 (j=1 endpoint consistency with paper boundary closure):
    STILL OPEN, unchanged.
U5 (periodicity class): STILL OPEN, now with the corrected
    fixed-normalisation telemetry spec (B4).
U6 (store-anomaly causality): STILL OPEN, now with the
    nonlinear reference-solution requirements (B5).
U7 (index-direction algebra): RESOLVED (B2). No inconsistency;
    code L and paper l share ascending-wavelength order.

New unresolved questions:

U8. Does the code's discretisation match Knop Paper I
    (A&A 496, 295-298, DOI 10.1051/0004-6361:200810290) Eqs.
    7 and 9 EXACTLY? The code cites those equations for
    chi_hat = chi + |GAM|/dLOG_NU and signed Q. Paper I's
    text was not retrievable in this audit; if Eq 7/9 define
    the logarithmic form, the code-vs-Paper-II difference is
    a paper-to-paper formulation choice, not an
    approximation error. REQUIRED before closing U2.

U9. Is a_l identically GAM (including the 3.33564e-6
    geometrical factor and the (1+sigma*mu^2) term), or does
    Paper I/II define a_l with additional factors (e.g. the
    lambda weighting of Eq. 1's d(lambda I)/d lambda)? The
    sign and zero-set correspond; the magnitude identity is
    unverified. Affects only coefficient normalisation, not
    the structural conclusions.

U10. Does the <=0.1% wavelength-coefficient discrepancy have
    any dynamical amplification through the iteration map
    (e.g. near a stability boundary)? Bounded a priori by the
    coefficient ratio, but the map's sensitivity is exactly
    what B3 measures; no claim either way until then.

------------------------------------------------------------------------
B7. CONSOLIDATED STATUS AFTER 5.3B
------------------------------------------------------------------------

[V] verified code behaviour (new this phase):
- chi_hat = chi + |GAM|/DLN for all L>1 coupling cells;
  Q = GAM/(DLN*chi_hat) > 0 for expanding cells.
- DLN = ln(lambda_L/lambda_{L-1}) > 0; code L and paper l
  share ascending-wavelength order; expanding donor L-1 <->
  l-1 consistent.
- Grid statistics: DLN in [1.66e-5, 1.94e-3], median 1.67e-4
  on the m15_2p82d grid; coefficient discrepancy <= 0.1%.

[M] measured numerical behaviour: unchanged from 5.3A -
persistent bounded oscillatory residual plateau ~0.6-0.7,
argmax j in {1,2} outer edge, OMEGA floor 0.45, M converges
in 2 passes.

[H] hypothesised cause: M1-M7 from 5.3A stand, with M4
    (wavelength-coupling discretisation) now bounded at <=0.1%
    coefficient error and explicitly noted as unable to
    explain the O(0.6) residual by itself.

[D] demonstrated defect: NONE. The wavelength-coefficient
    difference is a quantified small approximation
    (classification B), not a demonstrated violation of the
    intended discretisation; whether it violates Paper I's
    intended form awaits U8.

No solver change, test implementation, or production
parameter change is authorised or implied by this addendum.
The audit awaits ChatGPT's review of the B1-B7 corrections.
