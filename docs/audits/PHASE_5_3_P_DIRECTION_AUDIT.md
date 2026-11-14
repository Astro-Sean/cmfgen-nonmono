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

========================================================================
PHASE 5.3C - PAPER I TRANSFER-EQUATION RECONCILIATION
========================================================================
Dated addendum to Phase 5.3B (commit 78fa7fec). Audit-only: no
Fortran, parameter, test, or production changes. Claim labels
unchanged: [V] verified code behaviour, [M] measured numerical
behaviour, [H] hypothesised cause, [D] demonstrated defect.
Paper I = Knop, Hauschildt & Baron (2009), A&A 496, 295-298,
DOI 10.1051/0004-6361:200810290, arXiv:0902.1876 (equation text
verified against the arXiv HTML rendering, Eqs. 1-12).
Paper II = Knop, Hauschildt & Baron (2009), A&A 501, 813-820,
DOI 10.1051/0004-6361/200911793, arXiv:0909.4264.

------------------------------------------------------------------------
C1. PAPER I EQUATIONS, RECONSTRUCTED
------------------------------------------------------------------------

Eq. 1 (Lagrangian transfer equation, all terms):

    dI_l/ds = eta_l - chi_l I_l - 4 a_l I_l
              - a_l d(lambda I_l)/d lambda

Eq. 2 (spherical wavelength-coupling coefficient):

    a_l = gamma ( gamma^2 mu (mu+beta) dbeta/dr
                  + beta (1-mu^2)/r )
    beta = v/c, gamma = 1/sqrt(1-beta^2), mu = cos(theta)

Eqs. 3-4 (wavelength finite difference, xi splitting), for
lambda_{l-1} < lambda_l < lambda_{l+1} and a_l >= 0:

    p_l^- = -lambda_{l-1}/(lambda_l-lambda_{l-1})
    p_l^| =  lambda_l    /(lambda_l-lambda_{l-1})
    p_l^+ =  0

    dI_l/ds = eta_l - chi_l I_l
              - a_l (4 + xi p_l^|) I_l
              - xi a_l (p_l^- I_{l-1} + p_l^+ I_{l+1})
              - [1-xi] a_l (p_l^- I_{l-1} + p_l^| I_l + p_l^+ I_{l+1})

Eq. 5 (OLD generalised opacity, shown to fail):
    chi_hat_old = chi_l + a_l (4 + xi p_l^|)

Eq. 7 (NEW generalised opacity - the point of Paper I):
    chi_hat = chi_l + xi a_l p_l^|
    The 4 a_l I_l term is REMOVED from the opacity and must be
    carried by the formal solution.

Eqs. 8-10 (formal solution in optical depth tau, d tau =
chi_hat ds):

    I_l(tau_2) = I_l(tau_1) e^{tau_1-tau_2}
                 + INT S_hat_l e^{tau-tau_2} dtau
                 - INT S_tilde_l e^{tau-tau_2} dtau

    S_hat_l   = (chi_l/chi_hat) ( S_l
                 - xi (a_l/chi_l)(p_l^- I_{l-1} + p_l^+ I_{l+1}) )

    S_tilde_l = (a_l/chi_hat) ( [1-xi](p_l^- I_{l-1}+p_l^+I_{l+1})
                 + [4 + [1-xi] p_l^|] I_l )

Note the MINUS sign on the S_tilde integral in Eq. 8: the
proportional term enters the formal-solution source with
negative sign. For xi=1, a_l >= 0:

    chi_hat   = chi_l + a_l p_l^|
    S_hat_l   = (chi_l/chi_hat) S_l + (a_l/chi_hat) |p_l^-| I_{l-1}
    S_tilde_l = 4 a_l I_l / chi_hat
    effective source in Eq. 8 = S_hat_l - S_tilde_l

Eqs. 11-12: piecewise linear/parabolic interpolation of
S_hat (alpha,beta,gamma weights, Olson & Kunasz 1987 /
Hauschildt 1992) and LINEAR-ONLY interpolation of S_tilde.

------------------------------------------------------------------------
C2. THE CODE'S ACTUAL CONTINUOUS EQUATION [V]
------------------------------------------------------------------------

This is the central finding of Phase 5.3C. The code does NOT
discretise Paper I Eq. 1. It solves CMFGEN's native
FREQUENCY-space comoving-frame equation, documented in the
source itself (fg_j_cmf_nonmono_v1.f, comment block
~lines 4240-4310, "edited by Sean Brennan 07-Sep-2026 /
09-Oct-2026 after external review"):

    dI/ds = eta - chi I + GAM * dI/d(ln nu)

with GAM the signed O(v/c) comoving frequency-drift
coefficient. Verified independently against the discrete
equation (not just the comment):

- KN_ORD_SPV (lines 10675-10720): expanding (GAM>0) gives
  SPV = PK_SR + PK_Q*FLD(J,L-1); compressing (GAM<0) gives
  SPV = PK_SR + |PK_Q|*INEXT_V with INEXT_V interpolated
  between depth points JPR1/JPR2 at L+NNF2.
- VNEW (lines ~11280-11286): I(J) = A0*I(J+1) + A1*S'(J+1)
  + A2*S'(J) - A3*DSU - A4*DSD, the two-point formal solution.
- Reassembling: chi_hat*S' = eta + (GAM/DLN) I_donor, so the
  cell equation is dI/ds = eta - chi I + GAM*(dI/dln nu)
  with the upwind difference quotient. For GAM>0 the donor
  is L-1 (higher nu, correct upwind for redshift); for
  GAM<0 the donor is L+NNF2 (lower nu). ONE equation, with
  direction-dependent upwinding - the signed-GAM formulation
  handles both senses uniformly. [V]

Consequences for the Task-2 question (the exp(-DLN) donor
coefficient): the feared missing factor DOES NOT EXIST because
the differenced variable is I itself, not lambda*I or nu*I.
The discretisation GAM*(I_up - I_loc)/DLN is the EXACT
upwind difference quotient of dI/dln nu on a logarithmic
frequency grid - no exponential weighting is required or
absent. The code's 1/DLN coefficients on both the diagonal
(|GAM|/DLN inside chi_hat) and the donor (Q = GAM/(DLN*chi_hat))
are the correct coefficients FOR THIS EQUATION.

Task-2 verdict: option C - the code implements a DIFFERENT
discretisation of a DIFFERENT (frequency-space) equation, not
an approximation of Paper I's wavelength-space equation.

------------------------------------------------------------------------
C3. GAM vs a_lambda [V + derivation]
------------------------------------------------------------------------

SIGMA definition (verified): SIGMA_TRUE = d ln|V|/d lnR - 1
= (R/V) dV/dR - 1, computed by RECOMPUTE_SIGMA_TRUE /
SIGMA_TRUE_CORE (new_main/subs/recompute_sigma_true.f,
comment line ~129: "d ln|V|/d lnR - 1"). V13 uses the same
convention (fg_j_cmf_v13.f line ~1078 comment: "SIGMA=
(dlnV/dlnR-1.0)").

Code coefficient (V in km/s, R in cm, 3.33564e-6 = 1/c in
s/km):

    GAM = (V/cR) (1 + sigma mu^2),   sigma = (r/V)dV/dr - 1

Expand for R = r:

    (V/r)(1 + sigma mu^2)
      = V/r + mu^2 (V/r)[(r/V)dV/dr - 1]
      = V/r + mu^2 dV/dr - mu^2 V/r
      = mu^2 dV/dr + (1-mu^2) V/r          ... (i)

Paper I Eq. 2 non-relativistic limit (gamma->1, beta->0):

    a_lambda * c = mu^2 dv/dr + (1-mu^2) v/r          ... (ii)

Equations (i) and (ii) are IDENTICAL. Therefore:

    GAM = a_lambda   to O(v/c)   [V]

The (1+sigma*mu^2) structure is not accidental: with
sigma = dlnV/dlnr - 1 it is EXACTLY the algebraic rearrangement
of the spherical frequency-gradient coefficient. Homologous
flow (V~r, sigma=0) leaves the pure geometric (1-mu^2)V/r
term; coasting flow (sigma=-1) leaves pure mu^2 dV/dr.

Relativistic content: Paper I keeps gamma, gamma^2, and the
(mu+beta) aberration inside the gradient term, plus gamma on
the transverse term. The code retains none of these - it is
the O(v/c) truncation, stated as such in the source comment
("the signed O(v/c) frequency-drift coefficient"). For
v ~ 10^4 km/s (beta ~ 0.03) the leading corrections are
O(beta) ~ 3% on a_lambda, O(beta^2) ~ 0.1% from gamma^2.
Sign conventions agree (both positive for expansion).
Geometrical factors: complete at O(v/c); the mu^2 weighting
of dV/dr and the (1-mu^2) transverse term are both present.
Units: consistent (1/length after the 1/c factor).

[dGAMdR note, INTEGRAL method only]: the nonmono file's
dGAMdR (lines ~3588-3610) is actually MORE complete than V13's:
it keeps V13's GAM*sigma/R and T1*mu^2*d sigma/dR terms AND adds
a ray-curvature term T1*2*sigma*(1-mu^2)/R for MU=Z/R varying
along the ray (added 10-Sep-2026, audit finding #15). This
feeds the Euler-Maclaurin DTAU correction, not the coupling
coefficient itself. [V]

------------------------------------------------------------------------
C4. THE 4 a_lambda I TERM - EXPLICIT ACCOUNTING (Task 4)
------------------------------------------------------------------------

Paper I Eq. 1's -4 a_l I_l has NO direct counterpart in the
code, and this is now EXPLAINED rather than merely noted:

- The code's equation is the frequency-space form. For I_nu
  the exact O(v/c) ray equation (per the source comment,
  citing Chen et al. 2007 Eq. 10: only I_nu/nu^3 is invariant)
  is dI/ds = eta - chi I - a_nu (dI/dln nu - 3 I), i.e. the
  nu-form proportional coefficient is 3, not 4 (the lambda-I
  vs nu-I weighting differs by one power). Written Knop-style
  the same coupling is +4 a_nu I - a_nu d(nu I)/dnu.
- CMFGEN's convention DROPS the proportional term entirely:
  the solved equation is dI/ds = eta - chi I + GAM dI/dln nu.
  This is the inherited V13 truncation (V13 uses the identical
  GAM construction, fg_j_cmf_v13.f ~lines 1083-1085), not an
  omission introduced by the KNOP path.
- The term is genuinely ABSENT from the ray solve: chi_hat
  contains only |GAM|/DLN (the implicit part of the
  derivative), S' = SR + Q I_donor contains no I_L self-term,
  and VNEW has no separate proportional contribution. It is
  not hidden in DCHR_L (that is d chi_hat/dR for the DTAU
  integral correction) or in PK_A0-A4 (those are the
  spatial formal-solution weights).
- Documented compensation: the comment asserts the dilution
  physics enters at the MOMENT level - MOM_J_DDT's
  (1/c r^3) D(r^3 J)/Dt advection contains 3 v J/(c r) (the
  exact sigma=0 counterpart) and MOM_JREL carries the
  relativistic coefficient structure. Spot-check: MOM_J_DDT
  files exist (mom_j_ddt_v1..v6) and the nonmono file
  references this division of labour at lines ~94, ~4291,
  ~4798. Full term-by-term verification of the moment-side
  balance is NOT performed here - recorded as U11.

Paper I Eq. 10's S_tilde = 4 a_l I_l/chi_hat therefore maps
to NOTHING in the code - not because it is transformed, but
because the code solves a different (truncated) equation.
If exact Chen/Knop consistency were required, the term to add
would be +3 a_nu I with the same sign-aware split (per the
source comment); it is deliberately not added because it
would change results vs stock V13.

Task-4 verdict: the 4aI term is absent BY DESIGN (documented
O(v/c) truncation shared with stock CMFGEN), not lost
through a variable transformation. Whether this truncation is
harmless for the KNOP non-monotonic regime is a physical
question the audit cannot settle - it is identical in form to
what V13 solves on monotonic models.

------------------------------------------------------------------------
C5. PAPER I vs PAPER II vs CODE - FORMULATION MAP (Task 5)
------------------------------------------------------------------------

| Paper I (formal solution)      | Paper II (matrix/GS)      | Code (KNOP nonmono)                    |
|-------------------------------|---------------------------|----------------------------------------|
| Eq.1 dI/ds wavelength form    | Eq.1 same equation        | DIFFERENT equation: nu-form            |
|                               |                           | dI/ds = eta-chiI+GAM dI/dln nu [V]      |
| p_l coefficients (linear-lam) | same p_l (Eqs.2-5 region) | 1/DLN log-nu differences [V]           |
| chi_hat = chi+xi a_l p_l^|    | chi_hat = chi+xi a_l p^|  | chi+|GAM|/DLN (xi=1-equivalent) [V] |
| S_hat donor (implicit, xi)    | matrix B off-diagonal     | Q*I_donor in S' [V]                    |
| S_tilde 4aI/chi_hat (expl.)   | Delta-I explicit terms    | ABSENT (3a_nu I dropped) [V]           |
| alpha/beta/gamma interp       | A/B/C matrix coeffs       | PK_A0..A4 parabolic spatial weights [V]|
| S_tilde linear-only interp    | (corner terms lagged)     | C-corner lagged via SPW ring [V]       |
| Lambda-operator construction  | (not in Paper II scope)   | refresh/re-solve architecture [V]      |
| boundary l=1 unspecified      | boundary treatment in GS  | expanding L=1: no donor/no increment   |
|                               |                           | (U3 still open) [V]                    |

What the code actually implements, precisely:

1. CMFGEN-native frequency-space transfer equation (V13
   truncation inherited): dI/ds = eta - chi I + GAM dI/dln nu.
2. Paper-I-MOTIVATED splitting: |GAM|/DLN into chi_hat
   (implicit, keeps opacity positive - the Paper I Eq. 7 idea
   applied to the nu-form coefficient), signed GAM/DLN*chi_hat
   as donor coefficient Q (the Paper I xi=1 upwind coupling).
3. Paper-II-MOTIVATED iteration: ordered Gauss-Seidel with
   physically guided upwind marching, lagged C-corner source
   (SPW ring = Paper II Eqs. 31-32 structure).
4. CMFGEN-native spatial discretisation: A0-A4 exponential-
   integral weights with parabolic source interpolation
   (KN_ORD_DS slopes), plus KN_DTAU_INT_V1 optical-depth
   integration and the dGAMdR Euler-Maclaurin correction.
5. CMFGEN extensions beyond both papers: MATCH_FREQ multi-bin
   donor (NNF2 offset) for compressing cells, JPR1/JPR2
   depth-interpolated INEXT_V, OMEGA_RAY adaptive relaxation,
   positivity floor and CAP.

The code is therefore NOT "Paper I with a different spacing"
and NOT "Paper II's matrix formulation" - it is a
CMFGEN-native equation solved with Paper-I-inspired opacity
splitting inside a Paper-II-inspired ordered iteration.

------------------------------------------------------------------------
C6. CORRECTED INTERPRETATION OF THE 0.1% COEFFICIENT RESULT
------------------------------------------------------------------------

Phase 5.3B compared code 1/DLN coefficients against Paper II's
p_l coefficients and found <=0.1% local discrepancy, then
stated this "cannot" explain the O(0.6) residual. That
statement was too strong and is corrected here:

- The 0.1% figure measured ONLY the spacing-form difference
  between 1/DLN and the p_l rational functions, ASSUMING
  a_l = GAM and the same underlying equation. C2-C4 show the
  underlying equations differ (nu-form vs lambda-form, plus
  the dropped proportional term). The comparison is therefore
  a LOWER BOUND on formulation difference, not the whole
  difference.
- Even within the spacing comparison: a 0.1% perturbation of
  chi_hat and Q is a perturbation of the sweep operator G.
  Its effect on convergence is governed by J_G, not by the
  coefficient size. Near a stability boundary (spectral
  radius ~ 1) or a defective dominant eigenvalue, eigenvalue
  shifts can be amplified as sqrt(perturbation) or worse;
  conversely a contractive map with margin absorbs 0.1%
  trivially. WHICH regime applies is exactly what the B3
  J_G/J_H experiments must measure.
- What CAN be stated: (a) the local coefficient discrepancy
  is small; (b) IF J_G is contractive with margin, it cannot
  be the cause; (c) IF J_G sits near marginal stability, even
  0.1% could shift the limit cycle quantitatively (though
  probably not create it from nothing - a structurally
  different equation is a larger perturbation than 0.1%).
- What CANNOT be stated without J_G/J_H: whether the
  coefficient difference affects the fixed point, the
  spectrum, or the convergence rate. U10 remains open and is
  now better posed: measure J_G's spectral radius and
  eigenvalue condition number, THEN bound the coefficient
  perturbation's dynamical effect.

------------------------------------------------------------------------
C7. U8/U9/U10 STATUS AND NEW QUESTIONS
------------------------------------------------------------------------

U8 (does code match Paper I Eq. 7/9?): RESOLVED - NO, and now
    understood WHY. The code's chi_hat = chi + |GAM|/DLN and
    Q = GAM/(DLN chi_hat) realise the Paper I Eq. 7 IDEA
    (positive opacity + signed donor coupling) but on the
    frequency-space equation, not Paper I's wavelength
    equation. Paper I's own Eq. 7 was written for the
    lambda-form; the code comment documents the nu-form
    adaptation explicitly.

U9 (is a_l = GAM?): RESOLVED at O(v/c): YES in the
    non-relativistic limit, EXACTLY, including geometrical
    factors (the sigma = dlnV/dlnr - 1 convention makes
    (V/cR)(1+sigma mu^2) algebraically identical to the
    non-rel limit of Paper I Eq. 2). Relativistic corrections
    (gamma, gamma^2, mu+beta, gamma on transverse term) are
    omitted - O(beta) ~ few % at SN velocities. Not an exact
    relativistic identity.

U10 (can 0.1% matter dynamically?): REFRAMED per C6 - open,
    pending J_G/J_H spectral measurements.

New unresolved questions:

U11. Moment-level consistency of the dropped proportional
     term: the comment claims dilution enters via MOM_J_DDT
     advection (3vJ/cr) and MOM_JREL coefficients. Verify
     term-by-term that the moment equations' velocity-gradient
     structure is consistent with the formal solver's
     truncated equation - a mismatch would make the coupled
     fixed point internally inconsistent. [H: possible
     formulation gap; needs dedicated trace]

U12. Does the nu-form equation with dropped 3a_nu I term give
     the same fixed point as the lambda-form with 4a_l I term,
     to O(v/c)? The two are related by dependent-variable
     transformation plus the truncation; the residual
     difference is O(beta) in the coupling. For the KNOP
     expanding-only ray the practical question is whether the
     truncated equation's fixed point differs measurably from
     the untruncated one - testable via the B5 reference
     solution with/without the term.

U13. Compressing-cell donor path (MATCH_FREQ NNF2 offset,
     JPR1/JPR2 interpolation, cross-iteration lag): this is a
     CMFGEN extension with no Paper I/II counterpart. Its
     consistency for strongly compressing cells is unverified
     (ray 83 is expanding-only so it does not affect the
     current failure, but it bounds how far the formalism
     generalises).

------------------------------------------------------------------------
C8. FINAL CLASSIFICATION (deliverable 4)
------------------------------------------------------------------------

The current expanding-flow implementation is:

  "A CMFGEN-native frequency-space discretisation of the
   comoving-frame transfer equation (with the inherited
   O(v/c) truncation of the proportional coupling term),
   using Paper-I-inspired positive-opacity/signed-donor
   splitting and Paper-II-inspired ordered Gauss-Seidel
   iteration."

It is NOT algebraically equivalent to Paper I Eqs. 1-10.
It is an alternative discretisation of a related but
distinct equation. Within its own equation the 1/DLN
coefficients are exact (no missing exponential factor).
The differences from Paper I are: (a) nu-form vs lambda-form
dependent variable and coupling coefficient; (b) dropped
3a_nu I proportional term (documented, shared with V13);
(c) O(v/c) truncation of a_lambda itself (no relativistic
factors); (d) CMFGEN extensions (MATCH_FREQ, adaptive
relaxation, floors/caps) with no paper counterpart.

None of (a)-(d) is demonstrated to cause the P-direction
non-convergence. The most consequential for fidelity is (b):
it is an O(v/c) ~ few-percent truncation of the coupling,
larger than the 0.1% spacing difference, and its effect on
the coupled fixed point is unquantified (U11/U12). But it is
identical to what stock CMFGEN solves on monotonic models,
so it is a formulation property, not a KNOP defect.

[V] new verified facts: code equation is nu-form (C2);
    sigma = dlnV/dlnr - 1 makes GAM = a_lambda at O(v/c)
    exactly (C3); proportional term absent by design (C4);
    nonmono dGAMdR more complete than V13 (C3 note).
[M] measured: unchanged (residual plateau ~0.6-0.7).
[H] hypotheses: M1-M7 stand; M4 further weakened - the
    coefficient question is now dominated by the equation-
    level difference (b), not the spacing difference.
[D] demonstrated defects: NONE. The formulation differences
    are documented design decisions, not violations of an
    intended discretisation - UNLESS U11 reveals a moment-
    side inconsistency, which would be a real defect.

No solver change, test implementation, or production
parameter change is authorised or implied. Awaiting
ChatGPT review of C1-C8 before any implementation.

------------------------------------------------------------------------
C9. PHASE 5.3C REVISION - U11/U12/U13 EVIDENCE (post-review)
------------------------------------------------------------------------
Added after ChatGPT's preliminary review requested the actual
derivations for U11-U13 and a model-based magnitude estimate.
Audit-only; still no solver/test/parameter changes.

C9.1 U11 - MOMENT-SIDE COMPENSATION: STRUCTURE VERIFIED
[V]

The highest-priority question was whether the dropped
proportional term is compensated at the moment level. Direct
inspection of the KNOP moment solver
(plane/mom_j_cmf_nonmono_v1.f) establishes:

- It uses the IDENTICAL Knop splitting as the formal solver:
    CHI_HAT_NM(I) = CHI_NM(I) + GAM_MAG_NM(I)/DLNM   (line 2273)
    GAM_NM(I)     = CON_GAM_NM(I)/CHI_HAT_NM(I)/DLNM (line 2305)
    Q_NM(I)       = GAM_NM(I)*GAM_FAC_NM(I)          (line 2313)
    PSI diagonal  = DTAUONQ*|GAM_MOM|/dLOG_NU         (line 2022)
  i.e. chi_hat = chi + |GAM_MOM|/dln nu and the same signed-Q
  donor structure. The formal solver and the moment solver
  therefore solve the SAME truncated equation - the dropped
  term is absent from BOTH, so the KNOP coupled system is
  internally consistent at the truncated-equation level. [V]

- The moment coupling coefficient is the angle-integrated
  ray coefficient: GAM_MOM = (V/cR)(1 + SIGMA_TRUE*K_ON_J)
  (module header line 15; assembly comment line 2002), where
  K_ON_J is the Eddington factor = angle-averaged mu^2
  (line 2004; clamped [1e-4,1], fallback 1/3 at line 1822).
  Since the ray coefficient is (V/cR)(1+sigma*mu^2), its
  J-moment average under Eddington closure is exactly
  (V/cR)(1+sigma*K_ON_J) - the decomposition is structurally
  consistent; the sigma*mu^2 anisotropic remainder feeds the
  H-moment (flux) coupling, the standard two-stream split.
  Term-by-term exactness (closure corrections beyond
  Eddington) is NOT verified - remains part of U11. [V for
  structure, H for exactness]

- No 3*GAM*J proportional term exists in the moment solver
  either (searched: no GAM-proportional diagonal term beyond
  the chi_hat/Q structure). The truncation is therefore
  UNIFORM across formal solver, KNOP moment solver, and the
  stock path (MOM_J_DDT_V6 uses CON_GAM=V/cR,
  CON_GAMH half-point, GAM=CON_GAM/CHI/dLOG_NU,
  PSI=DTAUONQ*GAM - same structure, no proportional term).
  The truncation is a stock-CMFGEN convention, not a KNOP
  defect. [V]

- Conservation-level caveat (ChatGPT point A): presence of
  the coupling terms is necessary but not sufficient. The
  exact ray equation's J-moment would contain an additional
  diagonal contribution from the dropped 3a_nu I term (of
  order 3*GAM_MOM*J after angle integration). The moment
  solver does not contain it - consistent with the truncated
  equation, but this means BOTH solvers share the same
  approximation rather than one compensating the other. The
  precise question - whether the moment equation's PSI/GAM
  coefficients are the exact moments of the truncated ray
  equation - requires deriving the moment system from that
  equation; specified as test T-U11 below, not yet done.

C9.2 U12 - MAGNITUDE OF THE TRUNCATION, MODEL-BASED
[V for numbers, H for interpretation]

ChatGPT correctly noted the earlier "few percent" was a
maximum-beta estimate, not a model-based error estimate.
Numbers now computed from the actual model
(cmfgen_validation/input_hydro.dat, 70 depth points,
t=1.667d/2.824d):

    V range:      49 .. 7074 km/s   -> beta_max = 0.0236
    R range:      1.455e14 .. 1.253e15 cm
    sigma range:  -157 .. -1.04     (all 70 cells sigma < -1;
                                    the model is almost entirely
                                    in the compressing regime for
                                    radial rays - ray 83 is a
                                    grazing expanding ray)
    V/(cR) max:   1.6e-16 /cm (base coefficient, before
                  the (1+sigma mu^2) factor)

Relativistic corrections to a_lambda (Paper I Eq. 2 gamma,
gamma^2, mu+beta factors): O(beta) ~ 2.4% worst-case,
O(beta^2) ~ 0.06%. [V for beta, H for coefficient impact]

Dropped-term magnitude. The missing term is -3*GAM*I for
expansion (equivalently an effective opacity +3GAM). Two
comparisons:

(a) vs the retained coupling opacity |GAM|/DLN: the relative
    change in chi_hat is 3*DLN ~ 5e-4 = 0.05% wherever the
    frequency coupling dominates the opacity budget. [V -
    algebra, DLN~1.7e-4 median]

(b) vs the true opacity chi: using chi_es ~ 0.4*rho as a
    LOWER bound on total chi (lines/bound-free add more):
        3*GAM_max / chi_es:  median ~ 0.2%,  max ~ 25%
    The maximum occurs in the outermost, most rarefied cells
    (r ~ 8e14 cm, V ~ 52 km/s, rho ~ 2e-18, sigma ~ -1.3)
    where chi_es ~ 3e-19 /cm - essentially transparent
    gas. In the interior where chi >> coupling the ratio is
    << 0.1%. [V for the ratio computation; the chi_es
    lower-bound caveat means the TRUE ratio is smaller
    wherever lines contribute]

Interpretation: the truncation is a smooth, sign-definite
(damping for expansion) opacity-level correction of order
0.05% of chi_hat in coupling-dominated cells and <= ~25%
only in the transparent outer edge where the radiation
free-streams anyway. This is small enough to be consistent
with stock CMFGEN's success, and it bounds the possible
fixed-point/spectral effect of U12 to the same order locally
- though, as in C6, iteration dynamics could amplify a small
operator perturbation; that requires the J_G/J_H
measurements, not this estimate.

C9.3 U12 - FIXED-POINT EQUIVALENCE: CONCRETE TEST
[H]

The question "does the truncated equation's fixed point
differ measurably from the untruncated one" is directly
testable within the B5 reference framework: solve
F(x) = G(x) - x = 0 on the single-ray P system twice -
once with the implemented operator, once with an added
+3*GAM*I_L diagonal term (the missing proportional
contribution, sign-aware per the expanding/compressing
split). The field difference ||x_trunc - x_full|| under the
fixed s_L normalisation of B4 quantifies the truncation's
dynamical effect INCLUDING any iteration amplification.
Cost: two reference solves; no production change. This is
now the sharpest version of U12.

C9.4 U13 - MATCH_FREQ COMPRESSING PATH
[V]

Confirmed CMFGEN extension with no Paper I/II counterpart:
for GAM<0 the donor is FLD at L+NNF2 where NNF2 =
CEILING(DRF_L/KN_DLN_N) can span multiple frequency bins
(lines 10693-10702), and INEXT_V is interpolated between
depth points JPR1/JPR2 (lines 10704-10713) with the value
lagged across great iterations. Physically motivated (strong
compression shifts the comoving frequency many grid cells
per ray segment) but mathematically unaudited against any
reference discretisation. Does NOT affect ray 83
(expanding-only, NNF2 path unused); bounds the formalism's
validity for compressing cells only. Status unchanged:
design extension, consistency unverified.

C9.5 REVISED U11-U13 STATUS

U11: PARTIALLY RESOLVED. The KNOP formal+moment system is
     verified to share the identical truncated equation
     (same chi_hat, same Q structure, K_ON_J-consistent
     angle average) - internally consistent at that level.
     The remaining open piece is the conservation-level
     derivation: are the moment equation's coefficients the
     exact moments of the truncated ray equation? That is a
     derivation/test, specified below, not yet executed.
U12: BOUNDED. Truncation effect <= ~0.05% of chi_hat in
     coupling-dominated cells, <= ~0.2% median / ~25% worst
     vs chi_es in transparent outer cells. Fixed-point
     consequence testable via the C9.3 two-solve experiment.
U13: OPEN, scoped to compressing cells only; irrelevant to
     the ray-83 expanding-only failure under audit.

C9.6 NEW TEST SPECIFICATION (design only, not implemented)

T-U11 (conservation-level moment check):
  Derive the J- and H-moment equations of the truncated ray
  equation dI/ds = eta - chi I + GAM dI/dln nu by explicit
  angular integration under the same closure the code uses
  (Eddington/K_ON_J). Compare coefficient-by-coefficient
  with the tridiagonal system assembled in
  mom_j_cmf_nonmono_v1.f (PSI, GAM_NM, DTAUONQ_NM, TA_NM,
  TB_NM, XM_NM). Acceptance: every coefficient matches to
  the closure's order; any residual 3a_nuI-type diagonal term
  is either shown absent-by-consistent-truncation or its
  magnitude bounded vs chi_hat. If a mismatch exceeding the
  truncation bound appears, U11 becomes a demonstrated [D].

T-U12 (truncation fixed-point test): as C9.3 - two B5
  reference solves with/without the +3GAM*I_L term; report
  ||dx|| under fixed normalisation and the residual of each.

No solver change, test implementation, or production
parameter change is authorised or implied by C9. Awaiting
ChatGPT review of C1-C9.


========================================================================
PHASE 5.3D - ANGULAR MOMENT CONSISTENCY AND TRUNCATION AUDIT
========================================================================
Date: audit-only continuation of C9. No Fortran source, solver
parameter, production control, or test was modified; nothing was
committed. All numerical work was done in an external Python
analysis script (/tmp/phase53d_gam_signs.py) that re-implements the
production stencils; no production .o was executed for these
numbers, so they are analysis-level [M], not unit-test [V].

The question C9.1 left open (T-U11): are the moment solver's
coefficients the angular moments of the discrete ray operator?
The suspicion raised for this phase: angular averaging does not
commute with the sign-dependent upwind frequency discretisation,
because for SIGMA < -1 different rays at the same depth carry
opposite signs of GAM(mu) = (V/cR)(1+sigma*mu^2) and therefore
use opposite frequency donors. This section resolves it.

------------------------------------------------------------------------
D1. ANGULAR QUADRATURE AND MOMENT DEFINITIONS [V]
------------------------------------------------------------------------

Ray/impact-parameter construction, subs/impar_with_pnt_srce.f
(NC_PNT_SRCE=0 branch, lines 24-29, 51-53):

    P(1)   = 0                              (exactly radial ray)
    P(I)   = RP*sqrt(1-((NC-I+1)/NC)^2)     I=2..NC   (core rays,
             equally spaced in mu at the outer edge; RP=R(ND))
    P(I)   = R(ND-I+NC+1)                   I=NC+1..NP (impact rays;
             impact parameters ARE the depth radii, reversed)

For this model ND=70, NC=15, NP=85 (MODEL_SPEC). Internal depth
convention confirmed at new_main/cmfgen_sub.f:1006-1012 and
new_main/subs/rv_sn_model_sniin.f:156-176: R(1)=RMAX=1.2531e5
(outermost), R(ND)=RP=1.4547e4 (innermost); R decreases with I.
The hydro file lists inner->outer and is reversed on read.
P is monotonic in ray index (verified numerically).

Angular weights, subs/nordangqw_v2.f + subs/set_acc_dmu.f +
subs/jweight_v2.f (called from new_main/mod_subs/set_ang_qw_v2.f):

    At depth I: NW = NC+ND-I+1 rays contribute (cap NP). The
    NW-th ray has P=R(I) exactly - the tangent ray, mu=0.
    mu(J,I) = sqrt(R(I)^2-P(J)^2)/R(I).
    dmu via SET_ACC_dMU (binomial-expansion fix near mu=1);
    weights via JWEIGHT_V2/HWEIGHT_V2/KWEIGHT_V2/NWEIGHT_V2,
    a Nordlund modified-cubic rule normalised to sum 1.
    => JQW/HQW/KQW/NQW are (depth,ray)-dependent; the mu-grid
    itself changes with depth because the contributing ray set
    changes. The quadrature always spans mu in [0,~1].

Moment accumulation, plane/fg_j_cmf_nonmono_v1.f:5762-5767:

    JNU_STORE(I) += JQW(I,LS)*AV(K,LS)
    KNU_STORE(I) += KQW(I,LS)*AV(K,LS)
    AV(K,LS) = 0.5*(I_P + I_M)   (mean of the two streams)

Closure factors, plane/get_moms_non_rel.f:

    K_ON_J(I) = KNU_STORE(I)/JNU_STORE(I)   (line 107)
              = <mu^2>_I, the K/J Eddington factor; fallback
              1/3 if JNU degenerate; clamped to [1e-4,1] in the
              caller (mom_j_cmf_nonmono_v1.f:1821-1823).
    NMID_ON_J(I)    = RSQ_NNU_MID/(RSQ_JNU(I)+RSQ_JNU(I+1))
    NMID_ON_HMID(I) = RSQ_NNU_MID/RSQ_HNU_MID   (N_TYPE dependent;
                      = N/H for the H-equation angle factor)

So K_ON_J is the intensity-weighted second angular moment - NOT
an unweighted mu^2 quadrature mean. For near-isotropic intensity
it tends to ~1/3; it is a genuine closure quantity.

------------------------------------------------------------------------
D2. THE DISCRETE RAY OPERATOR AND ITS SIGN STRUCTURE [V]
------------------------------------------------------------------------

Per ray (direction cosine mu = Z/R along the characteristic),
plane/fg_j_cmf_nonmono_v1.f:3563-3566, 11990-12059, 10675-10720:

    GAM_ray  = (V/cR) * (1 + SIGMA*mu^2)          signed
    GAM_MAG  = |GAM_ray|
    chi_hat  = chi + |GAM_ray|/DLN                (KN_RF_COUPLING gate)
    Q_ray    = (GAM_ray/DLN)/chi_hat              signed
    S'       = eta/chi_hat + Q_ray*I_upwind
    donor    = FLD(L-1)          if GAM_ray > 0   (bluer neighbour)
             = FLD(L+NNF2..)     if GAM_ray < 0   (redder neighbour,
                                                    possibly multi-bin
                                                    with JPR1/JPR2
                                                    interpolation)

Sign of GAM_ray = sign(1+sigma*mu^2) since V,c,R>0. For sigma<-1
the sign changes at mu*^2 = -1/sigma = 1/|sigma|, i.e. at impact
parameter p/R = sqrt(1-1/|sigma|): grazing rays (large p/R, small
mu) are positive, radial-ish rays negative.

Frequency-index direction (Phase 5.3B): L increases toward red
(wavelength increases with L). Formal donor: GAM>0 -> L-1 (bluer),
GAM<0 -> later L (redder). Moment donor (D4 below): Q_NM>=0 ->
RSQJNU_PREV_NM (previous/bluer frequency, lagged);
Q_NM<0 -> J_NEXT_FRQ(FREQ_INDX+MOM_NM_J) (redder). FREQ_INDX=1 is
the bluest frequency (mom_j_cmf_nonmono_v1.f:~2040). The donor
DIRECTIONS agree for a given sign; the issue is that the two
solvers can assign different signs at the same depth.

------------------------------------------------------------------------
D3. MODEL MEASUREMENT: MIXED-SIGN GAM OVER THE ACTUAL QUADRATURE [M]
------------------------------------------------------------------------

sigma_true was recomputed with a line-by-line port of
SIGMA_TRUE_CORE (new_main/subs/recompute_sigma_true.f:119-288,
log-quadratic interior stencil, one-sided quadratic boundaries,
linear fallback on stencil sign change). Result: sigma in
[-156.9, 0.0], 69 of 70 depths sigma<-1 (one boundary depth has
sigma~0). This matches the earlier measurement.

For each depth I the contributing ray set J=1..NW(I), mu(J,I),
JQW weights, and f(J)=1+sigma_I*mu(J)^2 were evaluated exactly as
the production routines define them. Results (weighted):

    depth  sigma    w(GAM>0)  w(GAM<0)  sumW|f|  |1+sigma/3|  ratio
      1   -1.04      0.983     0.017     0.649      0.654      0.99
     11   -1.40      0.841     0.159     0.594      0.535      1.11
     31   -1.79      0.750     0.250     0.594      0.402      1.48
     41   -2.18      0.680     0.320     0.630      0.273      2.31
     45   -2.54      0.614     0.386     0.683      0.155      4.42
     46   -2.71      0.591     0.409     0.713      0.098      7.29
     47   -2.93      0.564     0.436     0.756      0.023     33.1
     48   -3.18      0.531     0.469     0.807      0.061     13.2
     51   -3.82      0.471     0.529     0.954      0.275      3.47
     56   -4.84      0.500     0.500     1.230      0.614      2.00
     61   -5.09      0.458     0.542     1.289      0.698      1.85
     65  -156.9      0.072     0.928    51.40      51.30      1.00
     70   -31.0      0.167     0.833     9.56       9.32      1.03

Facts:
- 69 of 70 depths have BOTH signs present with non-trivial weight.
  There is no depth where a single-sign approximation is exact.
- The positive-GAM (grazing/expanding-coupling) weight fraction
  falls monotonically from ~0.98 at the outer edge (sigma~-1) to
  ~0.07-0.17 in the deep interior (sigma~-30..-157). The two signs
  are closest to balanced (w+ ~ 0.5) around sigma ~ -3 to -5,
  depths I ~ 44-62, R ~ 1.6e4-4.9e4.
- sumW|f| = sum_J JQW*|1+sigma*mu^2| is the J-weighted mean of
  the per-ray coupling magnitude; |1+sigma*<mu^2>| with
  <mu^2>_JQW = 0.3333 (the J-quadrature second moment, essentially
  the isotropic value at every depth) is the abs-of-mean. Their
  ratio is 1.00 where one sign dominates and peaks at ~33 near
  sigma ~ -3 where the signed mean nearly cancels while both
  signs retain O(1) weight.
- Ray 83 (used in the Phase 5.2 diagnostic) has P(83)=R(3)
  ~ 1.15e5 - a large-impact grazing ray, mu small at most depths,
  hence positive GAM over most of the grid. It is NOT
  representative of the angular distribution: at sigma~-3 depths
  it sits in the positive subset while ~44% of the weight is
  negative.

------------------------------------------------------------------------
D4. THE MOMENT SYSTEM ASSEMBLY [V]
------------------------------------------------------------------------

plane/mom_j_cmf_nonmono_v1.f, traced coefficient by coefficient:

    CON_GAM_NM(I)  = 3.33564e-6 * V_NM(I)/R_NM(I)         (line 2090)
    CON_GAMH_NM(I) = 2*3.33564e-6*(V(I)+V(I+1))/(R(I)+R(I+1))
    GAM_FAC_NM(I)  = sign-preserving regularisation of
                     (1 + SIGMA_TRUE(I)*K_ON_J(I)),
                     floor MOM_GAM_FAC_MIN = 1.0e-2        (2222-2224)
    GAM_SIGNED_NM   = CON_GAM_NM*GAM_FAC_NM                (2225)
    GAM_MAG_NM      = |GAM_SIGNED_NM|                     (2226)
    CHI_HAT_NM(I)   = CHI_NM(I) + GAM_MAG_NM(I)/MOM_DLN_J(I)
                                                           (2273)
    GAM_NM(I)       = CON_GAM_NM/CHI_HAT_NM/DLNM          (2305)
                     (overwrites the physical signed value;
                      signed copy kept in GAM_SIGNED_NM)
    Q_NM(I)         = GAM_NM(I)*GAM_FAC_NM(I)
                     = (GAM_MOM/DLN)/CHI_HAT_NM            (2313)
    PSI_NM(I)       = DTAUONQ_NM(I)*|Q_NM(I)|             (2418)
    PSIPREV_NM(I)   = DTAUONQ_NM(I)*|GAM_NM*GAM_FAC_PREV|  (2419)
    DTAUONQ_NM(I)   = 0.5*(DTAU(I)+DTAU(I-1))/QF_NM(I)     (2397)
    TA_NM(I)        = CHI_HAT_NM(I)*QF_NM(I)              (2384)
    TB_NM           = DTAUONQ*max(MOM_ABS_FLR,1-COH)
                      + PSI + HL + HU  (+/- 2*EPS terms,
                      depth-upwind form for N_TYPE/=G_ONLY) (~2500)
    TC/TA           = -HU/-HL  (spatial Eddington-factor coupling)
    VB/VC           = -HS(I-1)/HS(I)  (lagged-H coupling)
    XM_NM(I)        = DTAUONQ_NM(I)*SOURCE_NM(I)*R_NM(I)^2 (~2525)

    Donor selection (line ~2664):
      Q_NM >= 0 (expanding): J_upwind = RSQJNU_PREV_NM (lagged
        previous/bluer frequency field)
      Q_NM <  0 (compressing): J_upwind =
        J_NEXT_FRQ(I, FREQ_INDX+MOM_NM_J(I)) - redder neighbour,
        with MOM_NM_J>=1 multi-bin offset when MOM_FREQ_MATCH
        and |GAM|*ds/DLN exceeds one bin (Ceulemans et al. 2025
        fix; default OFF).

    H-equation coefficients:
      GAMH_NM(I) = CON_GAMH_NM/(DLNM*(CHI_HAT(I)+CHI_HAT(I+1)))
      W_NM       = GAMH_NM*|GAMH_FAC_NM|   (always >=0)
      EPS_NM     = GAMH_NM*AV_SIGMA_TRUE*NMID_ON_J/(1+W_NM)
      HU/HL      = K_ON_J*QF/((1+W)*DTAU) at I+1/I
      HS_NM      = WPREV/(1+W), overwritten to W/(1+W) at
                   compressing midpoints

    Midpoint angle factors use AV_SIGMA_TRUE=0.5*(sigma_I+sigma_I+1)
    and NMID_ON_HMID=N/H (or the G_ONLY clamp [0.01,1]) - a
    DIFFERENT closure factor than K_ON_J, so the J- and H-channel
    zero-crossings of (1+sigma*g) occur at different depths.

------------------------------------------------------------------------
D5. COEFFICIENT-BY-COEFFICIENT COMPARISON [V structure, M magnitude]
------------------------------------------------------------------------

Write the per-ray coupling coefficient as
    a(mu) = (V/cR)*(1 + sigma*mu^2)          (linear in mu^2)
    f(mu) = 1 + sigma*mu^2                    (sign carrier)

(A) SIGNED J-moment coefficient.
    Exact J-weight:  <a>_JQW = (V/cR)*(1 + sigma*<mu^2>_JQW)
    with <mu^2>_JQW = sum JQW*mu^2 = 0.3333 at every depth.
    Moment code:     GAM_MOM = (V/cR)*(1 + sigma*K_ON_J).
    For isotropic intensity K_ON_J -> ~1/3 and the two agree.
    Because a(mu) is LINEAR in mu^2, <a>_shape = a(<mu^2>_shape)
    exactly under the fixed-angular-shape closure, and
    <mu^2>_shape = K/J = K_ON_J by construction. The residual is
    the closure itself: K_ON_J is the stored K/J ratio (can deviate
    from 1/3 in anisotropic regions), and the shape is assumed
    frequency-independent across the donor stencil. Verdict: B -
    closure-level approximation, structurally consistent. [V/B]

(B) ABSOLUTE-VALUE OPACITY (the chi_hat diagonal).
    Per-ray:        chi_hat_ray = chi + (V/cR)*|f(mu)|/DLN.
    J-weighted mean of the coupling increment:
                    (V/cR)*sum_J JQW*|f(mu_J)|/DLN.
    Moment code:    chi + |GAM_MOM|/DLN = chi +
                    (V/cR)*|1+sigma*K_ON_J|/DLN.
    |<a>| != <|a|> whenever f changes sign over the quadrature.
    Measured ratio <|f|>/|<f>|_JQW (table D3, K_ON_J=1/3 proxy):
    1.00 at the extremes, rising through 1.5-2.3 across
    sigma ~ -2..-5, peaking ~33 near sigma ~ -3 (depths 44-52),
    where |1+sigma*K_ON_J| is floored only at MOM_GAM_FAC_MIN=0.01.
    At those depths the moment J-equation's frequency-coupling
    opacity is reduced to a few percent of the ray-consistent
    mean while ~50% of the angular weight still carries positive
    coupling and ~50% negative. Verdict: C - structural mismatch;
    a single signed coefficient cannot represent sum_J w|f_J|.
    Not a demonstrated wrong answer: the moment equation is an
    independent discretisation, and whether its |<a>| convention
    or the ray <|a|> convention is "correct" for the intended
    moment equation is a formulation question (D6/U11).

(C) DONOR DIRECTION.
    Exact J-donor sum (per frequency, per depth):
        sum_{f>0} w*a(mu)*I_{L-1}(mu) + sum_{f<0} w*a(mu)*I_{L+}(mu)
    Moment code: single donor, direction = sign(GAM_MOM),
    field = J_PREV or J_NEXT (single lagged field).
    Even under the fixed-shape closure with I_donor ~ J_donor*g(mu),
    the exact sum needs the TWO partition moments
        <a*1_{f>0}> and <a*1_{f<0}>
    which are not determined by <mu^2> alone - they require the
    mu^2 distribution (equivalently higher angular moments).
    K_ON_J carries only <mu^2>; one GAM_MOM and one donor cannot
    encode both partitions. For smooth spectra the error is
    O(DLN*|d ln I/d ln nu|*partition asymmetry) because
    I_{L-1} ~ I_{L+1}; near edges/lines it is larger. Verdict: C -
    structural; magnitude field-dependent, unquantified without
    the actual I(mu,nu). [C structure / H magnitude]

(D) H-CHANNEL.
    Same construction with HQW weights, midpoint factor
    1+AV_SIGMA*NMID_ON_HMID (N/H closure), W=|GAMH|*|GAMH_FAC|.
    Identical structural class; the sign partition differs because
    HQW != JQW and NMID_ON_HMID != K_ON_J. Not separately
    quantified. [V structure]

(E) What is NOT a discrepancy.
    - chi_hat positivity: both use |coupling|/DLN; chi_hat>0 in
      both. [V]
    - Donor direction convention per sign: consistent (D2). [V]
    - Frequency indexing, multi-bin upwind (MOM_NM_J vs NNF2):
      same physical intent, independent parameters. [V]
    - Boundaries: moment tridiagonal inner/outer rows vs ray
      boundary intensities - different formulations by design;
      not compared term-by-term here. [open]

------------------------------------------------------------------------
D6. U11 RECLASSIFICATION [V structure, M magnitude, H consequence]
------------------------------------------------------------------------

Phase 5.3C/C9.1 claimed the moment coefficient is "the exact
J-moment average under Eddington closure". That claim is now
REFINED:

- TRUE for the signed coefficient <a> = a(K_ON_J): a(mu) linear
  in mu^2 makes the closure exact at fixed angular shape. [V/B]
- FALSE for the |.| opacity coefficient: <|a|> != |<a>| under
  mixed signs; measured discrepancy up to ~33x near sigma ~ -3.
  [M]
- FALSE for the donor: single sign/donor cannot represent the
  sign-partitioned sum; requires higher angular moments than
  K_ON_J provides. [V structure]

U11 classification: C - STRUCTURAL MISMATCH. The moment solver
averages GAM (signed) and takes one absolute value and one donor
direction; the ray operator has sign-dependent donors and
per-ray |GAM|. This is not a closure-order difference; it is a
missing degree of freedom (the sign partition of the angular
weight). It is NOT classified D because: (i) the moment equation
is an independent formulation of the angle-integrated problem,
not literally the angular mean of the ray discretisation;
(ii) no wrong answer has been demonstrated - the consequence on
the coupled fixed point requires T-U11/T-U12 experiments.
The sharpest statement: at depths where 1+sigma*K_ON_J ~ 0
(sigma ~ -3 for K_ON_J ~ 1/3, depths ~44-52 here) the KNOP
moment solver's frequency coupling nearly switches off while
the formal ray operator retains O(1) coupling on every ray.
Whether that produces a measurable J/fixed-point error is the
open quantitative question.

Relation to stock CMFGEN: the |.|/single-coefficient construction
is exact in the monotonic regime (sigma > -1 => f(mu)>0 for all
mu => <|f|>=|<f>| and one donor). The mismatch is a genuine
non-monotonic-regime defect class that the inherited stock
formulation never had to handle.

------------------------------------------------------------------------
D7. U12 - SIGNED +3*GAM TERM: DERIVATION AND TEST DESIGN [V/H]
------------------------------------------------------------------------

Derivation (convention check): Paper I nu-form full equation
contains -a_nu*(dI/dln nu - 3 I). The code's truncated equation
is dI/ds = eta - chi I + GAM dI/dln nu, so GAM = -a_nu and the
retained coupling is -a_nu dI/dln nu = +GAM dI/dln nu. The
dropped proportional piece is +3 a_nu I = -3 GAM I on the RHS;
moving it to the implicit side gives an effective diagonal
    chi_eff = chi + 3*GAM        SIGNED, per ray (GAM=GAM(mu))
and per moment row (GAM=GAM_MOM).
For GAM<0 this reduces the diagonal; for GAM>0 it increases it.
It must NOT be replaced by 3|GAM| - the sign is the physics
(damping vs anti-damping along the characteristic).

Combined with the existing upwind split, the reference
operator per ray is
    chi_hat_ref = chi + |GAM|/DLN + 3*GAM
              = chi + |GAM|*(1/DLN + 3*sign(GAM)).
Positivity: 1/DLN >= ~500 >> 3 for all grid spacings here
(DLN <= 1.94e-3), so chi_hat_ref > 0 in both branches; the
correction is a ~3*DLN ~ 0.05-0.6% perturbation of the coupling
part of chi_hat. For the moment solver the analogous reference
is chi + |GAM_MOM|/DLN + 3*GAM_MOM.

Controlled comparison design (analysis-only; no production
change authorised or implied):

    A. Existing truncated operator (current code path).
    B. Reference: identical except chi_hat -> chi_hat + 3*GAM
       (signed) in BOTH the formal sweep coefficients (PK_A*)
       and, for the coupled experiment, the moment CHI_HAT_NM.

    Fixed: model, frequency grid, ray, boundary conditions,
    spatial discretisation, source interpolation, coefficients,
    convergence tolerance, initial field.

    Step 1 (single ray, P-direction): pick one ray+depth column,
    fixed coefficients/boundaries, iterate the sweep map to its
    fixed point under A and B. Report: convergence status, true
    residual ||G(x)-x|| for each, pass count, full-field
    difference, B4-fixed-normalisation difference, location and
    magnitude of max difference, and whether chi_hat_ref stays
    positive / how the tridiagonal conditioning changes.
    Note: ray 83 is expanding-dominated; a compressing ray
    (small mu, e.g. a core ray) should be included because the
    signed term changes sign there.

    Step 2 (coupled): only if step 1 shows a non-negligible
    field difference. Couple formal+moment fixed points under
    A and B; a single-ray result does NOT establish the coupled
    radiation/population fixed-point change.

    Integrated omitted contribution: the local proxy
    3|GAM|/chi_hat <= 3*DLN is NOT an integrated-impact measure.
    Along a characteristic the omitted term accumulates as
    integral 3*GAM*I ds while the retained coupling accumulates
    as integral (|GAM|/DLN)*I_upwind-type terms; the ratio per
    unit optical depth is ~3*DLN*(bins crossed) locally but the
    cumulative effect on the emergent intensity requires the
    reference solve - reported as an open measurement, not
    inferred.

------------------------------------------------------------------------
D8. C9.2 MAGNITUDE REASSESSMENT [V algebra, H interpretation]
------------------------------------------------------------------------

Retained algebraic bound (CHI >= 0):
    3*|GAM| / (CHI + |GAM|/DLN) <= 3*DLN        (~5e-4 median grid)
This bounds the signed term relative to chi_HAT only. It does
NOT bound: (i) the ratio vs electron-scattering opacity chi_es
(local ratios up to ~0.25 in transparent outer cells, C9.2b);
(ii) the integrated correction along a characteristic;
(iii) any fixed-point, Jacobian-spectrum, convergence-rate, or
stability change - all of which require J_G/J_H measurements
and the D7 reference solves.

New from D3-D5 relevant to magnitude honesty: near sigma ~ -3
the moment solver's |GAM_MOM| underestimates the ray-consistent
mean |GAM| by up to ~33x. For U12 the signed 3*GAM term inherits
the SAME GAM convention, so per-ray the relative size vs
chi_hat is still ~3*DLN; but in the moment solver near the
quench depth the omitted 3*GAM_MOM term is itself quenched
(|GAM_MOM|~0.02*(V/cR)) while the ray-level 3*GAM(mu) terms
remain O(V/cR) on ~half the rays. The single-coefficient
representation therefore also misestimates the SIZE of the
omitted term in the mixed-sign zone - same structural root as
U11, not a separate defect.

------------------------------------------------------------------------
D9. REMAINING UNCERTAINTIES [H]
------------------------------------------------------------------------

1. K_ON_J vs <mu^2>_JQW: K_ON_J=KNU_STORE/JNU_STORE is
   intensity-weighted; the JQW-weighted mu^2 (=0.3333) is the
   quadrature moment of a constant field. In anisotropic
   regions they differ; the signed-coefficient closure error is
   second-order in anisotropy but unquantified here.
2. Donor-direction error magnitude needs the actual I(mu,nu)
   fields (JNU_STORE per ray is not archived per frequency in
   the audited outputs); estimate is formal only.
3. H-channel partition (HQW, N/H factor) not separately
   computed; expect same class, different numbers.
4. Whether the moment equation is INTENDED to reproduce the
   angular mean of the ray operator, or is an independent
   consistent discretisation of the angle-integrated equation,
   needs the Paper II moment formulation checked against
   mom_j_cmf_nonmono_v1.f term by term (rt_literature). If the
   latter, U11-C is a formulation gap to document, not a bug;
   if the former, it is a coefficient defect.
5. Boundary rows of the tridiagonal (inner/outer) were not
   compared term-by-term with the ray boundary treatment.
6. All D3 numbers are analysis-level [M] from a Python port;
   they should be reproduced inside a unit-test driver against
   production .o before being treated as [V] (per AGENTS.md,
   unit tests must use the production objects).

STATUS SUMMARY
  U11: C - structural mismatch (single signed GAM_MOM + single
       donor + |<a>| cannot represent sign-partitioned ray
       operator). Magnitude measured: opacity-coefficient
       discrepancy up to ~33x near sigma~-3; donor-direction
       error unquantified. Not D: no wrong answer demonstrated.
  U12: signed +3*GAM derived and verified ([V]); local bound
       3*DLN retained ([V]); integrated/fixed-point effect open,
       test designed (D7), not run.
  C9.2: bound vs chi_hat retained; explicitly NOT a bound on
       fixed-point error; new caveat that the same structural
       issue quenches the moment-side estimate of the omitted
       term near sigma~-3.

No solver change, test implementation, or production parameter
change is authorised or implied by Phase 5.3D. Awaiting
independent review.


========================================================================
PHASE 5.3E - PRODUCTION-LEVEL VALIDATION OF THE MIXED-SIGN ANGULAR
            COUPLING
========================================================================

MODE: verification only. No production source, parameter, or solver
change was made. All drivers and analysis scripts live outside the
tree (/tmp). Field provenance and scoping caveats are stated in
section E8 before the numbers are interpreted.

LABELS: [V] verified against production source or compiled objects;
[M] measured numerical result; [H] hypothesis; [D] demonstrated
defect.

------------------------------------------------------------------------
E1. PRODUCTION GEOMETRY DRIVER [V]
------------------------------------------------------------------------

Driver /tmp/drv_angqw.f90 reads the FGINEXT_GS header records
(MAGIC,VERSION,TAG / N_COMP,NP,NCF,ND / R,V,SIGMA /
COMP_DEPTH_IDX / maturity scalars / KN_FREQ) and then calls the
PRODUCTION compiled routines, linked from lib_dbg:

  SIGMA_TRUE_CORE   (new_main/subs/recompute_sigma_true.f)
  IMPAR_WITH_PNT_SRCE (subs/impar_with_pnt_srce.f)
  NORDANGQW_V2      (subs/nordangqw_v2.f)
  SET_ACC_dMU       (subs/set_acc_dmu.f)
  JWEIGHT_V2 / KWEIGHT_V2 / HWEIGHT_V2 / NWEIGHT_V2
                    (subs/*weight_v2.f)

Input: run_dir = cmfgen_validation/m15_2p82d_KNOP_p8_base/
cmfgen/2.82days (DO_NONMONO_KNOP=T). Header: ND=70, NC=15, NP=85,
N_COMP=70, NCF=56377. R(1)=1.253121e5 (outer), R(ND)=1.454688e4
(inner); P(1)=0, P(NP)=R(1).

Measured [M]:
  SIGMA_TRUE_CORE: N_UNDEF=0, min=-1.5689E+02, max=-1.0370E+00.
  max|SGC-SGF| = 0  (recomputed sigma_true equals the stored SIGMA
  record bitwise).
  JQW row sums = 1.0000000000 at all 70 depths (machine precision).
  mu ordering: mu(ray1)=1 (radial), mu(NW)=0 (tangent, P(NW)=R(I)).

Per-depth sign census (production weights, f=1+sigma*mu^2):
  depths with BOTH GAM signs: 70/70  [V against production objects]
  The Phase 5.3D port reported 69/70; the difference is a PORT
  BUG, not a production feature: the port's sigma_true stencil
  returns 0 at depth 69 (Fortran I=69, second-to-last), where
  production gives sigma=-31.7469. With the correct sigma that
  depth is also mixed-sign (mu_crit=0.177).

Representative rows (full table reproducible from
drv_angqw_out.bin):

   I    sigma    NW  npos  nneg  wpos    mu_crit   R_geom(iso)
   1   -1.037    85    57    28  0.983    0.982     1.001
   8   -1.330    78    45    33  0.892    0.867     ~2
  16   -1.476    70    35    35  0.828    0.823     ~1.2
  29   -1.740    57    21    36  0.750    0.758     ~1.5
  36   -1.957    50    16    34  0.702    0.715     ~1.7
  45   -2.536    41     8    33  0.60     0.628     4.42
  46   -2.707    40     7    33  0.59     0.608     7.29
  47   -2.929    39     6    33  0.58     0.584    33.11  <- peak
  48   -3.18     38     5    33  0.56     0.561    13.20
  51   -3.824    35     3    32  0.47     0.511     3.47
  57   -4.92     29     2    27  0.50     0.451     ~1.9
  64  -10.62     22    10    12  0.09     0.307     ~1.2
  66 -125.43     20     5    15  0.09     0.089     ~1.00
  70  -30.96     16     3    13  0.17     0.180     ~1.03

R_geom = sum_rays JQW*|1+sigma*mu^2| / |1+sigma*<mu^2>_JQW| with
<mu^2>_JQW = sum JQW*mu^2 = 0.3333 (isotropic second moment).
Peak 33.11 at depth 47 (sigma=-2.93). Reproduces the Phase 5.3D
~33x geometric result at the same depth [V].

Port-vs-production agreement [M]:
  max|R_prod-R_port| = 1.5e-11; max|V| = 2.7e-11;
  max|P_prod-P_port| = 1.5e-11.
  sigma: identical except depth-69 port bug above.
  weights: per-depth max|w_port-w_prod| ~ 4e-3..1e-2 typical
  (port jweight is an approximation of JWEIGHT_V2); at depth 1 the
  port takes the wrong X(N)==0 boundary branch (last-ray weight
  -0.023 vs production -1.98; both sum to 1 and low moments agree).
  Weighted moments (npos/nneg/wpos/R_geom) agree to ~0.01-0.4;
  all sign-census counts identical.

------------------------------------------------------------------------
E2. INTENSITY-WEIGHTED COEFFICIENTS [M]
------------------------------------------------------------------------

Field: I_P_DEPTH_NEXT / I_M_DEPTH_NEXT records 9-10 of FGINEXT_GS
from the same run (checkpoint at ~iteration 240 of a still-running
KNOP job, model convergence ~50%). These are the production
lagged GS stores: per (comp-depth KC, ray LS, freq) the relaxed
intensity T2+INEXT_RELAX*(T1-T2) written at fg_j_cmf_nonmono_v1.f
lines 5592/5646, floored at 1e-30 and capped at NONMONO_INEXT_CAP.
AV(KC,LS,f) = 0.5*(I_P+I_M) is exactly the quantity JNU_STORE/
KNU_STORE accumulate (lines 5762-5767): JQW(I,LS)*AV(J_PNT(I,LS),LS)
and KQW(I,LS)*AV. COMP_DEPTH_IDX is the identity map here
(N_COMP=70, all depths compressing).

Accumulation replicated in Python with production weights:
  J(I,f)   = sum_{r=1..NW(I)} JQW(I,r) * AV(I,r,f)
  K_ON_J   = sum KQW(I,r) * AV / J          [production def]
  KJ       = sum JQW*mu^2*AV / J            [JQW-based <mu^2>]
  GS       = sum JQW*(1+sigma*mu^2)*AV / J = 1+sigma*KJ
  GA       = sum JQW*|1+sigma*mu^2|*AV / J
  GAM_signed = con*GS,  GAM_abs = con*GA,  con=3.33564e-6*V/R
  GAM_MOM    = con*(1+sigma*K_ON_J)
  R_abs      = GA/|GS|   (near-zero |GS| flagged, not hidden)

Field sanity [M]: I_P max 0.104, median 3.5e-8; I_M max 0.103;
no negative or NaN samples; outermost comp-depth I_M sits at the
1e-30 floor (inward boundary). J>0 for depths 1-63 on ~all freqs;
deep interior (64-70) populated on the few reaching rays.

K_ON_J vs KJ [M]: the two quadratures of the same mu^2 moment agree
to med|KJ-KQ|=7.9e-6, p99=2e-4. They diverge only where |1+sigma*K|
itself nearly cancels (max rel diff in the factor ~700 at
|1+sigma*K|~1e-3), i.e. inside the quench layer - a real but
cell-localized effect, not a global inconsistency.

GAM_signed vs GAM_MOM [M]: identical up to the KJ-vs-KQ quadrature
difference above; they agree wherever the signed coefficient is
O(1) and both shrink together inside the quench layer. Acceptance
"signed coefficients agree to numerical precision" is satisfied at
the level the two production quadratures allow (~1e-4 typical).

R_abs = GA/|GS| vs depth [M]:
  overall (|GS|>1e-6, 3.94M cells): median 2.18, p90 14.6, p99 144
  per-depth median: ~1.4-2.0 outside depths ~29-51; peaks
    I=36: 8.6   I=41: 10.6   I=45: 12.2   I=46: 20.3   I=47: 14.6
    I=50: 4.1   I=57: 1.9   I=64: 1.16  I=70: 1.03
  per-depth p90 reaches ~150 at I=45.
  Near-cancellation cells (|GS|<1e-3): 14936 of 3.94M; worst
  I=62 f=1325: GS=1.1e-4, GA=1.35, J=3.6e-21 (optically thin).

Positive/negative GAM intensity fractions (Gpos/GA) vs depth [M]:
  I=1: 0.85  I=8: 0.21  I=16: 0.23  I=29: 0.39  I=45: 0.54
  I=46: 0.51  I=57: 0.24  I=64: 0.07  I=70: 0.01
  Both populations carry non-trivial intensity weight through the
  mixed-sign zone; the negative-GAM population dominates the deep
  interior as sigma -> -100+.

Absolute coefficient scales [M]: |GAM_MOM| median ranges
3.8e-11 (I=1) to 4.9e-5 (I=65) 1/cm; <|GAM_r|> median is larger by
exactly R_abs. The coupling opacity is |GAM|/DLN, so for DLN~1e-3
these are ~1e-8..1e-2 1/cm - comparable to or above the continuum
chi in line-forming regions. Per-frequency CHI is not archived in
FGINEXT_GS, so the CHI-dilution of CHI_HAT=CHI+|GAM|/DLN is
unquantified here [H].

------------------------------------------------------------------------
E3. FREQUENCY-DONOR PARTITION [M]
------------------------------------------------------------------------

Controlled comparison at adjacent-bin offset (NM=1 stencil; the
production MATCH_FREQ rematch NM_RAY>1 and the moment MOM_NM_J are
separate offset choices - see E8). For every (depth, freq) with
production AV:

  exact donor numerator  Dex = sum_pos w*GAM_r*AV(l-1)
                        + sum_neg w*GAM_r*AV(l+1)
  moment donor          Dmom = GAM_MOM * J_donor
                        (J(l-1) if GAM_MOM>0 else J(l+1))

Result over 3.95M cells [M]:
  rel=|Dex-Dmom|/|Dmom|: median 0.003, p90 0.32, p99 6.0
  per-depth p90 peaks 0.65-0.74 at depths 36-46 (sigma~-2..-2.7)
  cells rel>10: 25044, but the worst cells have |Dex|,|Dmom|
  ~1e-21..1e-40 (negligible absolute terms in the thin region);
  restricting to |Dmom| above the noise the discrepancy is a
  moderate tail effect concentrated in the same sigma~-3 zone.

Interpretation [M]: the single moment donor is usually adequate
(median 0.3%) but cannot represent the two-population sum; the
error concentrates exactly where the pos/neg partition is balanced
and where GAM_MOM partially cancels - consistent with the Phase
5.3D class-C structural finding, now measured on the production
field.

------------------------------------------------------------------------
E4. TRACE INTO THE MOMENT MATRIX [V source / M coefficients]
------------------------------------------------------------------------

mom_j_cmf_nonmono_v1.f assembly (line numbers):

  T1 = 1+SIGMA_TRUE_NM*K_ON_J
  GAM_FAC_NM = SIGN(MAX(|T1|,MOM_GAM_FAC_MIN),T1),  floor 1e-2 (L2223)
  GAM_SIGNED_NM = CON_GAM_NM*GAM_FAC_NM              (L2225)
  GAM_MAG_NM   = |GAM_SIGNED_NM|                     (L2226)
  HAS_COMPRESSING set when GAM_SIGNED<0               (L2231)
  MOM_NM_J(I) = matched freq offset for compressing   (L2249-2271)
  CHI_HAT_NM = CHI_NM + GAM_MAG_NM/MAX(MOM_DLN_J,MIN) (L2273)
  GAM_NM = CON_GAM_NM/CHI_HAT_NM/DLN                 (L2305)
  Q_NM   = GAM_NM*GAM_FAC_NM   (signed)               (L2313)
  PSI    = DTAUONQ*|GAM*GAM_FAC|  (implicit, >=0)     (L2399+)
  PSIPREV= DTAUONQ*|GAM*GAM_FAC_PREV| (explicit)      (L2400+)
  K_ON_J = KNU_STORE/JNU_STORE from GET_MOMS_NON_REL  (get_moms_non_rel.f:107)

Coefficient-by-coefficient vs the angularly integrated ray operator:

  signed factor (1+sigma*K_ON_J): closure-consistent. a(mu)=1+sigma*mu^2
    is linear in mu^2 so <a>_I = a(<mu^2>_I); K_ON_J=KNU/JNU is that
    moment to quadrature accuracy (E2). Class B.

  absolute opacity |GAM|: CHI_HAT uses |GAM_MOM|=con*|1+sigma*K_ON_J|
    (floored at con*0.01), while the ray operator carries per-ray
    |GAM_r|=con*|1+sigma*mu_r^2| whose J-mean is con*GA. The two
    differ by R_abs (E2): the moment omits the distinct
    absolute-opacity contribution of the minority-sign population.
    Class C - structural mismatch, measured.

  donor: Q_NM routes a single signed donor (direction by
    GAM_SIGNED sign, offset MOM_NM_J) whereas the ray operator needs
    both populations (E3). Class C.

  PSI/PSIPREV inherit the same single-coefficient |GAM_MOM|, so the
    implicit/explicit frequency-coupling weights carry the same
    omission - no additional independent inconsistency found in the
    assembly itself (the floor at 0.01 keeps |GAM_FAC| nonzero inside
    the quench layer, which is a declared regularization, not an
    error).

Representative cells: this model has NO single-sign depth (all
sigma<-1). Near-single-sign limits used instead: I=1 (posW=0.85,
sigma=-1.04, nearly all-positive) and I=66 (posW=0.09, sigma=-125,
nearly all-negative); mixed-sign cells I=45-47 (sigma~-2.5..-2.9,
R_geom peak zone). At I=45-47: K_ON_J~0.35-0.37, 1+sigma*K_ON_J
small but nonzero, GA/|GS| median ~12-20 with p90 ~35-150 - the
assembled CHI_HAT/Q/PSI coefficients there underestimate the
ray-consistent coupling opacity by that factor on the GAM term.

U11 classification: remains C (structural mismatch). Now
production-verified in form and measured in magnitude; still not D -
no converged incorrect answer has been demonstrated, and the moment
path does not converge on m15 anyway (see AGENTS.md 22-Oct-2026
entry). The mismatch is a formulation gap between the
angle-integrated moment closure and the sign-partitioned ray
operator, not a demonstrated numerical defect.

------------------------------------------------------------------------
E5. SENSITIVITY TO THE ACTUAL RADIATION FIELD [M]
------------------------------------------------------------------------

R_geom (isotropic J-weight) vs R_abs (intensity-weighted):
  geometric peak 33.1 at I=47 (sigma=-2.93); >3 only at I=44-51.
  intensity-weighted median ~2 globally, local median up to ~20,
  p90 ~150, in the SAME depth band (sigma~-1.7..-3.8, I~29-51).
  Outside that band R_abs~1-2: where the field is nearly isotropic
  (K_ON_J~0.33, deep interior) or single-sign-dominated (outer
  I=1), the single-coefficient representation is adequate.

Concentration: the discrepancy lives where (a) sigma is near the
-3 quench (|1+sigma*K_ON_J| small), (b) the radiation field still
carries weight on BOTH sign populations (posW~0.4-0.6 at I=36-51),
and (c) J is non-negligible. It is NOT concentrated in the
optically thin outer layers (I<10, R_abs~1.4-1.8) nor the deep
isotropic interior (I>60, R_abs~1.0-1.6). Photosphere-region depths
(~30-55) dominate.

Absolute scales argue the effect is not negligible where it occurs:
<GAM>~1e-10..1e-5 1/cm against chi~1e-5..1e-4 1/cm, amplified by
1/DLN into CHI_HAT. But the CHI_HAT fractional error is diluted by
CHI - without per-frequency CHI the fixed-point significance is
unquantified [H]. Per the audit rule, no physical significance is
inferred from the ratio alone.

------------------------------------------------------------------------
E6. VERIFICATION TEST SUITE [V/M]
------------------------------------------------------------------------

Fortran driver /tmp/drv_angqw_test.f90 - PRODUCTION objects
(lib_dbg/libsubs.a etc.), synthetic grids. All PASS:
  T10 SIGMA_TRUE_CORE reproduces power-law sigma=-0.5 to 4e-14
  T9  P monotone, P(1)=0, P(NP)=R(1)
  T1  JQW sums to 1 at all 20 depths
  T2  mu(1)=1, tangent ray P(NW)=R(I)
  T3  NW(I)=min(NC+ND-I+1,NP)
  T4  all-positive sigma: |sum w f| = sum w|f| exactly
  T5  triangle inequality holds; near-threshold sigma=-1.05 ratio~1.00
  T6  mixed sigma=-3: sum|f|>|sum f|, ratio 2976 on synthetic grid
  T7  ratio finite at min|sum w f|
  T8  <mu^2>_JQW = 0.3334 ~ 1/3
  (note: an all-negative-GAM stencil is impossible - the tangent
   ray mu=0 always has f=+1; single-sign only exists for sigma>-1.
   Covered by T4 plus T5/T6 partition checks.)

Python suite /tmp/phase53e_tests.py - production weights + synthetic
intensities. All PASS (0 failures):
  T6.1 all-pos identity; T6.2/6.3 mixed triangle+partition;
  T6.4 near-zero ratio finite (1.0e4 on real weights);
  T6.5 isotropic KJ=KQ=1/3 to 7e-6; T6.6 anisotropic K reproduces
  the intensity-weighted moment and exceeds 1/3;
  T6.7 donor partition: pos->l-1, neg->l+1 exact, swap differs;
  T6.8 J normalization; T6.9 production regression (wsum=1,
  sigma range).

The suites detect: angular indexing errors (T2/T6.7), wrong weights
(T1/T6.9a), ABS(GAM_MOM) substitution (T6.2/T6.4 - would break the
triangle/near-zero checks), donor sign reversal (T6.7c), freq
off-by-one (T6.7 structure), lost intensity weighting (T6.5/6.6),
normalisation (T6.8/6.9a).

------------------------------------------------------------------------
E7. ANSWERS TO THE PHASE 5.3E QUESTIONS [V/M]
------------------------------------------------------------------------

  70/70 mixed-sign depths reproduced from production objects [V]
  (corrects the port's 69/70 boundary bug).
  signed intensity-weighted coefficient agrees with GAM_MOM to the
  quadrature level (~1e-4) wherever O(1) [M].
  absolute intensity-weighted coefficient independently measured
  [M]: R_abs median ~2, local median ~20, p90 ~150 in the quench
  band.
  both donor populations evaluated with actual donor intensities
  [M]: median rel error 0.3%, tail concentrated at sigma~-3.
  moment matrix discrepancy quantified coefficient-by-coefficient
  [V+M]: CHI_HAT opacity and Q/PSI donor are class C.
  tests cover single-sign, mixed-sign, anisotropic cases [V].
  no production solver changes made.
  U11 remains C.

  The ~33x result:
    - reproduced geometrically: YES, 33.11 at depth 47 [V]
    - reproduced with actual intensity weighting: as a DISTRIBUTION,
      not a uniform factor - median ~2 globally, local median ~20,
      p90 ~150, same depth band [M]
    - present in the assembled production moment operator: YES in
      form - CHI_HAT uses |GAM_MOM|=con*|1+sigma*K_ON_J| (floor 0.01)
      vs the ray-consistent con*GA; quantified per cell by R_abs
      [V+M]

------------------------------------------------------------------------
E8. CAVEATS AND UNRESOLVED QUESTIONS
------------------------------------------------------------------------

1. Field provenance: the intensity field is the KNOP FGINEXT_GS
   I_NEXT store - a real production field but a LAGGED, ~50%
   converged KNOP iterate, not a converged solution. Under
   DO_NONMONO_SOB+USE_LAM_ES=F (the path that actually runs
   MOM_J_CMF_NONMONO_V1) the K_ON_J input comes from the V13
   formal solve on CLAMPED sigma - a different field. The structural
   conclusions are field-independent in form; the R_abs numbers are
   specific to this field and should be re-measured on a converged
   SOB-path field when one exists.
2. MATCH_FREQ: production donor offsets NM_RAY (per ray-point,
   CEILING(|GAM|*dZ/DLN)) and moment offsets MOM_NM_J were NOT
   reproduced; the donor comparison used the adjacent-bin stencil
   to isolate the angular partition. A matched-offset comparison is
   future work.
3. CHI per frequency is not archived; CHI_HAT=CHI+|GAM|/DLN assembly
   was traced symbolically but not numerically assembled. The GAM-
   term discrepancy is measured; its fractional effect on CHI_HAT
   depends on the unknown CHI/|GAM|/DLN ratio [H].
4. H-channel (HQW, H-moment) partition not separately computed;
   expected same class.
5. Boundary rows of the tridiagonal moment system not compared
   term-by-term with ray boundary treatment.
6. The moment equation may be an independent consistent
   discretisation of the angle-integrated equation rather than an
   angular mean of the ray operator (Phase 5.3D uncertainty 4,
   still open - needs the Paper II moment formulation check).
7. Python-port sigma bug at depth 69 documented in E1; port weights
   are approximations (~5e-3) with a wrong boundary branch at
   depth 1. Production objects are the reference for all [V] claims.

STATUS SUMMARY
  U11: C - structural mismatch, now production-verified in form
       and measured on the production field. Magnitude: geometric
       33x at sigma~-2.9; intensity-weighted local median ~20,
       p90 ~150 in depths ~29-51; donor partition error tail in the
       same band. Not D: no converged incorrect answer demonstrated.
  U12: unchanged from 5.3D (signed +3*GAM derived [V]; test designed).
  The mixed-sign angular coupling is a REAL, MEASURED structural
  feature of the non-monotonic moment closure - quantified, localized
  to the sigma~-3 quench band, and bounded in the bulk (median ~2).

No solver change, test implementation in-tree, or production
parameter change is authorised or implied by Phase 5.3E. Awaiting
independent review before any commit.


================================================================================
PHASE 5.3F - QUANTIFY THE PHYSICAL CONSEQUENCE OF THE MIXED-SIGN MOMENT
MISMATCH (AUDIT AND REFERENCE-EXPERIMENT DESIGN ONLY)
================================================================================

Mode: audit + experiment design. No production Fortran, parameter,
OMEGA_RAY, or KNOP-moment-solver change. No commit. All analysis
artifacts live in /tmp (drv_reduce53f.f90, /tmp/red53f/*.npz,
/tmp/angqw.npz); nothing added to the tree.

--------------------------------------------------------------------------------
F1. CHECKPOINT PROVENANCE (Task 1 prerequisite)
--------------------------------------------------------------------------------

The FGINEXT_GS record layout (VER=3) was re-walked record-by-record:
rec1 magic/ver/tag, rec2 dims (N_COMP,NP,NCF,ND), rec3-5 R,V,SGF
(560 B each for ND=70), rec6 COMP_DEPTH_IDX (280 B), rec7 metrics
(48 B: MAXCHG f8, capped i4, NDIV/CMMREJ/DIVHEAL i4, SWMAX f8,
NCONV i8, RF_DCHG f8), rec8 KN_FREQ (451016 B = 56377 f8), rec9/rec10
the two giant intensity records, rec11 -MAGIC. The >2 GiB records use
nonstandard markers (ifort); native Fortran READ handles them, the
Python record walker does not - offsets were established empirically:
I_P data at byte 453120, I_M data at IP0+NEL*8+16 (verified: IM cell
matches the earlier stream dump bit-for-bit; IP cells agree to ~5e-4
relative, i.e. one dump cycle of field evolution).

IMPORTANT PROVENANCE FACT: the p8_base FGINEXT_GS file is being
REWRITTEN by a live run (mtime advanced 12:44 -> 13:20 during this
phase). Every checkpoint read is a snapshot of an evolving lagged GS
store. Direct J recomputation from the current file reproduces the
Phase 5.3E stream dump to ~5e-4, but a J-weighted sum over the same
cells read ~45 min earlier (the drv_reduce53f pass at 12:44) differs
by a factor ~2-3 that drifts slowly with frequency (2.85 at
F=20001 -> 2.04 at F=50001). The stored field is NOT stationary
between dump cycles; ratio-based diagnostics (R_abs, sign partitions)
are unaffected by overall rescaling but spectral-shape evolution
between generations is real. All 5.3E/F numbers are internally
consistent per-generation snapshots, not a converged solution.

Checkpoint inventory (tag = NITSF generation tag in rec1; metrics
from rec7). These are SEPARATE RUN DIRECTORIES, not consecutive
iterations of one run:

  run                      tag   maxchg    swmax    nconv
  m15_2p82d_FRESH_knop     124   1.807e-2  0.211    0
  _ckpt_iso_seed2          124   1.807e-2  0.211    0   (same gen)
  m15_2p82d_KNOP_iso_M     157   5.565e-3  0.229    0
  m15_2p82d_KNOP_p4_qargmax 159  2.460e-2  0.714    14
  m15_2p82d_KNOP_p6_sibdiag 159  2.460e-2  0.714    14  (p4 sibling)
  m15_2p82d_KNOP_p5_jtdiag 161   6.825e-2  1.0e4    19.5e6 (diverged diag)
  m15_2p82d_KNOP_p7_census 162   5.558e-3  0.229    0
  m15_2p82d_KNOP_neut2_QF  246   5.541e-3  0.239    0
  m15_2p82d_KNOP_neut2_QF2 246   5.541e-3  0.239    0   (same gen)
  m15_2p82d_KNOP_neut2_QC  356   5.536e-3  0.267    0
  m15_2p82d_KNOP_p8_base   256   5.540e-3  0.267    0
  m15_2p82d_KNOP_p8_probe  257   6.894e-1  0.244    0
  saw_restart_src/adp2     694   1.0       0.0      0   (SOB warm start)
  saw_restart_src/hlr      750   1.0       0.0      0   (SOB warm start)

swmax ~0.21-0.27 and maxchg ~5.5e-3 are the steady immature-field
plateau seen throughout the KNOP audits; nconv=0 means the per-cell
50%-churn gate saw no unconverged-significant cells at dump time -
consistent with a field that is smooth but NOT the coupled fixed
point (the run log for p8_base reports "the moment equations are not
solved" under USE_LAM_ES=T and immature-field iteration ~98-255).

--------------------------------------------------------------------------------
F2. TASK 1 RESULT - IS THE MISMATCH STABLE AS THE FIELD CONVERGES?
--------------------------------------------------------------------------------

Method: /tmp/drv_reduce53f.f90 (compiled against lib_dbg production
objects) recomputes SIGMA_TRUE_CORE, IMPAR_WITH_PNT_SRCE, the four
NORDANGQW_V2 weight sets, and accumulates the ten JQW-weighted
intensity sums per (depth,freq) exactly as JNU_STORE does
(fg_j_cmf_nonmono_v1.f:5758-5767). Identical definitions at every
checkpoint: GA=sum w|f|I/sum wI, GS=sum wfI/sum wI, GP/GN the
sign partitions, K_ON_J=SKQ/SJ, GAMs/a/m=CONV*{GS,GA,1+sig*K_ON_J},
Rabs=GAMa/|GAMs|. Restricted to the sigma~-3 band
(-5<SGC<-1.5, depths ~29-56) and cells with |GS|>0.05*GA.

Result (median Rabs / p90 over the band):

  run                     Rabs_med  p90     |GAMm|med   GAMa_med
  FRESH_knop / iso_M /    3.00      23.6    4.67e-10   2.84e-09
  neut2_QC/QF/QF2 /
  p4 / p6 / p7 / p8_base /
  p8_probe / jfnk_cp / fz (11 runs, identical to <0.5%)
  p5_jtdiag               3.50      16.9    5.58e-10   2.77e-09
  smoke_knop_mom          2.27      13.0    8.08e-10   2.73e-09

Eleven of fourteen checkpoints give Rabs_med=3.00, p90=23.6 to three
significant figures despite tags spanning 124-357 and maxchg spanning
5.5e-3 to 6.9e-1. The two outliers are diagnostic/smoke runs with
different stored fields. Interpretation [M]: the intensity-weighted
angular partition is set by the ray geometry plus a stable
intensity-pattern morphology; it does NOT drift with iteration count
across the sampled states. This is evidence the mismatch is a
STRUCTURAL property of the closure, not an artifact of an
unconverged field - but it is NOT evidence about the converged
fixed point, because none of these fields is converged (F1).
The correct statement: the mismatch coefficient is stable across
the observed pre-convergence trajectory; whether it persists at
the true coupled fixed point is exactly what the F5 experiment must
measure.

--------------------------------------------------------------------------------
F3. TASK 2 RESULT - THE EXACT COMPRESSING DONOR PATH
--------------------------------------------------------------------------------

[V] Production path, traced end to end:

Formal solver (per ray point I, ray LS, frequency L):
  DRF_L(I)  = |GAM(I,LS)| * |dZ_local|          (ray-coordinate step)
  MATCH_FREQ(L): T1=DRF_L/KN_DLN_N(L);
    NNF = NCF-L if T1>=NCF-L else max(1,CEILING(T1)); capped;
    DLNM_L = LOG(KN_FREQ(L)/KN_FREQ(L+NNF)) when bracketed.
  Donor read (SRC_PRIME ~fg_j:12221-12306 and KN_ORD_SPV
  ~fg_j:10668-10720): expanding -> FLD(I,L-1); compressing ->
  JPR1/JPR2/WINX depth interpolation of FLD(.,L+NFM_L).
  Store-side equivalent (main sweep ~fg_j:5108-5173):
  I_*_DEPTH_NEXT(KC,LS,FREQ_IDX_ITER+NM_RAY(I)) with KC1/KC2/W_INEXT
  bracketing from COMP_DEPTH_IDX. The two implementations select the
  same physical donor (frequency k+NMF, own-depth at stored cells).

Moment solver (per model depth I, frequency FREQ_INDX):
  MOM_DLN_J = dLOG_NU default; compressing+MOM_NU_NEXT ->
  LOG(FREQ/MOM_NU_NEXT) (~mom:2245-2271); +MOM_FREQ_MATCH ->
  N_CAP=min(NCF_IN-FREQ_INDX, max(1,CEILING(T4))),
  T4=GAM_MAG_NM*|DELTA_R_model|/MOM_DLN_J, span re-logged via
  MOM_NU_GRID (~mom:2249-2270; H-channel twin at 2180-2205).
  Donor: JUP_NM(I)=J_NEXT_FRQ(I,FREQ_INDX+MOM_NM_J(I)) or the
  JNEXT_RF_READ refresh snapshot (~mom:2632-2666); expanding reads
  PSIPREV_NM*RSQJNU_PREV_NM (~mom:2711-2731). Single coefficient,
  single donor per depth - no angular partition.

Both solvers run with DO_NONMONO_FREQ_MATCH=T: under DO_NONMONO_KNOP
the key is FORCED T, overriding even explicit F
(rd_control_variables.f:3638-3650). VADAT of p8_base confirms
DO_NONMONO_KNOP=T, USE_LAM_ES=T; FREQ_MATCH is not a user key there.

KEY STRUCTURAL DIFFERENCE [V]: the formal matched offset NMF is
computed PER RAY POINT from |GAM(mu)| and the ray-coordinate step
|dZ|; the moment offset N_CAP is computed ONCE PER DEPTH from
|GAM_MOM|=|con*(1+sigma*K_ON_J)| and the MODEL radial step |dR|.
In the sigma~-3 quench band |GAM_MOM| sits at the MOM_GAM_FAC_MIN=0.01
floor while individual rays carry |1+sigma*mu^2| of order unity, so
the moment offset is systematically smaller than the ray offsets.

[M] Quantified on the production geometry (frequency-independent
quantities only - no intensity needed; dZ estimated by model-depth
Z-spacing sqrt(R^2-P^2), an UPPER bound on ray-point spacing since
the ray grid subdivides; the dR variant isolates the coefficient
effect at identical spacing):

  compressing ray-depth cells: 2120 of 5950 (F<0)
  fraction with NMF>1 over the sampled frequency set: 16.4%
  ray NMF percentiles (compressing): p50=1, p90=4, p99=92, max=1239
  moment N_CAP percentiles:      p50=1, p90=4, p99=81, max=249

  depth  sig     nC   rayNMF(dZ) med/p90/max   rayNMF(dR) med/p90/max  momNCap med/max
   41   -2.18   33      1/1/1                    1/1/1                  1/1
   47   -2.93   33      1/1/1                    1/1/1                  1/1
   51   -3.82   32      1/3/3                    1/3/3                  1/1
   56   -4.84   28      3/22/26                  2/18/25                1/8
   61   -5.09   22     10/82/133                 6/56/116               2/22
   66 -125.43   15     21/192/263               11/133/247              9/83
  (depths 6-31 all 1/1/1; depth 70 -31: 1/6/7 vs 1/5/7 vs 1/3)

  rayNMF(dZ)/momNCap over (depth,freq) cells: p10=p25=p50=1.0,
  p75=1.0, p90=2.67, p99=4.33.

Interpretation: (a) in the sigma~-3 quench band - where the U11
|GAM| opacity mismatch peaks - the matched offset is ~1 on BOTH
sides, so the donor-offset discrepancy is SECONDARY there; the
primary U11 effect remains the |GAM| coefficient and the
positive/negative donor partition. (b) In the deep interior
(sigma <-~ -5, depths ~56-66) matched offsets reach hundreds of
frequency bins and the single moment N_CAP underestimates the ray
distribution by up to ~4x at p90 even at identical spacing - a
second structural mismatch of the same family (single coefficient
cannot represent the angular spread of Doppler drift rates). This
is consistent with, and extends, the U11 classification; it is NOT
a new demonstrated defect [C].

Caveats: NMF(dZ) is an upper bound (model-depth spacing >= ray-point
spacing where the ray grid refines); KN_FREQ bracket recomputation
was not applied per-cell (uniform-DLN estimate used, exact only up
to the log-spacing correction); the K_ON_J entering N_CAP is the
stored-field value, which drifts between generations (F1).

--------------------------------------------------------------------------------
F4. TASK 3 RESULT - FREQUENCY-DEPENDENT EXTINCTION (CHI PROVENANCE)
--------------------------------------------------------------------------------

[V] Both solvers build opacity by the same log-log cubic machinery
(MON_INT_FUNS_V2) but on different grids:

Formal ray path (fg_j_cmf_nonmono_v1.f):
  CHI_EXT/ETA_EXT on the extended grid (lines 3752-3808): interior
  CHI_EXT(I)=CHI(I-ND_ADD) - the SAME model opacity array the moment
  solver sees - with power-law edge extension; CHI_COEF cubic
  coefficients in log(CHI) vs log(R).
  Per ray point (~3890-3910, ~4211-4221): CHI_RAY=CHI_EXT(K) at EXT
  nodes else EXP(cubic in LOG(R_RAY/R_EXT)) - opacity interpolated
  to the RAY-POINT radius.
  Compressing augmentation (~4465-4478): CHI_RAY += GAM_MAG/DLN_C
  for cells passing KN_COUPLING_OK - the |GAM(mu)|/matched-DLN term.

Moment path (mom_j_cmf_nonmono_v1.f):
  CHI_NM (~1905-1935): same cubic-log interpolation of CHI_SM to
  model-depth radii R_NM via R_PNT_NM; else-branch copies
  CHI_SM(1:ND_SM) directly.
  CHI_HAT_NM = CHI_NM + GAM_MAG_NM/MAX(MOM_DLN_J,MOM_DLN_MIN)
  (~2273-2275) - |GAM_MOM|/matched-DLN.

Decomposition of the assembled-opacity difference at a cell:
  CHI_HAT_RAY - CHI_HAT_MOM =
      [CHI_raypoint - CHI_depth]            (grid/interpolation only)
    + [|GAM(mu)|_I/DLN_C - |GAM_MOM|/DLN_J] (the U11 coefficient term)

The first bracket is a pure discretisation difference on the same
underlying opacity - generically small where opacity is smooth, but
not guaranteed near opacity edges (which coincide with small DLN_N,
where the second bracket dominates anyway). The second bracket is
the measured U11 quantity: |GAMa|/DLN_C vs |GAMm|/DLN_J with the
F3 offset mismatch folded in.

NOT MEASURED: per-frequency CHI/ETA/ESEC are not in FGINEXT_GS and
were not dumped for this run (WRITE_VEC 'CHI Vec'/'CHI_HAT Vec' to
unit 48 is debug-gated). The dilution factor
  eta_CHI = CHI / (CHI + |GAM|/DLN)
that converts the coefficient discrepancy into a fractional
CHI_HAT error is therefore UNQUANTIFIED per frequency [H]. Bound
argument: where |GAM|/DLN >> CHI the U11 term dominates CHI_HAT and
the opacity provenance is irrelevant; where CHI >> |GAM|/DLN the
relative perturbation is suppressed by ~|GAM|/(DLN*CHI). Which
regime holds at the sigma~-3 band frequencies is the single most
important missing number for Task 6 - it requires either a debug
dump of CHI_NM/CHI_HAT_NM (unit-48 path) or the paired solve below.
Electron-scattering vs bound-free/line split: ESEC_NM is available
in the moment path (~1919) but likewise not archived; same caveat.

--------------------------------------------------------------------------------
F5. TASKS 4-5 - REFERENCE OPERATORS AND CONTROLLED COMPARISON
--------------------------------------------------------------------------------

Reference A (= production moment operator, now fully traced):
  coefficients per (depth I, freq k):
    f_I  = sign-floor(1+sigma_I*K_ON_J(I,k)),  |f|>=0.01
    GAMs = CONV_I*f_I                      (signed moment GAM)
    GAMm = |GAMs|
    DLN_J = dLOG_NU, or LOG(nu_k/nu_{k+1}), or matched span
    CHI_HAT = CHI_NM + GAMm/DLN_J
    Q_NM = (CONV_I/CHI_HAT/DLN_J)*f_I       (signed coupling)
    PSI  = DTAUONQ*|Q_NM|                    (implicit diagonal)
    donor: expanding -> PSIPREV*J_{k-1};
           compressing -> PSI*J_{k+N_CAP(I,k)}

Reference B (angularly partitioned, proposed - NOT implemented):
  retain the two ray populations separately:
    A_+(I,k) = sum_{mu:f>0} w I |f(mu)| / J,   donor_{k-1}
    A_-(I,k) = sum_{mu:f<0} w I |f(mu)| / J,   donor_{k+NMF(mu)}
  coupling = CONV*[ A_+*(J_{k-1}-J_k)/DLN_P + A_-*(J_{k+NMF}-J_k)/DLN_N ]
  i.e. the exact angular integral of the discrete ray operator
  before collapse to one coefficient. Requires per-ray donor
  intensities and per-ray matched offsets - data the moment solver
  does not currently carry (it would need NW x NCF extra storage or
  an on-the-fly ray loop inside the tridiagonal assembly).

Controlled comparison - three feasible levels, in increasing cost:

  L1 (coefficient/residual level, no new solver): apply A and B to
     the SAME stored field and compare the implied frequency-
     coupling source term cell-by-cell. Extends the 5.3E donor
     comparison with exact MATCH_FREQ offsets (F3) and per-ray
     donors. Answers "how different is the coupling operator on a
     given field" - NOT "how different is the fixed point".
     Feasible now; partially done in 5.3E/F3.

  L2 (reduced solve): freeze populations/opacity at one checkpoint,
     solve the moment system on a reduced depth/frequency subset
     (or single depth column) with operator A vs a Python Reference
     B assembly, identical eta/chi/J boundaries/tolerances. Needs
     CHI_NM/ETA_NM per frequency - obtainable from the unit-48
     WRITE_VEC path or a one-off diagnostic run. Answers whether
     the operator difference moves J at all in a linearised step.

  L3 (full paired run): the decisive experiment. Requires a
     Reference-B solver variant (new code, outside audit scope) or,
     cheaper, an ablation: run the EXISTING moment solver with the
     quench-band coefficient floor MOM_GAM_FAC_MIN raised/lowered and
     DO_NONMONO_FREQ_MATCH off/on to bound sensitivity - these are
     parameter changes, not solver changes, and would quantify how
     much the coupling coefficient matters dynamically without
     implementing B. Measure: convergence status/iteration count,
     final moment+formal residuals, JNU/HNU/K_ON_J field diffs
     (max/RMS/argmax), emergent flux and spectrum vs the SOB
     baseline, and energy conservation (luminosity integral).
     Identical tolerances and normalisation; only compare
     converged-vs-converged.

Physical consequence (Task 6): NOT determinable from coefficient
diagnostics alone. Arguments on record: the discrepancy is
  - localized: sigma~-3 band depths ~29-56, where |GAM_MOM| is
    floored small while ~half the angular weight retains O(1)
    coupling of both signs;
  - bounded in the bulk: median Rabs ~2, p90 ~24 across the band,
    ~1 outside it;
  - amplified by 1/DLN: |GAM|/DLN can exceed CHI by orders of
    magnitude at bound-free-edge frequencies (DLN_N min 1.67e-5) -
    but the SAME amplification applies to CHI_HAT on both sides,
    so the fractional CHI_HAT error is what matters (F4, unknown);
  - accompanied by a donor-location mismatch (F3) that shifts the
    compressing read by up to hundreds of bins in the deep interior.
None of these establish a spectral or population consequence.
Classification (Task 7): U11 remains [C]. The matched-offset
spread is recorded as U13 [C] - same structural family (single
coefficient vs angular distribution), quantified, not a defect.

--------------------------------------------------------------------------------
F6. WHAT WOULD BE NEEDED BEFORE GENERALISING
--------------------------------------------------------------------------------

1. A converged (or near-converged, documented tolerance) KNOP
   coupled solution to re-measure R_abs - the current evidence is
   pre-convergence only.
2. The CHI-dilution number eta_CHI per (depth,freq) - unit-48
   diagnostics or a diagnostic dump run.
3. The L2/L3 controlled comparison - at minimum the parameter
   ablation (MOM_GAM_FAC_MIN, FREQ_MATCH) on a converged run,
   ideally a Reference-B operator.
4. Per-frequency donor-intensity weighting in the compressing band
   (the 5.3E donor comparison used adjacent-bin stencils; F3 gives
   the offsets but the full matched-donor intensity comparison is
   still pending).
5. Whether the deep-interior offset mismatch (U13) matters depends
   on whether those depths/frequencies carry significant flux -
   an emergent-flux weighting of the NMF discrepancy is needed.

STATUS SUMMARY
  U11: C - structural mismatch, production-verified, measured on
       stored fields, STABLE across 11 diverse checkpoints
       (Rabs_med 3.00, p90 23.6 in the sigma~-3 band). Not D:
       no converged incorrect answer; consequence unproven pending
       the controlled comparison.
  U13 (new): C - matched-frequency offset is per-ray in the formal
       solver vs single-valued per-depth in the moment solver;
       ray offsets exceed the moment's by up to ~4x (p90) in the
       deep interior; ~1 everywhere in the quench band.
  Checkpoint caveat: FGINEXT_GS is a live, evolving store; all
  snapshot analyses carry generation drift (~2-3x in J between
  reads ~45 min apart); ratio diagnostics are robust, absolute J
  comparisons across generations are not.

No solver change, parameter change, or commit is authorised or
implied by Phase 5.3F.


================================================================================
PHASE 5.3G - OPACITY PROVENANCE AND CONTROLLED U11/U13 REFERENCE EXPERIMENT
================================================================================
Mode: diagnostic / reference implementation only.  No production Fortran,
parameters, OMEGA_RAY, or the KNOP moment solver were modified.  Nothing
committed.  All artifacts live in /tmp (knop53g_*).

--------------------------------------------------------------------------------
G1. IMMUTABLE DATA PROVENANCE [V]
--------------------------------------------------------------------------------
The 14 FGINEXT_GS files from Phase 5.3F are separate runs, not
sequential checkpoints.  The p8_base file is rewritten by the live run,
so an immutable copy was taken before any 5.3G calculation:

  snapshot:   /tmp/knop53g_snapshot/FGINEXT_GS_snap.bin
  size:       5,367,543,560 bytes
  md5:        1c5be38bc50f9b41e124c80fb2420a82
  VER=3, tag(NITSF)=275, dims N_COMP=70 NP=85 NCF=56377 ND=70
  metrics:    MAXCHG=5.5383e-03, SWMAX=2.6686e-01, NCONV=0, NDIV=0
  end-magic verified; record walk clean through rec8; I_P data at
  byte 453120, I_M at 453120+NEL*8+16 (NEL=N_COMP*NP*NCF).
  source commit: c682aeb6ae50fb102aefdd31897cd5a00ef6d7ef
  exe: exe_dbg/cmfgen_dev.exe md5 36493c21c471a0a4f08b3852de82e610
       ("Phase-7 exe 0153ec58" per config_resolved.yaml)
  inputs: SN_HYDRO_DATA md5 ea6c5d5bd2d1c0ae785d70d4c9a56053
          MODEL_SPEC    md5 b8b94bcc1d3b3434508950ffcce3fd8e
  model time 2.82 days; extraction 2026-11-14 ~16:05 (+0100).
  Full manifest: /tmp/knop53g_provenance.txt.

Angular geometry/weights: /tmp/angqw.npz regenerated from production
objects (SIGMA_TRUE_CORE, IMPAR_WITH_PNT_SRCE, NORDANGQW_V2,
J/K/H/NWEIGHT_V2) - geometry is static hydro, generation independent.
Field coefficients: computed from the SNAPSHOT (tag 275), not the older
ifield bins, so coefficients and provenance share one generation.

--------------------------------------------------------------------------------
G2. OPACITY PROVENANCE [V]
--------------------------------------------------------------------------------
Both solvers interpolate the SAME model opacity with the SAME routine:

  subs/mon_int_funs_v2.f  MON_INT_FUNS_V2(COEF,CHI,R,ND)
    monotone log-log cubic; called on LOG_CHI_EXT/LOG_ETA_EXT
    (fg_j_cmf_v10.f:965-967 pattern) producing CHI_COEF/ETA_COEF(ND,4).
  plane/fg_j_cmf_nonmono_v1.f
    CHI_EXT built on an extended grid (ND_ADD power-law/ESEC-anchored
    extrapolation, :3760-3800); at ray points equal to model radii uses
    CHI_EXT(K), between points evaluates the cubic in
    T1=LOG(R_RAY/R_EXT(K)) -> CHI_RAY=EXP(T2); then adds
    |GAM|/DLN_C for compressing cells.
  plane/mom_j_cmf_nonmono_v1.f
    CHI_NM interpolated at model-depth radii; CHI_HAT_NM=
    CHI_NM+GAM_MAG_NM/MAX(MOM_DLN_J,MOM_DLN_MIN).
  new_main/mod_subs/comp_j_blank.f:1107
    moment solver receives CHI_CLUMP, CHI_SCAT_CLUMP (=CHI*CLUMP_FAC,
    CHI_SCAT*CLUMP_FAC) - i.e. clumping-corrected total opacity.
  new_main/mod_subs/comp_opac.f
    CHI assembled as: GENOPAETA_V10 (free-free + bound-free) +
    TWO_PHOT_OPAC_V3 (optional) + ESOPAC (electron scattering -> ESEC)
    + Rayleigh (CHI_RAY) ; CHI_NOSCAT keeps the pre-scattering part;
    CHI_CONT = total CHI.  Spectral lines are NOT in CHI - they are
    handled by the line/Sobolev machinery through source terms, the
    standard CMFGEN continuum/line split.

Consequence [V]: the U11/U13 terms are ADDITIVE corrections to a shared
base opacity; the only asymmetry is WHERE the opacity is sampled
(ray-point radius vs model-depth radius) and WHICH |GAM| coefficient
multiplies 1/DLN.  Per-frequency raw CHI is not archived (unit-48
WRITE_VEC is debug-gated), so the CHI-dilution factor F could not be
measured from the checkpoint; it is bracketed by manufactured CHI in
the reduced experiment below.

--------------------------------------------------------------------------------
G3. FOUR REFERENCE ARMS - REDUCED REFERENCE PROBLEM [M]
--------------------------------------------------------------------------------
Driver: /tmp/knop53g_arms.py (Python/scipy, snapshot field, production
angular coefficients).  Solves the fixed-point linear system for J(I,k)
directly (donors implicit) - identical fixed point to the lagged GS
iteration, without iteration-path noise.

Grid: all ND=70 depths; freq window W=1600 bins at k0=17994
(DLN 1.67e-5..1.67e-4, the small-spacing region where offsets are
largest).  Angular coefficients per (I,k) from the snapshot field:
K_ON_J, GA=<|f_mu|>_I, GS=<f_mu>_I with production JQW/KQW weights.

Equation per (I,k) (documented reduction of the production G_ONLY
coupling sector; radial skeleton shared identically by all arms):
  TA*J(I-1)+TB*J(I)+TC*J(I+1) - sum_g w_g*PSI*J(I,k+off_g) = XM0
  TB = DTAUONQ + PSI + HL + HU ;  PSI = DTAUONQ*|g_eff|/(CHI_HAT*DLN)
  CHI_HAT = CHI_manufact + |g_eff|/DLN ;  QF=R^2 approx;  COH=0;
  H/EPS lagged terms dropped; zero-gradient radial BCs; PSIPREV=PSI
  (MOM_GAM_PIN assumption).  These reductions are IDENTICAL across arms
  so they cancel in every DX - they limit absolute fidelity to
  production, not the inter-arm comparison.

Arms (exactly the task definitions):
  A: g_eff = g_mom = CONV*f_mom, f_mom=1+SGC*K_ON_J (|f|>=0.01 floor);
     donor: sign>=0 -> J(k-1); else J(k+N_CAP), N_CAP=ceil(|g_mom|*dR/DLN)
  B: |g_eff| -> g_abs = CONV*GA (intensity-weighted |f_mu|); sign and
     donor construction unchanged from A (U11 only).
  C: coefficients = A; donor ray-partitioned: each ray picks its own
     donor by sign(f_mu): f_mu>0 -> k-1, f_mu<0 -> k+ceil(|GAM_ray|*dZ_ray/
     DLN_k); group weights = JQW*I*|f_mu| fractions (U13 only).
  D: B coefficients + C donor.

CHI manufactured at three dilution scales (3.0, 0.3, 0.03) bracketing
chi >>, ~, << |GAM|/DLN since real per-frequency CHI is unavailable.
Dark cells (J < 1e-10*max per depth) get neutral closure - radiation
is zero there so arm differences vanish anyway.  Donor offsets outside
the window are dropped and counted (139 for A/B, 773 for C/D).

Validation performed: all 12 systems solve to machine precision
(spsolve, direct factorisation - no iteration tolerance ambiguity);
J positive in the interior; ARM A coefficient formulas verified
term-by-term against mom_j_cmf_nonmono_v1.f lines traced in 5.3E/F;
single-sign limit checked analytically (C->A when all f_mu>0).
NOT done: bit-level reproduction of the production radial operator
(simplified skeleton above) - ARM A reproduces the production
COUPLING-SECTOR coefficients, not the full matrix.

--------------------------------------------------------------------------------
G4. CONTROLLED REDUCED SOLVES - RESULTS [M]
--------------------------------------------------------------------------------
DX_XA = X_X - X_A, normalised by RMS(J_A) over interior freqs
(margin 560 bins each side; max donor offset ~700 so kl>900 in the
W=1600 grid is edge-affected - see caveat).

scale=3.0 (chi dominates):
  DX_BA  rms=3.57e-04  max=1.57e-02
  DX_CA  rms=1.14e-03  max=2.33e-02
  DX_DA  rms=4.05e-04  max=1.61e-02
  DX_int rms=1.02e-03
scale=0.3:
  DX_BA  rms=5.39e-04  max=1.92e-02
  DX_CA  rms=4.24e-03  max=5.11e-02
  DX_DA  rms=1.24e-03  max=2.31e-02
  DX_int rms=3.38e-03
scale=0.03 (GAM dominates):
  DX_BA  rms=4.91e-04  max=1.06e-02
  DX_CA  rms=9.75e-03  max=6.04e-02
  DX_DA  rms=4.52e-03  max=4.88e-02
  DX_int rms=6.32e-03

DEPTH PROFILES (rms DX per depth, scale=0.3 / 0.03):
  DX_BA peaks ONLY at I=64 (sigma=-156.9, deepest compressing point):
    4.5e-3 / 4.0e-3; all other depths <= 4e-4.
  DX_CA peaks at I=66-69 (sigma=-50.9,-39.9,-31.8,-31.0) plus I=64:
    up to 1.9e-2 / 4.1e-2 per-depth rms - the deep interior.
  The sigma~-3 band (I~44-52) shows NO local maximum in any arm.

INTERPRETATION:
  A. U11 has a measurable but SMALL independent effect (rms ~4e-4,
     ~0.1-0.5% typical, max ~1-2%), concentrated at the single deepest
     compressing depth.  It does not scale with the chi dilution
     bracket the way U13 does.
  B. U13 has a LARGER independent effect (rms 1e-3..1e-2, max ~2-6%),
     growing as |GAM|/DLN dominates over chi - exactly the regime where
     the frequency coupling matters.  Concentrated at the inner
     compressing depths, NOT the sigma~-3 quench band.
  C. Interaction is non-negligible (same order as DX_CA at low chi) -
     the two effects partially reinforce rather than cancel.
  D. Effects are NOT confined to sigma~-3; they extend into the deep
     compressing region (sigma ~ -30..-157).
  E. Inner-boundary caveat [H]: I=66-69 sit adjacent to the inner
     zero-gradient BC; part of the concentration may be a boundary
     artifact of the reduced skeleton.  A production-BC variant is
     needed before quantifying the physical magnitude.

--------------------------------------------------------------------------------
G5. CLASSIFICATION REVISION
--------------------------------------------------------------------------------
  U11: remains [C] (structural mismatch).  The reduced solve shows a
       real but small solution-level effect (~4e-4 rms), localized to
       the deepest compressing depth.  Still not [D]: the reduced
       problem is not the production operator and no conserved-
       quantity violation was demonstrated.
  U13: remains [C], strengthened to [C+] evidence-wise: it is the
       DOMINANT of the two approximations in the reduced experiment
       (up to ~30x DX_BA in rms at chi_scale=0.03) and grows in the
       GAM-dominated regime.  Its concentration near the inner
       boundary needs the production-BC check before physical claims.

--------------------------------------------------------------------------------
G6. FULL-MODEL EXPERIMENT - DESIGN ONLY (not run)
--------------------------------------------------------------------------------
Justified: YES as a next step, because the reduced problem shows both
approximations produce measurable (percent-level local, sub-percent
global) solution differences whose sign/magnitude depends on the
unmeasured CHI dilution - exactly the quantity the full model has.
Design: identical initial model/inputs/BCs/freq grid/convergence;
ARM A = production solution (already exists as reference);
B/C/D = same run with the two coefficient constructions swapped in a
DIAGNOSTIC copy of mom_j_cmf_nonmono_v1 (not production).  Compare
JNU, HNU, KONJ, formal+moment residuals, emergent flux, populations,
line profiles.  Hard gate: every arm must reach the same convergence
criterion before comparison (the KNOP non-convergence problem of
5.3A-F applies to all arms equally and must be solved first - the
full-model experiment is BLOCKED on KNOP convergence, not on U11/U13).

LIMITATIONS: reduced problem uses manufactured CHI (real per-frequency
CHI unarchived), simplified radial skeleton, zero-gradient BCs, frozen
tag-275 angular coefficients, single freq window.  The W=2400-bin
confirmation run (edge-clean margin) was in progress at write time;
W=1600 numbers above are edge-checked at the reported max-diff
locations but rms includes ~140 edge-affected freqs - expect rms to
drop ~10-30% in the clean rerun, not to change sign or ranking.

No production correction implemented.  No commit.


================================================================================
PHASE 5.3H - KNOP P-DIRECTION CONVERGENCE FAILURE + U13 BOUNDARY ROBUSTNESS
================================================================================
Mode: audit/investigation only.  No production changes, no commit.

--------------------------------------------------------------------------------
H1. P-DIRECTION MAP AND RESIDUAL DEFINITION [V]
--------------------------------------------------------------------------------
plane/fg_j_cmf_nonmono_v1.f, SWEEP_FLD_ORD P-pass (ASC=.TRUE., L=1..NCF):
  per ray point J:  KN_ORD_SPW_FILL(FLD,J+1) builds a PREVIEW of row j+1
    from its pre-march state (the lagged C-corner); KN_ORD_SPV gives the
    source SPJ; KN_ORD_DS(j-1,L),KN_ORD_DS(j,L) give the two short-
    characteristic coefficients; VNEW = FLD(J-1,L)*A0 + SPW*A1 + SPJ*A2
    + DSD*A3 + DSU*A4  (fresh radial upwind FLD(J-1,L); same-row donor
    handled inside SPV/DS via KN_FREQ ordering; SPW = lagged preview).
  Then VNEW is floored (1e-30) and capped (NONMONO_INEXT_CAP) BEFORE the
  residual is taken, and the commit is FLD = VOLD + OMEGA_RAY*(VNEW-VOLD)
  applied INSIDE the recurrence (per-cell, compounding along the ray).

RESIDUAL [V]:  RES = max over cells of |VNEW - VOLD| /
  MAX(VOLD, 1e-4*VMAXF(L)), evaluated on the UNRELAXED increment
  (VNEW vs pre-update VOLD) BEFORE OMEGA_RAY is applied.  It is a
  per-cell MAX norm (not RMS), gated to VOLD>0 cells, per ray then
  maxed into RES_P; boundary row J=1 is zeroed and also residual-scored.
  Positivity floors/caps DO enter the norm (they shape VNEW).

CONSEQUENCE [V]: the reported plateau (~0.6-0.7) measures the UNDAMPED
map increment ||G(x)-x||, not the damped step ||H(x)-x|| =
OMEGA_RAY*||G(x)-x||.  OMEGA_RAY cannot reduce the reported residual
directly; a plateau therefore means the undamped map genuinely fails
to contract at the worst cell - oscillation, slow contraction, or a
near-unit eigenvalue - not a damping-reporting artifact.

--------------------------------------------------------------------------------
H2. OMEGA_RAY CONTROLLER [V]
--------------------------------------------------------------------------------
Adaptation (fg_j_cmf_nonmono_v1.f:10416-10475), per direction IDIR,
same-direction history RES_PREV(IDIR) (pass N vs N-2 because M/P
alternate):
  RATIO_OM = RES/RES_PREV(IDIR)
  RATIO>1.05 and NPO>=5  -> OMEGA_RAY *= 0.7   (floor 0.45)
  RATIO<0.5              -> OMEGA_RAY *= 1.15  (cap 1.0)
  sign-flip detector: two consecutive flips of the signed update sum
    SGN_ACC*SGN_PRV<0 -> OMEGA_RAY *= 0.7 (period-2 damping, W6)
  stall detector: NSTALL consecutive same-direction ratios >0.98 ->
    early exit (default <2%/pass improvement = stagnation).
OMEGA_RAY is never >1 (comment: it is applied inside the
short-characteristics recurrence where amplification would compound
per cell).  The controller CAN produce artificial plateaus: repeated
damping to the 0.45 floor slows the committed field without changing
G, so ||G(x)-x|| can stay large while the iterated field crawls -
a damped-plateau is consistent with an undamped near-unit or
negative eigenvalue.  Distinguishing requires the field-level
iteration distances (Task 4), which are NOT stored: FGINEXT_GS is
written once per model iteration, not per ray-pass.  Per-pass d1..d4
requires a diagnostic driver dumping FLD between passes - designed,
not run this phase.

--------------------------------------------------------------------------------
H3. C-CORNER (LAGGED PREVIEW) [V]
--------------------------------------------------------------------------------
The P-sweep reads SPW(L,MOD(J-1,4)+1): a preview of row j+1 built by
KN_ORD_SPW_FILL from the PRE-MARCH state of that row.  This is a
block Gauss-Seidel splitting: the downwind neighbour's coefficients
are evaluated one step ahead but its field values are lagged.  It is
the standard short-characteristics ordering (the j+1 row cannot be
fresh because the sweep has not reached it); the lag is O(one cell)
and direction-consistent.  It is NOT a simultaneous-solve defect per
se, but it IS a delayed-feedback term that can sustain period-2
cycles together with the frequency donor lag.  Reference test
(lagged vs consistent C-corner) designed; not run - needs a driver.

--------------------------------------------------------------------------------
H4. U13 BOUNDARY ROBUSTNESS - RESULT [M]
--------------------------------------------------------------------------------
Repeated the 5.3G ARM A vs ARM C comparison (scale=0.03, W=2400) with
a 15-depth ghost extension replicating the innermost coefficients so
the inner boundary moved 15 cells deeper (/tmp/knop53g_arms53h.py,
N_GHOST=15; coefficients con_vr-rescaled in ghosts).

Per-depth rms DX_CA (real depths, interior freqs):
  I<=59 : 8.5e-4   (bulk, small)
  I=62-65: up to 2.1e-2
  I=66-69: ~4.1-4.3e-2  (formerly adjacent to inner BC)
  I=70-84 (ghosts, 15+ cells from ANY boundary): 4.1e-2 -> 4.5e-2,
    a flat plateau, NOT decaying with distance from the boundary.

CONCLUSION [M]: the U13 solution difference is a BULK property of the
deep compressing region (sigma ~ -30..-157), not an artificial
inner-boundary artifact.  The earlier concentration at I=66-69 was the
innermost real depths where |GAM|*dZ offsets are largest; extending
the domain shows the effect persists uniformly through the ghost
region.  U13 confidence upgraded: the reduced-solution effect is
physically located, not numerically boundary-generated.
(Max compressing donor offset actually used in this window: 168 bins,
not the 700 cap - earlier edge concerns were conservative.)

--------------------------------------------------------------------------------
H5. SPECTRAL ANALYSIS - DESIGNED, NOT RUN
--------------------------------------------------------------------------------
For the reduced system the lagged-donor GS map is E = T^{-1} D where
T = block-tridiagonal radial+diagonal part (factorisable once via
splu) and D = donor off-diagonal coupling.  Planned: splu(T) then
scipy eigs on the LinearOperator matvec = T^{-1}(D x) for the
dominant eigenvalues (magnitude + phase -> distinguishes
contraction vs period-2 vs slow mode), plus OMEGA_RAY sweep of the
damped map (1-w)I + w*T^{-1}D.  Not executed this phase (context).

--------------------------------------------------------------------------------
H6. CLASSIFICATIONS
--------------------------------------------------------------------------------
  P-direction plateau: [B/C] - the residual measures the undamped
    map increment; the plateau is consistent with a genuine
    non-contractive or near-unit mode of the discrete map (mixed GS
    ordering + lagged C-corner + frequency-donor lag), NOT shown to
    be a damping artifact.  Whether it is a true fixed-point defect
    [D] remains unproven: needs per-pass field distances and the
    spectral estimate.
  OMEGA_RAY controller: [C] convergence-control weakness candidate -
    it can mask an undamped plateau by slowing commits; cannot be
    exonerated or blamed without d1..d4 dynamics.  Do NOT tune it
    until the undamped map's spectrum is measured.
  C-corner lag: [B] discretisation-consistent block-GS splitting;
    possible period-2 contributor, unquantified.
  U13: upgraded confidence - bulk deep-interior effect, boundary-
    robust [M].  Still [C] for production (reduced model only).
  U11: unchanged [C], small effect.

NEXT CONTROLLED EXPERIMENT (before any solver change):
  1. Diagnostic driver dumping FLD per P-pass for one representative
     ray/freq window -> compute d1..d4, classify oscillation vs stall.
  2. splu-factorised spectral estimate of T^{-1}D on the reduced
     operator (ARM A) -> dominant eigenvalue magnitude/phase.
  3. C-corner consistency variant in the same driver.
  4. Only then decide whether the plateau is a discretisation
     property [B] needing a formulation change or a control weakness
     [C] needing a controller change.

No production solver, parameter, OMEGA_RAY, MATCH_FREQ, or BC change
made or proposed for implementation.


================================================================================
PHASE 5.3I - P-SWEEP MAP, FIELD-LEVEL CONVERGENCE, SPECTRAL ANALYSIS
================================================================================
Mode: diagnostic/reference only.  No production changes, no commit.

--------------------------------------------------------------------------------
I1. WHAT WAS MEASURED [M]
--------------------------------------------------------------------------------
Per-pass FLD states are NOT archived (FGINEXT_GS is per model
iteration), so Tasks 1-3 (complete P-pass field distances, independent
undamped-map evaluation) require a diagnostic build - designed, not
run.  What WAS executed: the exact matrix splitting and spectral
analysis on the Phase 5.3G reduced ARM-A operator, which isolates the
frequency-coupling sector of the moment J-equation.

Splitting derived from the assembled stencil (not assumed):
  A = T + E,  T = block-diagonal radial+diagonal part (per-frequency
  70x70 blocks: TB diag, -HL/-HU off-diag, BC rows), E = donor
  off-diagonals (-w*PSI at (I,k)->(I,k+off)).  Lagged-donor iteration
  T x_{n+1} = -E x_n + b gives G = I - T^{-1} A.  T is block-diagonal
  -> splu exact in 0.14 s (N=168000, W=2400, ND=70).

Spectral/transient results (ARM A, chi_scale=0.3, undamped G):
  - scipy eigs/ARPACK: NO convergence in 301 iterations - the operator
    is strongly non-normal; Krylov eigen estimation fails outright.
  - Power iteration ||G^k x||, 3 random starts (identical):
      k=1: 0.070   k=2: 0.020   k=3: 0.007   k>=4: <0.002
    -> the undamped map contracts ~3x per pass; NO transient growth,
    NO period-2 signature in the reduced operator.
  - Fixed-omega spectral radii (60 power iterations each):
      omega=1.00 -> ~0.33 ; 0.80 -> ~0.48 ; 0.60 -> ~0.62 ; 0.45 -> ~0.72
    Under-relaxation SLOWS convergence here: H=(1-w)I+wG inherits a
    non-contracting identity component (eigenvalue 1-w at lambda=0),
    so damping only helps against rho(G)>1 or oscillatory modes -
    neither is present in this sector.

--------------------------------------------------------------------------------
I2. INTERPRETATION [M]/[H]
--------------------------------------------------------------------------------
The reduced J-sector (frequency-coupling + radial skeleton) is
strongly contractive with rho ~ 0.3 and no non-normal transient.
THEREFORE the production P-direction plateau (~0.6-0.7) is NOT
explained by the moment J frequency-coupling discretisation itself.
The residual must originate in physics the reduction excluded:
  - the H/EPS (second-moment) coupling dropped from the skeleton;
  - production boundary conditions (diffusion inner, thin/thick outer)
    vs zero-gradient;
  - the formal-solver ray structure (short characteristics along rays,
    per-ray MATCH_FREQ offsets, C-corner) rather than the moment
    tridiagonal;
  - the OUTER model iteration - opacity/population lag means the J-map
    is re-linearised every pass; the plateau may be a coupled
    fixed-point of the whole model state, not the J-sector alone.

This is a valuable negative localisation [M]: it rules out the J-moment
frequency-coupling discretisation as the plateau source and directs
the search to the formal-solver sweep itself or the outer iteration.

C-corner experiment (Task 7): N/A in the reduced moment operator - the
C-corner/preview exists only in the formal ray sweep (SPW_FILL), not
in the moment tridiagonal.  Its test requires the formal-solver
driver (Task 1-3 instrumentation), not yet built.

--------------------------------------------------------------------------------
I3. CLASSIFICATIONS
--------------------------------------------------------------------------------
  J-moment frequency-coupling sector: [V-negative] strongly contractive
    reduced map (rho~0.3, no transient, no period-2) - exonerated as
    the plateau source within the reduced model's scope.
  P-direction plateau: [E] unresolved - narrowed to formal-solver
    sweep structure, boundary conditions, H/EPS coupling, or the
    outer opacity/population iteration.  NOT yet [B] or [D].
  OMEGA_RAY damping: [C-candidate] on contractive sectors damping only
    adds slow identity modes - but production plateau may sit in a
    genuinely non-contractive sector where damping IS appropriate;
    cannot conclude until the true undamped map is measured.
  Non-normality: [M] ARPACK non-convergence is direct evidence the
    reduced operator is far from normal - spectral radius alone will
    not describe production dynamics; pseudospectral analysis needed.

NEXT (blocking, in order):
  1. Diagnostic build dumping FLD per P-pass (Tasks 1-3) - the only
     way to get true complete-map residuals and d1..d4 on production.
  2. Extend the spectral splitting to the FORMAL ray sweep structure
     (per-ray short characteristics), not just the moment tridiagonal.
  3. Outer-iteration coupling test: freeze vs live opacities/populations
     to isolate whether the plateau is J-sector or model-state.

No production solver, OMEGA_RAY, coefficient, BC, or parameter change
made or proposed.


================================================================================
PHASE 5.3J - PRODUCTION P-SWEEP DIAGNOSTICS + PART-A CLAIM VALIDATION
================================================================================
Mode: diagnostic/reference only.  No production changes, no commit.

--------------------------------------------------------------------------------
J0. PART A - CORRECTIONS TO PHASE 5.3I CLAIMS [V]
--------------------------------------------------------------------------------
The 5.3I spectral measurements are RELABELLED as follows:

A1. The sequence 0.070, 0.020, 0.007, <0.002 was a SINGLE-VECTOR power
    measurement: Gmv applied repeatedly to one random unit vector,
    recording ||G^k x||/||x||-style norms (script /tmp/knop53i_spec.py,
    trials 0-2, rng seed 1).  It is NOT an induced operator norm and
    NOT a power-method estimate of sigma_max(G^k).  Correct label:
    "single-vector transient amplification along random directions".
    All three random starts agreed to ~3 digits (0.070/0.020/0.007),
    consistent with a dominant-contracting direction being generic.

A2. ARPACK eigs (k=8, which='LM', tol=1e-4, maxiter=300) returned
    "No convergence (301 iterations, 0/8 eigenvectors)".  This is NOT
    evidence of non-normality - it only shows the Krylov projection
    did not isolate eigenvalues, expected for a non-symmetric
    168k-unknown operator with clustered spectrum.  A follow-up
    sigma_max/ rho diagnostic (/tmp/knop53j_nonnorm.py,
    /tmp/knop53j_trans.py) gave INCONCLUSIVE results: power iteration
    on G^T G (80 iters) estimated sigma_max~1.67 but the returned
    vector only achieves ||Gx||=0.66 - the iteration did not converge
    to the true top singular direction, so the 1.67 figure is
    unreliable and is WITHDRAWN.  Along the approximate top singular
    vector ||G^k|| decays 0.66->0.19->0.05->... - no amplification.
    Honest status: non-normality of the reduced operator is
    UNRESOLVED; all measured transient sequences decay.

A3. The omega spectral radii (0.45->0.72, 0.6->0.62, 0.8->0.48,
    1.0->0.33) were 60-iteration power estimates on random vectors -
    estimates of the dominant-eigenvalue magnitude along generic
    directions, NOT certified spectral radii.  They correctly show
    the ordering (damping slows a contractive map) but the absolute
    values carry ~10-20% uncertainty.  Same G (I-T^{-1}A) and state
    ordering used for all omega - that part verified.

A4. Splitting validated: A = T + E with T block-diagonal per
    frequency (2400 blocks of 70x70), E = donor off-diagonals.
    Fixed point of x_{n+1}=T^{-1}(-E x_n+b) satisfies (T+E)x=b=Ax.
    Boundary rows (zero-gradient, I=0 and I=ND-1) are inside T and
    carry no donors - consistent with production interior-only donor
    logic.  Frequency-endpoint rows: donors outside the window are
    dropped (139-773 cells) - a window-edge truncation, documented.
    Matrix-free Gmv = x - LU.solve(A@x) uses the SAME assembled A -
    self-consistent but NOT an independent matrix-free check (the
    production Fortran sweep was not replicated; that is Part C,
    not done).

--------------------------------------------------------------------------------
J1. PARTS B-G - STATUS: DESIGNED, NOT EXECUTED
--------------------------------------------------------------------------------
The diagnostic build (per-pass FLD capture, R_MAP full-pass defect,
OMEGA_RAY/pass history, MATCH_FREQ donor census, C-corner variant,
freeze-vs-live opacity arms) requires Fortran instrumentation of
plane/fg_j_cmf_nonmono_v1.f behind a disabled-by-default flag plus a
production run on the interactive partition.  This is a multi-hour
implementation+run task that was not completed in this session.
Design is fully specified in the task text; the blocking item is a
dedicated build/run, not analysis.

--------------------------------------------------------------------------------
J2. WHAT IS ESTABLISHED vs UNRESOLVED
--------------------------------------------------------------------------------
[V] Reduced ARM-A J-sector: undamped GS map strongly contractive on
    every measured direction (single-vector transients 0.07->0.02->
    0.007->...; all trials identical); fixed-omega ordering correct.
[M] sigma_max/rho non-normality diagnostic: inconclusive (non-
    converged power iteration); withdrawn 1.67 figure.
[C] The plateau therefore remains unexplained by the reduced
    J-frequency-coupling operator; candidate sources unchanged:
    formal ray sweep (short characteristics, MATCH_FREQ, C-corner),
    H/EPS coupling, boundary conditions, or outer model-state
    feedback (opacity/population lag).
[E] Whether OMEGA_RAY helps or hurts the PRODUCTION map: unknown -
    the reduced operator is contractive so damping only adds slow
    modes there, but the production plateau may sit in a different,
    genuinely non-contractive sector.

BLOCKER for Part F answers 1-8: the per-pass diagnostic build +
one production run.  Until then no plateau explanation is claimed.

No production solver, OMEGA_RAY, coefficient, BC, parameter, or
commit-policy change made.


================================================================================
PHASE 5.3K - P-SWEEP DIAGNOSTIC INSTRUMENTATION (implementation in progress)
================================================================================
Mode: diagnostic-only instrumentation, disabled by default.  No
production numerical changes, no commit.

--------------------------------------------------------------------------------
K1. PART A - PASS-LOOP STRUCTURE (source-traced) [V]
--------------------------------------------------------------------------------
SWEEP_FLD_ORD (contained, plane/fg_j_cmf_nonmono_v1.f ~line 10204):
  called per ray per direction from KNOP_RAY_SOLVE_V1:
    :9628  SWEEP_FLD_ORD(FLD,  .TRUE.,  RES_M, ...)   M channel
    :9697  SWEEP_FLD_ORD(FLD_P,.FALSE., RES_P, ...)   P channel
  Two channels hold SEPARATE field arrays (FLD vs FLD_P) - no shared
  mutable field state between directions.

Pass loop: `DO ; NPO=NPO+1` (~10274) ... `END DO` (~10607).
  ASC = MOD(NPO,2)==1 : odd passes ascend L (LB=1,LE=NCF,LDI=1);
  even passes descend (LB=NCF,LE=1,LDI=-1) - the alternating-direction
  Gauss-Seidel sweep (Knop Paper II ordering).
  One complete pass = one (J,L) sweep: `DO L=LB,LE,LDI` outer over
  frequencies; inner `DO J=2,NI` depths (~11047).  FLD(1,L) zeroed per
  pass before the depth march (thin outer BC write).

Iterated field = FLD(NI,NCF), in-place per cell.
  Fresh: FLD(J-1,L) radial upwind (just relaxed), SPJ source.
  Lagged: SPW(L,MOD(J-1,4)+1) preview of row j+1 built by
    KN_ORD_SPW_FILL from pre-march state (the C-corner).
  Retained from previous pass: all FLD cells not yet visited this
  pass (the sweep is in-place).
  VMAXF(L) = MAX(MAXVAL(FLD(1:NI,L)),1.0e-25) computed once before the
  pass loop (~11017); the production residual normalisation is
  MAX(VOLD, 1.0e-4*VMAXF(L)).

Exit classes inside the pass loop:
  IXIT=0 HAS_CP false (exact single pass, RES forced 0)
  IXIT=1 RES < KN_RAY_TOL (clean)
  IXIT=2 bulk bad-fraction gate (NPO>=6 + contraction/RES floor opts)
  IXIT=3 stall (NSTALL>=2, no hot window)
  IXIT=4 pass cap KN_RAY_MAXPASS / hot-pass budget exhausted
  IXIT=5 diverged (inside sweep, ~10337)
  IXIT=6 plateau (L_PLAT, through the stall branch)

CAPTURE POINTS (for the diagnostic):
  X at pass head: after the ASC/RES=0 initialisation (~10302), before
    `DO L` - snapshot FLD(NI,NCF).
  Y at pass end: immediately before `IF(L_NOPROBE)CALL NOP_PASS`
    (~10535) - reached only when the full L-sweep completed.
  Incomplete pass: any EXIT inside the L-sweep (IXIT=5 path) or an
    exit fired before Y capture -> flag partial, never mislabel.
  Storage: LOCAL allocatable buffer inside SWEEP_FLD_ORD (thread-
  private by construction; the ray loop is the OMP parallel region).
  Output: buffered per-ray records, flushed once under
  !$OMP CRITICAL to a diagnostic unit; snapshots optional bounded.

STATUS: Part A tracing complete and documented; Parts B-G
(implementation, tests, controlled run) pending - see K2+.

======================================================================
PHASE 5.3K - PRODUCTION P-SWEEP DIAGNOSTIC INSTRUMENTATION
======================================================================

K3. Implementation summary
--------------------------

Source changes (all diagnostic-only; no production numerical
behaviour touched):

- plane/fg_j_cmf_nonmono_v1.f (+~210 lines):
  * Module flags (default F): KN_PSWEEP_DIAG, KN_PSWEEP_SNAP;
    ints KN_PSWEEP_LS (=0 all rays), KN_PSWEEP_PLO (=1),
    KN_PSWEEP_PHI (=9999), KN_PSWEEP_SNHI (=0). Setter
    SET_FG_PSWEEP_DIAG mirrors SET_FG_NOCP_PROBE; PHI is clamped
    >= PLO so a reversed window fails closed.
  * SWEEP_FLD_ORD local state: PSD_X (immutable pass input),
    PSD_REC(NMAX,18) record buffer, all ALLOCATABLE local - no
    SAVE, no shared state between rays/iterations.
  * Arm: PSD_ARM = KN_PSWEEP_DIAG .AND. .NOT. L_QP .AND.
    (LS<=0 .OR. LS==KN_PSWEEP_LS). The .NOT. L_QP guard excludes
    the KN_Q_RESID verify passes, which run the same routine on a
    scratch field and are not production iterations.
  * X capture: after the per-pass VMAXF rebuild, before the
    IF(INWARD) branch - the last point where FLD is unmodified
    this pass. VMAXF itself is the production column-peak scale
    (MAX(MAXVAL(FLD(1:NI,L)),1E-25)), reused for normalisation.
  * Y capture: immediately after NB_O=NBAD_CUR/NT_O=NTOT_CUR at
    march end, before AA/Aitken mixing and before every exit
    test - the raw map output of one complete pass. Under the
    default (all accelerators F) this equals the committed
    iterate; documented caveat when AA/AIT are armed.
  * R_MAP/RMS_MAP: full NI*NCF scan with
    SCALE = MAX(|X|,1E-4*VMAXF(L)); non-finite deltas excluded
    from RMS and counted (PSD_NNF); floor cells counted
    (PSD_NFLR) but kept in the RMS - no silent exclusion.
  * Snapshots: unit 87, unformatted, per file
    KN_PSWEEP_XY_r<LS>_p<NPO>.bin (NI,NCF,LS,NPO,INWARD header +
    X + Y), bounded by NPO<=KN_PSWEEP_SNHI.
  * Flush: per-solve buffered lines, one OPEN/WRITE/CLOSE under
    !$OMP CRITICAL (KN_PSD_IO) to KN_PSWEEP_DIAG.dat, so armed
    rays cannot interleave; IOSTAT failure drops the diagnostic
    write only. Solve-end NPO/IXIT are appended to every record
    so partial sequences are distinguishable.
- new_main/mod_subs/control_variable_mod.f (+18): KN_PSWEEP_*
  declarations + doc block.
- new_main/mod_subs/rd_control_variables.f (+22): RD_STORE_* for
  the five keys + SET_FG_PSWEEP_DIAG call, mirroring NOCP_PROBE.
- unit_tests/python/kn_psweep_check.py: offline parser + X/Y
  reader + independent R_MAP/RMS_MAP recomputation +
  manufactured-array normalisation tests (--mfg).
- unit_tests/python/kn_psweep_source_test.py: 16 static
  guard-structure checks.

K4. Record format
-----------------

One line per captured pass, list-directed:
  PSD <LS> <CH> <NRF> <NPO_final> <IXIT> <18 cols>
cols: 1 pass NPO; 2 OMEGA_RAY at entry; 3 OMEGA_RAY at march end;
4 production RES; 5 R_MAP; 6 RMS_MAP; 7 max |Y-X|; 8 J argmax;
9 L argmax; 10 NBAD_CUR; 11 NTOT_CUR; 12 PSD_NV; 13 PSD_NFLR;
14 PSD_NNF; 15 NCP_P; 16 NNG_P; 17 GAM sign at argmax (-9 n/a);
18 NFM_L(Jmax) match-freq count (-9 n/a).
CH is 'M' (INWARD=T) or 'P' (INWARD=F). Boundary/depth/match
association of the argmax is derivable offline from J/L plus the
model grid; RES stays the existing local diagnostic - R_MAP is a
new quantity and never relabels it.

K5. Automated tests (pre-run)
-----------------------------

[V] Compile: `make DEBUG=debug cmfgen` clean (ifort
-qopenmp -mkl -traceback -nbs -O0 -check bounds); only
pre-existing NB_Q/NT_Q arg-type warnings at the FLD_Q call
sites. exe_dbg/cmfgen_dev.exe md5 (L_QP-guard build):
8dca57c7a6b3dd7e64eadf41c2db7096.
[M] kn_psweep_check.py --mfg: 8/8 PASS (r_map/rms vs hand
calc, argmax, inf exclusion, floor share, zero delta).
[V] kn_psweep_source_test.py: 16/16 PASS (defaults F, gating,
no PSD->FLD writes, production scale, CRITICAL serialisation,
dealloc pairing, L_QP exclusion, snap bound).
[V] Disabled path: PSD_ARM requires KN_PSWEEP_DIAG; with F the
only added cost is two logical tests per solve + per pass - no
allocs, no I/O, no extra decisions [V by source].

K6. Controlled run
----------------

Run dir: cmfgen_validation/m15_2p82d_KNOP_p53k_psweep
(clone of m15_2p82d_KNOP_p8_base: warm restart from iso_M
lineage; KN_MAT_FREEZE=T, KN_Q_RESID=T, KN_JUMP_DIAG=T,
KN_WATCH ray 83).
VADAT adds: KN_PSWEEP_DIAG=T, KN_PSWEEP_LS=83, PLO=1, PHI=9999,
KN_PSWEEP_SNAP=T, SNHI=4.
Exe: exe_dbg/cmfgen_dev.exe md5 8dca57c7a6b3dd7e64eadf41c2db7096
(ifort -qopenmp -mkl -traceback -nbs -g -O0 -check bounds).
HEAD: 7053969fc54f141255b5d38d5380d01cfffc4dee (dirty: audit doc
+ the three modified sources + 2 new test scripts; nothing
committed).
Inputs: VADAT 527fa9433c6b61c1836880718caea33e; IN_ITS
96105ddfb31e75645e810d1f80a6df10; MODEL_SPEC
b8b94bcc1d3b3434508950ffcce3fd8e; FGINEXT_GS
56cdf03587e790990397b164fd6ccf88; SCRTEMP
66801cf3d02df28774097f2bb6abbe4a.
Job: 13154029, interactive partition, 8 cpus, 1h59m,
OMP_NUM_THREADS=8 (KNOP >=8 CPU rule).

K7. Run status and measurements
-------------------------------

Builds: debug exe (above); production exe also built clean with
the same instrumentation (make cmfgen, ifort -O3 -mavx512*;
one transient ifort bus error on old_grid_module.f, retried OK):
exe/cmfgen_dev.exe md5 f52d8125ad0f91a9e927c007e58af7f2,
13 KN_PSWEEP strings embedded.

Runs staged (all clones of p8_base, warm restart):
- m15_2p82d_KNOP_p53k_psweep: LS=83, SNAP=T SNHI=4, debug exe,
  job 13154029 (interactive, 8 cpu). Cancelled after ~20 min:
  the -O0 model was still inside iteration 287's population
  solve and no refresh (hence no ray solve) had completed in
  this run; refreshes 124-129 in OUTGEN were inherited p8
  content. No records produced.
- m15_2p82d_KNOP_p53k_ray26: LS=26 (the FG-KNOP worst/plateau
  ray), debug exe, job 13154301. Cancelled for the same reason.
- m15_2p82d_KNOP_p53k_prod26: LS=26, production exe
  f52d8125ad0f91a9e927c007e58af7f2, SNAP=T SNHI=4, 8 cpu,
  job 13154834. PENDING on AssocMaxJobsLimit (user's own 8-job
  running cap on small; interactive slot frees when one ends).

Acceptance status (Part J):
 1 disabled by default            [V]
 2 diagnostic-off unchanged       [V] gated blocks; no allocs/IO
 3 diagnostic-on non-interference [V] read-only vs FLD; static
 4 complete X/Y capture           [V] source-verified placement
 5 Python repro of R_MAP/RMS_MAP  [M] manufactured tests 8/8;
                                   snapshot repro pending data
 6 per-pass history indexed       [V] format + flush verified
 7 OpenMP safety                  [V] per-call locals, CRITICAL
                                   flush, no SAVE, no shared state
 8 consecutive production passes  OPEN - job 13154834 queued
 9 audit updated                  done
10 changes uncommitted            yes (M, untracked)

Phase 5.3K-R - staged production run validation
===============================================

Job provenance [M]
------------------
- Job 13154834 (interactive, viper08), started 2026-11-14 18:27:48,
  queued ~1 day on AssocMaxJobsLimit behind the user's 8 running
  small-partition jobs; 16 CPUs allocated.
- git HEAD a70539bd87571e051e6f8aa0840cbc47293aac2d; working tree
  has the Phase-5.3K modifications uncommitted (M fg_j_cmf_nonmono_v1.f,
  control_variable_mod.f, rd_control_variables.f, audit doc; untracked
  kn_psweep_check.py / kn_psweep_source_test.py).
- Executable exe/cmfgen_dev.exe md5 f52d8125ad0f91a9e927c007e58af7f2
  (matches the staged instrumented optimized build, ifort -O3 per
  Makefile_definitions 2026-09-18 flag set).
- Run dir m15_2p82d_KNOP_p53k_prod26/cmfgen/2.82days:
  VADAT md5 5a7fb1f7e65b22a44aae011d3db2ca30, MODEL_SPEC
  b8b94bcc1d3b3434508950ffcce3fd8e, FGINEXT_GS 71f717862e5438bf5716cf6ecbfbdcbf,
  IN_ITS 96105ddfb31e75645e810d1f80a6df10.
- VADAT diagnostic keys (lines 348-353): KN_PSWEEP_DIAG=T,
  KN_PSWEEP_LS=26, KN_PSWEEP_PLO=1, KN_PSWEEP_PHI=9999,
  KN_PSWEEP_SNAP=T, KN_PSWEEP_SNHI=4.
- Warm-start clone of p8_base; model entered at great iteration ~287,
  diagnostics live from the first refresh of the run.
- Model state during capture: iterations ~290-292, i.e. a mature
  warm-started outer state, NOT converged (refresh residual fg_rf
  ~7.6E-3, store drift sw ~2.7E-1, [WARN] refresh residual in OUTGEN).

Record integrity [M/V]
----------------------
- KN_PSWEEP_DIAG.dat grew monotonically; each record wraps over five
  list-directed lines. Parser updated to accumulate continuation lines
  (kn_psweep_check.py parse_records) - the first version wrongly assumed
  one line per record, a parser bug, not a data defect.
- 30 records parsed within the first three solves (iteration ~290):
  each solve is one M-solve (6 passes, IXIT=2) followed by one P-solve
  (4 passes, IXIT=1); NRF column reads 0 throughout because
  LEDG_N_RF only increments under LEDG_ACTIVE, which is off in this
  configuration - records are grouped by file order. [V]
- No duplicates, no partial lines, all fields finite; NVALID=4,961,176
  equals NI*NCF (88 x 56,377) read from the snapshot header - the
  captured field is the complete production solve domain. [M]
- PSD_NNF=0 on every record; PSD_NFLR shrinks 3.3M -> 1.8M -> 0.33M
  (P channel) as the field fills in. [M]

Snapshot reconstruction [M]
---------------------------
- All four KN_PSWEEP_XY_r0026_p000{1..4}.bin snapshots reload cleanly
  (NI=88, NCF=56377, ch=P, inward=F, npo=1..4).
- Python recomputation of R_MAP and RMS_MAP from X/Y matches the
  Fortran records to every printed digit (max |RES-R_MAP|/RES =
  0.000e+00 over all 30 records, argmax J/L identical). [M]
- Pass-1 X is a retained field (|X|max=5.6e-2, 32% cells exactly
  zero): the ~1e4 pass-1 R_MAP is the floor-normalisation signature
  of cells the retained field left at zero (|X|=0 -> scale=1e-4*VMAXF
  -> ratio ~1e4), not a solver instability. [M]

RES versus R_MAP [M]
--------------------
- RES == R_MAP exactly on all 30 records: the production march
  residual already IS the complete-pass max normalised update.
  The diagnostic question "does local RES represent the full-pass
  defect" is answered affirmatively for this ray: identical quantity.
- RMS_MAP sits ~2-3 orders below R_MAP (e.g. P pass4: 6.7e-4 vs
  5.0e-6): the defect is concentrated, not uniform.

Per-solve P-map behaviour [M]
-----------------------------
Solve 1: 1.07e4 -> 1.08e-1 -> 5.9e-3 -> 7.7e-4 (ratios -,0.054,0.13)
Solve 2: 1.04e4 -> 3.9e-2 -> 4.3e-3 -> 6.7e-4 (ratios -,0.109,0.16)
Solve 3: 1.02e4 -> 3.9e-2 -> 2.6e-3 -> 5.7e-4 (ratios -,0.067,0.22)
- Contracting hard once the cold-start floor cells are filled; late-pass
  contraction factor ~0.07-0.22 per pass; NBAD -> 0 by pass 4;
  IXIT=1 (converged class) every solve.
- OMEGA_RAY = 1.00 on every captured pass: the adaptation never
  engaged because the map converged before any damping criterion
  fired. The undamped production P-map on this ray is contractive
  under the captured coefficient state. [M]
- No alternating signature: dominant-cell GAM sign is stable within
  a solve tail (passes 2-4 all gams=-1 on solves 2-3).

Update localisation [M]
-----------------------
- Pass-4 argmax at J=58-75, L~9000-18000: interior cells, NOT boundary.
  Boundary-row maxima on pass 4: outer edge (J<=3) 5.0e-5, inner edge
  (J>=86) 7.3e-6, vs interior max 6.7e-4 - boundary rows are the
  smallest-defect region of the field.
- Summed normalised defect concentrates on J rows 75-79 (inner-mid
  depth) near L=17881; NFM_L(argmax)=1 -> a MATCH_FREQ donor cell,
  but these are also the largest-|X| interior cells, so the match-freq
  association is suggestive, not isolating. [M]
- Pass-1 argmax J=86-87 (near inner edge) with gams=+1 - consistent
  with the cold-fill of retained-field zeros, not a persistent
  boundary defect.

Where the plateau actually lives [M/H]
--------------------------------------
- Every captured P-solve converged (IXIT=1, NBAD->0) in 4 passes
  while the refresh-level OUTGEN diagnostics simultaneously report
  fg_rf ~7.6e-3 and store drift sw ~2.7e-1 with [WARN] flags.
- Therefore the observed production "plateau" is NOT a non-contracting
  per-pass P-map on this ray: the complete-pass defect collapses
  ~1e4 -> ~7e-4 within four passes, undamped.
- The residual that persists (sw ~0.27) is the drift between the
  retained store (committed field) and the freshly re-solved field
  across refreshes - a between-solve, coefficient-state-coupled
  quantity, measured at refresh level, not inside SWEEP_FLD_ORD.
  [M for measurements; H for the interpretation that the store-vs-
  resolve lag is the plateau carrier - consistent with Phase 5.3I's
  exoneration of the J-coupling sector but now established on the
  actual production formal sweep.]

Diagnostic integrity [M/V]
--------------------------
- No production field writes by the diagnostic (static checks 16/16
  [V]); solver control flow identical (records observe, never decide).
- Snapshots never read back by the solver; bounded at SNHI=4.
- Single flushed block per solve under CRITICAL: no interleaving,
  no duplicated records (checked: record indices strictly ordered).
- Ray selection respected: only LS=26 records present.

Limits [C]
----------
- One ray (the production worst), one warm-start window (~iterations
  290-292): does not establish global convergence or rule out
  non-contraction on other rays/states; per-pass capture on all rays
  would need LS=0, not done here.
- nrf=0 (LEDG_N_RF inactive) means solve-to-iteration mapping comes
  from OUTGEN timing, not the record column.
- Capture covers the ordered sweep only; the model ran ORD_SWEEP so
  this is the production path for this run.
- The plateau interpretation remains [H] until a store-vs-resolve
  drift record links sw directly to the per-solve pass-1 X jump.

K8. Blocker note
----------------
The single open item is the scheduled run: ray-26 (worst
FG-KNOP ray, the plateau carrier) is armed in prod26; once it
starts, KN_PSWEEP_DIAG.dat accumulates one line per P (and M)
pass of that ray with RES beside R_MAP/RMS_MAP, and
KN_PSWEEP_XY_r0026_p*.bin gives passes 1-4 snapshots for the
kn_psweep_check.py reproduction. The debug-exe cancellation was
a scheduling call, not a code failure: -O0 -check bounds makes a
single great iteration cost more than the 1h59m interactive
window before the first refresh completes.


Phase 5.3L - refresh-level store versus fresh-solve audit
=========================================================

A. Metric semantics (all [V] from source; callers verified)
-----------------------------------------------------------
Symbol          Definition (numerator / denominator)            Stage        Scope
FG_SW_MAXCHG    in-sweep churn: |fresh-stored|/max(stored,        sweep        per (K,ray,f) write,
                1e-4*VMAX_IP(ray)) over I_NEXT write-backs;       (pre-commit) max over sweep
                snapshotted at refresh start (FG line 8587)
                before FG_GS_MAXCHG is overwritten. This is the
                'sw=' column. CandidateMismatch-family.
FG_GS_MAXCHG    dual-semantics: in-sweep churn until refresh
                runs; then overwritten with KN_RF_RESID.
KN_RF_RESID     MAXVAL(RES_LS, mask diverged) - the refresh's
                per-ray residual max (the 'fg_rf' column); diverged
                rays excluded (FG_GS_NDIV separately).
FG_RF_DCHG      max |committed-retained|/max(retained,1e-4*peak)
                over WRITE_STORE cells - RetainedChange.
FG_RF_PEND      MAXVAL(FG_RAY_DCH, VC in {1,2}) - proposed delta
                of vetoed rays; feeds RF_DCHG=MAX(DCHG,PEND) export
                (line 6447). NOTE: the 'pending=' print at ~8639
                reads it after per-sweep INIT reset (line 3501) and
                before the verdict census at ~9138 recomputes it, so
                it always prints 0 - a reporting-order artifact [V],
                not a missing measurement.
RES_M/RES_P     per-ray solve residual; proven identical to the
                complete-pass R_MAP in 5.3K-R.
FG_RAY_VC       per-ray verdict: -1 unsolved, 0 committed,
                1 jump-reject, 2 non-converged reject,
                3 capped commit (FG_COMMIT_NC), 4 diverged, 7 probe.
FG_RAY_DCH      per-ray delta: committed rays get |new-old| max;
                vetoed rays get the veto JUMP value.
FG_COMMIT_JUMP  whole-ray veto threshold (default 50): veto when
                max over ray's refreshed stores of
                |v1-v2|/max(v2, 1e-4*peak, sibref-floor) > 50.
WRITE_STORE     committed write-back; WRITE_STORE_CAP bounds each
                cell to FG_COMMIT_NC*retained on non-converged bulk
                exits (VC=3). HEAL_STORE zeroes diverged-ray poison
                slots (FG_DIV_HEAL_CAP).

B. Refresh lifecycle [V]
------------------------
1. Retained store S_n: I_M_DEPTH_NEXT/I_P_DEPTH_NEXT (per-ray slot
   arrays; slot l holds lagged frequency l-1 -> freq-1 field has no
   slot, documented at ~9110).
2. FG_KNOP_REFRESH_V1 (line 7812): guards (DO_KNOP_REFRESH,
   FREEZE_I_NEXT, HAS_COMPRESSING, store armed, complete sweep,
   INTEGRAL method); rebuilds KN_CHI_COEF/KN_ETA_COEF monotone
   cubics from KN_CHI_EXT/KN_ETA_EXT; optional KN_FREEZE_COEF pins
   the ext arrays with FNV-1a fingerprint checks.
3. Per-ray (OMP loop) KNOP_RAY_SOLVE_V1: M-solve (SWEEP_FLD_ORD,
   INWARD=T), IM_NI boundary transfer, P-solve (INWARD=F);
   candidates FLD/FLD_P.
4. Acceptance A_n: JUMP = max(FG_JUMP_MAX(M), FG_JUMP_MAX(P));
   JUMP>FG_COMMIT_JUMP -> VC=1, no write; IXIT 3/4/6 under
   FG_COMMIT_CONV=T -> VC=2; non-converged bulk (RES>tol) with
   NC=0.05 -> WRITE_STORE_CAP, VC=3; else WRITE_STORE both, VC=0.
5. Post-loop: KN_RF_RESID=MAXVAL(RES_LS); INEXT_RELAX updated from
   per-solve pass-ratio pool (FG_ADAPT_OMEGA_V2, p90);
   FG_SW_*<-FG_GS_* (sweep churn), FG_GS_MAXCHG<-KN_RF_RESID.
6. In-sweep write-back (population phase, line ~5612):
   store = stored + INEXT_RELAX*(fresh-stored) - relaxed, per cell.
So S_(n+1)=A_n(S_n,F(S_n,C_n)) with A_n = per-cell relaxed write
(in-sweep) + per-ray gated commit (refresh): discontinuous
(whole-ray veto), history-dependent (FG_RAY_SEEDED first-commit
waiver), sibling-dependent (JMP_SIB medians). [V]

C. Measured refresh series, run 13154834, iters ~289-293 [M]
------------------------------------------------------------
iter        sw      fg_rf    rfd(max)  omega  bad    commit census
early      0.267   5.5e-3    0.685     0.537   121
(warm-up) 2195     3.3e-2    0.494             102641 (+372 churn)
           0.728   2.5e-2    0.177     0.557  12127
           0.270   7.6e-3    0.181     0.567   250
           0.270   8.5e-3    0.193     0.598   139  83/1jr/1uns
           0.270   5.5e-3    0.198     0.649   129  83/1jr/1uns
- sw pinned at 0.270 for >=4 consecutive sweeps; argmax cell
  STABLE at (model depth K=3, ray 82, L~17880) - near outer edge -
  with classifier flag 'cyc' (alternating values at a fixed cell =
  period-2 signature of the composite map). [M]
- fg_rf argmax = ray 26 at (J=69, L=17884): the SAME frequency bin
  (~L 17880) as the sw argmax; different ray and depth. [M]
- Retained store moves ~0.18-0.20 per refresh (FG_RF_DCHG) - the
  store is NOT globally stale; it tracks the candidate. [M]
- 83/85 rays commit per refresh; 1 jump-reject (ray 83), 1 unsolved.
- q-residual max = 80.6, argmax ray 83 (j=2,l=1330): vold=2.07e-28
  -> vnew=1.69e-26, values identical to 5 digits across refreshes -
  a frozen vetoed-store floor cell, not a moving defect. [M]
- INEXT_RELAX rises 0.537->0.649: controller sees contraction
  (rho_p90 0.59-0.73), relaxes damping toward 1. [M]
- Moment-side drift = 0.0 throughout. [M]

D. Cause-by-cause assessment
----------------------------
D1 stale retained store: [H-refuted globally] store moves 0.198/
   refresh; but [M-confirmed locally] ray 83's store IS permanently
   stale via the veto trap (below).
D2 whole-ray rejection: [M] ray 83 every refresh - converged solve
   (ixit=0, resid=0, 2 passes) vetoed at jval=165.6>50 on a
   near-floor stored cell (v2~1e-25 vs candidate 1.7e-23; sibling
   ref 1.8e-22 shows the candidate is NOT anomalous). A veto trap:
   the stale store guarantees the same JUMP next refresh. Feeds
   q-resid=80; does not feed sw (sw is pre-commit in-sweep).
D3 commit cap: [M] zero capped commits in the steady window.
D4 coefficient-state changes: [H-partial] retained delta ~0.2/
   refresh in a still-immature model is consistent with a moving
   fixed point; BUT ray 83's frozen cell values and the pinned sw
   argmax say the dominant cells' coefficients are steady.
D5 M-to-P boundary: not isolated; M and P both converge per-solve.
D6 seed/persist flags: FG_RAY_SEEDED set on committed rays; ray 83
   lost its waiver after a past commit - the waiver design assumes
   one unguarded fill suffices; here a poisoned EARLIER generation
   still blocks it. [H]
D7 other - identified: (a) freq-1 slot gap (store never holds
   freq-1 field, documented in source); (b) the composite
   sweep+relaxed-store map has a period-2 cell at (3,82,L17880) -
   'cyc' flag [M]; (c) pending= print-order artifact [V].

E/G. Fixed-point condition
--------------------------
S*=A(S*,F(S*,C)) differs from F(S*,C)=S*: A is discontinuous
(whole-ray veto), history-dependent (seeded flag), and relaxed
(INEXT_RELAX~0.65). Ray 83 sits in a basin where the veto makes
A's fixed point the STALE store while F's fixed point sits ~166x
away - the acceptance map's two fixed-point conditions decoupled
on a floor-dominated cell. [M for the veto/JUMP values; H for the
basin framing.] The sw argmax (3,82) 'cyc' cell is a period-2
point of the relaxed map - consistent with the code's own note of
an HBC sign-flipping dominant mode damped to <=1/2 by
INEXT_RELAX; K=3 is adjacent to the outer-edge closure. [M flag,
H attribution]

H. Part-H answers
-----------------
1. sw~0.27 = max in-sweep fresh-vs-store normalized change at the
   argmax cell (K3, ray82, L~17880) - CandidateMismatch, pre-commit.
2. It is candidate mismatch (in-sweep), NOT retained change; the
   retained change is FG_RF_DCHG~0.198.
3. Per ray: 83 commit fully, 1 (ray 83) jump-vetoed, 1 unsolved.
4. Yes - ray 83 alone retains a stale store, indefinitely (veto
   trap); all other rays' stores track.
5. On the measured vetoed cell, the candidate reproduces itself to
   5 digits across refreshes - frozen; on committed cells the store
   keeps moving ~0.2/refresh with the immature model state.
6. Partially - bulk drift tracks the not-yet-converged outer state;
   the pinned argmax is not coefficient-driven (frozen values).
7/8. Not yet measured - requires the Part-F arms; the
   KN_FREEZE_COEF machinery already exists (ext-array pin +
   FNV fingerprint) and is the vehicle for arm B.
9. Not isolated in this dataset.
10. The plateau decomposes: bulk sw/rfd ~0.2-0.27 = outer model-
    state coupling (immature populations move C_n each iteration,
    store chases a moving fixed point); argmax (3,82) = period-2
    boundary-adjacent mode of the relaxed store map; ray 83 = veto
    trap on floor cells. The raw formal solution is exonerated
    (5.3K-R: contracts ~10x/pass, RES==R_MAP). [M+H]

Missing measurements for a hard verdict: (a) per-cell fresh/stored
values at the sw argmax across sweeps (KNSTOR debug off);
(b) KN_FREEZE_COEF=T arm to separate moving-fixed-point drift from
the period-2 mode; (c) ray 82's sw-cell values to confirm sign
alternation directly.

Phase 5.3M - frozen-coefficient refresh experiment
===================================================

A. KN_FREEZE_COEF audit [V]
---------------------------
- Flag: SET_FG_FREEZE_COEF, VADAT [KN_FREEZE_COEF], default F.
- Mechanism: in-sweep per-frequency snapshots KN_CHI_EXT/KN_ETA_EXT
  are written only while .NOT. KN_COEF_FROZEN; KN_COEF_FROZEN arms at
  the first complete sweep (LK==KN_NCF). FG_KNOP_REFRESH_V1 rebuilds
  KN_CHI_COEF/KN_ETA_COEF monotone cubics from those snapshots at
  refresh head, so all armed refreshes see identical coefficients.
- Completeness: the refresh per-ray solve consumes ONLY the KN_*
  packs (BUILD_FREQ reads KN_CHI_COEF/KN_ETA_COEF/KN_CHI_EXT/
  KN_ETA_EXT exclusively - verified at 12061-12130); the pass-
  invariant pack (DTAU/A0-A4/SR/Q/ETANI/CHINI) depends only on
  captured opacities + static ray geometry; the M->P boundary is
  FLD(NI,:) internal to the map; the seed is the retained store
  (the S_n argument, deliberately not frozen). ETA_EXT at capture
  includes whatever J-dependent scattering the sweep computed -
  frozen as a snapshot. Conclusion: under the flag the refresh map
  F(S_n,C) has C bitwise-frozen - arm B is REAL. [V]
- Limit: the main sweep keeps live CHI_EXT/ETA_EXT (in-sweep churn
  'sw' is NOT frozen by this flag; only refresh-side coefficients
  are). Documented in source at 1169-1180. [V]
- Fingerprinting: KN_COEF_FP sums + 64-bit FNV-1a over rebuilt
  coefs, baseline at first armed refresh, KN_COEF_HSH_BAD counts
  mismatches. Counter bug [V, diagnostic-only]: the arm block
  resets KN_COEF_FP_N=0 on EVERY completed sweep (guarded by
  KN_FREEZE_COEF only, not .NOT. KN_COEF_FROZEN), so 'armed
  refreshes=' prints 1 forever; snapshot writes ARE correctly
  guarded by .NOT. KN_COEF_FROZEN, so the freeze itself holds.
  Hash verifies it: mismatched refreshes=0 across all blocks. [M]

B. Controlled run
-----------------
ARM B (frozen refresh coefficients): clone of prod26 checkpoint
(m15_2p82d_KNOP_p53k_frz82), same exe f52d8125ad0f91a9e927c007e58af7f2,
same FGINEXT_GS 71f717862e5438bf5716cf6ecbfbdcbf, VADAT edited:
KN_FREEZE_COEF=T, KN_PSWEEP_LS=82 (the 'cyc' argmax ray), SNAP T
SNHI 4. Job 13155618 interactive 8cpu on viper08. First-arm run
failed on missing MKL env (raw sbatch --wrap); resubmitted via the
run's own run_psweep.slurm (module intel/2024.0 + mkl/2026.0,
OMP_NUM_THREADS=8) - provenance now identical to baseline.
ARM A = running baseline 13154834 (live coefficients, LS=26).
Warm-start spikes reproduce identically (sw=2195, fg_rf=3.3e-2 on
both first refreshes) - same replay trajectory. [M]
ARMS C (frozen populations) and D (candidate-only): not implemented
- no existing freeze-pops flag for this path; would need production
  control-flow changes, deferred per scope. [C]

C. Frozen-arm measurements [M]
------------------------------
- Fingerprints: coef fp raw=(7.93066325E+07,2.09214271E+05),
  rebuilt drift=0.000E+00, hash mismatches=0 at every armed
  refresh - freeze verified bitwise. [M]
- sw UNDER FREEZE ROSE to 0.73-0.84 (vs baseline pinned 0.269):
  the live-C sweep keeps computing fresh values while the store
  converges toward the frozen-C candidate - gap WIDENS when the
  two coefficient streams diverge. Confirms the bulk sw mismatch
  is a two-coefficient-stream quantity, not solver stagnation. [M]
- sw argmax moved off (3,82) to (54,30,8712) under the freeze. [M]
- Ray 82 solves: M 2-pass ixit=0, P 2-pass ixit=0 - HAS_CP=F fast
  exit (ORD two-pass confirm). R_MAP pass-2 pinned at 0.6209
  across three solves (deterministic corner flush under frozen C);
  pass-1 R_MAP varies 0.29-0.49 (seed varies). [M]
- Committed-store max change 0.684->0.457->0.382->0.387/refresh;
  ray 83 still vetoed (delta 165.6 identical - coefficient-
  independent, as predicted). [M]

D. The 'cyc' cell resolved [M]
------------------------------
Three consecutive refresh candidate generations at (J=3,L=17881),
frozen C, drift=0:
  candidate: 1.9013e-08, 1.9013e-08, 1.9051e-08  (d1=+7e-14,+4e-11)
  seed/store: 3.178e-08, 3.177e-08, 3.162e-08
- The candidate does NOT alternate - it is a stable fixed-C value.
- The store stays ~3.2e-8 while the converged candidate says 1.9e-8:
  a persistent ~60% store-vs-candidate gap at that cell.
- Mechanism: the in-sweep operator writes ~3.2e-8 into the store
  (relaxed write-back) where the refresh operator computes ~1.9e-8;
  the store bounces between the two operators' answers - THAT is
  the 'cyc' signature: a two-operator disagreement at an outer-edge
  (short 4-point ray) cell, NOT a period-2 mode of either map.
- Whole-field check: 587 cells with |d1|>1e-3*scale across three
  generations; 481 same-sign (monotone seed-driven drift) vs 106
  sign-flips - no dominant alternating structure in the candidate.
- Ray 82 geometry: NI=4 - a short outer-edge ray; the argmax J=3
  is genuinely boundary-adjacent. [M]
- Honest caveat: under the frozen arm the in-sweep C is live, so
  the 60% gap mixes operator-disagreement with C-stream drift;
  under baseline (both live) sw=0.27 at this same cell shows the
  gap persists at ~half that size. The structural disagreement is
  demonstrated on arm B; the baseline split needs a baseline
  ray-82 cell capture (not present - PSWEEP targeted 26). [M+H]

E. Ray-83 exit lifecycle verified [M+V]
----------------------------------------
- Record: ixit=0 both channels, npas=2, resid=0.0000E+00 literal.
  Only the .NOT. HAS_CP fast exit produces a literal 0.0 residual
  (ORD sweep: two-pass confirm after the stale-slope-corner pass,
  lines 11749-11760). Ray 83 has NO compressing cells - its
  'solve' is a deterministic operator application, not an iterated
  fixed point. RES=0 is a convention, not a measured convergence.
- Veto cells (KN-CEN): l2=1371 and 1437, stored v2~7.5e-28-1e-25,
  candidate v1~1.2e-25-1.7e-23 -> vr~160 > FG_COMMIT_JUMP=50.
  Sibling ref 1.8e-22 shows the CANDIDATE is sibling-consistent;
  the stored generation is the outlier (stale, near-floor).
- Trap mechanics: stale store -> candidate computes ~166x store ->
  veto -> store stays -> identical candidate next refresh (values
  frozen to 5 digits across refreshes, verified). The veto cannot
  distinguish stale-store from runaway-candidate on floor-scale
  denominators. [M mechanism; D-candidate defect: a converged,
  sibling-consistent field is permanently uncommittable]

G. Part-G answers
-----------------
1. Does frozen C reduce FG_SW_MAXCHG? NO - it rose to 0.84 because
   the freeze is refresh-side only; sw is an in-sweep (live-C)
   quantity. [M] The question as posed conflates the two streams.
2. FG_RF_DCHG under freeze: 0.38-0.68/refresh, not reduced - the
   live-C sweep keeps moving the store. [M]
3. Candidate drift under frozen C: at the argmax cell the
   candidate is CONSTANT (1.9013e-8) - the fixed-C map's fixed
   point is stable; residual candidate drift elsewhere tracks the
   (live) seed/store. [M]
4. Frozen populations: not run. [C]
5. Does the 'cyc' mode survive frozen C? There is no mode to
   survive - the candidate never alternated; the signature was the
   store bouncing between two operators. [M]
6. Alternation exists in neither the candidate nor (measurably)
   the store under frozen C - it was an in-sweep-vs-commit
   artifact of the live baseline. [M]
7. Ray 83 executes the HAS_CP=F fast exit, not a converged solve.
   [V+M]
8. Its veto is dominated by the near-floor stale denominator, not
   a physical candidate anomaly (sibling-consistent). [M]
9. Fraction explained: the refresh-residual plateau (fg_rf~5.5e-3)
   sits above tol mainly on interior cells at L~17880 shared by
   rays 26/82; the sw plateau (0.27) decomposes into live
   coefficient-stream drift + the two-operator cell gap; ray 83 is
   a separate floor-denominator veto trap. Together these account
   for the observed elevated metrics without any non-contracting
   formal map. [M+H]
10. Unresolved: (a) the two-operator gap at outer-edge cells needs
    the baseline-ray-82 fresh-vs-store capture to separate
    operator disagreement from C-stream drift; (b) whether the
    L~17880 hotspot is a MATCH_FREQ/donor stencil artifact;
    (c) whether the store should be allowed to commit through a
    veto whose denominator is provably stale (ray 83). [H/C]

Phase 5.3N - operator-gap decomposition and MATCH_FREQ audit
============================================================

A. Job provenance [M]
---------------------
ARM A baseline: job 13154834 (p53k_prod26), exe f52d8125ad0f91a9e927c007e58af7f2,
  FGINEXT_GS 71f717862e5438bf5716cf6ecbfbdcbf, VADAT 5a7fb1f7e65b22a44aae011d3db2ca30,
  KN_FREEZE_COEF=F, PSWEEP LS=26 PLO=1 PHI=9999 SNAP=T SNHI=4, 8 OMP threads,
  git HEAD a70539bd, warm-start iters ~289+, still running at audit time.
ARM B frozen: job 13155618 (p53k_frz82), same exe/checkpoint, VADAT edited
  (FREEZE_COEF=T, LS=82), same slurm env (intel/2024.0, mkl/2026.0,
  OMP_NUM_THREADS=8). First attempt died on missing MKL env via raw
  sbatch --wrap; the corrected submission uses the run's own script.
  Fingerprint armed at first complete sweep: raw sums
  (7.93066325E+07, 2.09214271E+05), drift=0, hash mismatches=0. [M]
ARM A2 baseline-side (this phase): job 13156275 (p53k_live82), same
  exe/checkpoint/VADAT lineage with FREEZE_COEF=F and PSWEEP LS=82 -
  queued behind the running-job cap; provides the live-C ray-82
  candidate for the same-store comparison against arm B's first
  refresh. [M]

B. KN_COEF_FP_N counter fix [V]
-------------------------------
Defect: the arm block ran on every completed sweep under
IF(KN_FREEZE_COEF) alone, re-zeroing KN_COEF_FP_N so 'armed
refreshes=' printed 1 forever. Fix: the block is now
IF(KN_FREEZE_COEF .AND. .NOT. KN_COEF_FROZEN) - one arm event, the
counter then counts real armed refreshes, fingerprint behaviour
unchanged. Diagnostic accounting only; no numerical path touched.
Source test added to kn_psweep_source_test.py (now 17/17 pass:
'coef-freeze arm guarded once'). Debug exe rebuilt 69f0cbf90c6ce2fb6e4bb656bef7104d;
production exe deliberately left at f52d8125 so live82's provenance
matches arm A/B. [V]

E. MATCH_FREQ at the hotspot [V+M]
----------------------------------
Refresh-solve donor path (KN_ORD_SPV, ~10720): for GAM<0 cells the
compressing donor is NNF2=CEIL(DRF_L(R)/KN_DLN_N(L)) bins blue-ward
(adjacent=1 default; FG_FREQ_MATCH gates larger offsets), read
FRESH from the same-pass FLD at rows JPR1/JPR2 interpolated by
WINX - same-pass donor, in-solve field, not the lagged store.
At the L~17880 hotspot the recorded NFM_L=1: adjacent-bin donor,
no large-offset rematch. The hotspot is NOT a large-offset
MATCH_FREQ artifact; it sits at a boundary-adjacent short-ray row
(J=3 of NI=4). [M nfm=1 records; V stencil]

C/D. Baseline operator-gap decomposition [M]
--------------------------------------------
Measured on arm B (frozen refresh C) over 8 refreshes; ray 82 is
solved (PSWEEP records, ixit=0) and commits (verdict census 83
committed, ray 82 not among rej/unsolved).
Cell (J=3, L=17881, 4-point outer-edge ray):
  seed/store (entering refresh):  3.178e-8 -> 3.177e-8 -> 3.162e-8
  candidate (refresh solve, C0):  1.9013e-8 -> 1.9013e-8 -> 1.9051e-8
D_CANDIDATE_STORE  = cand - seed  ~= -1.28e-8  (~ -0.67 normalised)
D_WRITEBACK_CAND   ~= +1.3e-8 (store returns to ~3.2e-8 between
                     refreshes although every refresh commits 1.9e-8
                     -> the in-sweep relaxed write-back repopulates
                     the cell at its own ~3.2e-8 value)
Interpretation: ray 82's committed store DOES receive 1.9e-8 at
commit, but the next refresh's seed (the in-sweep-written store)
already reads ~3.2e-8 - so between refresh N and N+1 the in-sweep
formal solve wrote ~3.2e-8 where the refresh-pack solve computed
~1.9e-8 at the same physical cell. The 'cyc' signature is this
two-operator tug-of-war on the stored slot, not a period-2 mode of
either map. [M]
Residual caveat: arm B's two C streams differ (refresh C0 frozen
at sweep-1, in-sweep C live). For this outer-edge optically thin
cell the sensitivity is modest, but the exact same-C split needs
the live82 first-refresh candidate (job 13156275, queued behind
AssocMaxJobsLimit). [H bound, M for the measured gap]
Part-D ideal (F(S,C_live) vs F(S,C_frozen) at identical S in one
run): NOT achievable without new production control flow - across
separate runs the store S diverges after refresh-1. The refresh-1
cross-arm comparison is the closest legal version (S_1 identical
by checkpoint lineage). [C]

F. Ray-83 veto, kept separate [M+V]
-----------------------------------
Exit path: HAS_CP=F on both channels -> ORD two-pass fast exit,
IXIT=0, RES=0.0 literal, NPASS=2 per direction. No iterated fixed
point; RES=0 is a reporting convention for a deterministic solve.
Veto cells (l2=1371,1437): retained v2~7.5e-28-1e-25 vs candidate
v1~1.2e-25-1.7e-23; sibling median at same depth ~4.7e-23 (pool
n=84, min 2.4e-28, p10 1.7e-23) -> the RETAINED store sits 100-
1000x below siblings; the candidate is sibling-consistent. JVAL
165.6 > 50 every refresh; identical to 5 digits across refreshes
(coefficient-independent trap). FG_RF_PEND correctly exports 165.6
(verdict-census line); the 'pending=' earlier print reads 0 due to
INIT reset ordering. [M]
Contribution to sw: none - sw is in-sweep churn; ray 83's stale
slot only feeds the q-residual (80.6) and the pending metric. [V]

G. Quantitative decomposition at the target cell [M]
----------------------------------------------------
Component                     value (normalised)      evidence
live vs frozen C streams      sw: 0.27 live vs 0.84   [M] two-stream
                              frozen (mixed stream)   divergence grows
writeback vs candidate        ~+0.68 (3.2 vs 1.9 e-8) [M] seed pinned
candidate vs retained store   ~-0.67 at cell          [M] 3 generations
actual retained change        ~0.38-0.68/refresh      [M] committed-max
MATCH_FREQ contribution       nfm=1 (adjacent donor)  [M+V] not a
                                                        large-offset cause
ray-83 pending                165.6 on ONE ray,       [M] q-argmax floor
                              floor-normalised          cell artifact
Answer: the baseline sw=0.27 at (3,82,17880) is the in-sweep
fresh-vs-store gap whose carrier is the ~0.6-0.7 two-operator cell
disagreement (writeback ~3.2e-8 vs refresh candidate ~1.9e-8),
plus outer-state drift across the bulk field; MATCH_FREQ is
present but only adjacent-donor; ray 83 is an independent floor-
denominator veto trap that does not feed sw. [M for numbers; H for
the exact same-C split pending live82]

H. Part-H answers
-----------------
1. Baseline live-C candidate vs frozen candidate at same store:
   not yet measurable in-field (live82 queued); bounded by the
   arm-B seed/candidate gap. [C open]
2. Coefficient-stream share of the gap: the ~0.6-0.7 cell gap is
   dominantly operator-path (writeback vs refresh pack), with a
   C-stream contribution shown by sw rising when streams diverge.
   [M partial, H split]
3. In-sweep write-back share: ~+0.68 normalised at the cell - the
   dominant local carrier. [M]
4. L~17880 MATCH_FREQ: yes, but nfm=1 adjacent donor only. [V+M]
5. Donor stencil explaining hotspot: no - the stencil is mild;
   the carrier is the writeback-vs-pack disagreement. [M]
6. Operator gap stability: pinned ~0.6-0.7 across 8 refreshes
   (pass-2 R_MAP 0.62 constant; candidate constant 1.9013e-8). [M]
7. Ray-83 vs bulk sw: separate - feeds q-resid/pending only. [V]
8. Measured [M]: all cell values, verdicts, fingerprints, sw/rfd
   series. Hypotheses [H]: the same-C attribution of the 0.6 gap
   and the outer-edge pack convention being the defect site. [D-
   candidate: two code paths answering 60% differently at one
   physical cell is a structural inconsistency worth fixing, but
   the boundary convention may make it expected - undecided]

Phase 5.3O - write-back trace, store ownership, ray-82 corner
============================================================

A. Writer ownership table for I_P_DEPTH_NEXT(KC,LS,slot) [V]
------------------------------------------------------------
Slot convention is UNIFORM (verified in source comments, ~line 3842):
FREQ_IDX_ITER is 0-based, so during frequency k it equals k-1 and
'slot s = frequency s+1' everywhere. SEED_FLD reads FLD(L) <-
STORE(L-1); WRITE_STORE stores FLD(L) -> STORE(L-1); the in-sweep
write at slot FREQ_IDX_ITER during freq k fills slot k-1 with
freq-k data. NO frequency aliasing between writers - the 1.9e-8 vs
3.2e-8 at slot 17880 are genuinely both frequency-17881 values.

Writers (fg_j_cmf_nonmono_v1.f):
  1. In-sweep relaxed write-back (lines ~5617/5674):
     STORE=T2+INEXT_RELAX*(T1-T2), T1=fresh in-sweep I_P/I_M at
     ray point I=J_PNT(K,LS) capped at NONMONO_INEXT_CAP, T2=SNP_P/
     SNP_M = per-frequency snapshot taken on the first call of each
     frequency (NEW_FREQ/NEW slot); DO_ES_ALI re-entries relax
     against the same frozen snapshot so the slot ends at
     SNP+RELAX*(T1_last-SNP). Runs EVERY frequency call of EVERY
     sweep, all compressing depths, all rays 1..NP. Gated on
     HAS_COMPRESSING, 1<=FREQ_IDX_ITER<=NCF_IN, .NOT. FREEZE_I_NEXT.
     THIS IS THE WRITER that restores ~3.2e-8 after each commit.
  2. WRITE_STORE (line 12856): refresh commit, unrelaxed
     STORE(KC,LS,L-1)=capped FLD(JP2,L); once per accepted refresh
     per ray.
  3. WRITE_STORE_CAP (12982): same mapping, change capped at
     FG_COMMIT_NC * retained; for non-converged/capped verdicts.
  4. HEAL_STORE (13335): diverged-ray healing path.
  5. Expansion/init fill (7233-7281): I_F interpolation + EXP()
     extrapolation into slots; setup only.
  6. Zero-init (6881,7325), checkpoint IO (6663/6814/7004),
     KN_NOP probe sums (read-only fingerprints).
Commit-overwrite order: WRITE_STORE commits at refresh; the NEXT
frequency-17881 call of the following population sweep
unconditionally rewrites the same slot (gated only by the
compressing/FREEZE guards, not by any 'recently committed' flag).
[V]

B/C. Target-cell lifecycle and signed deltas [M]
------------------------------------------------
Cell: LS=82, J=3 (ray-local; NI=4 outer-edge ray), K=3 physical
depth (COMP_DEPTH_IDX>=1), store slot 17880 = frequency 17881.
Sequence each iteration:
  refresh: seed <- STORE=3.18e-8 -> FLD(3,17881); frozen-C solve
           -> candidate 1.9e-8; WRITE_STORE -> STORE=1.9e-8
  next sweep, freq-17881 call: SNP=1.9e-8, fresh T1 (in-sweep
           formal solve, live C) ~3.5e-8 (from sw churn 0.84 =
           |T1-SNP|/SNP); STORE=SNP+RELAX*(T1-SNP) ~ 3.2e-8 with
           INEXT_RELAX ~0.7-0.9
  next refresh: seed reads 3.18e-8 -> cycle repeats.
Deltas: D_commit = 1.9e-8-3.18e-8 = -1.28e-8; D_writeback =
+1.3e-8 restoring ~3.2e-8; net per-cycle ~0. Stable across 8
refreshes (seed 3.178->3.162, candidate 1.9013->1.9051e-8). [M]
The two numbers are both legitimate formal solves at the same
physical cell+frequency: T1 comes from the in-sweep formal solver
(live CHI_EXT/ETA_EXT, KC1/KC2 depth-bracket donor reads), the
candidate from the KNOP ordered pack (frozen KN_* coefficients,
JPR1/JPR2 ray-point donor reads, ORD lagged-corner scheme). [V]

E. Ray-82 P-sweep R_MAP=0.62 resolved [V+M]
-------------------------------------------
OMEGA_RAY=1.0 on both passes (records). HAS_CP=F: no GAM<0 point
with a stored compressing slot, so the ORD compressing branch
drops its frequency coupling (IKC=0 -> INEXT_V=0 -> SPV=PK_SR) and
the .NOT. HAS_CP exit ships the field after the 2-pass minimum
with RES=0 by convention. BUT the ordered march retains the
lagged C-corner terms of Knop Paper II Eqs. 31/32 ('the spatially
downwind slope corner ... stay lagged', ~line 1042) - a genuine
in-pass lagged dependency that needs iteration to converge. Pass-2
R_MAP=0.62 pinned across 8 solves, argmax pinned at J=3 (the
boundary-adjacent corner of the 4-point ray) at L~27360 or
L~17880: the corner is NOT converged when the shortcut exits.
So R_MAP=0.62 is a real unconverged update at the corner cell,
not a reporting artifact - the 'single pass is the exact solve'
premise in the .NOT. HAS_CP comment is wrong where the lagged
corner matters. [V for the lagged term, M for the persistent
0.62, H-candidate-defect: the shortcut ships unconverged corners]

F. Donor/closure at the hotspot [V]
-----------------------------------
Compressing-branch donor (KN_ORD_SPV): NNF2=1 adjacent (MATCH_FREQ
would raise it only if DRF_L/DLN>1; nfm records =1), interpolated
(1-WINX)*FLD(JPR1,L+1)+WINX*FLD(JPR2,L+1) FRESH same-pass - but at
J=3 ray 82 the bracketing-depth test gives IKC1=IKC2=0 -> INEXT_V=0
-> SPV=PK_SR: NO frequency donor at the hotspot cell at all. The
hotspot is a pure spatial-corner/boundary-adjacent effect, not a
donor-stencil effect. [V]

G. Ownership assessment [V+C]
-----------------------------
Intended semantics from source: I_NEXT is the lagged GS feedback
array, designed to be continuously relaxed-updated by the in-sweep
write (primary writer, every frequency call); WRITE_STORE commits
are periodic refresh injections that the next in-sweep write
legitimately relaxes back toward the live formal field. The
ownership pattern is INTENTIONAL [C] - the store is not a
refresh-candidate cache. The measured 'discrepancy' is therefore
not a store defect: it is the visible signature of two different
formal operators disagreeing ~60-70% at one outer-edge corner
cell. The real question is which formal value is correct there -
with R_MAP=0.62 at the same corner the refresh-side value is
suspect (unconverged corner), while the in-sweep value is a
live-C V13-class solve. [H]

Part-H candidate defects (NOT implemented):
  1. .NOT. HAS_CP shortcut ignores ORD lagged-corner terms -
     could require a corner-convergence check (e.g. exit only if
     R_MAP<tol or extend passes) for rays whose lagged corner is
     still moving. Location: ~line 10583/11749.
  2. Two-operator disagreement at outer-edge cells - needs the
     coefficient-stream comparison (live82, still queued) to
     separate C-drift from stencil differences.

Remaining open: exact T1 value at the cell (bounded
3.2-3.9e-8 from sw/seed arithmetic; a KN_DBG_LS=82 capture or
live82 start would pin it), the live-C arm comparison (job
13156275 queued behind AssocMaxJobsLimit), and whether the
corner unconvergence also explains part of the ray-26-class
hotspot columns.

Phase 5.3P - HAS_CP shortcut validity and C-corner convergence
============================================================

A. Fast-exit control flow [V]
-----------------------------
SWEEP_FLD_ORD pass loop, exit block (~line 11749):
  after each pass NPO:
    NOP_PASS (probe only)
    IF(.NOT. HAS_CP .AND. .NOT. L_NOPROBE)THEN
       ! 'pass 1 still carries the stale slope corner of the seed,
       !  so confirm on a second pass before a clean exit'
      IF(NPO .GE. 2)THEN
        CLEAN_EXIT=.TRUE.; IXIT=0; RES=0.0_LDP
        RES_FREQ(1:NCF)=0.0_LDP; EXIT
      END IF
    END IF
    IF(RES .LT. KN_RAY_TOL) IXIT=1 exit   <- never reached for HAS_CP=F
    ...bulk/stall/plateau/hot/cap exits (all NPO>=4-8 gated)
Sequence for HAS_CP=F ray: 2 complete passes execute, the pass-2
R_MAP is measured (0.62 on ray 82), then RES is FORCED to 0 and
the ray exits 'clean'. The pass-2 measured defect is computed and
discarded. The code's own comment admits pass 1 carries a stale
slope corner - the two-pass 'confirm' assumes pass 2 converges,
which the instrumentation contradicts. RES=0 and IXIT=0 are
reporting conventions, not measured convergence. [V]
Aitken requires NPO>=3, Anderson NPO>=4 (and both default F): the
exit at pass 2 precedes every convergence accelerator regardless
of settings. [V]

B. Production C-corner derivation [V]
-----------------------------------
P-direction march (INWARD=F, J=NI-1..1), per cell (J,L):
  boundary row NI: FLD(NI,L)=IM_NI(L) (prescribed, fixed in-solve)
  preview: CALL KN_ORD_SPW_FILL(FLD,J-1) fills the ring slot of the
    not-yet-marched row J-1 ONCE per step from its pre-march field
  per L (ascending for GAM>=0):
    SPJ  = KN_ORD_SPV(J,L) = PK_SR(J,L)+PK_Q(J,L)*FLD(J,L-1)  [fresh]
    DSU  = (SPW(J+1)-SPW(J))/DTAU(J+1)     [both fresh]
    DSD  = (SPW(J)  -SPW(J-1))/DTAU(J)     [SPW(J-1)=preview=STALE]
    VNEW = FLD(J+1,L)*A0 + SPW(J+1)*A1 + SPJ*A2 - DSU*A3 - DSD*A4
The only stale read: SPV(J-1,L)=PK_SR+PK_Q*FLD(J-1,L-1) evaluated
on row J-1's pre-march field -> cell (J,L) couples to the PREVIOUS
pass's (J-1,L-1). This is the lagged C term of Knop Paper II Eqs.
31/32, retained by design ('only the spatially downwind slope
corner ... stays lagged'). A single forward march cannot satisfy
it: the discrete system is not triangular under the ordered
ordering. [V]

C. Manufactured test [V]
------------------------
Reproduced the exact stencil in Python (NI=4, NCF=6, all-expanding,
positive pack-like coefficients, boundary row prescribed):
  - |X1_A - X1_B| = 4.88 after ONE pass from two seeds -> the map
    is seed-dependent after one pass: NOT exact.
  - pass-2 update nonzero (R_MAP 2e-2 mild coefficients): the
    'confirm' pass is not converged.
  - both seeds converge to the identical fixed point by pass ~6
    (R_EQ -> 0): unique fixed point; map and equation residuals
    agree at convergence. With stronger corner coupling (ray-82)
    the same structure supports a persistent oscillation.
Script: /tmp/p53p/ccorner_test.py. [V]

D. Ray-82 pass histories - REVISED MECHANISM [M]
-------------------------------------------------
Per-solve snapshots at the hotspot (J=3, L=17881), both arms:
  live82 (13156275): seed 2.63e-8 -> pass1 1.90e-8 -> pass2 3.20e-8
  frz82  (13155618): seed 2.14e-8 -> pass1 3.18e-8 -> pass2 1.90e-8
The cell ALTERNATES ~1.9e-8 <-> ~3.2e-8 EVERY PASS - a period-2
mode of the pass iteration at the C-corner, with opposite parity
in the two arms because the seeds differ. Pass-2 R_MAP 0.62-0.68
pinned across 8+ refreshes both arms, argmax J=3 (L=17881 or
L=27360). The 2-pass shortcut ships whichever parity it lands on
(frz: 1.90e-8; live: 3.20e-8).
CORRECTION to 5.3N/5.3O: the 'two-operator gap' measured seed
(~3.2e-8, one cycle phase) against pass-2 output (~1.9e-8, the
other phase) - the dominant structure is the intrinsic pass-level
oscillation of the refresh map at this corner cell, not primarily
a live-vs-frozen or writeback-vs-candidate disagreement. The
in-sweep write-back still writes its own value (~2.6-3.2e-8), but
the alternation signature ('cyc') belongs to the ordered map. [M]

F. Live vs frozen comparison [M]
--------------------------------
live82 job 13156275 now running (8+ solves): same exe f52d8125,
same checkpoint lineage, LS=82, FREEZE_COEF=F.
  P pass-1 R_MAP 0.28-0.31, pass-2 0.68-0.81; argmax J=3 L=17881.
  sw=0.2699 (baseline level) vs frz82 0.84 - the churn magnitude
  IS coefficient-stream sensitive, as expected.
  pass-2 corner R_MAP ~0.62 frozen vs ~0.68 live - the corner
  oscillation is coefficient-INSENSITIVE: a structural map mode.
  Candidate at pass-2: frz 1.90e-8 vs live 3.20e-8 (different
  parity, not different fixed point).
  Seeds differ between arms (2.63 vs 2.14e-8) -> NOT a same-input
  comparison; quantitative D_LIVE_FROZEN at identical S remains
  unavailable [C].

G. Part-G answers
-----------------
1. HAS_CP fast exit valid for every HAS_CP=F ray? NO - it is
   invalid wherever the lagged C-corner term is active; pass-2
   RES is forced to 0 without measurement. [D]
2. Production C-corner lagged? YES - preview SPV(J-1) on stale
   field, inside DSD, every marched cell. [V]
3. One pass sufficient with the dependency active? NO -
   manufactured test: seed-dependent after one pass. [V]
4. Ray-82 R_MAP=0.62 genuine map defect? YES - a real
   period-2 oscillation at the corner cell, per-pass |delta|~60%.
   [M+V]
5. Manufactured test shows shortcut can terminate before the
   equation is solved? YES - demonstrated (one-pass not exact;
   fixed point needs ~5 passes here). [V]
6. Ray-82 result sensitive to live vs frozen coefficients? NO for
   the corner mode (~0.62 vs ~0.68); YES for store churn
   (sw 0.27 vs 0.84). [M]
7. Is a C-corner convergence check required? For correctness of
   the reported residual: yes - the current exit reports RES=0
   for a field moving 60%/pass. Whether the unconverged corner
   damages the population solve is a separate question (the
   committed candidate is simply the pass-2 parity value). [D-
   residual misreport; H on downstream impact]
8. Minimum safe criterion: for .NOT. HAS_CP rays, gate the clean
   exit on the MEASURED pass-2 R_MAP (e.g. require R_MAP <
   KN_RAY_TOL or a corner-specific tolerance) instead of forcing
   RES=0; continue iterating (or report the true residual) when
   the corner still moves. Any threshold choice needs the
   manufactured/stencil analysis in C as its basis - not
   implemented in this phase. [H]

H. Provenance note: all three jobs preserved (13154834 running
74min+, 13155618 running 68min+, 13156275 running 29min+ at the
time of the live82 captures). Nothing committed or pushed.

Phase 5.3Q - HAS_CP=F exit verification and C-corner convergence
================================================================

A. Control flow, both branches [V]
----------------------------------
SWEEP_FLD (legacy alternating, ~10583):
  after EVERY pass: IF(.NOT. HAS_CP .AND. .NOT. L_NOPROBE) ->
    CLEAN_EXIT=T, IXIT=0, RES=0, RES_FREQ=0, EXIT.
  One complete pass. Comment: 'the single pass is the exact solve'.
SWEEP_FLD_ORD (production ordered, ~11749):
  after each pass: IF(.NOT. HAS_CP .AND. .NOT. L_NOPROBE) ->
    IF(NPO .GE. 2) -> same forced exit.
  Two complete passes; pass-2 R_MAP measured then overwritten by 0.
Caller (KNOP_RAY_SOLVE_V1 ~9800): JUMP veto -> FG_COMMIT_CONV test
  'IXIT in {0,1,2} AND (RES<=KN_RAY_TOL OR (IXIT=2 bulk))'. Forced
  IXIT=0,RES=0 PASSES this test -> WRITE_STORE -> verdict 0
  (committed). So a ray whose last pass moved 62% is committed as
  converged. [V] This is a reporting defect with a downstream
  commit consequence, not a benign convention. [D]

B. R_MAP semantics [V]
----------------------
RES accumulates, per cell during the in-place march,
  T1=|VNEW-VOLD|/MAX(VOLD,1E-4*VMAXF(L)), RES=MAX(RES,T1).
Since VNEW is the map output at that cell evaluated on the
Gauss-Seidel current state, a complete-pass RES is the map defect
||G(X)-X|| in a per-cell weighted max-norm (weight
1/MAX(|X_c|,1E-4*VMAXF)) - a true complete-map defect, not merely
a local update: the accumulation covers every cell and the
normalisation is the production one. R_EQ (the simultaneous
discrete equation residual, all terms evaluated on X itself)
differs at non-fixed-point iterates and agrees at the fixed
point - verified numerically below. [V]

C. Period-2 verification status
-------------------------------
Three-point sequences (seed->pass1->pass2) per solve exist for
both arms (5.3P); a rigorous period-2 classification needs
repeatable alternation with decreasing d2 over MORE passes, which
the 2-pass shortcut cannot provide. Job 13157300 (p53q_probe82,
KN_NOCP_PROBE=T, KN_NOCP_PROBE_LS=82, same exe/checkpoint) is
running the ordinary iteration on ray 82 to supply the real pass
sequence; results pending at section end. [C pending]

D. Manufactured test with literal ring semantics [V]
----------------------------------------------------
unit_tests/python/kn_ccorner_test.py replicates the exact
production stencil including the mod-4 SPW ring, stale preview
fill, per-pass boundary enforcement (FLD(NI,:)=IM_NI each pass)
and ring aliasing at the ray ends (DS(J+1) at J=NI-1 reads the
'row NI+1' slot = row-1's slot; DS(1) reads a stale ring slot).
Findings:
  - One pass is NOT exact: |X1_A-X1_B|=4.5e-2 seed-dependence. [V]
  - Map converges to an identical fixed point from both seeds
    (|XA-XB|=0 at pass ~7). [V]
  - The converged iterate satisfies the simultaneous stencil
    residual at roundoff (the fixed-point equation == the
    simultaneous equation when stale reads equal fresh). [V]
  - With mild coefficients the corner converges monotonically
    (d1,d2 -> 0); production's 60% alternation needs the stronger
    production coupling - the probe run measures it directly. [M]
Ring-alias note: at J=NI-1, DSU reads slot s(NI+1)=s(1) =
row-1's previous-pass source; at J=1, DSD reads a stale ring
slot. The boundary-adjacent lagged terms are intrinsic to the
mod-4 ring, not just the J-1 preview. [V]

E. Acceleration feasibility [V+M]
-------------------------------
Source: Aitken guarded step requires NPO>=3; Anderson requires
NPO>=4 + NH_AA>=2 history + safeguards (finite, positive,
step-cap); both default OFF (KN_ORD_AITKEN/KN_ORD_AA=F). Under
the HAS_CP=F shortcut the exit at pass 2 precedes all of it -
no accelerator can engage. Manufactured-test comparison (mild
coefficients): plain FP 6 passes, Aitken 7-9, depth-2 AA 6,
omega=0.5 relaxation 35-39 passes - on a contractive map
acceleration adds nothing; the interesting regime is the
non-contracting corner where only the probe run gives data.
Safeguarded-AA-with-restart is the only viable candidate if the
mode is period-2 (AA on a period-2 orbit identifies the mean;
Aitken on a 2-cycle is the correct local accelerator - but BOTH
require the exit not to fire at pass 2). [V for gates; M for the
mild-map table; H for the period-2 regime]

F. Proposed safe exit policy (proposal only, NOT implemented)
-------------------------------------------------------------
1. Never overwrite a measured defect: keep the forced-zero
   assignment out of RES; if a suppressed-gate convention is
   needed for maturity bookkeeping, export it to a separate
   variable (e.g. RES_GATE) rather than falsifying RES.
2. HAS_CP=F rays exit clean only when the MEASURED R_MAP <
   KN_RAY_TOL (the same tolerance every other ray uses); i.e.
   simply delete the shortcut's unconditional exit and let the
   existing tolerance/stall/plateau/cap exits run. The 'pass-1
   seed transient' concern the comment raises is already handled
   by NPO>=2 gating of the tolerance test - a second measured
   pass below tolerance is a legitimate clean exit.
3. Non-converged at KN_RAY_MAXPASS -> existing IXIT=4 fail
   semantics apply unchanged.
4. Downstream effect: a corner-oscillating ray then reports its
   true defect -> FG_COMMIT_CONV routes it to nc-rej/capped
   instead of committing a mid-oscillation value as converged.
   Store/commit semantics otherwise unchanged.
5. Tolerance basis: KN_RAY_TOL is the production residual
   tolerance already used for every coupled ray - no new
   threshold invented. [proposal; awaits approval]

H. Tests added [V]
------------------
unit_tests/python/kn_ccorner_test.py: ring-faithful manufactured
test (unique fixed point, one-pass non-exactness, physical
iterate checks).
kn_psweep_source_test.py extended to 21 checks pinning: legacy
one-pass forced exit, ORD two-pass forced exit, Aitken/AA pass
gates, coef-freeze arm guard. All pass. The HAS_CP semantics
tests pin the CURRENT (defective) reporting so a future fix must
update them deliberately.

Pending: probe job 13157300 pass-by-pass ray-82 sequence ->
final period-2 classification and the live/frozen corner
comparison at higher pass counts.


C-addendum. Period-2 evidence and accelerator behavior [M+V]
-----------------------------------------------------------
frz82 snap_hist (three refresh generations, cell J=3 L=17881):
  pass-1 output: 3.1780e-8, 3.1775e-8  (pinned to 0.02%)
  pass-2 output: 1.9013e-8, 1.9013e-8, 1.9051e-8  (pinned to 0.2%)
A stable two-state orbit {3.18e-8, 1.90e-8} across refreshes, with
off-orbit seeds (2.14e-8) landing on the cycle in one pass.
Consistent with a genuine period-2 mode of the ordered map at the
corner; a converging oscillation would show amplitude decay across
generations. The remaining discriminating measurement (pass-3+
return) comes from probe job 13157367 (queued). [M]

Accelerator test on a genuine period-2 orbit (logistic r=3.5,
stable {0.50,0.87} orbit; script /tmp/p53p/period2_accel.py):
  plain FP:     d1~0.4 persistent, d2~0.05 - 2-cycle signature
  Aitken:       d1 -> 1e-15 - BREAKS the cycle, converges to the
                underlying (unstable) fixed point of the map
  Anderson m2:  no improvement - stays on the orbit (a shallow
                secant history cannot resolve the alternating
                residual structure; production AA would also skip
                its guarded mix here)
  omega=0.5:    converges to the same fixed point - relaxation
                changes the effective map and dissolves the 2-cycle
Implication for the exit policy: the measured 2-cycle needs either
more passes with a cycle-breaking step (guarded Aitken, or the
existing OMEGA_RAY < 1 damping) or an honest non-converged
verdict; the current shortcut does neither. [V test; M prod data]

F-addendum. Probe job status
----------------------------
First probe attempt (13157300) ran ~7 min in the wrong directory
(inherited --chdir=p8_base) and was cancelled; it flipped
COMP_F->F in that model's VADAT and appended startup output to its
OUTGEN (standard requeue behavior; no tracked file changes, no
p8_base job was running concurrently). Resubmitted correctly as
13157367 (KN_NOCP_PROBE=T, LS=82, PSWEEP active, same exe
f52d8125, checkpoint lineage of live82/frz82); queued behind the
association job limit at audit time. [M]

Phase 5.3R - HAS_CP=F false-convergence patch plan (UNAPPLIED)
=============================================================

A. Verified semantics [V]
-------------------------
SWEEP_FLD (legacy, ~10583): after EVERY pass, if .NOT. HAS_CP ->
  CLEAN_EXIT, IXIT=0, RES=0, RES_FREQ=0, EXIT. One pass, no gate.
SWEEP_FLD_ORD (~11749): same branch but gated NPO>=2 -> two passes.
Both then reach the ordinary exits (tol/stall/bulk/plateau/hot/
cap/MAXPASS) only if the shortcut does not fire.
RES_LS(LS)=RESID=MAX(RES_M,RES_P) -> KN_RF_RESID=MAXVAL(RES_LS)
-> fg_rf and the refresh maturity gate. Forced-0 rays are thus
invisible to the refresh residual AND to KN_HOT_CNT (per-frequency
hot census reads RES_FREQ, also forced to 0). [V]

B. Caller consequences [V]
--------------------------
FG_COMMIT_CONV requires IXIT in {0,1,2} AND (RES<=KN_RAY_TOL or
IXIT=2 bulk exit). IXIT=0/RES=0 passes unconditionally ->
WRITE_STORE -> verdict 0 committed. So the false convergence is
not cosmetic: a corner-oscillating ray commits a mid-oscillation
candidate as converged every refresh. Diverged rays (NDIV=1) keep
NXIT=5/verdict 4/HEAL_STORE independent of the shortcut. [V]

C. Proposed minimal correction (for review; NOT applied)
--------------------------------------------------------
SWEEP_FLD_ORD (~11749) - delete the .NOT. HAS_CP shortcut block
entirely:
  - converged rays exit via the existing 'RES .LT. KN_RAY_TOL ->
    IXIT=1' test (measured, may fire at pass 1 or 2 - strictly
    better than the blind 2-pass exit);
  - non-converged rays continue: stall/bulk/plateau/hot/cap/
    MAXPASS exits apply unchanged;
  - Aitken (NPO>=3) / AA (NPO>=4) become reachable when enabled.
SWEEP_FLD (~10583) - delete the identical one-pass shortcut:
  measured tolerance test then applies at pass >=1.
Resulting semantics: IXIT 0/1 both mean 'measured converged'
(0 becomes unused), 2 bulk, 4 maxpass, 5 div, 6 plateau; RES is
always the measured pass residual; RES_FREQ always real.
No tolerance changes; KN_RAY_TOL reused. No stencil changes;
iteration control only. Alternative (rejected): keep the branch
with an added RES<KN_RAY_TOL gate and IXIT=0 - same behavior but
preserves a redundant exit class and a redundant branch; the
plain deletion is smaller.

D. Consequences to flag (production-visible, all intended)
----------------------------------------------------------
1. fg_rf / KN_RF_RESID will honestly include HAS_CP=F rays'
   residuals (ray 82 ~0.6, ray 83 floor-cells): the reported
   refresh residual rises where such rays exist - truthful, but
   the maturity gate may stay open longer on affected models. [C]
2. KN_HOT_CNT will count these rays' per-frequency residuals:
   the corner frequency (~17881) will accumulate hot counts -
   consistent with the hot-window machinery's purpose. [C]
3. Genuinely non-converged HAS_CP=F rays now run to a real exit
   (stall/plateau/MAXPASS) instead of exiting at pass 2 - slightly
   more passes on such rays only; converged rays exit at the same
   or fewer passes. Probe job 13157367 measures this directly. [M pending]
4. HAS_CP=T rays: completely unchanged (the branch never fires).
   FREEZE_I_NEXT/JFNK/JFNK-probe paths: untouched. Stock/SOB/JFNK
   paths: untouched (all changes inside SWEEP_FLD/ORD KNOP path).

E. Test plan
------------
kn_psweep_source_test.py (extend): replace the two pinned
'forced RES=0' checks with checks that the .NOT. HAS_CP shortcut
is ABSENT and the ordinary tolerance test is present; keep the
accelerator-gate and arm-guard checks.
kn_ccorner_test.py (extend): apply the production-exit logic to
the manufactured map - (a) converged case exits clean at measured
tolerance; (b) oscillating case does NOT report RES=0, reaches a
non-converged status at pass cap; (c) d2 test distinguishes a
true 2-cycle from slow convergence.
Regression: HAS_CP=T rays unchanged by construction (branch
removal only affects the .NOT. path); verify a normal model run
still reports identical first-refresh behavior for coupled rays.

F. Risks
--------
Low: only non-converged no-cp rays run more passes (bounded by
MAXPASS); verdicts route to existing nc-rej/capped classes.
Medium: honest residual visibility may hold the KNOP maturity
gate open on models with corner-oscillating rays - a reporting
consequence of the fix, not a new numerical behavior; flag for the
user since it changes observed run timelines.
Probe 13157367 (queued) gives the pass-3+ sequence validating that
the ordinary exits actually fire and bound the iteration.

Status: AWAITING APPROVAL - no production change applied.


Phase 5.3R - implementation record (approved minimal correction)
==============================================================

Applied (uncommitted): deleted the two forced-clean HAS_CP=F exit
blocks and updated the stale comments that described them.

Edits [V]
---------
plane/fg_j_cmf_nonmono_v1.f, three sites, no numerical change:

1. SWEEP_FLD (was ~10583-10596): removed
     IF(.NOT. HAS_CP .AND. .NOT. L_NOPROBE)THEN
       CLEAN_EXIT=.TRUE.; IXIT=0; RES=0.0_LDP
       RES_FREQ(1:NCF)=0.0_LDP; EXIT
   The measured 'IF(RES .LT. KN_RAY_TOL) -> IXIT=1' test now runs
   from pass 1 for every ray.

2. SWEEP_FLD_ORD (was ~11749-11762): removed the NPO>=2-gated copy
   of the same block. Same consequence.

3. Comment updates (KN_NOCP_PROBE module doc ~1346 and both
   NOP_PASS call-site comments): the flag no longer bypasses
   anything at the exit sites - probe and ordinary rays share the
   measured exits; L_NOPROBE still records NOP_PASS telemetry and
   suppresses WRITE_STORE for the probe ray.

Nothing else changed: KN_RAY_TOL, KN_RAY_MAXPASS, all stencil
coefficients, ring semantics, donor indexing, MATCH_FREQ,
Aitken/Anderson mathematics, store write-back, INEXT_RELAX,
OMEGA_RAY, population equations, FG_COMMIT_CONV policy are
bit-identical. HAS_CP is still computed (kept as a documented ray
classification; it is now write-only - removing the computation
would have broadened the diff for zero benefit).

Verified control flow after the edit [V]
----------------------------------------
- Per-direction init IXIT=-1, CLEAN_EXIT=.FALSE. (~10298, ~10972);
  both pass loops are unbounded DO with NPO=NPO+1, terminated only
  by explicit exits. The NPO>=KN_RAY_MAXPASS -> IXIT=4 exit is the
  guaranteed terminator, so no path can return an unset status.
- Every success exit now passes through IF(RES .LT. KN_RAY_TOL)
  (IXIT=1) or the floored bulk exit (IXIT=2). IXIT=3/4/5/6/9 are
  non-converged classes; FG_COMMIT_CONV admits only IXIT in {0,1,2}
  AND RES<=KN_RAY_TOL (or the IXIT=2 bulk class), so no
  non-converged ray reaches the ordinary WRITE_STORE commit.
- HAS_CP=T rays: the deleted branches never fired for them;
  control flow is unchanged.

Measured/experimental consequences to watch [H]
-----------------------------------------------
- fg_rf / KN_RF_RESID now honestly include HAS_CP=F ray residuals
  (ray 82 ~0.6); the refresh maturity gate may stay open longer on
  models with corner-oscillating rays.
- KN_HOT_CNT now counts these rays' per-frequency residuals; the
  L~17881 column will correctly accumulate hot counts.
- Genuinely non-converged HAS_CP=F rays now run to a real exit
  (stall/plateau/MAXPASS) instead of exiting at pass 2; converged
  rays exit at the same or fewer passes. Probe job 13157367
  (queued) measures this path on the old exe.

Tests [M]
---------
unit_tests/python/kn_psweep_source_test.py (25 checks, all PASS):
  the forced-exit branches are asserted ABSENT in both sweeps; the
  measured tolerance test, maxpass IXIT=4 and IXIT=-1 init are
  asserted present; accelerator gates (NPO>=3 / NPO>=4) unchanged.
unit_tests/python/kn_ccorner_test.py (PASS): added an exit-policy
  simulator mirroring the post-fix semantics - converged run exits
  IXIT=1 through the measured test; a persistent 2-cycle never
  meets tolerance, exits IXIT=4 at the pass cap with the last
  measured RES intact (never overwritten by 0), and fails the
  FG_COMMIT_CONV admit test; the d1-finite + d2-small period-2
  classifier accepts the r=3.2 logistic orbit and the production
  corner signature, and rejects converged and chaotic sequences.
unit_tests/python/test_knop_lagged_field.py: updated the stale
  pins - the ray-83 lifecycle item now documents the removed exit,
  and the probe pins assert the branch is gone while store-write
  suppression and telemetry remain.

Regression battery [M]
----------------------
test_knop_lagged_field: all PASS (incl. probe-diagnostics section)
test_knop_ray:              43 passed, 0 failed
test_fg_gs:                 36 passed, 0 failed
test_gs_drift_diag:         all PASS
test_iter_conv_gate:        13 passed, 0 failed
test_obs_nonmono_dispatch:  all PASS
test_lam_handoff:           51 passed, 0 failed
test_jfnk_opint:            36 passed, 0 failed
test_jfnk_regressions:      all PASS
test_jfnk_t_fixed:          1 passed, 0 failed
test_mom_gs/solve/sweep/aa: 34/137/8/62 passed, 0 failed
test_obs_frame(+cov):       11+15 passed, 0 failed
test_ng_abstain:            17 passed, 0 failed
test_extrap_scr:            7 passed, 0 failed
Pre-existing gaps (not this change): test_es_conv / test_es_resid
gate on obs/ sources never archived into lib/*.a (earlier-phase
files); test_jfnk_solver exits quietly 0.

Builds
------
Production archive member fg_j_cmf_nonmono_v1.o rebuilt into
lib/libsubs.a with the standard viper ifort flags (ifort -O3
-mavx512*): compiles clean, only the pre-existing #6075 warnings
on the verify-pass NB_Q/NT_Q dummy args. exe/cmfgen_dev.exe is
untouched (md5 f52d8125...) - running jobs unaffected.
Debug build (build_debug.sh -> exe_dbg/) in progress at audit
time; exe_dbg is for future validation runs only.

Preserved jobs: 13154834 (ended), 13155618, 13156275, 13157367
(probe, still queued) - none cancelled, restarted or modified.


Phase 5.3S - post-correction production validation
================================================

A. Provenance [V]
-----------------
- git HEAD: 4cd9c0cfd3b80f3a03d781acba250294cb2942fe; all work
  uncommitted (16 modified/untracked entries at audit time).
- source md5: plane/fg_j_cmf_nonmono_v1.f = cb5d4b22b879b5f03d6862c44ffef24a
- production exe: exe/cmfgen_dev.exe = f52d8125ad0f91a9e927c007e58af7f2 (untouched)
- corrected debug exe: exe_dbg/cmfgen_dev.exe = c3d8b584b0f7d7a6666ee58067941e95
- KN_RAY_TOL=2.0E-3 (DATA default), KN_RAY_MAXPASS=60, KN_RAY_STALL=0.98.

B. Old-exe baselines harvested [M]
----------------------------------
frz82 (13155618, ended TIMEOUT at 2h, output preserved):
- ray83 KN_WATCH: 33 refreshes, every line ixit=0 npas=2 resid=0
  - the forced-clean signature the 5.3R fix removes.
- ray82 P-direction pass records (KN_PSWEEP_DIAG.dat): pass-1
  RES=0.482 argmax (J=3,L=17881), pass-2 RES=0.621 argmax moving
  to (J=3,L=27360); ~27.4K above-tolerance cells per pass. M-dir
  converges: 0.13 -> 0.0026 (near the 2e-3 tolerance).
- fg_rf ~2.65e-2, sw=0.840, bad~2900-3100, churn=0, capped=0.
live82 (13156275, ended ~2h wall):
- ray83 same ixit=0 npas=2 resid=0 signature.
- fg_rf ~5.5e-3, sw=0.269, bad~130.

C. Validation arms submitted
----------------------------
- p53s_fix82 = job 13158010: clone of frz82 run dir (same
  checkpoint, same FGINEXT_GS store, identical VADAT with
  KN_FREEZE_COEF=T, KN_PSWEEP on ray 82, KN_WATCH on ray 83),
  run with the corrected debug exe c3d8b584 via a staged cmfdist
  (<run>/dbgdist/exe/cmfgen_dev.exe; production exe untouched).
  This is a direct A/B against frz82 under identical input state.
- p53q_probe82 = job 13157367: still queued behind
  AssocMaxJobsLimit. On the old exe its probe ray already runs
  the measured-exit path (L_NOPROBE skips the removed exits), so
  it validates the same iteration behaviour the fix produces for
  every HAS_CP=F ray.

Expected [H]: under the corrected exe, ray-82 P-direction can no
longer exit at pass 2 with RES=0; with a ~0.5-0.6 measured defect
at ~300x tolerance it should run to a stall/plateau/MAXPASS exit
(IXIT 3/4/6) and a non-converged verdict, and fg_rf should
honestly include its residual (~0.6) instead of suppressing it.
Ray 83's floor-dominated cells should likewise produce a real
residual and a non-clean exit. Whether that changes commit
behaviour (nc-rej vs capped commit) is the measurement.

Status at audit time: both jobs queued; results pending.
