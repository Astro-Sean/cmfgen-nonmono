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
