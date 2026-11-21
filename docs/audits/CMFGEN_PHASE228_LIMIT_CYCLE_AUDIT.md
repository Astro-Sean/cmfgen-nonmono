# Phase 2.28 - ray-81 stall identified as a period-2 limit cycle

Date: 21-Nov-2026. Model m15_2p82d_KNOP (all-compressing sigma, mixed
ray topology). Question: why is production ray 81 the persistent worst
ray when it has zero compressing cells, and does its difficulty come
from projected-velocity physics, direction handling, or the ordered
KNOP fixed-point iteration itself. Evidence: the phase226 retained
replay histories (KN_REPLAY_HST_L081/L082/L084.bin), the phase226
production deck summaries, OUTGEN telemetry, the decoded FGINEXT_GS
store, and the production solver driven directly through
unit-test-style drivers. No production source, solver setting or
acceptance policy was changed.

Verdict up front: the stall is a nonlinear period-2 limit cycle of the
ordered fixed-point map, proven bit-exact in the replay histories.
It is not a projected-gradient, multi-resonance, direction-selection,
storage or seeding defect. The physical velocity question is answered
in the negative for ray 81: its geometry is the simplest in the model.

---------------------------------------------------------------
1. THE RAY IS PROVABLY ALL-EXPANDING
---------------------------------------------------------------

Geometry. Ray 81 is a tangential ray (P(LS)=R(86-LS) for LS>KN_NC=15,
fort.188 echo; comp_j_blank.f caller). NI_RAY(81)=5, so its ray points
are model grid nodes R(1)..R(5) - no interpolation ambiguity.

SIGMA. SIGMA_TRUE reproduced from RVTJ V data matches the OUTGEN
echo exactly (min -155.61, max -1.037; all 70 depths sigma<-1).
Production GAM uses monotone-cubic interpolation of sigma onto ray
points (no overshoot; fg_j_cmf_nonmono_v1.f ~4170-4300), then

  GAM(I,LS) = 3.33564E-6 * (V/R) * (1 + sigma*mu^2)     (l.4451-4458)

which is the exact signed projected gradient
  (1/c) d(v_parallel)/ds = (v/R)(1 + sigma mu^2)
with sigma = dlnv/dlnr - 1, mu = Z/R. Analytic check:
  dv_parallel/dz = mu^2 dv/dr + (1-mu^2) v/r
  (1/c) dv_parallel/dz = (v/R)(mu^2 (r/v) dv/dr + 1 - mu^2)
                       = (v/R)(1 + sigma mu^2).              Verified.

For ray 81, sigma_true in [-1.04,-1.0] at depths 1-5 and mu^2 <= 0.14
at the outer ray points, so T2 = 1+sigma*mu^2 >= 0.85 everywhere:
GAM > 0 at every point, purely expanding along the ray, zero
compressing cells, no T2 sign bracket. The multi-root census (35 rays
with multiple turning points) does not include rays 62-85.

Consequence: none of the compressing-cell machinery (MATCH_FREQ slot
offsets, GAM<0 donor reads at L+NNF2, spatial interpolation through
JPR1/JPR2/WINX) is ever engaged on this ray. Hypotheses A (genuine
multiple-resonance complexity) and C (donor-direction error at sign
changes) are refuted for ray 81 by construction.

---------------------------------------------------------------
2. THE ORDERED MARCH ON AN ALL-EXPANDING RAY IS NEAR-DIRECT - YET
   THE PRODUCTION SOLVE DOES NOT CONVERGE
---------------------------------------------------------------

A faithful Python replica of SWEEP_FLD_ORD (same coefficient limits,
same lagged SPW preview, same DS clamp, same relaxed commits) was
driven with ray-81-faithful geometry (NI=5, all GAM>0), uniform
opacity, a Lyman-edge opacity step, a 20-rung edge picket, and
pathological seeds: it converges in 1-11 passes in every regime, and
the production KNOP_RAY_SOLVE_V1 driven with the same class of inputs
exits IXIT=1 at pass 11. The linear parts of the all-expanding map
are a direct solve; the only genuinely lagged terms are

  (i)  the C-corner source slope: DSU=KN_ORD_DS(J,L) /
       DSD=KN_ORD_DS(J+1,L), which reads SPW(J+1)/SPW(J-1) built by
       KN_ORD_SPW_FILL from the not-yet-marched row's pre-march state
       (l.13040, 13260; slope at l.12657-12694), and

  (ii) the tangent boundary feedback: for LS>KN_NC the P-field row NI
       is not set to IM_NI but relaxed toward it,
       FLD(NI,L) = VOLD + OMEGA_RAY*(IM_NI(L)-VOLD)   (l.13192,13255)
       while IM_NI is the M-march's last cell, whose update reads
       DSU=KN_ORD_DS(NI,L) - the endpoint-extrapolated slope built on
       the lagged SPW(NI) preview (l.13035-13064, 12671-12681).

Production telemetry nonetheless shows ray 81 stalling on every
refresh for 43+ refreshes: resid ~0.024-0.10, ixit=3, omega damped to
the 0.45 floor, worst cell persistently i=4 l=17881 (Lyman-series
region, a 5 km/s fine bin adjacent to a ~6x spacing step, ~6 bins
blueward of the Lyman edge at l~17887).

The phase226 production deck replay reproduces the stall class
(prod 3.14E-2 ixit=3 vs replay 3.28E-2 ixit=4 under a pass cap;
ray 82: 5.58E-2 vs 6.22E-2; ray 84 control converges in both). The
captured operator itself fails to contract - the stall is intrinsic,
not a stale-residual, refresh or exit-detector artifact.

---------------------------------------------------------------
3. DIRECT PER-CELL EVIDENCE: AN EXACT PERIOD-2 LIMIT CYCLE
---------------------------------------------------------------

Parsing KN_REPLAY_HST_L081.bin (31 pass snapshots of the P field,
5 x 56377 doubles each) and comparing pass k vs k-1:

  * The residual falls 2.1E-1 -> ~1.2E-3 by pass ~16, then pins.
  * The pin is an EXACT two-state orbit: cell (I=5, l=55273)
    alternates bit-identically between 5.4015557019E-12 and
    5.3952569463E-12 across all remaining passes. This is not noise
    and not slow contraction - a strict period-2 attractor.
  * 55,416 of 56,377 cells on row I=5 and ~50,000 cells on each
    marched row satisfy an exact 2-cycle test over the last 8 passes;
    the whole field rides a low-amplitude orbit (row I=4 is bit-static,
    row I=2 wobbles up to 6E-12, row I=5 up to 1.7E-14 absolute).
  * Only ~2 cells exceed res 1E-3 at the end: the pin is set by a
    handful of cells whose field (~5.4E-12) sits just above the
    significance gate 1E-4*VMAXF (VMAXF=3.6E-8 at l=55273), so their
    residual measures relative wobble, not absolute error.
  * The max-amplitude cycling cells cluster in the far-IR
    (l 52852-52878 at nu 8.7E-3, l 55252-55277 at nu 4.5E-3,
    1e15-Hz units) - the same neighborhood where KN_NO_GAIN reports
    its deepest floored cells (max chi_scat/chi ~11.5 at nu 4.84E-3,
    OUTGEN). These are the highest-albedo cells on the ray, where the
    local map is nearest to non-contracting.

Cross-ray:

  * Ray 82 replay: same orbit class, pinned res 7.7E-3 at I=1 l=24845.
  * Ray 84 (clean control): the SAME 2-cycle exists but its amplitude
    sits at field ~1E-283 - far under every floor - so res ~1E-270
    and the ray certifies IXIT=1.

So the orbit is a generic property of the map; what makes a ray
"worst" is where the orbit's amplitude lands relative to the
significance floor, not whether the orbit exists. Ray 81 is worst in
this model state because its orbit carries the largest relative
amplitude on cells that are still counted as significant.

Omega independence: the phase227 omega sweep (omsweep logs) shows
ray 82 stalling at res ~6E-2, ixit=3, at both omega=0.45 and
omega=1.00. A linear -1 eigenmode would be damped by under-relaxation;
an omega-persistent exact orbit is the signature of a nonlinear
(border-collision) limit cycle on a piecewise-smooth map.

---------------------------------------------------------------
4. MECHANISM
---------------------------------------------------------------

The ordered map is linear in the field except for two discontinuous
nonlinearities that can feed back with negative gain:

  * the sign/min parabolic-slope clamps in KN_ORD_DS (l.12669-12692):
    (sign(SA)+sign(SB or DD)) * min(...) toggles the slope between
    +2*min, 0 and -2*min when a secant or extrapolated DD crosses a
    sign boundary;
  * donor-index switching (JPR1/JPR2, WINX) on compressing cells -
    not engaged on ray 81.

The lagged reads that close the loop on an all-expanding ray are the
C-corner SPW preview (built once per row from the pre-march state,
l.13040/13260) and, at the inner end, the M(NI) -> IM_NI -> relaxed
FLD_P(NI) -> next-pass SPW(NI) -> endpoint DS(NI) -> M(NI) circuit
(l.13174-13260). At high-albedo cells (the KN_NO_GAIN floored far-IR
region and sharp source-structure regions like the Lyman edge) the
linear part of the map is near non-contracting, so the clamp's
discontinuous kick is not damped out: the iteration settles onto a
two-state orbit straddling a clamp sign boundary, machine-exact.

Why the global machinery cannot see it: the NOSC oscillation detector
integrates KN_SGN_SUM over all ~282K cells (l.13599-13605) - the
orbit's net signed sum is dominated by static rows and does not
alternate sign cleanly; RATIO_OM~1 is real, so the stall classifier
correctly exits IXIT=3 and the adaptive controller drives omega to
the 0.45 floor - which cannot break a nonlinear orbit. The per-cell
P2 flip tracker (DO_P2, l.13121-13140, 13222-13241) is exactly the
right instrument for this class (the cycling cells flip every pass),
but KN_ORD_P2DAMP is disabled pending evidence; this audit supplies
that evidence.

Residual-metric note: RES counts every marched cell including the
relaxed boundary row (l.13077-13084, 13204-13211) and normalises by
max(VOLD, 1E-4*VMAXF). A limit cycle of amplitude eps on a cell at
field F reads eps/F - the pin measures the orbit's relative
amplitude, not a distance to the true fixed point. The converged
rows of ray 81 replay are bit-identical pass to pass; the formal
solution has converged, the iteration's tail has not.

---------------------------------------------------------------
5. HYPOTHESIS RANKING (A-F) FOR RAY 81
---------------------------------------------------------------

A. Genuine multiple-resonance physical complexity - REFUTED.
   All-expanding ray, zero compressing cells, min T2=0.85.
B. Very small projected gradients - CONTRIBUTORY CONTEXT only.
   Small |GAM| makes chi_hat ~ |GAM|/DLN dominate chi and stiffens
   the advection coupling, but synthetic rays with the same GAM
   converge; not the cause.
C. Incorrect frequency-direction handling at sign changes - REFUTED.
   No sign change exists on this ray; the direction logic is
   mathematically verified.
D. Boundary/interpolation errors - PARTIALLY CONFIRMED as the locus.
   The pinned cells live on/adjacent to the tangent-boundary feedback
   (relaxed IM_NI copy, endpoint DS extrapolation on a lagged
   preview). Whether any of those steps is a defect per se is open;
   the orbit is a property of the assembled map, not an obviously
   wrong line.
E. Poor convergence of an otherwise correct formal solution -
   CONFIRMED, and made specific: a nonlinear period-2 limit cycle.
   The replay proves the map itself cycles; the committed field is
   essentially the fixed point (seed->end drift ~4E-4 median).
F. Radiation-field storage/refresh problems - REFUTED as the primary
   cause. The production deck replays the stall offline; seeds and
   store precision only modulate the orbit's siting.

---------------------------------------------------------------
6. CONFIRMED / SUSPECTED / NON-BUGS
---------------------------------------------------------------

Confirmed findings (evidence-backed):

  F1. The ordered KNOP iteration can reach a machine-exact period-2
      limit cycle on all-expanding rays (ray 81 replay, 54K+ cycling
      cells, omega-independent; ray 82 same class; ray 84 same orbit
      at the underflow floor). The residual stall is a property of
      the iteration, not of the formal solution.
  F2. The certified ray residual conflates orbit amplitude with
      distance-to-solution: converged rows are bit-static while the
      pin sits on barely-significant cells (ray 81 replay) or
      underflow cells (ray 84). Same map, different certified class.
  F3. The global NOSC detector cannot observe a localized period-2
      orbit; the stall classifier is working as designed on the
      aggregate signal.
  F4. Ray 81's zero-compressing-cell geometry makes it the cleanest
      demonstration that the cycle needs no velocity sign changes.

Suspected (structurally localised, not yet isolated by a
per-cell-instrumented replay):

  S1. The orbit is sustained by the KN_ORD_DS sign/min clamp acting
      on lagged C-corner previews and the endpoint-extrapolated
      DS(NI) inside the M(NI) -> IM_NI -> P(NI) feedback. Prime
      candidate; the linearised gain is too weak to explain the
      orbit, so the cycle must live on the clamp discontinuity
      (border-collision type) - consistent with machine-exact
      two-state values.
  S2. Relaxing the boundary row FLD_P(NI) toward IM_NI instead of
      assigning it adds one avoidable lag stage to that feedback and
      injects IM_NI motion into the residual. Flagged for review,
      not asserted as a defect.
  S3. Hot cells at sharp source structures (Lyman edge, far-IR
      high-albedo floor region) select where the orbit's measurable
      amplitude lands; the siting wanders with state (production
      i=4/l=17881 vs replay far-IR rows).

Non-bugs (checked and exonerated):

  N1. The GAM formula and its sign conventions are correct.
  N2. Frequency upwind selection (expanding->L-1, compressing->donor
      slot) is correct.
  N3. Seeding, store precision and refresh bookkeeping do not cause
      the stall (replay of the deck reproduces it; seed->end drift
      is small).
  N4. The Lyman-edge grid seam at l~17881 is not a bug; it is where
      this state's orbit happens to be measurable.
  N5. KN_NO_GAIN's opacity floor is not the cause; high albedo is
      where the map is least contracting, so the orbit's measurable
      cells gather there.

---------------------------------------------------------------
7. MINIMAL REGRESSION SUITE (proposed, not implemented)
---------------------------------------------------------------

R1. Replay-history cell test: reparse KN_REPLAY_HST_L08*.bin and
    assert the pass-15..end field is an exact 2-cycle on the pinned
    rows and bit-static on converged rows. Locks the phenomenon to
    data; zero new instrumentation needed.
R2. Manufactured stencil: extend kn_ccorner_test.py with an
    all-expanding NI=5 geometry plus a discontinuous source-slope
    state; assert a 2-cycle exists and that the per-cell P2 flip
    counter reaches FG_ORD_P2N - proving the existing DO_P2 tracker
    detects this class when enabled.
R3. Omega sweep assertion: the manufactured orbit's amplitude should
    be omega-dependent but its existence omega-independent (matches
    omsweep ray-82 evidence).
R4. Boundary-row check: assert FLD_P(NI) tracks IM_NI with exactly
    the relaxed lag (documented behaviour) and that its residual is
    reported separately from marched-row residuals in any proposed
    metric change.
R5. Production telemetry gate: after any orbit-breaking change, the
    ray-81/82/83 refreshes must show the same exit classes or better
    with unchanged committed fields (bit-compare against deck
    ends).

---------------------------------------------------------------
8. PRODUCTION-READINESS RECOMMENDATION
---------------------------------------------------------------

The formal-solution correctness is NOT in doubt for ray 81: committed
fields are stable and reproducible, and the converged rows are
bit-identical across passes. What is not production-ready is the
iteration's certification: rays certify clean or stalled depending on
where a generic limit cycle's measurable amplitude happens to sit,
which is a numerical-iteration defect class, not a physics defect.

The evidence now supports, as the minimal next correction to
evaluate (describe-and-test-first, per standing policy):

  C1. Enable KN_ORD_P2DAMP on a diagnostic arm only. The tracker
      already counts per-cell flips and commits the geometric
      midpoint inside the flip bracket; for a true period-2 orbit
      that midpoint is the fixed point. The replay histories are the
      missing evidence that justified keeping it off. Success
      criterion: the L081/L082 replay decks reach their committed
      end fields with res < KN_RAY_TOL, field bit-identical to the
      undamped orbit midpoint to within orbit amplitude.
  C2. Independently, report the boundary-row residual separately
      (S2) in diagnostics before any decision to change the residual
      or the boundary commit semantics.
  C3. Do NOT relax KN_RAY_TOL, the stall exits, or the significance
      floor: the orbit is real and the residual honestly reports
      non-convergence; the defect to remove is the orbit itself.

The velocity-geometry question that motivated the audit is answered:
projected-gradient complexity is not the ray-81 problem. For the
compressing rays (e.g. ray 83's documented 0.69 two-state cycle) the
same orbit class is already documented; whether those rays share the
same micro-channel or add donor-switching cycles remains open.

Status: VERIFIED for the limit-cycle claim and the refuted
hypotheses; PARTIAL on the exact micro-channel inside the lagged
corner/boundary feedback (needs a per-cell VNEW/DS instrumented
replay or a v2 deck carrying the M-field history).
