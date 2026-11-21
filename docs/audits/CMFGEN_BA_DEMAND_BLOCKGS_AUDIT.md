# Phase 2.29 - population demand audit + block-GS assessment

Date: 21-Nov-2026. Companion to CMFGEN_PHASE228_LIMIT_CYCLE_AUDIT.md
(which owns the ray-side findings). Question set: (a) which
species/levels/depths dominate the remaining BA_EX demand the strict
relative gate still reports at ~1e17-1e67 percent, and are they
physically meaningful populations or quasi-dead rows; (b) do the
observed hot-frequency runs justify a frequency-block Gauss-Seidel
strategy (assessment only - no implementation authorized).

Evidence: live OUTGEN of m15_2p82d_VAL_knop_formal (job 13473182,
warm restart, iterations 684-707 observed), the RUNAWAY/MODCTL/BA-CTX
diagnostic blocks, and the phase228 limit-cycle audit for the
ray-side context. Configuration under audit: BA_TRUST=1e4,
BA_DEM_PFLOOR=1e-30, BA_EX_RELFRAC=1e-30, BA_SAT_RAW=2e3,
EPS_TERM=5e-2. No solver change was made or proposed in-place.

---------------------------------------------------------------
1. WHO CARRIES THE DEMAND
---------------------------------------------------------------

Non-exempted binding demand (the "largest exit-exempted demand"
line and MAXCH) is carried by two quasi-dead families:

  * C6+ ground state: `CSIX:lev 1`, demand 1.24E+14 .. 2.60E+17 %
    wandering across depths 1-17 over its 699-705 (siting: d16 -> d14
    -> d3 -> d1 -> d4 -> d4 -> d17). eq index 246.
  * C5+ excited levels: `CV:1s_3d_1De` (eq 211), `CV:1s_2s_3Se`
    (eq 204), `CV:1s_2p_3Po` (eq 206), `CV:1s_3p_1Po` (eq 212),
    `CV:1s_2p_1Po` (eq 207), `CV:1s_7s_1Se` (eq 239) at depths 4-22.

RUNAWAY top-10 (it 706) gives the populations directly:

    d=10 CV:1s_3d_1De  n=3.6E-91  ionf=2.4E-93  dem%=-2.97E+67
    d=7  CV:1s_3d_1De  n=2.8E-96  ionf=3.5E-98  dem%=-8.31E+58
    d=9  CV:1s_2p_3Po  n=1.1E-75  ionf=9.6E-78  dem%=-2.70E+49
    d=11 CV:1s_3p_1Po  n=1.7E-72  ionf=5.8E-75  dem%=-8.66E+48
    ...all top-10 rows are C V excited states at depths 4-22,
    populations 1E-62 .. 1E-96, ion fractions ~1E-93.

The applied correction sits at the trust cap every iteration
(applied 9.9E+03 = BA_TRUST 1e4 -1% margin, its 684-706).

Interpretation: the demanded corrections are -1e67..-1e13 percent on
levels whose populations are 1e-62..1e-96 and ion fractions ~1e-93.
At these abundances the rows carry no measurable opacity or
emissivity; the demanded "correction" is the SE solver asking a
numerically dead row to change by -100%. The siting wanders
(depth 1 -> 17 between consecutive iterations) - the demand is a
wandering runaway tail, not a stable physical.

Verdict for the audit question: the strict relative-demand gate is
dominated by QUASI-DEAD ROWS - C V/C VI levels at outer depths with
populations 1e-62..1e-96 - not by physically meaningful populations.
The deepest binding contribution is cascade-refilling into nearly
unpopulated excited levels of high ions that are poorly conditioned
SE rows (negative ~100% demands alternating with positive 1e4 caps).

Separately, the worst SE equation is `var 27 = H2:lev 1` at depth 65
(max SE residual 2.13E+13, excl-worst-columns 1.02E+13) - a real,
populated equation; the residual gate and the demand gate are
driven by different mechanisms and should not be conflated.

---------------------------------------------------------------
2. GATE ACCOUNTING (ties to the Patch-B classifier)
---------------------------------------------------------------

Reporter classification per the new CONV_MARGIN_EVAL semantics:

  * CSIX:lev 1 demand 1e13-1e17 % >= BA_SAT_RAW(2e3): SAT-released.
    This is the demand that printed as `worst margin 78.85 dex (exm)`
    on the old reporter; it never imposed a hold.
  * CV excited-level demands likewise saturate vs BA_SAT_RAW or are
    floor-exempted (populations << BA_DEM_PFLOOR-adjacent scales);
    the RUNAWAY list is exactly the exempted population.
  * The residual active binding channel is ex_maxch (non-exempted
    MAXCH): 1.24E+14 at it706 - still astronomical because the same
    quasi-dead family feeds it; it binds the exit honestly as long
    as the demand gate is relative.

Policy unchanged: BA_EX_MAXCH and BA_DEM_EXIT_OK semantics verified
in phase 1 are untouched; this audit only names WHO dominates, not
WHICH gate should move. Any gate change is a separate proposal.

---------------------------------------------------------------
3. BLOCK-GS ASSESSMENT (Phase 6 - proposal verdict)
---------------------------------------------------------------

Question: do the ~27,582 hot frequencies / ~495 contiguous runs
justify a frequency-block Gauss-Seidel (block-GS) strategy?

Evidence considered: phase228 audit sections 3-4. The stall class is
a machine-exact period-2 border-collision limit cycle sustained by
the lagged C-corner/boundary feedback (lagged SPW preview into the
sign/min DS clamp, endpoint DS(NI), relaxed IM_NI copy). The cycling
cells are isolated at the significance floor; converged neighbouring
rows are BIT-STATIC pass to pass.

Why block-GS does not match the defect:

  * Block-GS accelerates convergence when the residual is carried by
    LAGGED CROSS-FREQUENCY COUPLING - neighbouring frequencies
    solving against stale values within a pass. The observed orbit
    is not a cross-frequency lag defect: the same two-state
    alternation exists per cell under an already-ordered sweep, and
    converged rows are bit-identical - the frequency coupling at
    those cells HAS converged; only the clamp straddles a sign
    boundary and alternates.
  * Ordering the hot runs into blocks would change nothing about the
    DS sign/min clamp discontinuity that sustains the orbit
    (border-collision is a local nonlinearity, not an ordering
    defect).
  * The 495 contiguous runs correlate with the high-albedo floored
    region and sharp source structures (Lyman edge); they are where
    the orbit's amplitude is MEASURABLE, not a coherent lagged
    coupling block that would relax faster as a unit.
  * Omega-independence (ray 82 stalls identically at omega=0.45 and
    1.00) rules out any relaxation-schedule remedy of this class.

Verdict: block-GS is NOT justified by current evidence and no
proposal is drafted. The evidence-supported candidate correction
remains phase228-C1 (KN_ORD_P2DAMP on a diagnostic arm - the
existing per-cell flip tracker that already commits the geometric
midpoint inside the flip bracket). A future trigger that WOULD
justify revisiting block-GS: residual analysis showing neighbouring
frequencies lagging each other (a cell's residual dominated by
cross-frequency lag terms rather than the clamp discontinuity), or
hot runs whose cells all sit inside one frequency-coupled stencil
dependency chain - measurable with the existing per-cell residual
breakdown instrumentation, no block-GS needed to measure it.

---------------------------------------------------------------
4. CONFIRMED / OPEN
---------------------------------------------------------------

  * CONFIRMED: demand carriers are C V/C VI quasi-dead rows at
    depths 1-22 (n ~1e-62..1e-96, ionf ~1e-93); applied correction
    trust-capped for 20+ consecutive iterations; demand siting
    wanders (runaway tail, not stable structure).
  * CONFIRMED: the old reporter's `worst margin (exm)` was the
    SAT-released CSIX row; Patch B's active/released split matches
    the release predicates observed in OUTGEN.
  * CONFIRMED: worst SE eq is H2 lev1 at depth 65 - a separate,
    real mechanism from the demand tail.
  * VERDICT: block-GS not justified; orbit is local
    border-collision, not cross-frequency lag.
  * OPEN (not authorized): whether quasi-dead rows merit a
    population-floor or demand-gate policy change - any such change
    is a separate proposal requiring its own evidence (opacity/
    emissivity contribution census of the demanding rows).

Status: VERIFIED for the demand-carrier identification and the
block-GS rejection; PARTIAL on the quasi-dead-row remedy (no
remedy proposed or authorized here).
