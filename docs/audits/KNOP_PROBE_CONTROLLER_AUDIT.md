# KNOP probe/controller audit

Date: 20-Nov-2026. Companion to `KNOP_THREE_SEED_EXPERIMENT.md`.
Scope: `plane/fg_j_cmf_nonmono_v1.f` (probe ray, refresh, omega pool,
seed arrays) and `new_main/mod_subs/solve_for_pops.f` (maturity gate,
freeze/release state machine). Read-only audit - no code changed, no
running job touched. All line numbers are against the current tree
(commit d369771e + this file).

Contents:

1. Probe data-flow map and inclusion/exclusion table
2. GS_MAX_ALL exact expression and maturity conditions
3. Adaptive omega (KN_OMEGA_POOL_V1) sample audit
4. Freeze/release state machine
5. Seed binary array validity
6. Defect list (classified)
7. Recommended corrections (not implemented)
8. Proposed MODE1/MODE2 follow-up experiment
9. Test cases required before any implementation

----------------------------------------------------------------------
## 1. Probe data flow

Ray loop (per refresh): `KNOP_RAY_SOLVE_V1(LS, RES_LS(LS), NST_LS(LS),
..., NDV_LS(LS), ..., NXT_LS(LS), ..., RHO_M_LS(LS), RHO_P_LS(LS),
NCL_LS(LS), NNL_LS(LS))` called at :8566 for every `LS=1..KN_NP`,
including the probe. Inside the solve `L_NOPROBE =
KN_NOCP_PROBE_RAY(LS)` (:9796) suppresses store writes only; every
output argument is still written with the probe's own solve results.

Flow map:

    probe ray 83 solve (same code path as production rays)
      -> RES_LS(83)    :8813-8816  MAXVAL mask NDV==0 -> KN_RF_RESID
      -> KN_RF_RESID   :8884       -> FG_GS_MAXCHG
      -> FG_GS_MAXCHG  getter      -> GS_MAXCHG (solve_for_pops)
      -> GS_MAX_ALL    :1248-1250  MAX(GS_MAXCHG, GS_SW_MAXCHG,
                                       GS_RF_DCHG)
      -> GS_FG_IMMATURE:1285-1314  -> GS_RAW_IMMATURE :1586
      -> GS_FREEZE     :1617       -> N_FROZEN, hold/release, counters

    NBD_LS/NTT_LS/NST_LS/NCL/NNL/NDV (probe rows included)
      :8885,8897,8900-8904 -> FG_GS_NBAD/NTOT/NCONV, FG_RF_CAPPED,
      FG_GS_NDIV -> share tiers :1272-1283, :1331-1333, :1378-1381

    RHO_P_LS(83) (RHO_M_LS(83)=0, npM=2 < NPO>=5 sample floor :10869)
      :8851 -> KN_OMEGA_POOL_V1 -> INEXT_RELAX :2083-2085
      -> I_P_DEPTH_NEXT/I_M_DEPTH_NEXT update :5837,:5891
      -> every later sweep's retained field (all rays, all frequencies)

    WRITE_STORE / WRITE_STORE_CAP :10214-10218, :10237-10239
      -> SUPPRESSED under L_NOPROBE -> DCH_M=DCH_P=0 (:9793)
      -> FG_RF_DCHG=MAX(...,0,0) no-op :10223,:10249
      -> FG_RAY_VC overwritten to 7 :10264-10268 (NOP_VC keeps the
         would-be verdict)
    HEAL_STORE :10410-10412 suppressed for a diverged probe
      (probe then reports RESID=1.0 like any diverged ray)

Inclusion/exclusion table (probe = ray KN_NOCP_PROBE_LS=83, NDV=0,
always solves):

| Quantity | Probe contributes | L_NOPROBE suppresses | Diagnostic or control | Downstream consumers | Alters later iterations |
|----------|-------------------|----------------------|-----------------------|----------------------|-------------------------|
| RES_LS(83) | yes (solve residual) | no | control input | KN_RF_RESID :8814, MAXLOC :9407, census :8946 | yes |
| KN_RF_RESID / FG_GS_MAXCHG / GS_MAXCHG | yes - argmax every refresh | no | control | GS_MAX_ALL :1248 | yes |
| GS_NBAD / GS_NTOT (NBD_LS,NTT_LS) | yes | no | control | bad-share tier :1310, :1281 | yes |
| GS_NBAD_BLK | yes | no | control | MIN() at :1252 | yes |
| GS_NCONV (NST_LS) | yes | no | control | churn tier :1333 | yes |
| FG_RF_CAPPED (NCL/NNL) | yes | no | control | GS_CAP_NUM :1273 -> UNSAFE/PIN | yes |
| GS_NDIV (NDV_LS) | yes (0; would count) | no | control | diverged share :1280 | yes |
| NX_*/exit census (NXT_LS) | yes (ixit=3 row) | no | diagnostic | print :8944+ | no |
| FG_RAY_VC | yes, value 7 | partial (would-be verdict kept in NOP_VC) | diagnostic | verdict census :9452+, pending mask :9450 | no |
| FG_RAY_DCH | yes, value 0 | yes (DCH never computed) | diagnostic | excluded by VC masks anyway | no |
| DCH_M/DCH_P -> FG_RF_DCHG | NO contribution (stays 0.0) | yes - the write skip makes it 0 | control | GS_MAX_ALL :1250 | no |
| FG_RF_PEND | no (VC=7 fails mask :9450) | yes | diagnostic | print :8934 | no |
| FG_COMMIT_REJ | yes if probe JUMP>veto (:10143 ungated) | NO | diagnostic counter | print :6832, FGINEXT_GS :6878 | no |
| FG_NC_CAP / FG_NC_REJ | yes if branch taken (:10220,:10230 ungated) | NO | diagnostic counters | prints :8926 | no |
| FG_RAY_SEEDED | no (:10246 gated) | yes | bookkeeping | jump-veto waiver | n/a |
| RHO_P_LS(83) | yes (one pool sample) | NO | control input | KN_OMEGA_POOL_V1 :8851 | yes |
| INEXT_RELAX | yes (indirect, via pool) | no | control | store update :5837,:5891 for all rays | yes |
| KN_RES_PASS envelope | yes (dominated by probe) | no | diagnostic only under OMEGA_V2 | telemetry | no |
| KN_RF_RESID_PRV | yes (inherits) | no | legacy fallback path (:8867, unused under OMEGA_V2) | old omega | marginal |
| I_P_DEPTH_NEXT / I_M_DEPTH_NEXT (store) | NO - writes gated | yes | the protected state | all later sweeps/seeds | no |
| LEDG_NP_RF ledger | yes (npP counted :8912) | no | diagnostic | ITER_LEDGER | no |
| KN-NOP telemetry files | yes | n/a | diagnostic | OUTGEN | no |

Measured corroboration (live experiment): probe store sums identical
before/after every refresh in all arms (delta M/P = 0.000E+00);
FG_RF_DCHG differences between arms are therefore NOT probe-caused -
they ride through the omega/pool drift (section 3) and OMP
nondeterminism, not through the delta accumulator.

## 2. GS_MAX_ALL and the maturity gate (solve_for_pops.f)

Exact expression:

    GS_MAX_ALL = GS_MAXCHG                                    (:1248)
    IF(GS_NTOT .GT. 0)
        GS_MAX_ALL = MAX(GS_MAX_ALL, GS_SW_MAXCHG, GS_RF_DCHG) (:1250)

with, on the FG side:

    GS_MAXCHG    = FG_GS_MAXCHG = KN_RF_RESID
                 = MAXVAL(RES_LS, MASK = NDV_LS .EQ. 0)      (:8814,:8884)
                 = post-refresh per-ray residual max - INCLUDES probe
    GS_SW_MAXCHG = FG_SW_MAXCHG = in-sweep store drift snapshot
                 taken before the refresh overwrites FG_GS_* (:8882).
                 Includes ray 83's slot updates - the in-sweep lagged
                 write at :5837/:5891 is production state (only the
                 probe's refresh commit is suppressed), so this term
                 legitimately contains the probe ray's share.
    GS_RF_DCHG   = FG_RF_DCHG = MAX over commit writes of
                 DCH_M/DCH_P (:10223,:10249) - probe contributes 0.0

Gate decisions (GS_NTOT > 0 branch, i.e. refresh ran):

    GS_NBAD_EFF = MIN(GS_NBAD, GS_NBAD_BLK)                  (:1251-1255)
    GS_NCONV_ALL = MAX(GS_NCONV, GS_SW_NCONV)                (:1256)
    GS_CAP_NUM = GS_RF_CAPPED, GS_CAP_DEN = GS_NTOT          (:1272-1274)

    GS_UNSAFE = (CAP_NUM > MAX(50, CAP_DEN/10))              (:1279-1281)
             OR (GS_NDIV > MAX(1, GS_NRAY/100))
             OR (GS_NBAD > MAX(50, GS_NTOT/10))
    GS_PIN = GS_UNSAFE OR (CAP_NUM > MAX(50, CAP_DEN/1000))  (:1282-1283)
    GS_UNSAFE = GS_UNSAFE OR
             (GS_NBAD_PLAT > MAX(50, GS_NTOT/100))           (:1378-1381)

    GS_FG_IMMATURE = GS_MAX_ALL > GS_MATURE_TOL              (:1285,
                     GS_MATURE_TOL = 0.05 at :213)
    IF GS_FG_IMMATURE: recompute as                          (:1310-1312)
        GS_FG_IMMATURE = (GS_NBAD_EFF > MAX(50, GS_NTOT/1000))
                      OR (GS_MAX_ALL > GS_UNSAFE_TOL)   (0.5, :221)
    churn OR (unless L_SETTLE):                              (:1331-1333)
        GS_FG_IMMATURE |= (GS_NCONV_ALL > MAX(50, GS_NTOT/1000))

    GS_RAW_IMMATURE = GS_PIN OR GS_FG_IMMATURE
                   OR (MOM_RAN AND MOM_FREEZE_B)             (:1586-1588)

So the max-norm veto inside the gate is `GS_MAX_ALL > 0.5`
(GS_UNSAFE_TOL); the 0.05 level only decides whether the refinement
runs at all. A single ray can hold the gate immature through term 1
of GS_MAX_ALL whenever its residual exceeds 0.5 - no share test
applies to that term. RF_PEND is NOT in GS_MAX_ALL; it is a print
only. The current/post-refresh residual (GS_MAXCHG) is used, not
the sweep value, with the sweep value kept as a separate term.

Answers to the posed questions:

- Every term: refresh residual max (probe-inclusive), in-sweep store
  drift, retained-field delta (probe-free by accident of DCH=0).
- GS_UNSAFE_TOL = 0.5 (parameter, :221). GS_MATURE_TOL = 0.05 (:213).
- Maturity requires: not GS_PIN, GS_NBAD_EFF <= 0.1% census,
  GS_MAX_ALL <= 0.5, GS_NCONV_ALL <= 0.1% census (plus settle/plateau
  variants). Immaturity = any of those failing.
- The gate uses the POST-REFRESH residual (probe's current solve).
- RF_PEND is not read by the gate at all.
- Bad-cell fraction: >0.1% bad cells immatures; >10% bad cells or
  >1% plateau-class bad cells is UNSAFE. Diverged rays >1% of solved
  is UNSAFE.
- Can a single probe ray prevent maturity for the whole model? YES -
  through the max-norm term alone (resP > 0.5). Observed: A and A/A
  pinned at 0.69614 > 0.5 on every refresh, B at 2.2799, so the gate
  has been immature on every call of those arms.

## 3. Adaptive omega audit (KN_OMEGA_POOL_V1, :2041-2089)

Pool construction: for every ray with NDV==0, append RHOM(I) if >0,
RHOP(I) if >0 (:2059-2067). Sort ascending; RHO_MX = last,
RHO_MED = middle, RHO_P90 = element MIN(NSAM, (9*NSAM+9)/10) clamped
to [0.05,20]; needs NSAM >= KN_RHO_MIN_SAM=10 (:565). Update:
`omega_new = 0.5*omega_prev + 0.5*clip(1/(1+rho_p90), 0.05, 1)`
(:2083-2085). Consumed at :5837/:5891 where the in-sweep lagged store
update `STORE = old + INEXT_RELAX*(new - old)` runs for every
compressing slot of every ray on every frequency call.

Per-ray sample provenance: RHO_M/RHO_P are the median of a 6-element
ring of same-direction pass ratios RATIO_OM = RES/RES_PREV(idir),
sampled only when NPO >= 5 and not hot-mode (:10857-10873,
:11060-11081, :12040+). The probe ray runs the identical code;
L_NOPROBE does not touch RHO_WIN, so the probe contributes its own
ratio like any other ray.

A. Probe -> KN_RF_RESID: yes, it IS the max (:8814).
B. Probe -> pool samples: yes, exactly ONE sample - RHO_P.
   RHO_M = 0 because npM = 2 < 5 never fills RHO_WIN. The probe's
   'rat' column in the KN-NOP pass report is RATIO_OM, so RHO_P_LS(83)
   can be reconstructed offline as the median of the last six
   P-direction rats. Reconstructed values (last refresh blocks):
   A ~0.83, B ~0.875, C ~0.98-1.03.
C. Probe -> INEXT_RELAX: indirect and arm-dependent.

   - rho_max in the print is the pool MAX, NOT the probe residual:
     A 0.9895 > probe 0.83; the earlier interim report's claim that
     rho_max "is the probe residual" was wrong. rho_max is some other
     ray's sample; the probe does not set it on A or B.
   - The applied omega reads only RHO_P90. With pool sizes ~58-59 the
     p90 index is (9*NSAM+9)/10 ~ 53-54; removing a sample that sits
     BELOW the p90 position shifts every element above it down by one
     while the index also drops by one, so the p90 element is
     IDENTICAL - a below-p90 probe sample contributes exactly nothing
     to the applied omega (A: probe ~0.83 < p90 ~0.947, B: ~0.875 <
     ~0.95 - inert).
   - Arm C's probe rho ~1.0 sits ABOVE the printed p90 ~0.95 - it
     directly raises the P90 element, lowering omega by ~one sorted
     position. This is a real controller input: the chained seed's
     stall behaviour feeds the relaxation the whole model uses next
     sweep (C applied omega prints 0.5131-0.5137 vs A 0.5137-0.5140,
     consistent with a one-position upward shift of rho_p90).
   - The measured cross-arm omega differences (mean |d| ~6E-4, max
     ~1.3E-3) are consistent with that one-index shift plus
     pool-level nondeterminism; they are NOT explained by the probe
     being rho_max (it is not, on A/B). Because per-ray pool samples
     are not logged, the exact counterfactual pool without ray 83
     cannot be replayed from OUTGEN alone; the bound above stands
     (<= one sorted position, sub-1% on p90, and only when the probe
     sample lands in the top decile).
D. Probe -> retained store update: indirect via INEXT_RELAX when (C)
   holds; zero via FG_RF_DCHG (DCH=0 always for the probe).

## 4. Freeze/release state machine (solve_for_pops.f)

State variables: GS_FREEZE (:1617), N_MATURE (:1595-1599),
N_FROZEN/N_FROZEN_TOT (:1720-1728), GS_JUST_REL/GS_REL_DEM
(:2150-2154), N_REL_NOPROG (:1156-1177), N_PRECOND/N_PRECOND_TOT
(:1643-1651), GS_WAS_FROZEN (:2161), GS_HOLD_POPS/GS_POPS_HELD
(:2186-2217).

Transition logic:

    GS_RAW_IMMATURE = GS_PIN OR GS_FG_IMMATURE OR MOM_FREEZE_B (:1586)
    N_MATURE: 0 on immature call, else +1                    (:1595-1599)
    GS_FREEZE = GS_RAW_IMMATURE OR (N_MATURE < GS_REL_STABLE)(:1617;
                GS_REL_STABLE = 2 in this run)
    frozen call: N_FROZEN++, N_FROZEN_TOT++                (:1721-1726)
    release (GS_FREEZE false after frozen): GS_JUST_REL=T,
        GS_REL_DEM = entry MAXCH                           (:2150-2154)
    next call judges release: no-progress if
        MAXCH >= 0.9*GS_REL_DEM  OR  GS_SW_MAXCHG > 0.5      (:1156-1160)
        -> N_REL_NOPROG++; else N_FROZEN_TOT=N_PRECOND_TOT=0 (:1174-1175)

Stops:

    N_FROZEN_TOT  >= GS_FREEZE_LIM(60)   -> STOP (:2087-2094)
    N_PRECOND_TOT >= GS_PRECOND_LIM(50)  -> STOP (:2105-2113)
    N_REL_NOPROG  >= GS_REL_NOPROG_LIM(10)-> STOP (:2122-2130)

Critical counter subtlety: on EVERY call where the precondition gate
is not tripped and MAXCH <= GS_PRECOND_MAX (= 1.0E+03, default
rd_control_variables.f:1154), both N_FROZEN_TOT and N_PRECOND_TOT
are cleared (:1685-1688). This run sits at demand ~82 << 1000, so
N_FROZEN_TOT is reset to 0 on every single call - GS_FREEZE_LIM is
UNREACHABLE in the low-demand held regime. The only live stop is
N_REL_NOPROG, and it can only arm if the gate releases.

What applied=0 and conv=50% mean: under GS_HOLD_POPS=T (VADAT), a
frozen call still runs SOLVEBA, prints `GS_HOLD_POPS: populations
held; attempted max change=...`, then restores POPS (:2216+) - the
applied correction MAXCH_APP is 0 by construction. CONV_NPCT =
NINT(100*CONV_NPASS/8) (:4496): 4 of 8 convergence gates pass, the
hold keeps the rest unmet. This is state B of the posed options: an
intended temporary hold pending field maturity (Knop Paper-II nested
iteration), here held permanently because the gate never matures.

Why arm C died and the others cannot:

- A/B/A-A: probe resP pinned above 0.5 -> GS_MAX_ALL > 0.5 every
  call -> GS_RAW_IMMATURE always true -> N_MATURE never reaches 2 ->
  GS_FREEZE never clears -> no releases -> N_REL_NOPROG stuck at 0;
  N_FROZEN_TOT cleared every call by the precondition-deny branch ->
  GS_FREEZE_LIM also unreachable. Absorbing state: the only exits
  are field maturity (probe would need resP < 0.5 twice running -
  it cannot, its seed is fixed) or the 24h wall/requeue cycle.
- C: chained seed produced resP < 0.5 at 53 of 292 refreshes. Run
  telemetry shows exactly TEN sub-0.5 streaks of length >= 2 - refs
  46-47, 152-153, 182-184, 202-203, 207-209, 239-240, 244-245,
  251-252, 258-259, 290-292 - one per judged no-progress release,
  matching GS_REL_STABLE=2 1:1: each 2-mature streak released, the
  release applied a real SOLVEBA step (hold lifts with GS_FREEZE
  false), demand stayed ~82 (held pops cannot have moved the field
  while frozen), so MAXCH >= 0.9*GS_REL_DEM judged no-progress on
  the next call and the re-poisoned gate refroze. The tenth streak
  (refs 290-292) reached the limit -> STOP at it449.

Correlated-or-independent: the counter counts individual
frozen->live transitions, each requiring its own 2-mature streak -
not a per-episode unique key. The episodes are strongly correlated
(same oscillating probe residual opens the same gate, same doomed
release outcome), so N_REL_NOPROG measures "number of gate-openings
the seed produced", not "number of independent non-convergence
trials". The no-progress test itself (demand must fall >10% or the
field must not re-poison) is reasonable where releases are expected
to make progress; it is degenerate in a regime where the held pops
guarantee no progress is possible - every release is judged by a
test it was designed to fail.

## 5. Seed binary array validity

Record writer :14426-14442 (NI=3 x NCF=56377, verified: 2,706,216 B
records, consistent markers, header fields match filename).

Cell classification for FLD_P/KN_DK_XREF/KN_DK_X0 (NI x NCF):

1. Valid transported cells: (JP2,LS2) for JP2=1..3, LS2=2..NCF where
   STORE holds finite physical values. ~88.5K cells.
2. Valid-but-inactive (dead) cells: FLD(:,1) - FLD is zeroed at
   :9995/:10064, SEED_FLD (:10627-10641) writes only LS2>=2, so
   column 1 stays 0.0 -> KN_DK_W=0 -> KN_DK_NDEAD=3. These ARE read
   by MARCH (boundary column) as zero - safe.
3. Dead cells by weight: KN_DK_W=0 wherever the whole XREF column
   is zero - same 3 cells only.
4. Unused workspace cells: none - SEED_FLD covers every (JP2,LS2)
   mapped cell and NI x NCF is exactly the mapped set (169128 live +
   3 column-1 dead = 169131).
5. Uninitialised/non-finite cells INSIDE the transported set:
   ~80.5K cells hold |x|>1E15 (uninitialized-scale) and ~70-95 hold
   NaN. Origin: they come straight out of
   STORE(KC,LS,LS2-1) -> SEED_FLD -> FLD - they are retained-store
   contents, identical for every ray and every arm at refresh 1, not
   a seed-arm artifact. Positions identical between A and A/A
   (identical store trajectory), differing for B/C (their stores
   drifted via omega/nondeterminism).

Do invalid cells reach the formal solve? YES - MARCH (:12967+)
iterates the full NI x NCF workspace; `FLD(I,L)` enters the upwind
recurrence `FLD(I,L)*A0L(I)` directly. Three containment mechanisms
bound the damage, all pre-existing solver behaviour, not seed-arm
code:

- value guards: VNEW < 0 or NaN -> 1.0E-30 (NNG_P++), VNEW > CAP ->
  NONMONO_INEXT_CAP (NCP_P++) (:13000-13009, :13087-13096).
- residual masking: VMAX_L = column max -> a 1E305 cell makes the
  whole column's significance floor ~1E301 -> every physical residual
  in that column is masked out of NTOT/NBAD and the max-norm
  (:12979-12985, :13021-13036).
- NaN cells: |VNEW - NaN| = NaN -> DCELL comparisons are all FALSE ->
  NaN cells are invisible to residual AND census silently.

Does KN_SEED_SIB write cells the solve reads? Yes, all of its
169128 writes land inside the marched workspace - that is its
purpose. Its sample filter (:14516-14518: finite, >0) means every
write is a finite positive median; it REPLACES garbage cells with
live values rather than injecting garbage. Consequence: arm B's
"perturbation" is two different things at once - changed values on
~46K live cells AND resurrection of ~42K dead columns (which then
re-enter the residual census because their VMAX_L becomes physical).
DAB ~8.26E4 is dominated by dead->live transitions, not by the
live-cell distance it was designed to measure.

Does KN_SEED_SETUP_P compute metrics on invalid cells? YES: the
DAB/RMAX/SCHG loops (:14395-14414) iterate all NI x NCF cells with
no finite-value mask; KN_DK_W only weights by magnitude (1/xref),
so dead->live transitions enter with weight ~1/|1E-4*colmax|.
Documented caveat: the reported distances are inflated by
non-physical cells; the live-cell subset needs a finite mask.

Does KN_SEED_SAVE_P store only valid cells? It copies FLD verbatim
(:14468) - the saved candidate includes the same garbage/NaN cells,
which arm C then re-injects next refresh (:14376). Consistent
chaining, including of the dead columns (they stay dead through the
chain).

FLD(:,1) read by MARCH: yes - but it is zeroed, not uninitialised;
SEED_FLD skipping LS2=1 is by design (store index LS2-1 starts at 1).
Not a defect; the column-1 zero behaves as the boundary value.

## 6. Defect list (classified)

CRITICAL - none found that corrupts committed results in production.
In production there is no probe ray; the contaminated channels below
only exist because the diagnostic probe shares the same accounting
paths as real rays. For the EXPERIMENT's validity, however, finding
H1 is effectively critical: the assumption "the probe is an
observer" does not hold.

HIGH
- H1. Probe residual enters the maturity gate unfiltered.
  RES_LS(83) -> KN_RF_RESID (:8814) -> FG_GS_MAXCHG (:8884) ->
  GS_MAX_ALL term 1 (:1248) -> gate :1312 -> GS_FREEZE :1617.
  A non-committing diagnostic ray is the argmax on every refresh of
  every arm and single-handedly keeps the model frozen (A/B/AA) or
  drives release episodes that end in a designed STOP (C). Earliest
  source location: the NDV==0 mask at :8813-8816 - the probe passes
  it because it legitimately converges/diverges like a ray, but it
  is not a ray of the model's field. Conflates "a real ray that has
  not converged" with "a diagnostic instrument that cannot fix the
  field" - exactly the distinction section 2 asks about.
- H2. Freeze state has no reachable exit in the low-demand held
  regime. N_FROZEN_TOT is cleared every call by the precondition-deny
  branch (:1685-1688) whenever MAXCH <= GS_PRECOND_MAX=1E3; at
  demand ~82 the only live stop is N_REL_NOPROG, which requires the
  probe-controlled gate to open. A/B/AA are in an absorbing state:
  they can neither converge populations (held) nor die (no stop path)
  until external budget limits. Design gap, not a code bug per se.
- H3. Probe rho sample enters the omega pool (:8851 via :2060-2066).
  Arm C's ~1.0 sample sits above p90 and measurably shifts the applied
  INEXT_RELAX (cross-arm omega differences ~6E-4 mean); that omega is
  then applied to every ray's store update (:5837,:5891) - a global,
  persistent feedback channel from the seed choice into the model
  field. Sub-1% per refresh but it does not decay; it is why B/C's
  retained seed drifts ~20x above the A/A nondeterminism floor.
- H4. Probe cells enter all share-tier censuses: NBD/NTT -> bad-share
  (:1310), NST -> churn (:1333), NCL/NNL -> cap-share (:1273). Small
  (~1/59 of cells) but nonzero; combined with H1 the gate's entire
  input set is probe-contaminated.

MEDIUM
- M1. Arm B conflates two perturbations: changed live-cell values AND
  resurrection of ~42K dead columns (garbage -> finite). DAB ~8E4 is
  dominated by the dead-column resurrection; it does not measure the
  intended "sibling-median live seed" distance. The resurrected
  columns also change the census denominators.
- M2. KN_SEED metrics computed over non-physical cells
  (:14395-14414, no finite mask). DAB/RMAX/SCHG inflated; the
  reported values should be re-derived offline with the documented
  physical mask (done in the experiment doc).
- M3. FG_COMMIT_REJ (:10143), FG_NC_CAP (:10220), FG_NC_REJ (:10230)
  increment ungated on the probe's would-be path; the counters feed
  prints and the FGINEXT_GS checkpoint state. Harmless to control but
  pollutes run-state diagnostics.
- M4. NaN columns are invisible to the residual AND the bad-cell
  census (NaN comparisons all false). Pre-existing solver condition,
  not probe-caused, but it means "converged" reports can silently
  exclude NaN-corrupted frequency columns. Worth flagging since the
  seed experiment makes those columns live again on arm B.

LOW
- L1. FLD(:,1) never written by SEED_FLD (zeroed at :9995; store
  index LS2-1 -> col 1 has no store). Harmless - it is the boundary -
  but it is why KN_DK_NDEAD=3: the "dead" count is a bookkeeping
  artifact, not pathology.
- L2. RHO_WIN ring is direction-agnostic until the last pass assigns
  its whole median to whichever direction ran last (:11077-11081);
  with M=2 passes RHO_M stays 0 and is skipped - fine here, but the
  ring can attribute mixed samples on rays where the last direction
  ran <6 passes.
- L3. 'rho_max' in the omega-est print is easily misread as the probe
  residual (it was, in the interim experiment doc). Naming/reporting
  clarity only.

## 7. Recommended corrections (NOT implemented)

In order of leverage; all gated on KN_NOCP_PROBE so production is
bit-identical when the diagnostic is off:

1. Exclude the probe from the residual pool: at :8813-8817 compute
   the max over `NDV==0 .AND. .NOT. KN_NOCP_PROBE_RAY(LS)` (needs a
   small mask or an explicit skip of LS==KN_NOCP_PROBE_LS). Keep the
   probe's own RES_LS for telemetry.
2. Exclude the probe's RHO_P from KN_OMEGA_POOL_V1 (mask array or a
   `LS .NE. KN_NOCP_PROBE_LS` guard in the pool loop :2059-2067).
3. Exclude probe cells from the share tiers: subtract its
   NBD/NTT/NST/NCL/NNL before the SUM() calls at :8897-8904, or gate
   the accumulations on .NOT. probe.
4. Gate FG_COMMIT_REJ/FG_NC_CAP/FG_NC_REJ on .NOT. L_NOPROBE (their
   would-be counts are already in NOP_VC telemetry).
5. Do NOT weaken GS_REL_NOPROG_LIM - it is doing its job. The fix is
   upstream: either the probe stops feeding the gate (corrections
   1-4) or the hold regime gets its own watchdog for "no releases
   possible" (e.g. count consecutive held calls without a release
   and treat N_REL_NOPROG as satisfied only after at least one
   release; that is a design decision, not a quick patch).
6. KN_SEED metrics: add the finite/physical mask to the DAB/RMAX/SCHG
   loops, or report both variants.
7. If a future arm wants the probe to also test commit behaviour,
   MODE2 below; do not relax the store gates in the same experiment.

## 8. Proposed follow-up experiment (two modes, do not mix)

Same checkpoint, executable, grids, controls, and A/A repeatability
arm in BOTH modes; only the probe's accounting changes.

MODE 1 - PROBE-ISOLATED:
- Apply corrections 1-4: the probe contributes its own telemetry
  (KN-NOP passes, KNSEED_* binaries, seed metrics) and NOTHING to
  RES pools, share censuses, omega pool, or commit counters.
- Arms A/B/C/A-A re-run identically otherwise.
- Interpretation: whatever differences remain between arms are pure
  local seed sensitivity of the ray solver. The background model
  trajectory should be identical across arms up to nondeterminism -
  a direct measurable check that isolation is complete.
- Expected outcome prediction: all arms converge/freeze identically;
  C's abort should NOT recur (its death was gate-feedback-driven).

MODE 2 - PRODUCTION-CONTROL:
- Remove the L_NOPROBE store gates for the seed ray (it becomes a
  normal committing ray) while keeping the same seed arms. Probe
  telemetry stays for comparison.
- Explicitly define and test: jump veto applies normally (no waiver
  past first commit), FG_RAY_SEEDED set on first write, NC-cap path
  behaves on a seeded ray.
- Interpretation: does a seed strategy change a REAL ray's commit
  quality, downstream field, and global convergence - the production
  question the probe was built to proxy.
- Run A/A first: on a committing ray the seed difference enters the
  store, so the repeatability floor must be re-derived before B/C are
  interpreted.

## 9. Tests required before any implementation

- Bit-identity: with KN_NOCP_PROBE off, production runs must be
  bit-identical to current (corrections gated on the probe flag).
- With probe on, MODE1 arms must show: identical KN_RF_RESID* (prod-
  only variant), identical applied omega to print precision across
  arms, identical GS gate decisions across arms, store sums still 0.
- Replay check: OUTGEN omega-est line must be reproducible offline
  from a logged per-ray rho table (requires adding the per-ray sample
  print behind a diagnostic flag, or the pool itself).
- Counter test: synthetic ray with forced JUMP>veto verifies
  FG_COMMIT_REJ no longer counts probe rejects.
- Census test: GS_NBAD/NTOT decrease by exactly the probe's NBAD/NTT
  when isolation is on.

## 10. Provenance

Executable under test: exe/cmfgen_dev.exe, md5
3f37f22bd133468fb2854758d19abf12 (commit e1df4324 + 3 diagnostic
JRESC-guard modifications; see KNOP_THREE_SEED_EXPERIMENT.md).
Audited tree state: commit d369771e + this document.
Analysis scripts: scripts/parse_knop_seed.py;
telemetry: docs/audits/knop_seed_telemetry.csv.
Run evidence: cmfgen_validation/m15_2p82d_KNOP_seed_{a,b,c,aa}/
cmfgen/2.82days/OUTGEN + KNSEED binaries (paths in the experiment doc).
