# USE_LAM_ES=F moment route - checklist-driven diagnosis

Date: 26-Nov-2026. Read-only audit against
`checklists/use_lam_es_eq_f.md`, executed alongside the
`MOMENT_LITERATURE_AUDIT_2026-11-25.md` and
`MOM_J_CMF_NONMONO_STABILITY_AUDIT_2026-11-24.md` evidence base.
Line numbers cite the working tree (mom file now 9169 lines).

Live evidence: r01 arm `m15_2p82d_VAL_sob_lamesF_blkgmres_r01`
(CMFGEN 13671873, USE_LAM_ES=F + MOM_BLKGMRES=T +
JNEXT_RELAX=0.1 + KN_ADAPT_RELAX=F, cold start, FIX_T=T,
VEL_LAW=14, ND=70) - FAILED at it~349 after 10 h wall with
Conv 67% plateau, every refresh non-converged, outer
luminosity oscillating 5E11-3E12 against inner 4.27E8.

## 0. Route structure (verified in code)

Per frequency, inside the INACCURATE loop of `comp_j_blank.f`:

1. `TA = ETA_CLUMP + CHI_SCAT_CLUMP*RJ` (coherent ES, :847) -
   the scattering emissivity is lagged at the current iterate.
2. `FG_J_CMF_V13` runs (:928) - the stock formal solver,
   documented in-code at :891 as "expansion-biased at sigma<-1".
   It supplies the Eddington closures `K_ON_J`, `NMID_ON_J`,
   `NMID_ON_HMID`, boundary HBC/NBC via MOD_RAY_MOM_STORE and
   `GET_MOMS_NON_REL` (mom_j_cmf_nonmono_v1.f:2002).
3. `MOM_J_CMF_NONMONO_V1` runs (:1122) - Knop moment solve.
4. `ES_ALI_MOM_ON` OAB update (:1132-1156) - active under
   auto-set `DO_ES_ALI=T` + `COHERENT_ES` on this route.
5. Up to 20 passes/frequency; RJ -> EDDFACTOR -> populations.
   `MOM_KNOP_REFRESH_V1` re-drives the solve on captured KN_M_*
   coefficients each great iteration (cmfgen_sub.f:3367).

## 1. Phase 1 - temperature correction & preconditioning

The checklist's root-cause frame (temperature-correction
Jacobian loses diagonal dominance when the ES lambda operator
is dropped) does not match this route: `FIX_T=T` on the
failing decks - the temperature correction is never exercised.
The m15 FAIL signature cannot come from a T Jacobian that does
not run. The items still map to real mechanism analogues:

### 1a. "Add dS_es/dT to the Jacobian diagonal" -> CONFIRMED mechanism, wrong location named

In the moment rows the scattering J term is already implicit
in the analytic sense: the diagonal carries
`DTAUONQ*MAX(MOM_ABS_FLR,1-COH_VEC)` with
`COH_VEC=ESEC_NM/CHI_NM` (:2676-2703) - i.e. the absorption
fraction chi_a/chi. For coherent scattering the chi_s*J
cancellation is exact: at epsilon->0 the local J-destroying
term vanishes by construction (photon conservation). The
scattering emission enters only lagged, through
`SOURCE_NM=ETA_NM/CHI_NM` (:2145) with ETA holding
`ETA+chi_s*RJ` from :847. So at ES-dominated cells the row's
contraction must come entirely from the spatial and frequency
couplings - and at compressing cells the frequency coupling is
itself lagged (`PSI*J_next` in XM, :2969). That is the
checklist's predicted stiffness, confirmed structurally - the
defect is not a missing derivative term, it is the lagged
splitting of the dominant couplings.

The code's compensation for the lagged scattering iterate is
`ES_ALI_MOM_ON`: an Olson-Auer-Buchler step on the moment
diagonal `dJ/dETA = DTAUONQ/(TB*CHI)` (:4307-4336) with
`DEN = max(1 - DJDETA*chi_s, 0.01)`. VERIFIED installed and
VERIFIED unstable on r01 OUTGEN: every iteration reports
`ES_ALI: freqs unconverged= ~6139-6665` of 8331 (about 80%),
with `max final dJ/J` growing 5E3 -> 1E84 and
`J-floor hits run` accumulating ~3E5-4E5 per iteration. At
epsilon->0, DJDETA*chi_s -> 1, DEN hits the 0.01 floor, and
the update amplifies the lagged map's pass-to-pass noise by up
to 100x - then feeds it back into TA on the next pass. The
prescribed remedy is present and is itself a divergence
driver on a non-contracting map.

### 1b. "epsilon denominator must not drop kappa_es" -> VERIFIED correct

`COH_VEC = ESEC_NM/CHI_NM` (:2149); the absorption term is
`MAX(MOM_ABS_FLR, 1-COH)` with MOM_ABS_FLR default 0.0
(:1157). kappa_es is not zeroed; the epsilon->0 limit is the
physically correct one. The ESEC>CHI interpolation overshoot
is clipped and documented (:2661-2670). Boundary row 1
carries the symmetric MIN form (:3320). No defect.

### 1c. "Enforce Ng acceleration, order 5-7" -> misdirected

NG_ABSTAIN=T (auto-set under SOB) is a guard on NG
extrapolation (abstain on zero-diff or non-positive
extrapolation), not an NG shutdown. The moment route has its
own accelerators: inner/outer Anderson acceleration in the
refresh, same-direction contraction estimator, adaptive omega.
Measured on r01: AA applied 0 times per refresh -
`||corr||=1.7E2 exceeds ||f||=2.3E1` skip, history restarts
on newly damped cells. Acceleration cannot contract a rho ~ 1
map; the map is the defect, not its accelerator depth.

## 2. Phase 2 - formal solver & source function

### 2a. "Smooth S_es at dv/dr=0 turning points" -> REJECTED as stated

No source smoothing exists and none is derived: the
turning-point treatment is (i) sign-keyed local upwinding
(verified against Paper I in the 25-Nov audit) and (ii)
`MOM_GAM_FAC_MIN=1E-2` - the `|1+sigma*K|` floor at the
zero-GAM resonance (:2291-2298, :2417-2419). On r01 the floor
binds at 10,144 cells every pass (gam-floor census). A
Savitzky-Golay or harmonic-mean smoothing of S_es is not the
cited literature's mechanism (Rybicki & Hummer 1991 is an ALI
method) and would be an equation modification per root
AGENTS.md section 6. The real turning-point defect is
measured: `nu-census amax 9.947E-01` - the contraction factor
of the compressing map sits at rho ~ 0.995, not noise ringing.

### 2b. "CMF opacity transformations, dropped (1-v/c)" -> VERIFIED not defective

The moment equations carry the O(v/c) coupling through
`CON_GAM = (v/r)/c` (:2285), `GAM = CON_GAM/chi_hat/dlnnu`
(:2500), and the gradient factors `1+sigma_true*K` /
`1+sigma_true*NMID` - the same construction as stock
MOM_J_CMF_V11, which assumes expansion. The generalized
opacity `chi_hat = chi + |GAM|/dlnnu` (:2468) is the Knop
Paper I Eq 7 positive split. No frame-transformation factor
was dropped relative to the V11 baseline.

### 2c. "dtau on total optical depth" -> VERIFIED

`DTAU` is built on `chi_hat*QF` (:2578-2582) - absorption +
scattering + the frequency-advection opacity, i.e. more than
tau_total. There is no short-characteristics step limiter on
this route (the solve is a depth-tridiagonal PDE per
frequency). The frequency-side analogue is MOM_FREQ_MATCH
(see 4a).

## 3. Phase 3 - SE & mean intensity

### 3a. "SE preconditioner diagonal / continuum rates" -> VERIFIED not route-specific

No route-specific SE preconditioner exists. Populations
consume J via the standard EDDFACTOR/RJ channel; the
moment-store rollback `RESTORE_MOM_JNEXT_LATEST` on rejected
BA steps is wired at solveba_v13.f:801,920, and
`GET_MOM_GS_MATURITY` gates the warm-up check
(solve_for_pops.f:1149). USE_LAM_ES=F adds nothing
off-diagonal to the rate machinery. What the populations see
is the committed J_NEXT store - honest state, but the store
stops moving when refreshes stall (the bounded stall).

### 3b. "J floor J >= min(B(Tmin),1e-30)" -> REJECTED as a remedy

The implemented policy is stricter and more honest:
`MOM_XM_POS` floors the assembled RHS at 0 (a negative RHS
means lagged debits outweighed the local source - counted,
~16-17k rows/pass on r01); negative-J post-solve rejection
restores the previous iterate (KEEPNEG=F default), ~11k cells
per blk-gmres call on r01; a quarantine path exists. A floor
at B(Tmin) would inject a thermal radiation floor - an
equation modification that hides the real deficit (lagged
upwind debits at compressing rows). Forbidden by the clamp
policy unless derived; not indicated.

## 4. Phase 4 - non-monotonic checks

### 4a. "Resonance point search must not skip" -> N/A literally; analogue VERIFIED with caveat

Resonance search is a Sobolev concept; the moment analogue is
the frequency-coupling graph. It is implemented: per-row
sign-keyed upwind (`Q_NM<0` reads lagged `J_next`,
:2863-2874) plus `MOM_FREQ_MATCH` characteristic-foot
rematching across `N_CAP` bins (:2379-2398, :2444-2465) with
CEILING keeping the foot inside the stencil. Measured on r01:
386,787 of ~1.6M rows frequency-matched, largest bin offset
8166 of 8331 - at the CDS the characteristic foot reaches the
red edge. 184 off-grid reads get `JUP=0`: a free-outflow
frequency boundary (literature audit F8 - physically
unverified). Integer-stride matching is not the Ceulemans
boundary machinery (F9).

### 4b. "Integration weights sum to 1" -> N/A literally; two real analogues

No short-characteristics quadrature exists on this route.
Analogue (i): the Eddington closures come from V13's ray
moments via GET_MOMS_NON_REL - V13 is expansion-biased at
sigma<-1 (comp_j_blank.f:891), so the moment equations'
closures are measured by a solver not built for this regime
(literature audit F10; degenerate fills SIGN(1E10,N) in
get_moms_non_rel.f:216, clamps at mom :2016-2026). Analogue
(ii): energy bookkeeping defects live in the equation-
modifying guards - LAG_BND (xs<=0 on ~244-250k rows/pass),
EPS_LIM (~67-100k cells/pass), XM_POS (~17k rows/pass) -
counted, not silent.

## 5. The actual dominant defects (diagnosis)

Ordered by measured residual damage (augres arm 13685499
plus r01 refresh detail):

1. Lagged frequency coupling is rho ~ 1 by construction at
   compressing cells (`nu-census amax 9.947E-01`; compressing
   passes stall at residual ~4E3; every r01 refresh exits
   "both directions stalled"; store retained). The compressing
   map is a near-pure frequency translation - a property of
   the upwind discretization, not a bug. The only
   contraction-capable lever identified is the augmented
   implicit J+H solve (FIX 4; prototype measured rho 2.55
   lagged -> exact solve unmassed on the captured hcon case).
2. H_NEXT_FRQ is not re-derived after the block J solve
   (MOM_BLK_SOLVE :8129-8270 updates J only); the committed
   J/H pair does not satisfy the H rows - augres measured
   s2H ~ 1E14-1E56, orders of magnitude over the J residual.
3. ES_ALI_MOM_ON is installed but counterproductive on this
   route (see 1a): ~80% of frequencies never satisfy
   ES_ALI_TOL; the DEN-floor amplification injects the map's
   oscillation into the next pass's source.
4. Equation-modifying guards bind at scale (LAG_BND ~245k
   rows xs<=0, EPS_LIM ~1E5, XM_POS ~17k, GAM floor ~10k) -
   a converged field there is the bounded map's fixed point,
   not the target equations'. Per the 25-Nov audit these are
   continuation devices, not convergence evidence.
5. Closures measured by expansion-biased V13 at sigma<-1 -
   consistent input, wrong-regime provenance (F10).

Checklist-frame verdict: the checklist's physical intuition
(epsilon->0 scattering dominance produces a stiff map and the
missing ES implicitness needs compensation) is correct in
spirit, but every named mechanism either does not exist on
this route (temperature Jacobian, SC weights, resonance
search, dtau limiter), is already correctly handled
(epsilon bookkeeping, dtau, ES-ALI present), or is an
unjustified equation modification (source smoothing, J
floor). The true failure is the lagged H/EPS coupling's
rho ~ 1 map plus the non-verified committed field - the same
verdict the 24/25-Nov audits reached, now re-measured on the
r01 arm.

## 6. Recommended next steps (no code changed)

1. Do not implement the checklist's literal prescriptions;
   they target a different architecture. The actionable items
   map onto the existing FIX-4/augmented-solve plan in
   `docs/plans/MOMENT_H_EPS_IMPLICIT_DERIVATION_2026-11-25.md`
   and the damage ranking already recorded in CHECKLIST.md.
2. Cheap controlled arm worth considering: `DO_ES_ALI=F` on
   the r01 deck to remove the OAB amplifier and isolate how
   much of the pass noise is self-inflicted (one flag, bounded
   question, no code change).
3. H reconstruction after MOM_BLK_SOLVE (audit F1/FIX 4
   stage 2) remains the largest measured residual damage.

## 7. Dated verdict

26-Nov-2026: DIAGNOSIS COMPLETE - the checklist's mechanism
(epsilon->0 scattering stiffness without an implicit ES
operator) is confirmed as the correct class of defect, but
the named prescriptions are N/A or rejected for this
architecture. The verified failure chain on m15 r01 is:
lagged rho~0.995 compressing map + unreconstructed H +
counterproductive ES-ALI OAB + equation-modifying guards at
scale + V13-sourced closures. Route stays FAIL/opt-in; the
contraction-capable path is the augmented J+H implicit solve
already queued as FIX 4. No source files were modified.
