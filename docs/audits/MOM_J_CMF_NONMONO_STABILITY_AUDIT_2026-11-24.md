# MOM_J_CMF_NONMONO_V1 stability audit (USE_LAM_ES=F moment route)

Date: 24-Nov-2026. Read-only audit of `plane/mom_j_cmf_nonmono_v1.f`
(7048 lines) plus its dispatch and VADAT wiring. Working tree was
dirty during the audit (uncommitted edits to
`control_variable_mod.f`, `rd_control_variables.f`,
`do_ng_accel_band_v4.f`, docs, tests); all line numbers cite the
working tree.

Question: where does the numerical stability of the non-monotonic
moment solve come from, where does it end, and what does each
stabilization cost the equations actually solved?

Route audited: `DO_NONMONO_SOB` + explicit `USE_LAM_ES=F` (opt-in;
absent-key default is T since 10-Nov-2026) and the moment half of
`DO_NONMONO_KNOP`. Dispatch: `new_main/mod_subs/comp_j_blank.f`
~:1119 (else-branch after the LAM_ES path); with the auto-sets
`DO_KNOP_REFRESH=T` and `DO_NONMONO_FREQ_MATCH=T` armed
(`rd_control_variables.f:3930-3942`). Status ledger already records
this route as FAIL on m15 (job 12316397, bounded non-contraction,
~183K persistently frozen cells) - this audit is the code-level
explanation of that signature, not a new run.

## 1. Executive summary

The solve has three layers with different stability verdicts:

1. Per-frequency tridiagonal solve (depth operator): VERIFIED
   stable by construction. The Knop generalized-opacity split
   produces an M-matrix with sign-definite off-diagonals, the
   boundary rows are guarded against the two degenerate closures
   that were measured breaking them, and every failure path
   (DGTSV INFO, non-finite, runaway, negative J, quarantine) ends
   in a counted repair with a sane fallback - nothing propagates
   silently.
2. Lagged frequency map (the fixed point the sweep and refresh
   iterate): PARTIAL. At compressing cells the row degenerates
   toward a pure frequency translation (same rho -> 1 stiffness
   the FG audit measured, rho_p90 ~ 0.97), and the code stabilizes
   it with a layered stack: LAG_BND row bound, XM_POS RHS
   positivity bound, EPS_LIM donor cap, GAM_FAC_MIN resonance
   floor, JNEXT_RELAX + adaptive omega, per-cell damp schedule,
   inner and outer AA, directional scheduling with write
   ownership, and the refresh's best-pass retention + divergence
   abort + jump veto. Each mechanism is documented and counted.
   Two of them (LAG_BND scaling, XM_POS) change the equation at
   the rows where they bind, so "converged" there means the
   bounded map's fixed point, not the unmodified equation's.
3. Sweep vs refresh consistency: one real defect and one latent
   guard gap (F1, F2 below). The refresh re-solve does not apply
   the same NMID_ON_HMID regularization as the sweep, so the
   refresh iterates a different map at sign-transition midpoints -
   a structural reason the retained-vs-sweep delta cannot close
   there.

Verdict: the local solve is sound; the fixed-point machinery is a
defensible, heavily-instrumented bounded-iteration scheme whose
stability is purchased with equation modifications at the hard
cells. The m15 signature is what this architecture predicts in the
translation-dominated regime. Fix list in section 6.

## 2. Layer 1 - the tridiagonal operator: stable by construction

Discretization per frequency slot (`mom_j_cmf_nonmono_v1.f`,
assembly ~:2900-3120, solve :3232):

- `chi_hat = chi + |GAM|/dLOG_NU` is always > 0; the advection
  split takes `PSI = DTAUONQ*|GAM|/dLOG_NU >= 0` on the diagonal
  and puts the signed coupling in the RHS as an upwind source.
  GAM >= 0 (expanding) reads `J_prev` from the fresh store;
  GAM < 0 (compressing) reads lagged `J_next`. Upwind on both
  axes.
- Interior rows are M-matrix rows: positive diagonal
  (`TB = chi_hat*DTAUONQ + PSI + |HU| + |HL|`), non-positive
  off-diagonals. The EPS midpoint scatter terms are signed but
  enter the RHS only (lagged), with MOM_EPS_LIM capping the
  debit at the donor row's positive balance so the pair stays
  conservative (:2750-2760, census MOM_EPS_CAP).
- DGTSV (partial pivoting) at :3232 with a fallback re-solve
  path at :3341; both failures restore the previous valid field
  (XM_FB) rather than the assembled RHS, and are counted
  (MOM_N_DGTSV).
- Boundary rows: inner BC verifies `DTAU(ND) > 0` before
  division (loud STOP on a zero optical-depth step, :3040s).
  Outer BC for OUT_BC_TYPE <= 1 keys the expanding/compressing
  branch on `A_BND` - the boundary row's own coupling
  coefficient - not on Q_NM(1); the comment at :2910-2928
  records the defect this fixed (sigma(1) ~ -1 misclassification
  injecting negative RHS). The HBC/NBC Eddington closure is
  floored to [0,1] only where the measured closure is negative
  or > 1 (`BND_HBC_FLOORED`, :2940-2960), counted in
  MOM_BND_HBC_FL and logged to MOM_BND_HBC_FLOOR - a physical
  closure replacement at a degenerate measurement, applied to
  the row only.
- Post-solve gates: non-finite rejection, runaway rejection
  (|dJ| > 1e3 relative), negative-J rejection with the
  fallback rebuild applying the same LAG_SCL and EPS cap to the
  same terms (:3359-3420). XM_BAD cells restore the previous
  iterate (KEEPNEG=F default) or quarantine (KEEPNEG=T).

This layer carries no found defects. The coefficient signs,
boundary keying, and failure containment all check out against
the Paper I discretization and the audit history embedded in
the comments.

## 3. Layer 2 - the lagged frequency map: stiff by construction, bounded by mechanism

At a compressing cell with `|GAM|/DLN >> chi`, the discrete row
is `J_k ~ J_{k+1} + O(chi*DLN/|GAM|)` - a near-pure translation
along frequency. The lagged map's spectral radius is then
`1 - O(chi*DLN/|GAM|)`, measured on the FG side at rho_p90 ~
0.97 (`CMFGEN_KNOP_LAGGED_MAP_NONCONTRACTION.md`, hypothesis H);
the moment map has the identical `chi_hat` structure, so the
same stiffness applies. No per-solve fix can remove it: it is a
property of the upwind discretization, not a bug.

The code's response is a stack of stabilizers, all gated,
counted, and documented in-module:

| Mechanism | What it does | Cost |
|---|---|---|
| MOM_LAG_BND (T, :2860-2884) | scales each row's whole lagged contribution by min(1, excess/weight), excess = TB-|TA|-|TC| | where it binds, the converged field solves the scaled map, not the equation; EPS implicit/explicit pair match becomes approximate (documented :2853-2858). Census: MOM_LAGB_CNT/MIN |
| MOM_XM_POS (T, :2897-2909) | clamps assembled RHS to >= 0 on interior rows | capped rows land on XM=0, suppressing the row's response; physically justified (nonneg source -> nonneg J), seeded by the production it57 failure. Census: MOM_XMNEG, MOM_XMWORST |
| MOM_EPS_LIM (T) | donor-cell cap on lagged EPS transfer | conservative pair, small bias at advection-dominated rows |
| MOM_GAM_FAC_MIN = 1e-2 | floors |1+sigma*K| at the resonance | at floored cells the row solves a regularized equation by construction; counted (PASS_NRES) |
| MOM_DLN_MIN = 1e-6 | floors the log-frequency spacing | at bound-free edge doublets (dLOG_NU ~ 2e-11) the advection is weakened ~2e4x and the row stays spatially coupled where the correct limit is J_k = J_upwind; bounded (unit test 'dbl' pins it) |
| MOM_GAM_PIN (T) | PSIPREV = PSI, explicit weight = implicit weight | makes the expanding-branch map non-expansive by construction |
| JNEXT_RELAX + adaptive omega (:6988-7004) | under-relaxed store writes; V2 estimator uses median same-map drift ratio, bounded [0.05,20], omega in [0.05,1] | feedback controller on the residual sequence - the FG audit flagged this controller class as a bounded-limit-cycle suspect (W5); the V2 median is more honest than the old cross-map ratio |
| per-cell damp schedule | slower writes where a cell failed | freezes hot cells toward their last good value |
| inner AA / outer AA | Anderson acceleration on the in-sweep composite and on committed refresh iterates; omega-change FHIST rescale (:7007+) | affine-map assumption degrades under an adapting omega; rescale applied, documented |
| MOM_RF_DIR='ALT' + MOM_RF_OWNW | alternating pass ordering; writes suppressed on the non-owned direction | suppressed cells hold refresh values; requires MOM_RF_ENABLED so the freeze is bounded by refresh cadence |
| MOM_KNOP_REFRESH_V1 (:5163-7048) | re-drives the solver on captured KN_M_* coefficients in alternating frequency order; best-pass retention, divergence abort (MOM_RF_DIV_BULK), 50x jump veto, bulk-exit (MOM_RF_BULK_EXIT) | converges the compressing lagged coupling within a great iteration; the retained store is the best pass, not the last |
| MOM_RF_ACTIVE (T) | active-set dilation over the hot rows | skips provably-cold frequencies; mask resets each refresh |
| checkpoint MOMJNXT_GS (:4364+) | magic/ver/tag/dims + full R_NM/R_PNT_NM fingerprint + finite-and-nonnegative value validation | rejects zero the store; a stale-mesh checkpoint cannot seed the map |

Assessment: the stack is honest engineering - every repair is
counted and every floor is justified in a comment with the
failure that seeded it. The structural caveat is that LAG_BND
and XM_POS are equation modifications at the cells where they
bind. At a shock row (PSI ~ 0 from the resonance floor, excess
~ 0), LAG_BND nulls nearly the whole lagged source and the
converged "fixed point" is the clipped map's. That is the right
engineering trade against a non-contractive map, and the code
says so in the comments (:2853-2858), but it means a converged
m15 field at flagged rows is not the unmodified moment
solution. The bias is bounded and reported via MOM_LAGB_*.

## 4. Findings

### F1 - Refresh solves a different map than the sweep at sign-transition midpoints (defect)

The sweep regularizes `NMID_ON_HMID` (g = N/H, the signed ratio
that feeds `|1+sigma*g|` and the H-row direction split
`GAMH_SIGNED = sign(1+sigma*g)`) as signed with a magnitude cap:
`SIGN(MIN(|g|,1e12), g)` (:1826-1831). The comment at :1815-1819
explains why the sign must be kept - the solver uses
`|1+sigma*g|` and the sign of `(1+sigma*g)` picks the H upwind
direction.

The refresh reload path applies a different regularization -
`MIN(1, MAX(0, g))` - at three sites:

- :5703 (current-pass reload from KN_MNMH)
- :5730-5734 (save-slot reload for the L-1 pass)
- :6908 (module-state restore before the next great iteration)

Consequences where |g| > 1 (the SIGN(1E10,N) degenerate-H fill
in `get_moms_non_rel.f:210-218` guarantees such cells exist):

- `sigma < 0, g > 1`: sweep sees `1+sigma*g` possibly < 0
  (compressing H row, reads redder store); refresh clamps g to
  1, sees `1+sigma` and may classify expanding. Direction flip.
- `sigma > 0, g < -1/sigma`: same flip in the other direction.
- Everywhere |g| > 1 the midpoint coefficients (W, HS) differ in
  magnitude, not just sign.

So the refresh converges a map that is not the sweep's map at
exactly the cells the direction split exists for. The best-pass
retention then stores a "fixed point" that is off-manifold for
the sweep equations, which is a permanent floor under the
retained-vs-sweep delta at sign-transition cells - a concrete
contribution to the observed bounded non-contraction that is a
defect, not stiffness.

Fix direction: apply the same `SIGN(MIN(|g|,1e12),g)`
regularization in the three refresh reload sites (the captured
KN_MNMH store is already sanitized; the clamp just needs to
stop changing it).

### F2 - JUP_NM lacks the finite guard the comment claims it has (latent)

The compressing-branch J upwind read fills JUP_NM from
`JNEXT_RF_READ` / `J_NEXT_FRQ` at :2653 and :2664 with no finite
test. HUP_NM gets the explicit non-finite -> 0 + MOM_INEXT_CAPPED
treatment (:2697-2704), and the comment at :2694-2696 states
"J_NEXT_FRQ gets the same finite test" - it does not.

Within-run this is mostly latent: every write path into
J_NEXT_FRQ is guarded, and the checkpoint restore validates
finite-and-nonnegative. But the guard exists for HUP because a
stale or corrupted store entry would enter the tridiagonal RHS
unchecked - the same holds for JUP (it feeds both the PSI*J_up
term and the EPS<0 lagged scatter, :2716-2730). A NaN/Inf here
propagates through DGTSV into a row rejection - caught, but as
a cell failure with the wrong attribution, and only
post-solve. One-line fix matching the HUP pattern.

### F3 - Stale default-value comments (doc defects, 3x)

- :641 "MOM_RF_KEEPNEG (VADAT key, default TRUE)" - DATA and
  VADAT defaults are FALSE (:955, `rd_control_variables.f:1213`;
  flipped F 06-Nov-2026 because the m15 stall carried
  quarantined negative J into STEQ).
- :676 "MOM_RF_ACTIVE (VADAT key, default FALSE)" - default is
  TRUE (:958, :1249-1252).
- :697 "MOM_RF_DIR (VADAT key, default 'AUTO')" - DATA and
  VADAT defaults are 'ALT ' (:959, :1309). The AUTO census
  machinery (:2316, :4293-4303) exists but is not the default.

Per project rules a stale comment is a defect; these three all
misreport the active stabilization defaults.

### F4 - Red-edge compressing boundary is a silent free outflow in frequency (caveat)

JUP_NM is zeroed per call (:2618) and the compressing read is
gated on `FREQ_INDX+MOM_NM_J(I) <= NCF_IN`; for the reddest
slots with compressing rows the upwind contribution is silently
dropped - a zero-field boundary at nu_min on the advection
axis. The refresh counts it (PASS_NOFF); the sweep does not.
Bounded and only wrong if the grid's red edge carries real
field; worth a sweep-side counter for visibility.

### F5 - MOM_DLN_MIN doublet floor computes the wrong limit, bounded (documented)

`MOM_DLN_MIN=1e-6` keeps |dLOG_NU| away from zero. At the
bound-free edge doublets the continuum grid places (dLOG_NU ~
2e-11, ~18% of slots), the correct discretization degenerates
to J_k = J_upwind; the floor instead keeps the row
spatially coupled with the advection term weakened ~2e4x. The
'dbl' unit test pins the bounded deviation. This is a
documented numerical regularization, not a defect, but it is
a model-level bias at every edge doublet and should stay on
the clamp-policy inventory.

### F6 - The m15 FAIL signature is consistent with the architecture, not a hidden bug

~183K persistently frozen cells + bounded non-contraction maps
onto: (a) resonance-floor cells where PSI collapses and the
map is non-contractive by construction (LAG_BND rows); (b)
translation-dominated compressing rows at rho ~ 1 where no
bounded number of passes converges; (c) under the KEEPNEG=T
era, quarantined cells excluded from metrics that could never
move (the F flip on 06-Nov-2026 stores repaired fallbacks
instead, so cells keep drifting); (d) F1 holding a residual
floor at transition cells. Nothing in the audit contradicts
the recorded verdict: the route is stable in the sense that
it never emits an unguarded bad field, and it does not
contract in the hard regime.

## 5. Clamp/floor inventory vs project clamp policy

AGENTS root section 6 prohibits artificial clamps that hide
physical problems. The full inventory for this solver:

| Clamp | Value | Hides physics? | Justification on file |
|---|---|---|---|
| CHI_HAT floor | 1e-30 | no - void-cell guard | prevents 0/0 on the Eddington ratios |
| QF clamp | [1e-20,1e30] | no | overflow guard on R^2-weighted flux factor |
| LOG args | 1e-290 floor via MERGE | no | prevents -Inf on empty cells |
| K_ON_J | [1e-4, 1/3->1] | marginal | physical Eddington bound; can only over-flag compression, counted |
| HBC/NBC row floors | [0,1] | no - physical closure replacement | degenerate measured closure, logged + counted (:2940s) |
| NMID_ON_HMID cap | |g| <= 1e12 signed | marginal | regularizes degenerate-H fill; sign kept (the F1 defect is the refresh breaking this rule) |
| MOM_GAM_FAC_MIN | 1e-2 | yes - regularizes the resonance | documented zero-GAM regularization; counted (PASS_NRES) |
| MOM_DLN_MIN | 1e-6 | yes - see F5 | documented; bounded |
| MOM_XM_POS | RHS >= 0 | yes - changes the RHS | documented; physically-motivated positivity bound; counted |
| MOM_LAG_BND | min(1, excess/weight) | yes - changes the map | documented consistency cost; counted |
| MOM_EPS_LIM | donor cap | marginal | conservative pair cap; counted |
| runaway/negJ/50x veto | rejection thresholds | no - failure containment | restores previous valid field; counted |
| H clip | |H| <= 0.99999 J | marginal | causality clip on reconstructed H; V11-inherited |

Every clamp in the table is flag-gated or counted; none is
silent. The two equation-modifying ones (XM_POS, LAG_BND) are
the stability purchase itself - flag for the user, not for
removal: removing them re-opens the demonstrated non-contractive
regime.

## 6. Recommendations (ordered, all small)

1. F1: make the three refresh reload sites use the same signed
   1e12 regularization as the sweep (:5703, :5730-5734, :6908).
   Highest-value fix: removes a structural floor under the
   retained-vs-sweep residual at the cells that matter most.
2. F2: add the HUP-style finite test to the JUP_NM fills
   (:2653, :2664) and correct the :2694 comment.
3. F3: fix the three stale default comments (:641, :676, :697).
4. F4: extend the PASS_NOFF census to the sweep path (cheap
   visibility into how much boundary is exercised).
5. Report: surface MOM_LAGB_CNT/MIN and MOM_XMNEG per iteration
   alongside the existing GS drift line - the equation-modifying
   clamps are the interesting census on this route.
6. The route stays opt-in/FAIL until a contraction-capable map
   replaces the lagged compressing coupling (the FG audit's
   lever 1: solve the per-ray frequency chain directly rather
   than iterate at rho ~ 1). Nothing found here changes that.

## 7. Test coverage assessment

`unit_tests/f90/drv_mom_solve.f` + `python/test_mom_solve.py`
exercise the per-frequency solve hard: first-principles
coefficient recompute, tridiagonal residual, admissibility,
drift contraction, the doublet limit, the seeded it57 XM_POS
case, ND=48 refinement invariance, checkpoint round-trip, and
the thin-compressing-row regime pinned as bounded
non-contraction. `test_mom_gs.py`, `test_mom_aa_state.py`,
`test_mom_sweep_table.py`, `test_obs_nonmono_dispatch.py`,
`test_vadat_keys.py` cover the surrounding machinery.

Gap: no unit test drives MOM_KNOP_REFRESH_V1 itself (the
drivers call the solver directly per iteration), so the F1
map mismatch is invisible to the current battery. A driver
that captures KN_M_* state, runs one refresh pass, and
compares coefficients against a direct re-solve would pin it.

## 8. Dated verdict

24-Nov-2026: PARTIAL - stable per-solve, bounded-not-contracting
fixed point. The tridiagonal operator and failure containment
are VERIFIED sound. The lagged frequency map is stiff by
construction (rho -> 1 in compression) and is stabilized by a
documented, counted stack of mechanisms that includes two
equation-modifying bounds; "converged" fields at LAG_BND/XM_POS
rows are the bounded map's fixed point. One genuine defect found
(refresh/sweep NMID_ON_HMID regularization mismatch, F1) plus a
latent guard gap (F2) and three stale comments (F3). The m15
FAIL recorded in docs/CHECKLIST.md is consistent with this
architecture; the recommendations in section 6 are the fix
list, none of which converts the route to contraction in the
translation-dominated regime.
