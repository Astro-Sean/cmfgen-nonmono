# SOB slow-mode analysis - where the tail lives

Date: 11-Nov-2026. Basis: FRESH_sob (converged it722, DIE-on),
FRESH_sob_dieoff (running, it241+ at writing), SOB_l14_lte_adp2_tauto
(converged reference), plus the X1-X4 warm arms once they produce
data. Audit basis: ES_ALI_BA_AUDIT.md.

This document addresses Phase 4 of the acceleration task: whether
the remaining slow convergence mode sits in the continuum
electron-scattering channel (what ES_ALI/ES_ALI_BA attack) or in
the Sobolev line-rate / inter-ion coupling (outside their scope).

## 1. What the converged run actually did

FRESH_sob ITER_HISTORY tail (it689-722):

- maxapp falls 5.9e-2 -> 4.9e-2 over the last ~30 iterations with
  sawtooth spikes at it695-698 (maxapp 12.0) and it711-713 (se
  spikes to 1.7e12). The spikes coincide with BA rebuilds /
  extrapolation events, then the geometric tail resumes.
- se (BA_SE_RES, max |STEQ| over all species/depth) stays at
  6e10-3e12 even at "convergence". It is diagnostic-only: the
  convergence gate is maxapp < EPS_TERM, not the SE residual.
  The se floor is dominated by the degenerate inner-boundary
  equations noted in solveba_v13.f (the ~1E262 saturation comment;
  on this model the floor sits ~1e10-1e12).
- se2 (RES2, K-th largest per-depth column residual,
  K=min(ND,max(4,ND/10))=7 for ND=70) tracks se at roughly half
  its magnitude - the worst columns are not a tiny boundary set,
  the residual is distributed.
- nrej=0 throughout the tail: no BA_RES_ACC rejections (flag was
  off anyway; the column is the resacc counter).

Implication: "convergence" here means the applied population
correction fell below EPS_TERM. The raw SE residual never gets
small in absolute terms; what matters is whether the *correction*
driven by it contracts. The slow mode is therefore a property of
the fixed-point map's contraction rate, not of residual magnitude.

## 2. Where the largest corrections live (species x depth)

Extracted from every `largest increase/decrease % @(J,I)` line in
FRESH_sob OUTGEN (1340 records over 722 iterations):

Species histogram (times a level was the max correction):

| species | count | ion stage |
|---|---|---|
| AlV    | 298 | Al V  |
| OSIX   | 212 | O VII |
| OIII   | 131 | O III |
| CV     | 131 | C V   |
| OI     | 112 | O I   |
| Fe2    |  90 | Fe II |
| SIII   |  48 | S III |
| NeIII  |  36 | Ne III|
| ArIII  |  36 | Ar III|
| CrIII  |  32 | Cr III|
| CIV    |  29 | C IV  |
| MgV    |  25 | Mg V  |
| OIV    |  24 | O IV  |
| CaV    |  24 | Ca V  |
| He2    |  23 | He II |
| CrSIX  |  17 | Cr VI |

Depth histogram of the same records: depth 1 leads (218), then a
mid-outer cluster at depths 6-16 (114+105+78+72+69+61+60+56+52+45
= ~800 records across d6-d16) and a second cluster at depths
53-59 (~160 records). Two active zones: the inner boundary region
and the metal-line-forming shell around d6-d16 plus the outer
recombination/front region d53-d59.

Reading: the max corrections are dominated by *high ionisation
stages of metals* (Al V, O VII, C V, Cr VI, Ca V, Mg V) plus the
oxygen sequence O I/III/IV/VII and Fe II. This is the same
population sector the POPOXY audit flagged as over-ionised
(O I/II/III -> 0, O V/VI dominant). The slow tail is therefore
plausibly an ionisation-balance / line-blanketing problem among
trace metals, not a pure continuum-ES problem.

Caveat: "largest correction" measures where the applied step is
biggest, not where the residual is worst. RES_LOC (the
BA_SE_RES argmax, printed under DO_BA_RES_ACC + CONV_DIAG_LVL>=1)
gives the residual-side localisation; X4 produces it directly.

## 3. Continuum vs line-rate coupling - current evidence

For the continuum-ES hypothesis:

- The audit established ES_ALI/ES_ALI_BA only touch the
  continuum coherent-scattering channel. The Sobolev line
  Jacobian (SECTION='LINE' frequencies) uses the stock
  CONT_VEL=F path and is untouched by all three flags.
- The species histogram above is dominated by bound-bound /
  ionisation-balance adjustments of metals - the sector where
  Sobolev line rates and photoionisation balance couple most
  strongly.
- adp2_tauto runtime evidence: ES_ALI converged the J loop to
  dJ/J ~ 2.7e-5 with zero fallbacks - the continuum inner
  iteration was already fast; the outer tail persisted anyway.
  That run still took the same slow trajectory, which already
  weakly suggests the bottleneck is not the continuum J lag.

Against (i.e. continuum still suspect):

- ES_FAC=0.1 in the recipe means electron scattering is a major
  continuum opacity source in the dense shell; the radiation-
  pinned slow mode of Rybicki & Hummer is exactly this regime.
- se2 staying at ~50% of se means the residual is spatially
  broad, consistent with a diffuse radiation-field lag rather
  than a few bad line rates.

Decisive test = the X-matrix: X2 isolates ES_ALI (continuum J
loop only), X3 adds the BA-diagonal amplification, X4 adds the
resacc guard. If the tail contraction rate (maxapp ratio per
iteration in the post-transient regime) is unchanged across
X1->X4, the slow mode is outside the continuum-ES channel and
the Sobolev-line Jacobian / JFNK-CPL is the next target
(decision gate, ES_ALI_BA_TODO X4).

## 4. Diagnostics now available per arm

| quantity | source | gate |
|---|---|---|
| maxraw/maxapp/maxatt | ITER_HISTORY cols 7-9 | always |
| se / se2 (SE residual, RES2 rank) | ITER_HISTORY cols 10-11 | always (BA calls) |
| trust / nrej | ITER_HISTORY cols 12-13 | always |
| fg_rfd/badf/churnf/cap, mom_* | ITER_HISTORY cols 15-28 | SOB route |
| ES_ALI freqs/dJ/fallbacks/floor hits | OUTGEN `ES_ALI:` line | DO_ES_ALI=T (X2-X4) |
| ES_ALI_BA floor hits run | do_var_cont ALO banner | DO_ES_ALI_BA=T (X3,X4) |
| RES_LOC worst (var,depth) + EQ_DESC | OUTGEN BA_RES_ACC block | DO_BA_RES_ACC=T (X4) |
| BA_COL_RES per-depth profile | ITER_SUM export | CONV_DIAG_LVL>=1 (default 2) |
| largest correction species/depth | OUTGEN BA-TRACE blocks | always |

Phase 3 gap notes: particle/charge conservation are enforced by
the SE solve structure (conservation rows are equations in the
same system; RES_LOC would name them via EQ_SPECIES rows); a
dedicated conservation-error column does not exist - the
MASS_FRACTION_SUM_CHK file provides a posteriori inventory.
Energy balance: FIX_T=T in this recipe so the temperature
equation is not iterated; radiative equilibrium residual is not
in COL_RES. heating/cooling diagnostics exist per-iteration in
OUTGEN but are not part of the convergence gate. Recorded as a
known diagnostic gap, not a blocker: the experiment compares
arms at identical physics, so relative convergence is measurable
without absolute energy residuals.

## 5. Resolution (11-Nov-2026, X-arm data)

The X1-X4 arms answered the continuum-vs-line question:

- X2 (DO_ES_ALI only) converged clean in ~30 new
  iterations / 33 min wall. X1 (identical, all accel off)
  is still running at ~it774 with maxapp ~0.65 and
  repeated reject spikes - the unaccelerated warm arm
  reproduces the cold-run slow tail.
- X3 (ES_ALI+ES_ALI_BA) converged in ~28 new iterations /
  35 min - statistically identical to X2. The ES_ALI_BA
  floor engaged 50 times with no tail effect.
- VERDICT: the slow mode WAS the continuum electron-
  scattering J lag after all. The metal high-ion species
  histogram (section 2) marks where populations are stiff,
  but the radiation-field lag pinning them was the ES
  continuum channel - exactly what ES_ALI's OAB update
  converges. The Sobolev line Jacobian coupling did not
  need to be the bottleneck on this model.
- X4 (resacc) still running: 5 rejections, trust 1e4->312,
  descending with recovery cost - stability mechanism,
  not accelerator.

Implication for JFNK-CPL: the coupled tier remains the
right tool IF a future model shows the line-rate channel
dominating (e.g. lower ES_FAC, stronger line blanketing,
or KNOP route where ES_ALI is excluded). On THIS model
the cheaper ES_ALI fix suffices; JFNK-CPL stays a
proposal, not a commitment (see JFNK_CPL_PROPOSAL.md).
