# CMF_FLUX V12 vs KNOP formal-solver comparison (Phase 2.4)

Status: IN PROGRESS.  Controlled observer-frame comparison of
`FG_J_CMF_V12` (clamped-sigma formal solver) against
`FG_J_CMF_NONMONO_V1` (Knop et al. 2009 arbitrary-velocity solver)
on identical frozen model inputs.

Baseline source revision: `19bfdd22b9bbfc1652bff9f5b714c70ce8e4a7bc`
(two commits newer than the Phase 2.3 audit commit `f6d99e78`; the
newer commits record the KN-SEED telemetry fix and resubmission).

Frozen executables (this phase):

| role | freeze dir | cmf_flux sha256 |
|------|-----------|------------------|
| intermediate (superseded) | `8a5737badd61716c_60428671b9551a4d` | `60428671b9551a4d...` |
| key-fix rebuild (superseded) | `8a5737badd61716c_3638b740fd42c410` | `3638b740fd42c410...` |
| final comparison exe | `8a5737badd61716c_eacd7bc0525a1195` | `eacd7bc0525a1195...` |

The superseded hashes are recorded because the first-attempt render
matrix ran on them; see section 5 notes.  The final exe adds the
JFRQ/SCAN diagnostic records and the absent-key reader fix; its
numerics on the submitted PINITs are identical to the earlier build
(the difference is absent-key tolerance plus new record types).

## 1. Executive summary

On identical frozen populations, the Knop et al. 2009 formal solver
produces an inner-boundary radiation field that is systematically
below the V12 field everywhere the V12 field is inflated, agrees
with V12 identically where V12 is valid, and cuts the rescue-event
count by a factor of ~17 (4,529 vs 75,738 total JEV lines across
both ES sweeps).

Headline numbers (74,620 aligned continuum frequencies, mature
sweep):

| quantity | A: V12+res | B: V12 nores | C: KNOP nores | D: KNOP+res |
|----------|------------|--------------|---------------|--------------|
| med max\|J\| | 1.000e-3 (pinned) | 3.05e-3 | **2.68e-5** | 2.68e-5 |
| p95 max\|J\| | 1.000e-3 | 7.01e-2 | 7.64e-4 | 6.67e-4 |
| absolute max | 1.000e-3 | 9.48e-2 | 8.72e-2 | 1.000e-3 |
| fired triggers | 30,814 | 46,159 | **2,192** | 1,752 |
| JEV lines | 75,738 | 92,318 | 4,529 | 4,080 |
| med J/zeta | 289 | 850 | **9.6** | 9.4 |
| OBSFLUX sha256 | bd78a108... | fb4145e6... | 6adf9511... | 67a77641... |
| spectral-mass vs A | - | +3307% | +91% | +78% |

Interpretation:

- V12 without the restore produces boundary fields up to 95x the
  rescue cap and a spectrum distorted by a median factor of ~89.
  The restore is load-bearing in production: it holds the artifact
  at the cap rather than curing it.
- KNOP fields are ~5x below the V12 pinned field at the median and
  ~60x below the V12 unrescued field; KNOP sits near or under the
  electron-scattering ceiling (J/zeta ~ 435 at depth 70 where
  omega=0.9977) while V12 exceeds it.
- KNOP does not eliminate over-cap fields: 2,192 frequencies still
  trigger, almost all in the same depths 67-70 inversion column,
  and its absolute max (8.7e-2) approaches V12's unrescued max.
  The mechanism is strongly suppressed, not abolished.
- The KNOP lagged I_NEXT field is NOT converged after one
  inter-sweep refresh (residual 0.98 vs tol 2e-3; sweep-to-sweep
  drift med 18%, max 7.8x). KNOP's magnitudes are therefore an
  immature iterate; the direction of the correction (down) is
  robust, the absolute level is not established.
- The emergent CMF spectrum differs materially: KNOP flux ~0.2-3%
  of the V12-rescue flux in the flux-carrying region
  (med C/A = 0.002-0.024 depending on mask). This is the size of
  the correction the J artifact was injecting through the
  albedo-coupled source function.

Recommendation (section 13): the comparison supports the Phase 2.3
hypothesis that V12's clamped-sigma formal solve is invalid in the
velocity-inversion region and that its ~1e-3 boundary fields are an
artifact, not a diffusion-limit asymptote. KNOP is the correct
direction but cannot yet be adopted as a production replacement:
its Gauss-Seidel wavelength coupling does not converge within
NUM_ES=2. Until a mature KNOP render exists, retain
RESC_JMAX=1.0e-3 and the restore, and treat the resulting spectrum
as contaminated at the ~percent-to-tens-of-percent level in the
flux-carrying region.

## 2. Observer dispatch audit (Task 1)

`DO_NONMONO_KNOP` is not a formal-solver selector.  It appears at
~57 source sites in `obs/` and ~210 more across 24 population-solver
files; setting it re-routes the whole pipeline, not just the J call.

Observer-side `DO_NONMONO_KNOP` sites (`obs/cmf_flux_sub_v5.f`,
`obs/COMP_JCONT_V4.INC`, `obs/rd_cmf_flux_controls.f`,
`obs/knop_spec_gate.f`, `obs/cmf_flux_cntrl_var_mod.f`):

| site | role | needed by isolated selector? |
|------|------|------------------------------|
| sub:658-678 | nonmono-RVTJ admission gate (ALLOW_NONMONO_OBS overrides) | no - model already passes via DO_NONMONO_SOB/ALLOW_NONMONO_OBS |
| sub:718 | SIGMA_TRUE allocate/fill | shared gate already active via DO_NONMONO_SOB=T; OR'd for safety |
| sub:791 | VINF=MAXVAL(V) for grid refinement | shared gate already active via SOB |
| sub:1546 | luminosity-DTDR boundary correction (J_RESCUE-combined gate) | already active via J_RESCUE |
| sub:1605 | KNOP_SPEC_GATE population-provenance check | NO - stays KNOP-only (correctly: SOB populations were not KNOP-converged) |
| sub:1644 | ES_RESID_SEED arg | no - diagnostic seeding, SOB-semantics retained |
| sub:1749 | SET_FG_GS_OBS_GRID(NU,NCF) registration | YES - solver's lagged-field remap needs the OBS grid |
| sub:2361 | FG_KNOP_REFRESH_V1 between non-final ES sweeps | YES - lagged I_NEXT matures here (NUM_ES=2, one refresh) |
| sub:2481 | GET_FG_GS_MATURITY census | YES - diagnostics only |
| sub:2811/3241 | OBS_FRAME_SOBSUB_V* flag arg `(SOB or KNOP).AND.DIF` | shared gate already active via SOB |
| INC:404 | extended-grid KNOP guard (INC_GRID) | hardened to include OBS_FG_KNOP |
| INC:925-959 | THE formal-solver call | THE selector point |
| INC:117 | negative-DBB floor (SOB-or-KNOP-or-J_RESCUE gate) | already active via SOB/J_RESCUE |
| INC:623 | extended-density quick path exclusion | KNOP-only by design, unreachable (INC_GRID=F) |
| INC:786 | FG_SOL_NM naming (pre-dispatch) | cosmetic; actual call site re-names correctly |
| rd:521-796 | KNOP key-dependent defaults/enforcements | left untouched - OBS_FG_KNOP is not meant to inherit them |

Population-solver side (NOT entered by this diagnostic): 24 files,
~210 sites - solve_for_pops.f (47), rd_control_variables.f (47),
cmfgen_sub.f (27), comp_j_blank.f (17), plus RVTJ readers, NG
acceleration, wind-velocity, S_tilde/grey paths, JFNK hooks.
Setting DO_NONMONO_KNOP=T in the OBS param file alone cannot reach
these (different control file) but documents why reusing the flag
for the observer would mislabel provenance.

## 3. Isolated selector (Task 2)

New diagnostic-only keys, both read unconditionally with defaults
equal to production behaviour:

- `OBS_FG_KNOP` (default F): when T, `COMP_JCONT_V4.INC:925`
  dispatches the per-frequency formal solve to
  `FG_J_CMF_NONMONO_V1` with the OBS-built `SIGMA_TRUE`, and the
  three solver-machinery hooks arm (grid registration, refresh,
  maturity census).  Everything else stays on the production
  `DO_NONMONO_SOB` gates.  INC_GRID=T + OBS_FG_KNOP aborts loudly
  (same guard as KNOP).
- `J_RESCUE_RESTORE` (default T): when F, the rescue trigger,
  event counting, throttled prints, JEV/JDMP/JFRQ records all
  still fire, but the seven `*_PREV` restore copies and the
  snapshot-validity abort are skipped; the unrescued field reaches
  LU_EDD and chains into JNU_PREV naturally.  J_RESCUE stays T;
  the DBB/DTDR boundary corrections it also gates are untouched.

Call-argument differences (INC:944 vs :953): KNOP receives
`SIGMA_TRUE` and `NCF` instead of `DIF`; `SET_FG_NU_NEXT(NU(k+1))`
supplies the forward log-spacing compressing cells need.  TA,
CHI_CLUMP, CHI_SCAT_CLUMP, V, SIGMA, R, P, TC, boundary args,
quadrature weights, METHOD, THK_CONT are identical.

## 4. Restore/boundary separation (Task 3)

J_RESCUE gates three operations: (a) the prev-frequency restore
copies, (b) the negative-DBB floor (INC:117), (c) the luminosity-
derived DTDR replacement (sub:1546).  The bypass touches ONLY (a):
the `IF(J_RESCUE_RESTORE)` wrappers cover the snapshot-validity
abort and the seven copies; the trigger evaluation, NRESC++, JEV,
JDMP and the armed-but-inactive boundary corrections are outside
the wrappers.  Finite-field validation stays armed on every
frequency (a genuine non-finite field is still FATAL'd upstream of
the restore site, not masked by the bypass).

## 5. Input equivalence (Task 4)

All four cases run the Phase 2.3 frozen inputs (hashes identical
across manifests): VADAT `7a3ac0c8...`, RVTJ `2bb5214c...`,
MODEL `499efc9d...`, MODEL_SPEC `d6fd4c41...`.  PINITs differ only
in the two new keys:

| case | solver | restore | PINIT tail |
|------|--------|---------|------------|
| A | V12 | on | (diag init as-is) |
| B | V12 | bypassed | `F [J_RESCUE_RESTORE]` |
| C | KNOP | bypassed | `T [OBS_FG_KNOP]`, `F [J_RESCUE_RESTORE]` |
| D | KNOP | on | `T [OBS_FG_KNOP]` |

Render matrix (definitive, exe `eacd7bc0525a1195`): job 13340788.
Preview matrix (superseded exe `60428671b9551a4d`): job 13340484 -
caseA aborted on the `L_TRUE` reader bug (key required, absent in
the diag init); fixed to `L_FALSE`.  Preview caseB/C/D retain the
keys explicitly so their physics is unaffected; they are kept as
auxiliary evidence only.

## 6. Radial radiation fields (Task 5)

Source: `J_RESCUE_DUMP` cadence records (JSCN every 500th
frequency, 150 dumps/case on identical grids) and the per-frequency
`J_RESCUE_FREQ` records (74,620 aligned continuum frequencies;
`jfrq_compare.csv`).

At optically thick frequencies the two solvers agree
bit-for-bit: at JSCN idx 500 (chi~40 at depth 70) every depth-70
value matches to the printed precision (J=zeta=6.9e-10, the
thermalization limit). This is the control: where V12 is valid,
KNOP returns the same field.

In the inversion column (depths 59-70, sigma_true = -31 to -38):

| JSCN idx | depth 70 J (A) | depth 70 J (C) | A/C | zeta(70) | J/zeta (A) | J/zeta (C) |
|----------|----------------|----------------|-----|----------|------------|------------|
| 20,000 | 9.9e-4 (at cap) | 2.0e-4 | 4.9 | 8.5e-7 | 1165 | 237 |
| 50,000 | 7.7e-4 | 1.7e-5 | 46 | 2.3e-6 | 334 | 72 |

Profile shapes (`figs/phase24_radial.png`): V12's |J| climbs
monotonically inward toward the cap, decoupled from the zeta
profile. KNOP's |J| follows the zeta shape (dip at depth ~61-62,
rise to the boundary) and stays 5-50x below V12's. V12-unrescued
(B) has the same shape as A but 10-100x higher (the chain
propagates the unrescued field). Argmax sits at depth 70 on ~99%
of frequencies for all cases.

## 7. Frequency dependence

`figs/phase24_rjmax.png` left panel: B sits above the cap nearly
everywhere (chain amplification), A is the dense band pinned at
the cap from idx ~3,000 to ~60,000, C drops below the A band from
idx ~5,000 upward, ending ~100x lower at the red end. The three
agree at idx < ~3,000 (pre-inversion onset in frequency order)
and re-merge at the far red tail.

Right panel (mean J/zeta per 500-freq bin): B sits at 1e3-1e5
between idx 10k-50k - far above the 435 scattering ceiling, so
the unrescued chain fields are unambiguously unphysical. A hovers
at ~100-300 (below the ceiling but inflated vs zeta). C sits at
~10-100, well under the ceiling.

Per-frequency ratio statistics (mature sweep): A/C med 5.42
(p95 62.5, p99 79.2, max 103.6); B/C med 60.5 (p95 204, max
1312).

## 8. Convergence analysis (Task 6)

- V12: F_DID_NOT_CONV=0 in the JFRQ records of A and B (the JEV
  event stream carried 36% non-converged flags in Phase 2.3; the
  JFRQ flag is sampled at a later call point and reads F
  throughout - record limitation, not evidence of convergence).
- KNOP sweep 1: 1,256 fdnc=T records; sweep 2: 0 - the refresh
  pass cleared the flagged non-converged solves.
- KNOP GS maturity (OUT_FLUX census): refresh residual 0.977 vs
  tol 2.0e-3, 49.1M unconverged cells of 728M (6.7%), sweep max
  rel change 1.0e4, 179M cells >50% change. The lagged I_NEXT
  field did not converge within NUM_ES=2.
- JFRQ sweep1->sweep2 drift: median 17.8%, p95 135%, max 7.8x.
  KNOP's reported fields are a mid-iteration estimate.
- J_RESCUE under KNOP: 2,192 mature-sweep triggers (4,529 JEV
  lines over both sweeps; 307 freqs once, 2,111 twice), all
  trig=3 (the |RJ| cap). The depth histogram is the same
  inversion-column signature as V12: 2,172 of 2,192 at depths
  67-70 (481 at depth 70), 20 scattered at depths 22-59. The
  mechanism persists at ~14x lower incidence - consistent with
  an immature I_NEXT occasionally overshooting, though a genuine
  residual cannot be excluded.

## 9. Scattering-source decomposition (Task 7)

Source structure (obs/COMP_JCONT_V4.INC): `SOURCE = ZETA +
CHI_SCAT*RJ/CHI = ZETA + THETA*RJ`, `ZETA=ETA/CHI`,
`ESEC = 6.65e-15*ED` (Thomson on electron density). At depth 70:
chi(es)=0.5131, chi(total)=0.5143, omega=0.9977 -> diffusion
ceiling J* ~ zeta/(1-omega) ~ 435*zeta.

Measured against this ceiling: V12's J/zeta med 289 (p99 2346,
max 5241) - the bulk of A sits below the ceiling, but ~2.4x of
its top tail exceeds it and its boundary field is pinned by the
cap, not by physics. B's unrescued chain (med 850, p95 3.9e4)
blows through the ceiling. KNOP's med 9.6, p95 343: inside the
scattering asymptote at essentially all frequencies.

This revises the Phase 2.3 open question: the ~1e-3 boundary J is
NOT a clean diffusion-limit asymptote. A correct scattering-limit
field would sit ~435x zeta and track zeta's profile; V12's field
climbs to the cap regardless of zeta, while KNOP's lower field
does track zeta.

## 10. Observer spectra (Task 8)

`OBSFLUX` = the CMF flux spectrum on the 45,541-point continuum
grid (the OBS_FRAME step aborts identically in all cases on the
pre-existing GLOBAL_CMF_GRID_UNDERCOVERAGE; equal footing, the
same abort the Phase 2.3 baseline renders ended with).

Spectral-mass metric = sum|F_case - F_A| / sum|F_A| over aligned
grids (the Phase 2.3 convention), plus ratio statistics under a
defined flux mask:

| pair | spectral-mass | med ratio (all pts) | med ratio (A>1e-3 max) | p95 rel |
|------|---------------|---------------------|------------------------|---------|
| B/A | 3307% | 89 | 2.08 | 425 |
| C/A | 91% | 0.024 | 0.002 | 2.26 |
| D/A | 78% | 0.99 | 0.002 | 3.90 |
| D/C | 31% | 1.00 | 1.00 | - |

Luminosity record (the 70-value cumulative tail) is identical in
all four cases (final 2.6499e9) - the model luminosity is
conserved; the differences are pure spectral redistribution.

Case D (KNOP + production restore) completes the matrix: its
per-frequency field is identical to C at the median (D/C rjmax
ratio med 1.000, p95 1.00 - the restore only touches the 1,752
fired frequencies), its max sits exactly at the cap
(9.998e-4), and its spectrum is 31% mass-closer to A than C's.
Notably D's p99 J/zeta = 435.6 - the restore pin lands almost
exactly on the electron-scattering ceiling, which is a
coincidence of this model's depth-70 albedo, not a designed
property of the cap.

The KNOP spectrum is ~40-500x BELOW the V12-rescue spectrum in
the flux-carrying region: the rescued V12 field, amplified by the
omega~0.998 albedo coupling into the emissivity, was injecting
tens of percent to tens of times excess flux. The unrescued V12
spectrum is distorted ~89x at the median - showing what the
restore suppresses.

`obsflux_spectra.npz` holds the aligned arrays; `jfrq_compare.csv`
the per-frequency merged records.

## 11. KNOP implementation audit (Task 9)

## 11. KNOP implementation audit (Task 9) - pre-render findings

Method actually used: short-characteristics INTEGRAL solve
(`SOLUTION_METHOD='INTEGRAL'`; the DIFFERENCE path is STOP-gated
dead code).  Structure vs Knop et al. 2009 (A&A 496, 295; A&A 501,
813):

- Upwind direction: per-point sign of the frequency-coupling
  coefficient `GAM=(V/cR)(1+sigma*mu^2)`; GAM>0 couples to the
  previous (bluer) frequency, GAM<0 to the next (redder) -
  Paper I Eq 3 / Paper II Eqs 25-28.  FREQ must strictly decrease
  per call (hard check).
- Generalized opacity `chi_hat = chi + |GAM|/DLN` is the
  positive-only velocity term of Paper I Eq 7; the upwind
  intensity enters as the explicit source `S' = eta/chi_hat +
  Q*I_upwind`, Q signed for the expanding side and ABS(Q) for the
  compressing side (documented signed-Q fix).
- Lagged Gauss-Seidel store `I_*_DEPTH_NEXT` on a sparse
  compressing-depth set (SIGMA_TRUE<-1 or V<0, +1 neighbor),
  remapped onto the OBS frequency grid by SET_FG_GS_OBS_GRID;
  refresh via FG_KNOP_REFRESH_V1 between ES sweeps (Paper II
  Sec 3.2 ordering, adaptive INEXT_RELAX).
- Source interpolation: monotone-cubic chi/eta in ln R; formal
  integral uses parabolic source + slope clamps (Olson & Kunasz
  coefficients E0-E3), linear-explicit side per Paper I Eqs 11-12.
- chi_hat stays strictly positive by construction; DTAU via
  Euler-Maclaurin corrected trapezoid (bounded between endpoint
  estimates), negative/NaN DTAU falls back to plain trapezoid.
- The lambda-form `4*a*I` term is NOT carried: the code solves
  the nu-form where the corresponding invariant-intensity piece
  is argued to cancel identically (dev audit comment); residual
  piece = the implicit upwind part absorbed into chi_hat.  This
  is the one structural deviation from the published equations -
  flagged for the correctness discussion, not treated as an error.
- Clamps present: intensities floored at 1e-30, INEXT store and
  ray intensities capped at NONMONO_INEXT_CAP=1e20 (VADAT),
  IN_HBC validated to [0,1], sigma floor -0.999 only inside the
  THK beta-extension (main grid uses unclamped SIGMA_TRUE).
- Boundary conditions: outer I_M=0 (no incident), inner supports
  DIFFUSION/ZERO_FLUX/PNT_SRCE/SCHUSTER with the same DBB/IC
  semantics as the V12 call.
- Moments: J/H/K/N accumulated with the caller's JQW/HQW/KQW/NQW
  quadrature; N_ON_J/RSQHNU assembled downstream identically for
  both solvers.

Post-render answers to the pre-render questions:

- GS maturity does NOT converge with NUM_ES=2: refresh residual
  0.977 vs tol 2e-3, 6.7% of cells unconverged, sweep drift to
  1e4. One refresh is not enough; the OBS-path I_NEXT needs more
  ES sweeps or a dedicated outer iteration to mature.
- The cold-start bias is visible: sweep-1 fields carry 1,256
  fdnc=T records (cleared in sweep 2) and sweep-to-sweep JFRQ
  drift is med 18%. KNOP magnitudes are mid-iteration values.
- Rescue events persist under KNOP (2,192 mature triggers, 481
  at depth 70): the boundary-value solver still produces
  occasional over-cap fields - consistent with a not-yet-mature
  lagged field, but a genuine residual cannot be excluded.

## 12. Remaining uncertainties

1. KNOP maturity is the dominant uncertainty. Its median field
   could still move (up or down) as I_NEXT converges; the med-18%
   inter-sweep drift bounds the current imprecision. A NUM_ES>=4
   run or a standalone I_NEXT iteration is needed for a mature
   value.
2. The caseC J_RESCUE_EVENTS file was shared briefly with a
   superseded preview job (killed mid-write). The file is
   internally consistent (no byte-identical duplicates; the
   307-once / 2,111-twice sweep structure matches the JFRQ counts
   exactly), and the JFRQ/JEV agreement corroborates it; the
   caveat is recorded here regardless.
3. The JFRQ record's fdnc column under-reports V12
   non-convergence (sampling point is after the JEV check); the
   Phase 2.3 JEV stream is the authoritative V12 fdnc source.
   Fixed in source (format + sampling) for future renders.
4. OBSFLUX is the CMF flux spectrum, not the observer-frame
   spectrum: all cases (and the Phase 2.3 baseline) terminate on
   the pre-existing GLOBAL_CMF_GRID_UNDERCOVERAGE abort before
   the observer-frame product is written. The comparison is at
   equal footing; an absolute observer-frame spectrum was not
   produced by any arm.
5. J/zeta interpretation assumes zeta=ETA/CHI is the true local
   source scale; at line frequencies the emissivity structure is
   more complex. The optically-thick control (JSCN 500: A=C
   exactly) bounds this concern.
6. The `S_tilde` (lambda-form 4aI term) omission noted in the
   code audit is a deliberate documented deviation, not verified
   against a reference implementation.

## 13. Recommendation for production use

1. Keep DO_NONMONO_SOB + USE_LAM_ES=T as the production route.
   NOTE (22-Nov-2026 addendum): the RESC_JMAX=1.0e-3 flat rail
   measured here was recalibrated to 1e2 after the m15
   interacting-model A/B proved it clipped verified-legitimate
   fields (~0.12) and collapsed that model's continuum ~100x.
   The flat rail now bounds only true runaway (>=1e2); the
   1e4*JNU_PREV relative bound is the working detector. This
   means the tauto-model boundary artifact quantified in this
   audit is NO LONGER clipped by the restore - its renders
   revert to pre-rescue stock fields. If the restore's
   ~89x-inflation guard is still wanted on that model, the
   per-cell J/J_eq bound (threshold-attribution audit sec. 6a)
   is the model-independent criterion to implement; no flat
   value separates its artifact (~9e-2) from m15's legitimate
   fields (~0.12).
2. Accept the measurement that the V12-rescue spectrum carries a
   solver artifact: the independent solver's field is ~5-60x
   lower at the same points, and the artifact's spectral imprint
   is tens of percent (mass) with ~40x local flux errors in the
   flux-carrying region. Report this as a systematic uncertainty
   on m15_2p82d-class non-monotonic models.
3. KNOP (OBS_FG_KNOP) is now a working diagnostic selector for
   the observer frame: physically structured, paper-faithful, and
   ~17x cleaner at the boundary. It is NOT production-ready: its
   GS-coupled field needs a convergence program (more ES sweeps
   or a dedicated I_NEXT outer loop) before its absolute levels
   can be quoted.
4. Do not treat KNOP's smaller J as automatically correct - its
   p95 sits below the scattering ceiling, consistent with an
   immature iterate as easily as with the true boundary-value
   solution. The direction of the correction is the robust
   result; the magnitude is pending.
5. Next experiment: NUM_ES=4 (or a KNOP-observer restart chain)
   to mature I_NEXT, then recompute this matrix. The JFRQ/JSCN
   records and the comparison scripts (obs_ensemble/
   phase24_analyze.py, phase24_figs.py, jfrq_lib.py) are reusable
   verbatim.

## Appendix - diagnostic flag record

- OBS_FG_KNOP: T in cases C,D only; F elsewhere (default).
- J_RESCUE_RESTORE: F in cases B,C; T in A,D (default).
- Both read with MUST_BE_PRES=L_FALSE (optional keys; absent =
  production default).
- RESC_JMAX=1.0E-3, J_RESCUE=T, DO_NONMONO_SOB=T, USE_LAM_ES=T,
  ALLOW_NONMONO_OBS=T, INC_GRID=F, NUM_ES=2 in all cases.
- Machine-readable: phase24_summary.json, jfrq_compare.csv
  (74,620 rows), obsflux_spectra.npz, per-case manifest.json,
  jfrq_lib.py + phase24_analyze.py + phase24_figs.py (all in
  obs_ensemble/).
- Case D (KNOP+restore) confirms the restore is solver-agnostic
  plumbing: it pins the same over-cap fields at 9.998e-4 and
  leaves the rest of the field untouched (D/C field ratio med
  1.000). Under KNOP only 1,752 frequencies need it - the rescue
  machinery is ~18x less busy with the correct solver.
- JFRQ record caveat: the shipped record has the field-misaligned
  format fixed post-matrix (commit 99fcc18e); all four case files
  use the shifted layout documented in phase24_analyze.py /
  jfrq_lib.py.
