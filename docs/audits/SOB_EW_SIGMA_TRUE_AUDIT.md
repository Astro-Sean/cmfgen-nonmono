# OBS Sobolev EW diagnostic: overflow on non-monotonic fields

- Date: 2026-11-13
- Tree: `cmfgen_jfnk`, branch `main`, HEAD at change: `eed3f328`
- Model: `m15_2p82d_FRESH_sob_dieoff`, epoch 2.82408d, run dir
  `/ptmp/sbrennan/SN2024pba/cmfgen_validation/m15_2p82d_FRESH_sob_dieoff/cmfgen/2.82days`
- Populations: `DO_NONMONO_SOB=T`, `GLOBAL_LINE='SOB'` (job 13084146,
  converged iter 756, warm restart). OBS: `GLOBAL_LINE='BLANK'`,
  `DO_SOB_LINES=T` + `WR_SOB_FORCE=T` enabled for the diagnostic runs
  only, then restored.
- Scope: `EWDATA`, `SOB_FORCE_DATA`, `SOB_FORCE_MULT` only. OBSFLUX,
  `obs_fin_10`, populations and the formal solver are not modified.

## 1. Failure mechanism

Production `EWDATA` carried 113 entries with |EW| > 1e30 (max
5.6e99; an earlier iteration of the same model reached 3.8e169 for
CIV(3p-3s) 5801). `SOB_FORCE_DATA` carried 100,666 values > 1e30
(max 3.0e198). Three ingredients combine:

1. **Real inversions.** C IV 3s -> 3p is pumped through the
   resonance doublet; `CIV_IN` shows 3p/3s inverted at 63/70 depths
   (up to ~100x the stimulated-emission threshold). `CHIL < 0` is a
   physically valid population state, not a defect.
2. **Weak-negative opacity passes the guard.** `CHECK_LINE_OPAC`
   crushes `CHIL` only for `T2 = CHIL*3.0D-10*R/V/FL < -0.5`.
   `NEG_OPAC` logs crushed C IV Rydberg lines (6d-6p, 7d-7p,
   10d-10p) but not 3p-3s, which stayed in `-0.5 <= T2 < 0`.
3. **Clamped SIGMA makes the angular denominator pathological.**
   `SOBEW_GRAD_V2` forms
   `T1 = 1 + mu^2*SIGMA`, `GAM = GAMH/T1`. Under `DO_NONMONO_SOB`
   the RVTJ `SIGMA` column is the work gradient re-clamped to
   >= -0.999 (`RECOMPUTE_SIGMA_TRUE`); in this model **all 70
   depths store -0.999**. For a radial ray, `T1 = max(0.001,
   -0.2*(1+sigma), 1e-4) = 0.001`, so `|GAM| = 1000*|GAMH|`.
   `EXPONX(GAM) = (1-exp(-GAM))/GAM` then gives unsaturated
   exponential amplification: GAM ~ -180 -> e^180 ~ 1e78.
   No single element is buggy; the EW estimator was fed a velocity
   gradient the model never had.

`EXPONX` verification (unchanged, matches `(1-e^-X)/X`):

| X | EXPONX(X) |
|---|---|
| -390 | 6.1e166 |
| -380 | 2.8e162 |
| -190 | 1.7e80 |
| -50 | 1.0e20 |
| 5 | 0.1987 |
| 5000 | 2.0e-4 |

IEEE double overflow of `e^|X|` needs |X| > ~709; the reachable
GAM here (~500 max) stays finite but astronomical. Positive-GAM
lines are also wrong: `EXPONX(5000)=2e-4` vs `EXPONX(5)=0.2`
suppresses saturated contributions ~1000x.

## 2. SIGMA vs SIGMA_TRUE and the projected gradient

- `SIGMA = dlnV/dlnR - 1` stored in RVTJ is the work gradient,
  re-clamped >= -0.999 under `DO_NONMONO_SOB` so stock routines
  that assume sigma >= -1 stay valid.
- `SIGMA_TRUE` (`SIGMA_TRUE_CORE`, `new_main/subs/
  recompute_sigma_true.f`) is the same signed quantity computed on
  the same ND grid from the unclamped R,V: log-space quadratic
  interpolant in ln R / ln|V|, exact for power laws, linear
  fallback near sign changes, undefined cells fail closed unless
  `ALLOW_SIGMA_UNDEF=T`. It is the gradient the population solver
  and the Knop path already use.
- The Sobolev line-of-sight velocity gradient along a ray of
  impact parameter p is `(V/R)*(1+sigma*mu^2)`, mu = cos(theta) =
  Z/R (standard result; Lamers et al. 1987 as cited in
  `sobjbar_sim.f`). `SOBEW_GRAD_V2` correctly regularizes
  `T1 = max(|1+sigma*mu^2|, -0.2*(1+sigma), 1e-4)` for sigma < 0
  (Dessart, Audit & Hillier 2015 Sec 3.1 footnote 6; the 1e-4 arm
  is a code-only extension, see
  `SOBOLEV_DAH2015_EQUATION_MAP.md`). The mathematics is correct;
  the input array was wrong for this model class.
- Recomputed sigma_true for this model: range -130.2 to -1.037
  (worst compression depth 62). Minimum |1+sigma_true| over depths
  is 0.037. The clamped path therefore inflated |GAM| by a factor
  37 to ~129,000 depending on depth, before angular quadrature.

## 3. The change (`obs/cmf_flux_sub_v5.f`)

Two hunks, both in the OBS driver; `SOBEW_GRAD_V2` itself is
untouched:

1. The `SIGMA_TRUE` derivation block now runs under
   `DO_NONMONO_SOB .OR. DO_NONMONO_KNOP .OR.
   (DO_SOBOLEV_LINES .AND. N_NONMONO_OBS .GT. 0)`. The fail-closed
   abort on undefined cells stays gated on the two nonmono flags;
   the EW-only extension warns instead (a stock model never enters
   the block).
2. The EW call selector changes from `ELSE IF(DO_NONMONO_SOB)` to
   `ELSE IF(ALLOCATED(SIGMA_TRUE))`, passing `SIGMA_TRUE` to
   `SOBEW_GRAD_V2` whenever it has been derived.

Why this shape:

- `DO_NONMONO_SOB=T` correctly requires `GLOBAL_LINE=SOB`
  (`check_cmf_flux_param_consis.f`), but production spectra need
  `GLOBAL_LINE='BLANK'` for full blanketing, so the flag could
  never gate the fix for this model.
- `N_NONMONO_OBS` is the existing field detector (V-reversal and
  stored sigma < -1 cells). Any sigma_true < -1 implies V rises
  inward at that cell, so the V test catches every compressing
  cell for V > 0.
- On a monotonic field `N_NONMONO_OBS=0`, `SIGMA_TRUE` stays
  unallocated, and the call is byte-identical to before - stock
  monotonic models are unaffected by construction.
- Under `DO_NONMONO_KNOP` the stored SIGMA is not clamped, so the
  selector change is a no-op there (same values, same call).
- Nothing is bypassed: `CHECK_LINE_OPAC`, `EXPONX`, CHIL, the
  populations, the formal solver and the ordered frequency sweep
  are all untouched.

## 4. Same-input A/B evidence

Binaries: `exe/cmf_flux.exe` md5 `da38547a82edace90752ecc0b569ab43`
(new, 13-Nov 07:15) vs `/tmp/cmf_flux_oldsel.exe` (same object list,
`cmf_flux_sub_v5.o` compiled from `git show HEAD:obs/cmf_flux_sub_v5.f`).
Identical inputs: iter-756 RVTJ, `DO_SOB_LINES=T`,
`WR_SOB_FORCE=T`, `DO_NONMONO_SOB=F`, `ALLOW_NONMONO_OBS=T`,
OMP_NUM_THREADS=4. Both ran `batobs.sh`; detector reported
`non-monotonic velocity signature in RVTJ (69 flagged cells)` in
both runs. Saved under `2.82days/sobfix_check/` as `.oldsel` /
`.sigtrue`.

| output | old selector (clamped SIGMA) | new (SIGMA_TRUE) |
|---|---|---|
| OBSFLUX | identical | identical |
| obs_fin_10 | identical | identical |
| OUT_FLUX | identical | identical |
| EWDATA rows | 10,935 | 11,410 (EW_CUT_OFF now admits lines the clamped path suppressed) |
| EWDATA \|EW\|>1e30 | 113 | 0 |
| EWDATA max \|EW\| | 5.6e99 | 1877 |
| EWDATA NaN/Inf | 0 / 0 | 0 / 0 |
| EWDATA negative EW | 3,926 | 5,164 (emission-convention; more rows admitted) |
| SOB_FORCE_DATA >1e30 | 100,666 | 0 |
| SOB_FORCE_DATA max | 3.0e198 | 200 |
| SOB_FORCE_DATA NaN/Inf | 0 / 0 | 0 / 0 |
| SOB_FORCE_MULT | per-depth values to -3e198 | finite (max ~90) |

C IV 3p-3s, same model state:

| line | old EW | new EW | old LFlux | new LFlux |
|---|---|---|---|---|
| 5802.92 (3p[3/2]-3s) | (under -0.5..0 guard window) | 53.1 | - | 3.3e-5 |
| 5813.58 (3p[1/2]-3s) | 6.1e82 | 26.1 | 3.8e76 | 1.6e-5 |

`obs_fin_10.sigtrue` is also bit-identical to the production OBS
job 13084147 output (`obs_fin_10_20261113`) run on the same RVTJ
with `DO_SOB_LINES=F` - the EW block is a terminal diagnostic, and
the spectrum is unchanged to machine precision.

Diagnostic-only proof: `FORCE_MULT` is a local array zeroed at the
top of the EW block, written direct-access to `SOB_FORCE_DATA`
(unit 82) and summarized to `SOB_FORCE_MULT`; EW goes to `EWDATA`
via `WRITE_OUT_EW`. Neither file is read by the population solver,
the formal transfer, the radiation-field update, or the
statistical-equilibrium equations. The whole block sits after
`IF(.NOT. DO_SOBOLEV_LINES)RETURN`.

## 5. Regression tests

- `unit_tests/python/test_sobew_short_ray.py`: 41/41 PASS
  (NI=1, NI=2, NI>2, multiple sigma, negative GAMH) - unchanged,
  the fix is above the routine.
- `unit_tests/python/test_sob_tau_floor.py`: 17/17 PASS
  (sigma = +2 identical, sigma = -0.5 identical, sigma <= -1 floor
  cases) - unchanged.
- Monotonic equivalence (Test A): by construction - a monotonic
  field never allocates `SIGMA_TRUE`, so the same object code takes
  the same `SIGMA` call. No monotonic production model exists in
  `cmfgen_validation` to re-verify end-to-end; the claim rests on
  the unreachable branch, the untouched `SIGMA` call, and the two
  untouched routines' unit suites.
- Near-zero projected gradient (Test D): sigma_true ~ -1 cells are
  the residual risk, not sigma = -0.999. Minimum floored T1 here is
  0.037 (DAH15 floor inactive since |1+sigma| > -0.2*(1+sigma));
  GAM stays bounded and the output is finite everywhere.
- Negative GAM / EXPONX (unchanged): verified numerically above;
  inversions are retained (negative CHIL is not touched).

## 6. Interpretation and limitations

- Negative Sobolev optical depth gives exponential amplification
  in the *unsaturated* formal solution (Goldreich & Kwan 1974;
  Gray 2012). A formal EW of 1e80 is not evidence of a physically
  saturated maser calculation - the estimator carries no
  saturation physics, and none was added.
- A near-zero projected gradient `|1+sigma*mu^2| -> 0` means a
  large Sobolev escape length. Where that length is not small
  against the relevant physical scale lengths (density,
  population-gradient scales), the local Sobolev estimate itself is
  questionable. Near a velocity extremum the true sigma < -1
  region still yields the smallest T1 (0.037 here); those cells
  dominate the corrected EW of inverted lines. The corrected EW is
  a diagnostic estimate, not a substitute for non-local transfer of
  a genuinely masing transition (Rybicki & Hummer 1978 non-local
  resonance treatment is the relevant machinery where needed).
- Negative EWs remain (emission convention) and inverted lines
  still report elevated |EW| - the inversion is preserved, only
  the wrong denominator is removed.
- `disp/maingen.f` has its own sigma-based EW display path and
  still reads the stored column; it is outside the OBS EW
  diagnostic and was intentionally not touched.
- `SOB_SIGMA_MAX` (default off) applies to `SIGMA_TRUE` on the new
  path exactly as under the nonmono flags - same treatment, same
  opt-in knob.

## 7. References

- Lamers, Cerruti-Sola & Perinotto 1987, ApJ 314, 726 - projected
  velocity-gradient factor `1+sigma*mu^2` used in the tau
  denominator.
- Dessart, Audit & Hillier 2015, MNRAS 449, 4304 (DOI
  10.1093/mnras/stv609), Sec 3.1 + footnote 6 - tau floor
  `max(|1+sigma*mu^2|, -0.2*(1+sigma))`; the 1e-4 arm is
  CMFGEN-specific and not from the paper.
- Rybicki & Hummer 1978, ApJ 219, 654 - non-local resonance
  coupling (relevant to genuinely masing lines).
- Gray 2012, "Maser Sources in Astrophysics", CUP - unsaturated
  exponential vs saturated maser regimes.
- Goldreich & Kwan 1974, ApJ 190, 27 - maser saturation physics.
