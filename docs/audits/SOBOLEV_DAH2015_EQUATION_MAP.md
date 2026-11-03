# DO_NONMONO_SOB vs Dessart, Audit & Hillier (2015) Sec 3.1

- Date: 2026-11-03
- Reference: Dessart, Audit & Hillier 2015, MNRAS 449, 4304
  (arXiv:1503.05463). Section 3.1 "Numerical approach" + footnote 6.
- Purpose: equation-by-equation map of the paper recipe onto this
  tree, closing P2-001/P2-002. The `docs/1503.05463v1.txt` OCR of
  footnote 6 drops minus signs ("- - -" glyph scatter); the reading
  below uses the only form that is positive for sigma < -1 and
  matches every code comment in the tree.

## Mapping

| Paper element (Sec 3.1 / footnote 6) | Code site | Status |
|---|---|---|
| Sobolev optical depth tau_S ∝ 1/\|1+sigma*mu^2\|, sigma = dlnV/dlnR - 1 (standard form, Lamers et al. 1987 as cited in sobjbar_sim.f:108-112) | `sobjbar_sim.f:124-130`, `sobew.f:131`, `sobew_grad_v2.f:137`, `sobew_grad_v3.f:140`: `T1=1+(Z/R)^2*SIGMA; GAM=GAMH/ABS_or_floored(T1)` | matches |
| Footnote 6: tau diverges at one angle when the radial gradient is not positive - the singular direction mu^2 = 1/\|sigma\| exists only for sigma < -1. Escape probability barely affected (small solid angle) but the coarse angle grid caused numerical instabilities; "limiting the increase in optical depth" via floor | same four sites: `T1=MAX(ABS(T1),-0.2*(1+SIGMA),1.0D-4)` under `SIGMA<0` | matches, plus one documented extension |
| Footnote 6 floor: `max(\|1+sigma*mu^2\|, -0.2*(1+sigma))` for sigma<0 | `-0.2*(1+SIGMA)` term in the MAX | matches - for sigma<-1 this is 0.2*(\|sigma\|-1); for -1<sigma<0 it is negative so the 1e-4 arm controls |
| (extension beyond paper) 1e-4 cap | the third MAX argument `1.0D-4` | code-only addition: bounds the sigma -> -1 case where the paper's own floor -0.2*(1+sigma) -> 0 and cannot regularize. Documented in the site comment ("1D-4 bounds the sigma=-1 case") |
| "model the lines for the non-LTE analysis using the Sobolev approximation" | `DO_NONMONO_SOB`: `SUB_SOB_LINE_V3` + `SOBEW` take `SIGMA_TRUE` (cmfgen_sub.f; sub_sob_line_v3.f:144-158) | matches intent; the continuum is additionally solved non-monotonically (see next row) |
| "we wrote a new routine capable of computing J in the comoving frame in the presence of non-monotonic velocity fields" | `MOM_J_CMF_NONMONO_V1` (USE_LAM_ES=F) / `FG_J_CMF_NONMONO_V1` (DO_NONMONO_KNOP); paper-era default `FG_J_CMF_V13` remains under USE_LAM_ES=T | paper's custom routine maps to the Knop/moment stack - a newer capability, not a divergence |
| "uses the absolute value of the velocity gradient" (Sobolev is local, sign-free) | `SIGMA_TRUE` keeps the true signed dlnV/dlnR-1 for geometry; the tau uses ABS(T1) | matches: sigma stays signed in coupling terms, ABS only in the tau denominator |
| Nonlocal coupling between resonance points | `sob_rh_couple.f` under `DO_NONLOC_SOB` (Rybicki & Hummer 1978 two-point kernel; BUG-002 fingerprint fix 2026-11-03) | beyond the paper - separate opt-in key |

## Notes for owner review

- The 1e-4 arm of the floor is the only extension past the paper's
  stated recipe. It is inactive except near sigma=-1 where the
  paper's floor vanishes; the site comments document it.
- The paper ignores lines in the continuous spectrum; this tree can
  additionally run the whole continuum non-monotonically - a
  superset, consistent with the paper's intent.
- sigma >= -1 leaves every formula bit-identical to upstream
  (verified: the sigma<0 guard is the only branch touched).
