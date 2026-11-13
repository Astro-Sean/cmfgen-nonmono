# P8 - Electron-scattering implementation plan

Companion to `CMF_FLUX_P8_ES_REDISTRIBUTION_AUDIT.md`. Proposed
work only - no production source is changed here. P8.1 is not
complete until the acceptance criteria below are met.

## 1. What the audit established

- Exactly-one-fold structure: `ESEC*RJ_ES` is the scattered
  emissivity; the observer-frame transfer receives it inside
  `ETA_RND_ST` and adds nothing itself.
- The OBS-side ES iteration is a fixed-count Lambda iteration
  with no measured convergence.
- The Rayleigh convention differs between `comp_j_blank`
  (`CHI_SCAT*RJ_ES`) and the obs side (`ESEC*RJ_ES + CHI_RAY*RJ`);
  inert while `INC_RAY=F`.
- The KNOP I_NEXT refresh already requires `NUM_ES>=2` and is
  fail-closed at input validation.

## 2. Manufactured tests (P8.1 numerical acceptance)

Extend `unit_tests/f90/drv_obs_frame_rel.f`-style drivers; new
driver `drv_es_fold.f` plus `test_es_fold.py`.

T1. **Exactly-one-fold census (algebraic).** Manufacture
    `ETA, CHI, ESEC, RJ_ES` fields; run the obs-side emissivity
    construction and assert
    `TA = ETA + ESEC*RJ_ES + (CHI_RAY+DUST)*RJ` elementwise
    (non-coherent) and `TA = ETA + CHI_SCAT*RJ` (coherent).
    A double fold shows up at O(ESEC*RJ_ES) - tolerance 1e-12.
T2. **Coherent limit.** With `COHERENT_ES=T`, `RJ_ES` is never
    read: assert output is bitwise `TA = ETA + CHI_SCAT*RJ` and
    that the `ES_J_CONV` unit is never opened.
T3. **Convolve identity oracle.** Feed `J = B_nu(T)` into
    `COMP_J_CONV_V2` in 'J PASSED VIA CALL' mode. The
    Planck-preserving scale should return `Jbar ~= B` at depth:
    acceptance `|Jbar - B|/B < 1e-6` for depths where the Wien
    shift applies, `< 1e-3` overall (the RH94 fit itself is not
    exact). A delta-function J input checks normalization
    `sum(R) = 1` within fit tolerance.
T4. **Render fold-in ordering.** Manufacture `RJ_CMF_ST`,
    `RAY_CMF_ST`, `DUST_CMF_ST`, `ESEC`: assert the rendered
    `ETA_RND = ETA + RAY*J + DUST*J + ESEC*conv(J)` - i.e. that
    Rayleigh/dust multiply the RAW stored J and only ESEC
    multiplies the convolved J. This pins the :2553-2579
    ordering against regression.
T5. **Lambda-iteration timing.** Drive two ES sweeps with a
    stubbed solver returning a scripted J sequence; assert sweep
    2's emissivity uses `Jbar_1` (the previous sweep), and the
    render uses `Jbar_2`. Catches a same-sweep RJ_ES feed
    (self-inconsistent) and an extra-refresh bug.
T6. **KNOP guard.** Existing `rd_cmf_flux_controls` STOP for
    `DO_NONMONO_KNOP + (NUM_ES<2 | COH_ES=T)` is exercised by a
    driver-mode negative test - assert the named STOP fires.
T7. **OBS_FRAME pass-through.** Feed `ETA_RND` containing a
    marked scattered-emissivity spike into the existing
    `obs_v9_adapter` harness; assert the render adds no second
    scattering term (compare against a hand-computed attenuation
    of the supplied emissivity only).
T8. **INC_RAY convention census.** With `INC_RAY=T` on a
    manufactured atmosphere, dump model-side TA and obs-side TA
    for the same fields; REPORT (do not fail) the
    `CHI_RAY*RJ_ES` vs `CHI_RAY*RJ` difference. Failing bound is
    set after the magnitude is known - a convention decision,
    not a test failure, pending the Q-RAY resolution.

Acceptance for P8.1: T1-T7 pass; T8 produces the reported
magnitude used to resolve Q-RAY. Exact pass counts reported as
with P5/P6/P7.

## 3. P8.2 maturity gate (proposal, not implemented)

Add a measured-residual option to `CMF_FLUX_SUB_V5`:

    ES_RESIDUAL = max over (I,ML) |Jbar_N - Jbar_{N-1}| / Jbar_N

computed inside the existing `COMP_J_CONV_V2` call (both arrays
already in memory there - zero new I/O). Print it after each
sweep; optionally STOP if the last sweep exceeds a VADAT-keyed
tolerance (default off). This reuses the Rybicki & Hummer
ALI-consistency idea at the level the code already implements:
the lag between consecutive convolved fields is the convergence
metric.

## 4. Open items carried into P8/P9

- Q-MATURITY: NUM_ES fixed count, no residual check (P8.2).
- Q-RAY: model-vs-obs Rayleigh convention under `INC_RAY=T` (T8).
- Q-STEQ-SIGN: sign convention of `ESEC*(RJ-RJ_ES)` in `STEQ_T`
  against the temperature-equation definition.
- Q-WASTE: last-sweep extra convolve write (harmless).
- P7-CARRYOVER: `3|H| < J` positivity diagnostic on production
  moments - evaluate during P8 or final validation; no clamp
  without evidence of violation and a derivation.
- Q-LEGACY: `disp/` legacy path not audited.
