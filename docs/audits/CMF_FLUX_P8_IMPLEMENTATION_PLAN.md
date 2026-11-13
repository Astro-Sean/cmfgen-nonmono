# P8 - Electron-scattering implementation plan

Companion to `CMF_FLUX_P8_ES_REDISTRIBUTION_AUDIT.md`. Proposed
work only - no production source is changed here. P8.1 is not
complete until the acceptance criteria below are met.

Status 15-Nov-2026: the manufactured battery is implemented and
passing - `unit_tests/f90/drv_es_conv.f`,
`unit_tests/f90/drv_es_render.f`, `unit_tests/f90/drv_es_seq.f`,
`unit_tests/f90/drv_knop_guard.f`,
`unit_tests/python/test_es_conv.py` (67/67). Measured results
recorded in section 2.1 and folded into the audit document.
T5 and T6 are executable production-linked tests, not source
inspection.

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
- The convolver is NOT a discrete kernel matrix. `COMP_J_CONV_V2`
  builds, per depth and per Rybicki & Hummer fit term L, a
  Helmholtz-type tridiagonal solve
  `(1 - (beta/BCOEF_L)^2 * d^2/dlnnu^2) D = J` and sums
  `J_ES = sum_L ACOEF_L * D_L` (`subs/comp_j_conv_v2.f:249-294`).
  `beta = 1.84E-3*sqrt(T)` (`:269`) gives an effective smearing
  width of beta/BCOEF ~ 3-4% in ln nu at T=1e4 K. No row-sum
  normalization test applies; normalization properties come from
  the operator structure plus the downstream Planck-preserving
  wavelength shift and rescale (`:302-354`).

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
T3. **Convolve production-operator battery.** Feed manufactured
    J fields into `COMP_J_CONV_V2` in 'J PASSED VIA CALL' mode.
    Replaces the earlier proposal of a delta-function
    `sum(R) = 1` check: the production operator is an implicit
    tridiagonal (Helmholtz) solve, not a stored kernel matrix,
    so no discrete kernel normalization exists to test.
    Required manufactured cases:
      - Planck preservation: `J = B_nu(T)` returns `Jbar ~= B`
        at depth (the routine is built to recover B). Metric on
        the significant-J support only; Wien-tail underflow
        points are masked since a relative error there measures
        denormal noise, not the operator.
      - Constant-field response: document, do not assume - the
        Planck-preserving rescale maps a constant to `B/B_es`
        where the rescale is active, and outputs zero where B
        underflows. Assert `Jbar >= 0` and finite.
      - Narrow spectral feature: a Gaussian line on a continuum;
        assert positivity, no overshoot below the continuum, and
        record the preserved contrast.
      - Grid boundaries: sharp features near NU(1) and NU(NCF);
        assert no boundary ringing beyond a recorded tolerance.
      - Temperature dependence: `B(T)` at several T; error may
        grow with beta (wider kernel) - record.
      - Grid refinement: repeat the Planck case at several NCF;
        record whether the error shrinks or saturates (the
        operator width is physical, not a discretization error,
        so saturation at ~1e-5 is expected).
    Record max and RMS error, positivity, and convergence
    behaviour for each case in the audit document.
T4. **Render fold-in ordering.** Manufacture `RJ_CMF_ST`,
    `RAY_CMF_ST`, `DUST_CMF_ST`, `ESEC`: assert the rendered
    `ETA_RND = ETA + RAY*J + DUST*J + ESEC*conv(J)` - i.e. that
    Rayleigh/dust multiply the RAW stored J and only ESEC
    multiplies the convolved J. This pins the :2553-2579
    ordering against regression. Discriminator: the coherent
    emissivity (no convolution anywhere) rendered through the
    real `OBS_FRAME_SUB_V9` must give a measurably different
    emergent flux, proving the convolved fold propagates.
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

## 2.1 Implemented coverage and measured results

`test_es_conv.py` (67/67 PASS; ifort production objects -
`compile_and_link_prod` links the exact `.o` files the
production executables use, plus `-qopenmp -mkl` per the
probed production link line; MKL sequential). Test-class
breakdown:

- **Executable production-linked**: T5 (sequencing harness on
  the real `COMP_J_CONV_V2`), T6 (the real
  `RD_OPTIONS_INTO_STORE` + `RD_CMF_FLUX_CONTROLS` path on a
  manufactured control file).
- **Real-convolver manufactured**: T3 battery (planck / const /
  linepair / bnd / tempu / edge / pos / grid), T4 render
  discriminator through the real `OBS_FRAME_SUB_V9`.
- **Algebraic reconstruction**: T1 emissivity rebuild, T8 Q-RAY
  census.
- **Source-inspection only**: the per-sweep call-order
  verification in audit section 6 (complements, not replaces,
  the T5 executable check).

T1 (algebraic): elementwise rebuild of both emissivity
conventions matches to <=1e-12 of the emissivity scale;
`DTA = CHI_RAY*CLUMP*(RJ_ES-RJ)` holds to machine noise at all
three Rayleigh amplitudes.

T3 battery, measured:

    planck  max 1.46e-5 / rms 1.28e-6 (significant support)
    const   non-negative; 3 zeros at the Wien-underflow edge;
            B/B_es rescale behaviour documented
    linepair narrow-feature battery (below)
    bnd     extended-domain boundary reference (below)
    tempu   uniform-T series (below)
    edge    max rel 9.1e-4, no ringing
    pos     Jbar >= 0 with a zeroed interior band
    grid    NCF=128/256/512 max err 1.32e-5/1.46e-5/1.66e-5 -
            saturates as expected (operator width is physical)

Narrow-feature smearing (linepair, baseline-subtracted):
`J_line = B + A*Gaussian` and `J_cont = B` convolved separately;
the excess is measured against its input, not the global max.
The production operator is a weak second-derivative smoother
((beta/B_L)^2 ~ 1e-6 in ln nu at 1e4 K), so FWHM growth is
unmeasurable; the asserted signature is the redistribution
pattern itself: peak drop, line-core depletion, wing lift,
recorded with the integrated-excess drift:

    NCF=256: peak 3.9986e-4 -> 3.9978e-4, core -7.5e-8,
             wing +2.1e-8, max|dexc|/pk 2.25e-4,
             int drift -3.1e-5, no negative ringing
    NCF=512: peak 3.9983e-4 -> 3.9975e-4, core -7.5e-8,
             wing +2.3e-8, max|dexc|/pk 2.26e-4,
             int drift -2.7e-5, no negative ringing

The response is resolution-consistent (2.25e-4 vs 2.26e-4),
so the measurement is not dominated by under-resolution
(feature: sigma_ln_nu = 0.2; 5.3 / 10.6 grid points per sigma).

Frequency boundaries (bnd/bndp): the extended domain appends
NCF/4 same-spacing nodes beyond BOTH production ends, so every
production node coincides exactly with an extended node - this
isolates the boundary treatment (one-sided -2/dx^2 stencil
rows plus coherent, non-interpolated endpoints per the source
comment) from resolution/interpolation error on steep slopes.
Result at NCF = 256 and 512:

    interior (middle 60%, significant mask): max 9.1e-5 / 8.5e-5
    upper edge (last 6 nodes):               0.0 / 0.0
    lower edge (first 6 nodes):              max 1.2e-3 / 2.4e-3

The upper edge is numerically identical because the field
there is Wien-dead (the extended content beyond NU=300 is
zero); the lower edge - the only end where J is alive at the
boundary - shows ~1-2e-3 local discrepancy, consistent with
the manufactured feature's ~0.4% tail lying beyond the
production boundary. Earlier versions of this test conflated
grid-resolution differences on steep slopes with boundary
distortion and are superseded.

Temperature dependence (tempu): uniform `T(k)=T_TEST`,
`J=B(T_TEST)` so the experiment isolates
beta(T)=1.84e-3*sqrt(T):

    T=5e3   beta=1.30e-3  beta/B=8.1e-4,6.0e-4  max 1.55e-4
    T=1e4   beta=1.84e-3  beta/B=1.14e-3,8.5e-4 max 1.82e-13
    T=5e4   beta=4.12e-3  beta/B=2.55e-3,1.91e-3 max 4.46e-14

RMS is computed on the same significant-J mask as the maximum.
The 5e3 case sits at the cold-temperature end of the fit's
intended range; the error there is still <=2e-4.

T4 (render fold-in): the coherent emissivity rendered through
the real `OBS_FRAME_SUB_V9` differs from the non-coherent
variant by 8.4e-3..1.4e-2 in emergent flux - the convolved
fold measurably propagates; the render adds nothing of its own
(same T7 evidence).

T5 (executable, `drv_es_seq.f`): the smallest sequencing
harness that mirrors the production ordering, NOT a link of
`CMF_FLUX_SUB_V5`. It reproduces the operations the audit
identified - per sweep: consume Jbar from the file-slot array,
produce J (scripted: `J_k = LAM*Jbar_{k-1} + NEWK*(1+PERT_k)`),
store ST=J_k, refresh slot (counter only, standing for
`FG_KNOP_REFRESH_V1`/`MOM_KNOP_REFRESH_V1`), then convolve the
STORED J through the real `COMP_J_CONV_V2` into the next
file-slot Jbar; the render folds ESEC*conv(J_N). NEWK is a
distinguishable spectral feature so the generations differ.
Broken-ordering modes fail as designed:

    'stale'   sweep 2 fed Jbar_0: perturbing NEW_1 leaves TA
              bitwise silent (the J_1 generation cannot
              propagate) while NEW_2 still responds - detected.
    'preconv' convolve before the store: same signature -
              detected.
    'oldrend' render folds Jbar_1: TA responds to NEW_1 but is
              bitwise silent to NEW_2 - detected.
    Correct mode: both scripted generations propagate;
              REFRESH_COUNT=1 between the sweeps.

T6 (executable, `drv_knop_guard.f`): the production store
reader + `RD_CMF_FLUX_CONTROLS` on a control file built from a
real model VADAT plus the OBS-required keys. Verified by the
named diagnostic, not the exit code (ifort STOP returns 0):

    A  KNOP=T NUM_ES=1 COH_ES=F -> STOP
       'DO_NONMONO_KNOP=T requires NUM_ES>=2 and
        COH_ERENT_ES=F' (emitted)
    B  KNOP=T NUM_ES=2 COH_ES=T -> same STOP (emitted)
    C  KNOP=T NUM_ES=2 COH_ES=F -> reader returns, GUARD_PASS

Q-RAY (T8) resolution, unchanged from the earlier battery:
scattering-dominated ETA with a narrow feature at the Rayleigh
band gives max|DTA/TA| = 1.5e-4/5.6e-4/8.6e-4 and
max|dF/F| = 1.6e-3/6.0e-3/9.2e-3 at chi_ray/sigma_es =
0.1/0.5/1.0; the smooth-continuum thermal-dominated case gave
<=1.8e-8 / <=7.2e-11. Quantitatively resolved; convention
decision deferred (not a complete production Jacobian
verification - see audit section 8).

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
- Q-RAY: RESOLVED quantitatively (section 2.1). Convention
  decision documented in the audit; no source change proposed
  until the solver-level implications on the model side
  (`comp_j_blank` feeding the Jacobian) are mapped.
- Q-STEQ-SIGN: sign convention of `ESEC*(RJ-RJ_ES)` in `STEQ_T`
  against the temperature-equation definition.
- Q-WASTE: last-sweep extra convolve write (harmless).
- P7-CARRYOVER: `3|H| < J` positivity diagnostic on production
  moments - evaluate during P8 or final validation; no clamp
  without evidence of violation and a derivation.
- Q-LEGACY: `disp/` legacy path not audited.
