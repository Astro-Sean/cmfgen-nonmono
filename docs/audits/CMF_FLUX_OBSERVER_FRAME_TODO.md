# CMF_FLUX Observer-Frame Implementation TODO

Companion to `CMF_FLUX_OBSERVER_FRAME_AUDIT.md` (HEAD `ce82a9f0`).
Work order follows the user's recommended priorities: confirmed
defects first, one pillar at a time. Status values: OPEN /
IN PROGRESS / DONE / BLOCKED. An item is DONE only after its
acceptance criteria are met numerically, not merely when code
compiles.

---

## P1. OpenMP index safety (Priority 1)

### P1.1 Exhaustive INDX write-before-read proof + hardening
- File/routine: `obs/obs_frame_sub_v9.f`, `OBS_FRAME_SUB_V9`.
- Problem: `INDX` is `PRIVATE` allocatable (:785-787); thread contents
  are undefined on region entry (:789-799 documents ifort behavior).
  Write-before-read appears to hold within each ML iteration but is
  unproven for `K=INDX(SM_NRAY)` at :1200 under TAU_MAX truncation.
- Solution: (a) add `INDX(1:NRAY)=1` (or FIRSTPRIVATE) per thread
  after the :804 conditional allocate, making the invariant explicit;
  (b) verify :1200 only executes when `SM_NRAY == NRAY` and
  `INDX(NRAY)` was written at :980 in the same iteration; document or
  guard.
- Tests: manufactured driver calling `OBS_FRAME_SUB_V9` with a
  truncated-tau atmosphere; `-check bounds` (ifort) build.
- Acceptance: no uninitialized-index path reachable; bounds-clean run.
- Status: DONE. Explicit per-thread `INDX=1` init added after the
  conditional allocate; INDX write-before-read verified across all
  branches. Committed `d1e09a31`. Bounds/uninit/fpe0 build clean at
  1/4/8 threads.

### P1.2 OpenMP reproducibility matrix
- Problem: no repeat-run evidence at OMP_NUM_THREADS = 1/2/4/8.
- Solution: run the OBSFRAME render (or the manufactured driver) x3
  per thread count on a fixed input; compare spectra.
- Tests: script comparing OBS_FLUX bitwise / within documented FP
  tolerance (target: identical, since each ML row is an independent
  reduction-free accumulate; the CLAMP_CENSUS criticals only count).
- Acceptance: NaN/Inf-free, bounds-clean, spectra equal within
  tolerance across thread counts.
- Status: DONE. 16-frequency case added to `drv_obs_frame_ib.f`
  (`8f4b19da`); all outputs bitwise identical (md5) at 1/2/4/8
  threads x2 reps, including the new ETAZ paths.

## P2. Inner-boundary transfer (Priority 2)

### P2.1 Derive ETAZ boundary term (defect D1)
- File/routine: `obs/obs_frame_sub_v9.f` STAU ~:1191-1210 vs ETAZ
  ~:1243-1263.
- Problem: STAU injects `I=J+3*mu*H` (*NU_ON_NU0^3) at the ray foot
  for `SM_NRAY==NRAY` under `DO_IB_DIFF_FLUX`; ETAZ starts `PAR_FLUX=0`.
- Solution: derive the boundary term in the emissivity-formal
  solution: `I_out = I_IB(mu_IB) * exp(-tau_total)` where tau_total is
  the ray's integrated optical depth foot-to-surface, I_IB built from
  the same interpolated IB_J/IB_H and the same mu aberration
  convention STAU uses. Implement in ETAZ branch; do not copy STAU's
  polynomial update - the attenuation factor differs.
- Tests (manufactured): A zero-eta positive-chi; B zero-eta thin;
  C thick envelope; D nonzero J,H; E negative-opacity fallback
  STAU->ETAZ; F monotonic V; G reversing V. Doubling IB_J,IB_H doubles
  boundary flux. STAU vs ETAZ agree in common regime at documented
  tolerance.
- Acceptance: tests A-G pass; STAU=ETAZ within tolerance on a shared
  atmosphere with DIF=T.
- Status: DONE. ETAZ branch now seeds `I_IB * exp(-TAU(IEND))` for
  `DO_IB_DIFF_FLUX .AND. SM_NRAY==NRAY .AND. LS<=NC` rays, built from
  the same interpolated IB_J/IB_H, foot mu, and NU_ON_NU0^3 factor as
  STAU. Manufactured tests in `test_obs_frame_ib.py`: thin STAU/ETAZ
  agree to 5e-9; doubled moments scale exactly 2x; thick model agrees
  at 9.6% (coarse ray) closing to 3.9% at 4x insertions - quadrature,
  not a missing term; negative-chi fallback and reversing-V cases
  finite.

### P2.2 IB_H provenance check (P3)
- Problem: `IB_H_CMF_ST = DBB/(3*CHI(ND))` analytic while
  `IB_J_CMF_ST = RJ(ND)` solved (`COMP_JCONT_V4.INC:114`).
- Solution: document the diffusion-limit convention; test consistency
  `J+3*mu*H` against solver-returned H where a solver provides one
  (USE_J_REL path) or against a manufactured diffusion-limit model.
- Acceptance: documented convention + test D above.
- Status: OPEN

### P2.3 Loc-1 blue-edge guard (defect D3)
- File/routine: `obs/obs_frame_sub_v9.f:824-866`.
- Problem: no blue-edge clamp at loc 1; `K=1` reads `SM_FREQ_CMF(0)`
  if `T1 >= SM_FREQ_CMF(1)`.
- Solution: add the symmetric blue-edge clamp/abort mirroring
  :899-924 (bounded guard + census + abort beyond).
- Tests: unit test with `T1 > SM_FREQ_CMF(1)` injected; bounds build.
- Acceptance: no element-0 read on any input; census counts the clamp.
- Status: DONE. Loc-1 blue-edge guard added mirroring the I>1 branch
  (K=2 + T1 clamp + census, abort beyond 2% margin). Committed
  `d1e09a31`; pin test updated for the fourth census site (`3c12e629`).

### P2.4 `K=INDX(SM_NRAY)` lower-bound guard (P2)
- File/routine: `obs_frame_sub_v9.f:1200-1203`.
- Problem: `SM_FREQ_CMF(K-1)` unguarded for K=1.
- Solution: same clamp pattern; or assert `2 <= K <= SM_NCF`.
- Acceptance: bounds-clean with K=1 manufactured input.
- Status: DONE. Resolved by construction: with the loc-1 blue-edge
  guard, every stored `INDX` value is in [2,SM_NCF], so
  `K=INDX(SM_NRAY) >= 2` and `J-1 = STRT_INDX_CMF+K-2 >= STRT >= 1`.
  No element-0 read reachable; bounds build clean.

## P3. Clumping consistency (Priority 3)

### P3.1 Convention proof
- Problem: establish that `CLUMP_FAC` multiplies effective opacity
  and emissivity exactly once across the pipeline; verify ED in
  `ESEC=6.65E-15*ED` is the smooth (unclumped) electron density and
  that the ES-as-rho^2 convention matches the population solve.
- Evidence so far: stored arrays unclumped; J solve uses *_CLUMP
  (`COMP_JCONT_V4.INC:123-128`); render folds once (driver ~:2558);
  production deposit unclumped (D2). Spot-check `comp_j_blank.f`
  uses `ESEC_CLUMP` identically.
- Solution: write the convention into the audit; add an assert or
  comment trail where each fold happens.
- Tests: CLUMP_FAC = 1/2/5 manufactured run; deposited + continuum
  opacity scale linearly with f in both stored and rendered arrays.
- Acceptance: exactly one fold verified numerically.
- Status: DONE. `test_sob_flux_obs` now runs the production routine
  on a depth-varying CLUMP_FAC=(1.0,2.0,0.5,1.5) model: the stored
  deposit normalizes identically at every depth (clump-independent),
  integrates to CHIL*1e-15, and static pins assert the .INC deposit
  carries zero CLUMP_FAC factors while the render fold
  (`*_RND_ST*CLUMP_FAC`, cmf_flux_sub_v5.f) exists exactly once.
  ED in ESEC=6.65E-15*ED is the smooth RVTJ electron density -
  consistent with the unclumped-storage convention.

### P3.2 Resolve ADD_SOB_LINE_PROF divergence (defect D2)
- Files: `obs/ADD_SOB_LINE_PROF.INC`, routine in `cmf_flux_sub_v5.f`.
- Problem: .INC deposits clumped values and lacks wings/floor args;
  unit tests exercise the stale version.
- Solution: single authoritative implementation - regenerate the
  .INC from the production routine (unclumped convention, current
  signature), update the unit-test driver accordingly; remove the
  dead `CLUMP_FAC` arg from production or use it in the DROP_SUM/
  INJ_SUM diagnostics whose comments claim clumped units.
- Tests: `test_sob_prof` equivalent on the production routine;
  integral of deposited opacity over frequency = CHIL_MAT (unclumped)
  within FP tolerance.
- Acceptance: one implementation; tests assert the production
  convention.
- Status: DONE. The .INC now holds the production routine verbatim
  and `cmf_flux_sub_v5.f` INCLUDEs it - one copy, compiled by both
  production and tests. `drv_sob_prof.f` signature updated;
  `test_sob_flux_obs` expectations corrected to the unclumped
  convention (integrals CHIL*1e-15, drop/inj sums unclumped). Stale
  "clumped deposit" comments fixed. Also renamed the .INC driver
  binary to drv_sob_prof_inc - it previously collided with the
  f90 driver's bin/drv_sob_prof and could run the wrong routine.

### P3.3 `sob_line_es_update.f` convention reconciliation (P5)
- Problem: orphan routine assumes clumped-ETA_ST storage and computes
  `IB_H_ST` differing by `CLUMP_FAC(ND)` from production.
- Solution: either align it with the production storage convention
  and wire it into the phase-2 ES rescatter (its stated purpose), or
  mark it test-only with a corrected convention comment.
- Acceptance: documented decision; no silent convention drift.
- Status: OPEN

## P4. Sobolev deposition accounting (Priority 4)

### P4.1 Transition accounting diagnostic
- Problem: no per-transition record of CMF-blanketed vs
  Sobolev-deposited vs omitted; a line could in principle be dropped
  from both or deposited twice.
- Solution: extend `TRANS_INFO` (or a new table) with treatment flag:
  BLANK / SOB-deposited / SOB-dropped(maser) / off-grid / EW-only.
- Tests: synthetic line list with one of each class; totals reconcile.
- Acceptance: every transition in the model has exactly one recorded
  treatment; off-grid lines reported (already logged).
- Status: DONE. A census now runs in cmf_flux_sub_v5.f right after
  VEC_TRANS_TYPE is finalized (before INS_LINE_V5): counts
  BLANK/SOB/other, names every transition whose type matches neither
  channel (invisible to both - e.g. a 'CMF' type, which no inserter
  in the tree blanketed, upstream included), and warns when SOB
  lines exist with DO_SOB_FLUX_OBS=F (EW estimates only, no observed
  flux - upstream semantics). Double-treatment is excluded by
  construction: INS_LINE_V5 inserts only 'BLANK', the deposit only
  'SOB', and each line belongs to one blend group. Static pins in
  test_obs_nonmono_dispatch.py:test_trans_type_census.

### P4.2 Deposition edge cases
- Tests: positive/negative opacity, maser, blends, off-grid,
  near-boundary, multiple resonances, nonuniform NU, non-unity
  CLUMP_FAC; verify frequency-integrated deposit against independent
  quadrature; check `DO_SOB_DEP_WINGS` Voigt normalization per depth.
- Acceptance: integrated deposited chi/eta match independent calc
  within tolerance; truncated-edge profiles not renormalized as if
  complete (report truncation fraction).
- Status: MOSTLY DONE. Covered by test_sob_flux_obs (nonuniform NU,
  integral = CHIL*1e-15 vs independent quadrature at FP tolerance,
  clump-independence, maser/blend-repair drops, off-grid and
  non-SOB skips, positive deposit) and test_sob_prof (off-grid
  naming, >10% masing-drop report, non-monotonic NU abort). STILL
  OPEN: DO_SOB_DEP_WINGS Voigt path (per-depth A_VOIGT and
  normalization) has no manufactured test; profile truncation at
  grid edges is not reported.

## P5. Non-monotonic frequency coverage (Priority 5)

### P5.1 Window completeness + strict validation mode
- File/routine: `obs_frame_sub_v9.f` STRT_INDX_CMF/END_INDX_CMF,
  MAX_VMU, GAM_RAY, Doppler extrema.
- Problem: window built from `MAXVAL(ABS(VMU_RAY))` + margin; under
  DO_RELO the extremal Doppler factor `D=gamma*(1-mu*v/c)` needs both
  extrema checked; three clamp sites substitute boundary opacity
  within 2% silently (census only).
- Solution: compute min/max D per ray analytically from VMU/GAM
  extrema; assert window covers it; add strict mode reporting
  total/blue/red clamps, max Doppler excursion, coverage margin,
  affected (ray,ML) indices; gate production acceptance on zero
  unexplained clamps.
- Tests: reversing-V model engineered to hit edges; non-rel vs rel.
- Acceptance: zero unexplained clamps on reference models; margin
  documented.
- Status: DONE. See docs/audits/CMF_FLUX_P5_FREQUENCY_COVERAGE_AUDIT.md
  for the full derivation. In short: the required CMF set for a band
  is exactly [nu_lo*X_min, nu_hi*X_max] with X=(1-VMU)*GAM; the hunts
  place the bounds inside FREQ_CMF and the +-10 pad keeps the SM
  edges strictly outside them, so LOCAL undercoverage is unreachable
  by construction - every clamp is GLOBAL undercoverage. Changes:
  (a) a per-band serial assertion verifies the two bound comparisons
  before the OpenMP loop and aborts with LOCAL_WINDOW_UNDERCOVERAGE
  or GLOBAL_CMF_GRID_UNDERCOVERAGE; (b) the four edge clamps now
  tolerate only the derived EDGE_TOL=1e-12 round-off band and abort
  with the same named diagnostics plus (LS,OUT_ML,ML,I,T1,SM edges,
  FREQ_CMF edges); (c) the census splits N_CLAMP_LOCAL/N_CLAMP_GLOBAL
  with per-class worst records, all under OMP CRITICAL; (d) both
  interpolation sites got a bracket-invariant guard (T1 finite, K in
  [2,SM_NCF], denominator >0, T2 in [-eps,1+eps], ETA/CHI finite);
  (e) X_MIN_RAY/X_MAX_RAY hoisted out of the band loop. Tests:
  drv_obs_frame_cov (mono/reversing/negative V, rel mapping + non-rel
  limit, NUM_BANDS 1-vs-2 bitwise equality, grid-extension bitwise
  stability, CHI feature detection, blue+red undercoverage aborts);
  OMP matrix bitwise at 1/2/4/8. Remaining: LOCAL coverage can only
  fire on a construction regression (by design).

## P6. Angular integration and flux normalization (Priority 6)

### P6.1 Quadrature + constant verification
- File/routine: `obs_frame_sub_v9.f` P/MU_AT_RMAX/HQW_AT_RMAX; final
  `6.59934*R(1)^2`.
- Problem: verify `F=2*pi*int(I*p*dp)/D^2` implementation and the Jy
  constant from first principles (R in 1e10 cm, D=1 kpc, Jy=1e-23).
- Solution: dimensional derivation using FUN_PI()/PARSEC() values;
  constant-surface-intensity test vs `pi*I*(R/D)^2`; luminosity check
  `L=4*pi*R^2*F`.
- Tests: manufactured constant-I atmosphere; COMP_OBS_V3 vs
  OBS_FRAME_SUB_V9 cross-check on the same model.
- Acceptance: numerical integral matches analytic within quadrature
  tolerance; constant reproduced to quoted digits or corrected.
- Status: DONE. See docs/audits/CMF_FLUX_P6_FLUX_NORMALISATION_AUDIT.md
  - Jy constant derived to 6.5993417128 (code uses 6.59934, upstream
    precision); PAR_FLUX identified as cgs specific intensity;
    HQW_AT_RMAX verified as the Nordlund mu dmu quadrature;
    constant-I and limb-law analytic tests pass; angular convergence
    demonstrated NP=9..129; flat-spectrum velocity invariance holds;
    pure D^2 boost scales as beta^2; tree V9 bitwise identical to
    frozen upstream V9 on all 24 manufactured outputs.
    No production defect found; no source change made.

## P7. Relativistic transformations (Priority 7)

### P7.1 GAM/NU_ON_NU0 construction + aberration audit
- File/routine: `obs_frame_sub_v9.f` VMU_RAY/GAM_RAY build (~:560-700),
  transform applications (:872-875, :987-990, ~:1209).
- Verified by inspection: eta*NU0^2, chi/NU0, I*NU0^3 powers.
- Open: whether GAM uses projected vs total v; whether boundary mu in
  `J+3*mu*H` should be aberrated; no `|v|>=c` input rejection seen.
- Solution: confirm GAM construction; aberrated-mu decision with
  error estimate; add `ABS(V)/C >= 1` hard stop (no radicand floor).
- Tests: zero-V limit = non-rel solution; uniform-V slab vs analytic
  Doppler transform.
- Acceptance: powers + aberration documented; superluminal input
  rejected.
- Status: DONE. See docs/audits/CMF_FLUX_P7_RELATIVISTIC_AUDIT.md
  - Doppler convention derived from the photon four-momentum and
    verified (NU_ON_NU0 = 1/GAM/(1-beta*mu), signed V_RAY kept);
    GAM uses total radial speed; eta*D^2, chi/D, I*D^3 applied
    exactly once at the local ray point.
  - Defect fixed: diffusion IB seed evaluated J+3*mu*H with mu_obs;
    now aberrated to mu_cmf=(mu_obs-beta)/(1-beta*mu_obs) at both
    STAU and ETAZ sites, gated on DO_REL_CORRECTIONS.
  - Boundary-condition physics verified (audit sec 4.3-4.4): the
    inner boundary is a fixed Eulerian surface, so entering
    characteristics are mu_obs>0 regardless of beta; the aberration
    maps them onto mu_cmf in (-beta,1].  For beta>0 and
    0<mu_obs<beta the photon is comoving-inward - interior radiation
    advected across the fixed boundary, not envelope incidence - and
    J+3*mu_cmf*H on the populated inward hemisphere is the correct
    seed.  Validity requires the diffusion premise at R(ND) and the
    linear closure; documented.
  - Defect fixed: PP boundary cosine was Z_RAY/R_RAY = 1/MU (a slant
    depth, not a cosine); now MU_AT_RMAX(LS) at both seed sites.
  - Defect fixed: T2<=0 warning + T2=1e-12 floor (GAM=1e6) replaced
    by fail-closed validation - NON_FINITE_VELOCITY on the V grid,
    SUPERLUMINAL_VELOCITY on grid and interpolated ray values
    (ulp-level overshoot inside V_RTOL=64*eps reprojected only),
    BAD_DOPPLER_DENOM on a non-positive 1-VMU.
  - test_obs_frame_rel.py: char battery - 80 single-ray seeds
    (beta 0/+-0.1/+-0.3 x mu 1.0/0.8/0.3/0.1 x STAU/ETAZ x sph/PP)
    vs independent forward-boosted diffusion-field oracle; pp
    battery - 8-ray PP grid, V=0 bitwise, +/-V vs oracle, OMP 1/8
    bitwise; full suite pass count in the audit doc.

## P8. Electron scattering + KNOP maturity (Priority 8)

### P8.1 ES redistribution verification
- Problem: verify ESEC*RJ(_ES) folded exactly once; coherent vs
  non-coherent coverage; RJ_ES provenance (COMP_J_CONV_V2 between
  NUM_ES_ITERATIONS).
- Tests: coherent vs non-coherent vs converged redistribution
  reference.
- Acceptance: single fold verified; spectra converge in NUM_ES.
- Status: AUDIT DONE + NUMERICAL TESTS DONE, acceptance pending
  independent review. Source audit in
  docs/audits/CMF_FLUX_P8_ES_REDISTRIBUTION_AUDIT.md; plan in
  CMF_FLUX_P8_IMPLEMENTATION_PLAN.md.
  - Single fold confirmed at every ESEC*RJ_ES site; OBS_FRAME adds
    no scattering term (ETA_RND arrives pre-folded).
  - RJ_ES provenance: EDDFACTOR convolved by COMP_J_CONV_V2 after
    each ES sweep (OBS side) and after each great-iteration (model
    side); render uses the freshest in-memory convolve.
  - test_es_conv.py 29/29 (drv_es_conv + drv_es_render; ifort
    -O0 -g -traceback -qopenmp -mkl=sequential): convolver
    battery (Planck 1.46e-5 max/1.28e-6 rms; temp 1.4e-3; edge
    9.1e-4; grid-saturated), one-fold elementwise census, T4
    coherent-vs-noncoherent flux discriminator (8e-3..1.4e-2),
    Q-RAY bound measured (DTA/TA <= 8.6e-4, dF/F <= 9.2e-3 at
    chi_ray=sigma_es, narrow-feature scattering-dominated case;
    <=1.8e-8 / <=7.2e-11 for smooth thermal-dominated fields).
  - Q-RAY resolved quantitatively (audit section 8 item 1):
    measurable only where Rayleigh is a large scattering
    fraction AND the field has feature-scale contrast. Inert
    under INC_RAY=F. No source change proposed - the model-side
    expression feeds the Jacobian, so a convention fix first
    needs the solver-implication map.
  - KNOP/ES call sequence verified at source: sweep k consumes
    Jbar_{k-1} from ES_J_CONV; FG/MOM refresh runs before the
    convolve-write inside each non-final sweep; render uses
    Jbar_N from the in-place RJ_CMF_ST convolve; no stale-RJ_ES
    path to the render in the covered chain.
  - Open items: no NUM_ES residual gate (P8.2); STEQ_T sign
    unverified; P7 3|H|<J positivity diagnostic deferred to
    P8/final validation.

### P8.2 KNOP maturity gate
- Problem: no OBS-side acceptance gate preventing a spectrum from an
  immature lagged I_NEXT field.
- Solution: configurable gate on the refresh maturity metric; derive
  threshold from convergence experiments, not invented values.
- Tests: spectra vs refresh count / NUM_ES_ITERATIONS / OBS grid /
  freq-match toggles.
- Acceptance: documented threshold + behavior below threshold.
- Status: OPEN

## P9. Spatial resolution convergence (Priority 9)

### P9.1 Insertion criterion + convergence study
- File/routine: `obs_frame_sub_v9.f` insertion on MAX_DEL_V_RES_ZONE,
  N_INS_OBS, OBS_ES_DTAU.
- Problem: insertion keyed on preceding-cell velocity width; reversal
  cells may be under-resolved.
- Tests: monotonic accel/decel, extrema, thin shocked shell, abrupt
  reversal, large gradients; N_INS_OBS and ES_DTAU convergence.
- Acceptance: spectrum converges under increased ray resolution on
  all test fields.
- Status: OPEN

## P10. Output files (Priority 10)

### P10.1 Product map + self-describing output
- Problem: OBSFLUX (COMP_OBS_V3) vs OBSFRAME (ray integral) are
  different methods sharing an array and similar names; negative
  flux is censused, not diagnosed; headers don't record mode/controls.
- Solution: document product map (done in audit); add header records
  (INT_METH, DO_RELO, DIF/IB_METH, clump, Sobolev phase, clamp
  census); investigate CHK_OBS_FLUX negative-bin causes before
  accepting warnings as benign.
- Acceptance: each output file's physical meaning and units written
  into its header; negative-bin report classified.
- Status: OPEN

## P11. Required numerical test suite (Pillar 13 of spec)

Tests 1-15 from the spec, implemented as production-source drivers:
1 static V=0 vs reference CMF; 2 monotonic vs stock; 3 non-monotonic
vs hi-res reference; 4 reversal coverage/finiteness; 5 STAU vs ETAZ
DIF boundary; 6 clumping scaling; 7 Sobolev deposit integrals;
8 ES coherent/non-coherent; 9 relativistic transforms; 10 OpenMP
1/2/4/8; 11 NCF convergence; 12 N_INS_OBS convergence; 13 flux
conservation vs analytic; 14 thin/thick limits; 15 negative-opacity
no-silent-corruption. All compile production objects.
- Status: OPEN (built per-pillar as items land)

## P12. Final reporting

- Record per pillar: files changed, exact tests run, results,
  remaining concerns, git SHA. Commit each verified source change
  individually, push to `origin/main` per tree rules, after
  `run_all_tests.py` + `./build.sh` are clean.
- Status: OPEN

---

## Completed in prior work (reference)

- VTURB extended-grid normalization (max|VEXT|) - fixed.
- MAX_VMU signed -> ABS for CMF window - fixed.
- MOM_J_CMF_NONMONO_V1 dispatch under DO_NONMONO_SOB+LAM_ES=F - added.
- ACCURATE/INC_GRID fail-closed on non-monotonic detection - added.
- Edge-clamp census (N_EDGE_CLAMP, worst indices) - added.
- CHK_OBS_FLUX finite check + negative census - added.
- DO_NONMONO_FREQ_MATCH OBS auto-set - added.
- Knop OBS dispatch + NUM_ES>=2/COH_ES=F gate - added.
