# Audit: conventional (non-JFNK) CMFGEN nonlinear solver

- Date: 2026-11-03
- Tree: `/viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk`
- Branch: `main` @ `22d63df97e2a7c5596987d4d6983c894c54cfa27` (clean at audit time)
- Outer repo: `main` @ `4b2b97c...`, unrelated dirty files preserved untouched
- Scope: the conventional nonlinear iteration (Lambda + BA population solve)
  on monotonic and non-monotonic velocity fields. JFNK internals excluded;
  JFNK code read only for dispatch/maturity boundaries.
- Baselines: upstream `cur_cmf_18jun25`; prior audits
  `NONMONO_COMPREHENSIVE_AUDIT_2026-10-30.md`,
  `NONMONO_NUMERICAL_ARCHITECTURE_AUDIT_2026-11-03.md`,
  `KNOP_AUDIT_2026-10-30.md`, `KNOP_STALL_DIAGNOSIS_2026-10-31.md`,
  `ROBUSTNESS_AUDIT_2026-10.md`. Those audits' source traces and unit-test
  evidence are cited, not repeated. This audit adds the full conventional
  call graph, ordering analysis, and stock-path defect assessment.

## 1. Executive summary

- The conventional iteration is a fixed-point outer loop: per-frequency
  formal solve with populations frozen, SE residuals and BA Jacobian
  accumulated during the sweep, then one linearized population/Ne/T
  correction per great iteration. Ordering is mathematically consistent
  for monotonic fields; it is upstream CMFGEN structure.
- On a non-monotonic input with all new flags off, the conventional
  solver does not solve the non-monotonic problem: `SIGMA` is clamped
  at `-0.999` before transfer (`recompute_sigma_true.f:50`,
  `fg_j_cmf_v13.f:1069`). Compression is silently flattened. No defect
  for expanding winds; a hard domain boundary for the stated purpose.
- `DO_NONMONO_SOB` defaults still route the *continuum* formal solve
  through the clamped `FG_J_CMF_V13`; only opt-in `DO_SOB_FG_KNOP`
  gives a fully non-monotonic continuum field. Under
  `DO_NONMONO_SOB + USE_LAM_ES=F` the moments are non-monotonic but
  the Eddington factors come from the clamped V13 ray solve — a
  coefficient inconsistency inside an opt-in path.
- Only `DO_NONMONO_KNOP` (or `SOB + DO_SOB_FG_KNOP`) solves the
  continuum formal transfer non-monotonically; that path carries the
  lagged-field maturity gates (freeze/refresh) that hold population
  corrections until the I_NEXT/J_NEXT field stops moving. Those gates
  implement the Knop Paper-II nested iteration correctly in principle;
  the previously documented stall/oscillation/NaN issues remain open.
- Multiple-resonance physics: a single line's repeated resonances along
  a ray are handled by the ray optical-depth integration; nonlocal
  coupling between distinct resonance depths is the two-point
  Rybicki-Hummer kernel (`sob_rh_couple.f`) which drops pairs straddling
  `V <= 0` and has no third-resonance extinction. Overlapping CMF lines
  are explicitly not installed (`cmfgen_sub.f` near the SIM_LINES sum).
  This is an approximation layer, not a full multi-resonance solver.

## 2. Conventional iteration call graph

Outer ("great") iteration, `new_main/cmfgen_sub.f` loop over MAIN_COUNTER:

```
[iteration start: read IN_ITS, SPECIFY_IT_CYCLE_V3 mode control]
  for FREQ_INDX = 1..NCF (continuum frequencies):
    opacities/emissivities from current POPS (cross sections, bound-free,
      free-free, e-scattering; resonance-zone line opacity via SIM_LINES;
      negative-opacity check; optional KN_NO_GAIN)
    ZETA = ETA/CHI, THETA = CHI_SCAT/CHI
    COMP_J_BLANK  -> radiation field this frequency (below)
    accumulate SE residuals: STEQ (pops), STEQ_ED (charge), STEQ_T (rad.eq.),
      STEQ_T_EHB (e- energy balance), BA Jacobian when COMPUTE_BA
      (EVALSE_QWVJ/PRRR/DO_VAR_CONT/UPDATE_BA_FOR_LINE/SUB_SOB_LINE_V3)
  [end freq loop]
  FG_KNOP_REFRESH_V1 / MOM_KNOP_REFRESH_V1  (nonmono only, gated)
  two-photon, charge exchange, advection (STEQ_ADVEC_V4) or
    comoving d/dt (STEQ_CO_MOV_DERIV_V3), adiabatic cooling,
    shock power, radioactive deposition -> STEQ/BA terms
  SOLVE_FOR_POPS -> mode controller -> SOLVEBA_V13 or JFNK_SOLVE_CORR_V1
  grid revision / hydro revision / RECOMPUTE_SIGMA_TRUE /
    grey-T iteration / pop rescaling / checkpoint
  convergence tests (MAXCH vs conv criteria; safeguards; exit logic)
```

COMP_J_BLANK per frequency (`new_main/mod_subs/comp_j_blank.f`):

```
BCs: THK_CONT outer, DBB diffusion inner (PNT_SRCE/ZERO_FLUX opts)
clumped opacity CHI_CLUMP = CHI*CLUMP_FAC
branch:
  CONT_VEL & ACCURATE   : extended grid, FG_J_CMF_V13 + MOM_J_CMF_V11
                          (nonmono solvers NOT dispatched here)
  CONT_VEL & !ACCURATE  : ND grid, f-iteration DO WHILE(INACCURATE):
     formal solve:  KNOP ? FG_J_CMF_NONMONO_V1
                   : (SOB & LAM_ES & DO_SOB_FG_KNOP) ? FG_J_CMF_NONMONO_V1
                   : else FG_J_CMF_V13                  (ll. 831-872)
     J update:      USE_LAM_ES   -> RJ = TC (+DO_ES_ALI OAB step)
                   : nonmono mom -> MOM_J_CMF_NONMONO_V1
                   : else MOM_J_CMF_V11                   (ll. 940-1112)
     test |FEDD-FOLD|max vs ACC_EDD_FAC; forced exit >20 passes (warn only)
  EDDINGTON / else      : static J (JFEAU_IBC / NEWJSOLD)
|RJ(1)| > 1e30 -> STOP (production)
SOURCE = ZETA + THETA*RJ ; store RJ in EDDFACTOR record
```

SOLVE_FOR_POPS (`new_main/mod_subs/solve_for_pops.f`):

```
maturity gates (nonmono): GS_FREEZE -> skip SE update entirely
  (transfer-only sweep); release needs N_MATURE >= GS_REL_STABLE
JFNK admission gates -> JFNK_SOLVE_CORR_V1  (out of scope)
else SOLVEBA_V13:
   CMF_BLKBAND_V3 banded solve (stored/rebuilt BA per COMPUTE_BA)
   optional DO_BA_RES_ACC residual acceptance (rollback of POPS and
     lagged FG/MOM fields, trust region, per-equation quarantine)
   apply: POPS *= (1 - STEQ*SCALE), per-element BIG_LIM/LIT_LIM clamps,
     T floor T_MIN, optional POPS(2)=POPS(1) (SET_POPS_D2_EQ_D1)
GS_POPS_HELD rollback; MAXCH/MAXCH_APP demand metrics
NG acceleration (scheduled), AV averaging, FIDDLE_POP_CORRECTIONS
```

## 3. Per-stage table

| Stage | Inputs | Outputs / state modified | Held fixed | Iter nesting | Convergence | On failure |
|---|---|---|---|---|---|---|
| freq sweep | POPS, CHI/ETA, V, SIGMA(_TRUE) | RJ per freq; STEQ, BA | POPS, T, ED | inner (per freq) | none (single pass) | RJ>1e30 -> STOP |
| f-iteration | TA=ETA+ESEC*RJ(_ES) | FEDD, RJ | opacity | inner (<=20/freq) | max\|dFEDD\|<ACC_EDD_FAC | exit at 20, warning only |
| FG vs MOM | FEDD from ray solve | J from moments (or RJ=TC) | — | same pass | implicit via f-test | divergent f -> warning |
| refresh (KNOP) | lagged I_NEXT/J_NEXT | matured ray/mom stores | POPS | within iteration | RF residual < tol | stall counters, freeze |
| SE accumulation | RJ, rates, POPS | STEQ*, BA matrix | POPS | inner | residual evaluated | — |
| SOLVEBA | STEQ, BA (or stored) | SOL corrections | J (linearized) | outer | residual-accepted opt-in | rollback/demote/hold |
| pop apply | SOL, POPS | POPS, T, ED | J until next sweep | outer | MAXCH vs conv crit | clamps cap step |
| NG accel | clean step history | extrapolated POPS | — | outer | N clean records | deferred |
| T/EHB update | STEQ_T_EHB, BA_T_PAR_EHB | T correction | — | outer | part of MAXCH | FIXED_T holds |

## 4. Ordering and state consistency

Verified:

- The SE residual for each frequency uses the radiation field computed
  *this* sweep with the *current* populations — F(x, J(x)) with x frozen
  inside the sweep. The BA step then moves x; the next sweep recomputes
  J. Standard, consistent outer iteration.
- Radiative terms accumulated in the same sweep the J was computed —
  no cross-iteration mixing inside the residual build.
- Under DO_NONMONO_KNOP the lagged coupling field (I_NEXT for the
  compressing branch, J_NEXT/H_NEXT for the moment path) means the
  first sweeps' J is not converged in the wavelength boundary-value
  problem. The GS maturity machinery (GS_MATURE_TOL=0.05 bulk,
  GS_UNSAFE_TOL=0.5 single-cell, N_MATURE release streak, refresh
  routines, freeze = transfer-only iteration) is designed so the BA
  correction is only committed once the lagged field stops moving.
  This is the right ordering principle (Knop Paper II). Where the
  population update is held, MAXCH_APP is zeroed so convergence
  bookkeeping cannot read a hold as progress (solve_for_pops.f:2654-2675).
- BA linearization uses solver-matched diagonals: KNOP+LAM_ES reads
  the FG_J_CMF_NONMONO_V1 diagonal, KNOP+!LAM_ES the MOM diagonal,
  SOB+LAM_ES the V13 ALO (opt-in DO_FG13_ALO), SOB+!LAM_ES the MOM
  diagonal (`do_var_cont.f:341-525`). Verified consistent pairing.

Lagged/stale state inventory (by design unless noted):

1. EDDFACTOR-stored RJ with COMPUTE_EDDFAC=F — J lags one great
   iteration; standard Lambda-iteration structure.
2. RJ_ES (non-coherent ES) read from ES_J_CONV — previous iteration's
   J (comp_j_blank.f:739-757). One-iteration lag in the scattering
   source; interacts with frozen sweeps under nonmono (scattering
   field can lag the matured field an extra iteration).
3. Stored BA matrix when COMPUTE_BA=F (FIXBA / judged reuse) — stale
   Jacobian, residual still evaluated each call; with DO_BA_RES_ACC=F
   acceptance is not residual-checked (standard quasi-Newton).
4. Knop lagged stores I_NEXT/J_NEXT — refreshed within the iteration;
   maturity gates cover them.
5. FOLD warm start across sweeps — benign.

## 5. Conventional solver on monotonic fields

Verified facts:

- Stock path = `FG_J_CMF_V13` (ray formal solve, Eddington factors)
  + `MOM_J_CMF_V11` (moment J) or USE_LAM_ES Lambda update. Both assume
  expansion: V13 clamps `SIGMA_RAY >= -0.999` (`fg_j_cmf_v13.f:1069`),
  MOM_V11 uses `W = GAMH*(1+AV_SIGMA*NMID_ON_HMID)` with `1+W`
  denominators (`mom_j_cmf_v11.f:548-580`) — singular for sigma<-1.
  For a monotonic expanding wind these are inert. No defect demonstrated
  for the monotonic domain itself.
- Pre-existing defect already fixed: INNER_BND_METH='ZERO_FLUX' branch
  tested 'ZERO FLUX' (dead code; diffusion flux always ran ->
  spurious HFLUX_AT_IB). Fixed at comp_j_blank.f:1003-1006 and 1094-1097.
  Models run before the fix with ZERO_FLUX were silently wrong.
- f-iteration convergence tests only the Eddington factor; an
  unconverged pass count >20 exits with a LUER warning, no flag in the
  convergence accounting (comp_j_blank.f:1147-1156). For ES-dominated
  continua the Lambda map can need more than 20 passes; the residual
  passes silently into rates/SE. Medium-severity latent issue, upstream.
- RJ(1)>1e30 STOP is a correct hard failure, not false convergence.

## 6. DO_NONMONO (KNOP) formal transfer

The solver internals were verified term-by-term in
`NONMONO_COMPREHENSIVE_AUDIT_2026-10-30.md` and
`NONMONO_NUMERICAL_ARCHITECTURE_AUDIT_2026-11-03.md`:
`chi_hat = chi + |GAM|/dLOG_NU`, signed `Q = GAM/(dLOG_NU*chi_hat)`,
sign-chosen upwind (J_PREV for GAM>=0, lagged J_NEXT/I_NEXT for GAM<0),
forward-spacing fix for compressing depths, refresh routines for the
lagged field, phantom/zero-outflow slot at NCF for the lagged red edge.

This audit confirms at dispatch level:

- `SET_FG_NU_NEXT`/`SET_MOM_NU_NEXT` pass nu(k+1) for asymmetric
  spacing (comp_j_blank.f:851, 1048).
- FIRST_FREQ dLOG_NU one-sided fix is gated to nonmono
  (comp_j_blank.f:663-669); stock path untouched.
- Frequency-space transport direction follows the sign of GAM via the
  upwind branch — verified against the solver code and prior audits.
- Under KNOP the formal solver also supplies the Eddington factors via
  GET_MOMS_NON_REL when USE_LAM_ES=T (comp_j_blank.f:994-997) and the
  ray-measured HNU_AT_OB for ZERO_FLUX outer flux (1013-1015).

Open items carried from prior audits (not resolved here):

- Strongly compressing production models stall in the lagged moment map.
- Lagged-field oscillation prevents population release on a warm run.
- A cold run produced NaN CHI_RAY/ETA_RAY near a near-zero-GAM boundary
  cell (GAM=-4.96e-11, SIGMA=-1.04).
- Moment zero-GAM regularization has documented O(1e-2) bias and a
  discontinuity at exact zero.
- Minimal NCF=2/3 runtime tests and end-to-end moment-refresh runs
  remain pending; population-coupled validation is model-level only.

## 7. DO_NONMONO_SOB pathway

Verified:

- Line operators take `SIGMA_TRUE` (unclamped) unless the diagnostic
  `DO_SOB_CLAMPED_LINES` is set (`sub_sob_line_v3.f:144-158`);
  SOBEW gets the same gradient (193-201).
- Sobolev tau uses the Dessart floor for the singular denominator:
  `MAX(|1+sigma*mu^2|, -0.2*(1+sigma), 1D-4)` when sigma<0
  (`sobjbar_sim.f:124-126`). Physically motivated (tangent-point
  resonance), inactive for sigma>=0, but it bounds tau exactly in the
  compressing regime — activation is not counted in the file section
  audited; recommend a census per AGENTS 5.2.
- Nonlocal RH78 coupling (`sob_rh_couple.f`): two-point kernel only —
  third-resonance extinction omitted (documented limitation,
  bounded overestimate); pairs straddling V<=0 or a stagnation depth
  are dropped (lines 129-140). For STELLA-type profiles with
  sign-crossing velocity, some physical couplings are excluded.
- The remote term is applied to STEQ only, kept out of the BA Jacobian
  — correct: a lagged remote-pumping term has no local derivative
  (documented, verified on TypeIIn 10d).
- Overlapping CMF lines are not installed (`cmfgen_sub.f` SIM_LINES
  region): CHIL/ETAL are summed over simulated lines; the CMF path does
  not treat line-line overlap as coupled. Single-line repeated
  resonances are covered by the ray quadrature; line-line blending at
  one depth/frequency is an approximation.
- Continuum under SOB+LAM_ES defaults to V13 (clamped sigma);
  `DO_SOB_FG_KNOP` opts into the Knop formal solver. Under
  SOB+!LAM_ES, `MOM_J_CMF_NONMONO_V1` computes J but V13 supplies
  Eddington/ray moments — mixed-fidelity coefficient path (opt-in).
- `STEQ_ADVEC_V4` refuses INC_ADV with any V<=0 under nonmono flags
  (steq_advec_v4.f:66-69): advection terms unavailable for profiles
  containing infalling cells — correct refusal, but a physics limit for
  interacting-SN structures.
- Negative-opacity handling: VB_SIM/VB_2 derivatives zeroed at
  NEG_OPACITY depths (sub_sob_line_v3.f:208-215) while the net rate
  still enters STEQ — semi-lagged masing treatment; KN_NO_GAIN can
  remove gain entirely.

## 8. Population / SE / T / Ne audit

Verified:

- STEQ construction (`steq_multi_v10.f`): collisional rates via
  SUBCOL_MULTI_V6 (CNM + dlnT derivative), SE equations 1..N_S,
  ionization-conservation EQION, number/charge conservation
  EQ_NUM_CONV; STEQ_T_EHB accumulates collisional heating-cooling;
  BA blocks for pops/ED/T/ion are assembled when COMPUTE_BA.
  Structure is standard CMFGEN; no sign or index defect found in the
  audited equations.
- SOLVEBA_V13: CMF_BLKBAND_V3 banded solve; correction applied as
  `POPS *= (1 - STEQ*SCALE)` with per-element caps, T floor T_MIN,
  optional depth-2=depth-1 closure; optional residual acceptance
  (DO_BA_RES_ACC) rolls back POPS and the lagged FG/MOM field stores.
- T and ED are unknowns in the same banded system (VT=NT, VED=NT-1);
  FIXED_T zeroes the T correction; FIXED_NE likewise via the ED row.
  STEQ_T_EHB vs radiative-equilibrium STEQ_T is selected by mode.
- LAMBDA_ITERATION holds JBAR fixed (Sobolev) and uses only the local
  derivative VB_SIM=-JBAR/ETAL (guarded against ETAL<=0,
  sub_sob_line_v3.f:282-304) — consistent with fixed-point semantics.
- Population corrections use the same-sweep radiation field; under
  nonmono the maturity gates can hold the correction entirely rather
  than commit it to a moving field (section 4). This is the correct
  ordering.

## 9. Clamp / floor / guard inventory (conventional path)

| Site | Bound | Trigger | Physically justified? | Counted? | Notes |
|---|---|---|---|---|---|
| fg_j_cmf_v13.f:1069 | SIGMA_RAY>=-0.999 | every compressing ray pt | No — rewrites dv/dr | No | stock; inert for expanding winds; hides compression |
| recompute_sigma_true.f:50 | SIGMA>=-0.999 | nonmono modes | infrastructure compat | No | SIGMA_TRUE kept unclamped |
| sobjbar_sim.f:124-126 | tau denom floor | sigma<0 | Partial — Dessart recipe | Not seen | bounds tau in compression |
| comp_j_blank.f:1147 | f-iter cap 20 | per freq | solver cap | Warn only | unconverged ES proceeds |
| comp_j_blank.f:1419 | RJ(1)>1e30 STOP | blowup | yes (hard fail) | STOP | correct |
| solveba apply | POPS*(1-STEQ*SCALE), caps | per element | solver safeguard | via MAXCH_APP | standard CMFGEN |
| solveba T floor | T>=T_MIN | per depth | solver safeguard | partial | |
| solveba D2=D1 | boundary closure | option | approx | option | |
| DO_BA_RES_ACC | residual rollback/trust/quarantine | opt-in | globalization | NREJ etc. | must not read as clean conv |
| sub_sob_line_v3.f:289 | ETAL_MAT>0 guard | lambda deriv | yes | branch | prevents wrong-sign deriv |
| steq_advec_v4.f:66 | refuse V<=0 | nonmono+INC_ADV | refusal, not clamp | message | physics limit |
| KN_NO_GAIN | zero gain | opt-in | physics choice | — | masing control |
| KNOP internals | GAM->0 regularizations | near-zero grad | documented bias | census files | prior audits |
| NEG_OPACITY | zero VB derivs | neg. opacity | semi-lagged | check loop | |

## 10. Controls and dispatch

- Flags default F: DO_NONMONO_SOB, DO_NONMONO_KNOP, KN_NO_GAIN,
  DO_NONLOC_SOB, DO_SOB_FG_KNOP, DO_SOB_CLAMPED_LINES, DO_ES_ALI,
  DO_JFNK. Mutual exclusion KNOP/SOB enforced; prerequisites checked
  (`rd_control_variables.f` ~lines 818-1478): SOB requires GLOBAL_LINE
  SOB and USE_LAM_ES=T for the intended path; DO_JFNK requires a
  nonmono flag; DO_NONLOC_SOB requires DO_NONMONO_SOB.
- USE_LAM_ES default handling updated Oct-2026: SOB default is F
  ("un-forced"); ordering dependence of derived defaults was
  specifically addressed (rd_control_variables.f:14-30, 1102).
- check_param_consistency.f:42-48 exempts nonmono modes from the
  SN-model USE_J_REL/USE_DJDT/USE_LAM_ES requirement so honoured F is
  not rejected — correct.
- Aliases: USE_JFNK deprecated for DO_JFNK with warning.
- Stock dispatch preserved when flags off: V13 + MOM_V11, no nonmono
  code executes — verified at the comp_j_blank dispatch.

## 11. Findings by severity

Numerical defects / capability limits (non-JFNK):

- HIGH — silent sigma clamp: on non-monotonic input the fully
  conventional path (all flags off) cannot represent compression; it
  produces a flattened-gradient solution that may converge and look
  plausible. A domain boundary users must not misread as capability.
- HIGH (opt-in path) — SOB+LAM_ES=T default continuum uses clamped
  V13 while lines use SIGMA_TRUE: inconsistent velocity-gradient
  treatment within one iteration unless DO_SOB_FG_KNOP=T.
- MEDIUM — SOB+!LAM_ES pairs nonmono moments with V13-derived
  Eddington factors (clamped ray solve): coefficient inconsistency.
- MEDIUM — f-iteration exit at >20 passes reports only a warning; an
  unconverged ES source function flows into rates/SE (upstream).
- MEDIUM — Sobolev tau floor bounds optical depth in compressing
  regions; no activation census observed in the audited section.
- MEDIUM — RH coupling drops sign-crossing pairs and lacks
  third-resonance extinction; interacting-SN geometries are partially
  covered only.
- MEDIUM — overlapping CMF lines not installed; blended lines are
  summed, not coupled.
- LOW — non-coherent ES (RJ_ES) lags one extra iteration relative to
  the matured field under frozen sweeps.
- FIXED pre-existing defect — ZERO_FLUX dead branch (see 5).

Physical-correctness boundary (not code defects):

- Convergence of the outer iteration does not by itself establish
  physical correctness under nonmono; the maturity/rollback metrics and
  the unresolved stall/NaN issues (section 6) govern that question.

## 12. Hypotheses (not verified this session)

- Under DO_NONMONO_KNOP + ACCURATE=T the extended-grid branch
  (comp_j_blank.f ~325) dispatches V13/MOM_V11 only — if ACCURATE can
  remain T under KNOP, the Knop solver is bypassed there. Likely
  guarded in practice (SN pipelines run ACCURATE=F); not traced to the
  VADAT guard in this session. **Verify before running ACCURATE=T
  nonmono models.**
- HFLUX_AT_OB under stock USE_LAM_ES inherits a lagged HBC_CMF(1)
  (standard); whether that lag matters for RLUMST diagnostics on the
  first iterations is untested.
- Whether the f-iteration's 20-pass cap has ever actually been hit in
  production SN models is unknown; grep OUTGEN for
  'Possible error converging f'.

## 13. Unresolved risks

- KNOP production stalls (warm oscillation, fresh ~75% stall),
  cold-run NaN in CHI_RAY/ETA_RAY near GAM~0 — carried from prior
  audits; not retested.
- Moment zero-GAM regularization bias ~1e-2 near threshold.
- No frequency-integrated conservation test in this audit.
- NCF=2/3 minimal runs and moment-refresh end-to-end runs still pending.
- No model-level run performed in this audit session; no claim of
  physical correctness is made.

## 14. Tests

Run this session: none — this was a source-tracing audit; no source
edits, no builds, no model submissions.

Relied upon (prior audits, same or nearby revisions):
unit_tests for KNOP formal solver (expansion/mono limits, asymmetric
spacing), MOM_NONMONO coefficient checks, SIGMA_TRUE derivative tests,
SOB path tests, JFNK harnesses.

Not run (required before any correctness claim):
monotonic stock regression vs cur_cmf reference; mild reversal, strong
compression, multiple-reversal, shocked STELLA profiles; difficult
ionization; conservation checks (flux, frequency-integrated J/H,
boundary flux); spectrum + population + temperature comparison;
safeguard census and exit-reason verification per AGENTS 7-9;
`run_model_grid.py` end-to-end.

## 15. Recommended next steps

1. Resolve the ACCURATE=T dispatch question for KNOP (section 12).
2. Add an activation census for the sobjbar_sim tau floor and the V13
   SIGMA_RAY clamp; report counts per great iteration.
3. Decide whether SOB+LAM_ES should default DO_SOB_FG_KNOP=T for
   nonmono inputs, or warn loudly when the mixed-fidelity path runs.
4. Run the pending minimal NCF=2/3 tests and the moment-refresh
   end-to-end case; then the staged profile suite (mono, reversal,
   compression, multi-reversal, STELLA).
5. Treat the f-iteration cap exit as a countable non-convergence event
   in OUTGEN rather than a warning line.
6. Extend SOB_RH_COUPLE for sign-crossing pairs or document the bound
   of the dropped coupling for interacting-SN profiles.
7. Re-verify ZERO_FLUX-era models predating the typo fix.
