# Comprehensive non-monotonic CMFGEN audit -- 30-Oct-2026

Branch `coupled-sob-jfnk` @ 7ea11323.  Scope per the 12-phase brief:
correctness of the radiative-transfer, moment and non-LTE population
equations under non-monotonic velocity fields -- not merely iterative
stopping.  Builds on and re-verifies the earlier audits in
`docs/audits/KNOP_AUDIT_2026-10-30.md`, `KNOP_FORMAL_SOLVER_AUDIT.md`,
and `cmfgen_obsidian/audits/KNOP_MATH_AUDIT.md`.  New this session: the
first runtime driver for `MOM_J_CMF_NONMONO_V1`
(`unit_tests/f90/drv_mom_solve.f` + `python/test_mom_solve.py`), giving
coefficient-level verification of the moment discrete operator and the
mandatory zero-GAM runtime evidence.

## Phase 1 -- execution path

Per great iteration `cmfgen_sub.f` DO 10000 ML=1,NCF:

1. `comp_j_blank.f` assembles TA (scattering-lagged emissivity under
   COHERENT_ES), then one of three formal paths:
   - `FG_J_CMF_NONMONO_V1` when `DO_NONMONO_KNOP` or
     (`DO_NONMONO_SOB`+`USE_LAM_ES`+`DO_SOB_FG_KNOP`) (comp_j_blank:844;
     `SET_FG_NU_NEXT` at :843 passes nu_{k+1} for the forward spacing).
   - `FG_J_CMF_V13` otherwise (stock path, :857).
   - Plane-parallel / USE_FORMAL_REL branches bypass all of this.
2. Moment solve: `MOM_J_CMF_NONMONO_V1` when
   `(DO_NONMONO_KNOP .OR. DO_NONMONO_SOB) .AND. .NOT. USE_LAM_ES`
   (comp_j_blank:1032; `SET_MOM_NU_NEXT` :1030); `MOM_J_CMF_V11` stock
   otherwise.  Under `USE_LAM_ES` the Eddington factors instead come
   from `GET_MOMS_NON_REL` on the ray store (:978).
3. Post-sweep (cmfgen_sub.f): `FG_KNOP_REFRESH_V1` (:3282,
   `DO_KNOP_REFRESH` default T under KNOP/SOB_FG_KNOP) converges the
   lagged I_NEXT store on captured opacities; `MOM_KNOP_REFRESH_V1`
   (:3293) the lagged J/H_NEXT stores.
4. `SOLVE_FOR_POPS` (:4544) reads `GET_FG_GS_MATURITY` /
   `GET_MOM_GS_MATURITY` for the GS warm-up/maturity gate before any
   BA/Newton update, then `SOLVEBA_V13` / `JFNK_SOLVE_CORR_V1` (:2632)
   when DO_JFNK is armed.
5. `RECOMPUTE_SIGMA_TRUE` at model setup (cmfgen_sub:1134), after
   `DO_WIND_VEL_V1`/hydro revision (:4638), after R-grid revision
   (solve_for_pops:3667) and after NG acceleration (:3784).
   `ADJUST_R_GRID_V4` rejects VEL_SPACE interpolation under the flags
   and leaves the recompute to the caller.
6. Checkpoints: `DUMP_FG_INEXT_GS`/`DUMP_MOM_JNEXT_GS` at every
   `SCR_RITE_V2`, tagged with NITSF; `RESTORE_*_GS` validates tag +
   expanded-grid fingerprint before acceptance.

Silent-fallback audit: every non-mono dispatch is flag-gated; with
flags off the binary is the stock path.  Within the non-mono paths the
identified fallback-like behaviours are all explicit and counted:
`KN_NOSLOT` (compressing row with no bracket -> keeps chi_hat
absorption, drops the |Q|I_next source; physically the correct "no
stored upwind" closure, counted per ray), `XM_FB` previous-valid-field
rejection on failed solves, `GS_FREEZE`-style maturity rollback.
No path silently reverts to the monotonic operator inside a non-mono
run.

## Phase 2 -- transfer equation derivation

Independently re-derived.  Stock CMFGEN solves the nu-form O(v/c) CMF
equation `dI/ds = eta - chi*I + a_nu*dI/dln nu` with
`a_nu = GAM = (V/R)(1+sigma*mu^2)/c`.  The Knop Paper-I split in nu
form, verified against the code:

- `chi_hat = chi + |GAM|/DLN` (>0 always; `a*|p|` is positive for both
  signs of a) -- fg BUILD_FREQ and mom CHI_HAT_NM both implement this.
- Signed coupling `Q = GAM/(DLN*chi_hat)`; source
  `S' = eta/chi_hat + |Q|*I_upwind` with upwind = bluer k-1 (GAM>=0)
  or redder k+1 lagged store (GAM<0).  Sign convention: expanding cells
  see J_{k-1} live (Gauss-Seidel), compressing cells see J_{k+1} from
  the lagged store (Jacobi, converged by the refresh) -- matches Paper
  II Sect. 3.2.
- Frequency derivative: `DLN_C = dLOG_NU` expanding,
  `ln(nu_k/nu_{k+1})` compressing (asymmetric-spacing fix) --
  verified.
- Ray/moment consistency: the moment solver uses the angle-integrated
  analogue `a_J = CON_GAM*(1+sigma_true*K_ON_J)` at depth points,
  `a_H = CON_GAMH*(1+av_sigma_true*NMID_ON_HMID)` at midpoints;
  `chi_hat(midpoint) = chi_avg + |a_H|/DLN`.  Same Eq-7 form per
  equation.
- Zero-GAM condition and boundary closures analysed in Phases 4/6.
- The `3*a_nu*I` aberration term is absent in stock V13 and absent
  here -- inherited upstream, deliberately not re-added.

## Phase 3 -- formal solver (fg_j_cmf_nonmono_v1.f)

State after the 30-Oct static audit (committed fixes in 4f3966c0):
sweep vs refresh solve the SAME discrete operator -- verified term by
term (NI=2 end-slope, INSERT/MON DTAU clamp parity, identical A0-A4
coefficients, chi_hat, |Q| upwind source, monotone-parabola dS
clamps, boundary set).  Runtime coverage:
`unit_tests/f90/drv_knop_ray.f90` + `test_knop_ray.py` T1-T8
(manufactured compressing ray): alternating+ordered sweeps converge
(resid ~6e-7), ord/alt committed stores agree to 2.6e-8,
FREQ_MATCH multi-bin slots fire under strong drift, forced stall/cap
exits visibly, FGINEXT checkpoint round-trip bit-exact with
fingerprint/tag/foreign-seed rejection, FG_KNOP_REFRESH_V1 end-to-end
converges and is OpenMP bit-identical 1 vs 4 threads.

## Phase 4 -- moment solver (mom_j_cmf_nonmono_v1.f)

New runtime verification this session (`test_mom_solve.py`, 85 checks,
all PASS).  On manufactured atmospheres (ND_SM=24, sigma_true-
integrated V profiles), every production coefficient was recomputed
independently and compared:

- `GAM_SIGNED_NM = CON_GAM*GAM_FAC_NM` -- bit-exact (0 diff).
- `CHI_HAT_NM = chi + |GAM|/max(DLN_J,1e-6)` -- bit-exact.
- `Q_NM = GAM_NM*GAM_FAC_NM`, `GAM_NM = CON_GAM/chi_hat/DLN` -- exact.
- `PSI_NM = DTAUONQ*|Q|` interior -- exact; `PSIPREV=PSI` under
  MOM_GAM_PIN -- exact.
- `W = GAMH*|1+av_sigma*g|`, `GAMH = CON_GAMH/DLN/(chihat_i+chihat_{i+1})`
  -- exact; `EPS = GAMH*av_sigma*NMID_ON_J/(1+W)` -- exact.
- `HU/HL = K_ON_J*QF/((1+W)*DTAU)`, `HS = W/(1+W)` (pinned form in
  both direction branches) -- exact.
- Interior TA/TB/TC incl the EPS depth-upwind scatter (+-2E on
  diagonal vs off-diagonal by sign) -- exact.
- Tridiagonal residual |Ax-b| ~1e-14; interior rows M-matrix
  (TB>0, TA,TC<=0) on every manufactured profile including the
  sigma=-60 jump.
- Compressing rows confirmed to read the lagged J_NEXT slot
  (Q<0 rows exist; store populated; drift contracts over sweeps).

Zero-GAM (mandatory item).  `MOM_GAM_FAC_MIN=1e-2` floor +
exact-zero exception analysed (prior audit's limit derivation stands):
T1=1+sigma*K -> 0 gives GAM->0, chi_hat->chi, Q->0, W->0, EPS->1;
the unregularised operator is smooth through the resonance.  Runtime
evidence ('zero' case, 36 exact-resonance cells): Q=0 and
chi_hat=chi hold exactly at T1==0; cells at |T1|<1e-15 but !=0 take
the floored +/-1e-2 factor (a ~1e-2*CON_GAM local coupling bias);
R^2*J varies smoothly through the band (max cell-to-cell step 3.9%,
pure slope, no discontinuity).  Verdict: the floor is a defensible
regularization -- it biases only |T1|<1e-2 cells by ~1% of nominal
coupling; the exact-zero exception is measure-zero-inconsistent but
harmless (36/36 cells verified).  Recommendation stands as
documented, not adopted: a C1 sign-preserving map
(e.g. `T1*sqrt(1+eps^2/(T1^2+eps^2))`) would remove the residual
jump but the measured field continuity is already at the
discretization-noise level; no evidence it changes a converged
solution.  Classification: robustness/documentation, not a physics
error.

## Phase 5 -- velocity gradient (RECOMPUTE_SIGMA_TRUE)

`SIGMA_TRUE_CORE` uses the hybrid log/linear non-uniform quadratic
stencil (commit 7cd11cb0): `sigma = d ln|V|/d lnR - 1`, exact for any
power law V~R^b on any grid (homologous -> 0 to roundoff), with
signed-V stagnation handling and linear fallback near sign crossings;
RTOL-relative spacing guard.  `test_sigma_true.py` 16/16 incl.
analytic power-law, exponential, and zero-crossing cases with bounded
near-zone error.  The 1e6 cap is a numerical bound only; the audit
notes the classification distinction (physical shock vs gradient
artefact) is carried by sigma magnitude into chi_hat, so a
finite-but-large sigma is treated as a resolved steep gradient --
physically conservative, not wrong.  Grid-refinement convergence of
GAM is implied by the exact-power-law property; a dedicated
refinement sweep remains a test gap (see findings).

## Phase 6 -- frequency boundaries

Audited + runtime: freq-1 lagged slot intentionally absent (no reader
needs it); red edge keeps `|GAM|/DLN` in chi_hat without a source
partner -- photons advected below the grid are absorbed, the
documented zero-incident BC (the one place absorption lacks its
source mate; a BC choice, consistent between sweep and refresh).
Blue edge keeps the V13 INIT convention.  Matched-frequency reads are
bounds-checked (`FREQ_INDX+MOM_NM_J <= NCF_IN` at all four sites).
Driver runs exercised NCF=8 grids; smallest-NCF edge cases verified by
index audit (slot k+1 guard) rather than by 2-bin run -- noted as a
residual test gap.

## Phase 7 -- physical admissibility

Every repair path catalogued with reason and consequence:

- DGTSV failure -> whole-frequency quarantine to XM_FB (previous
  valid field), XM_BAD mask, GS metric forced unconverged.  Correct:
  no partial solves accepted.
- Significant negative J (>1e-6 col peak or >1e-3 local prev) ->
  reject to XM_FB; roundoff-level -> |XM|.  Non-finite/runaway
  (>1e3*prev) -> XM_FB.  H |H|<=J causality clip (V11 stock option).
- Non-finite lagged H/J reads -> zeroed + counted (MOM_INEXT_CAPPED).
- Boundary HBC/NBC floored to [0,1] when the measured closure is
  unphysical (negative flux frequency) -- replaces an unmeasurable
  Eddington with H(1)=0, logged to MOM_BND_HBC_FLOOR.
- MOM_EPS_LIM donor-cell limiter + MOM_LAG_BND row bound scale the
  lagged transfer so the fixed-point map cannot amplify past the
  diagonal excess -- both counted; they regularize the ITERATION, the
  converged fixed point is unchanged where the unbounded map already
  contracts.
- MOM_ABS_FLR (default 0) floors the absorption fraction so
  ESEC>CHI interpolation overshoot cannot cancel the diagonal.

Repaired cells are excluded from convergence metrics (XM_BAD -> a
>50% event, not a converged value): the acceptance invariant "a
repaired intensity is not a converged physical value" holds.
Manufactured-model runs: zero repairs, zero DGTSV failures, J>=0 on
all 9 velocity profiles.

## Phase 8 -- lagged-field fixed point

Residual `RES=max|VNEW-VOLD|/scale` is the unrelaxed fixed-point
residual |g(x)-x| of the discrete operator (verified in the refresh
path; bulk exit also gated on contraction + floor).  Driver evidence:
drift contracts on homo/comp1/dec/sigx/zero (monotone-ish ~x0.5/sweep
tail); strong/jump contract slower with non-monotone residuals but
still decreasing by it8 (35-46 from ~1e3 initial).  Alternating vs
ordered refresh agree to 2.6e-8 (test_knop_ray T3) -- different
acceleration strategies converge to the same field.  MOM_DAMP_NM and
JNEXT_RELAX alter the map multiplier, not the fixed point.  Period-2
oscillations are detected via the argmax 'cycle' classifier
(test_mom_sweep_table).  Restart: FG and MOM checkpoints bit-exact
with fingerprint+tag+foreign-seed rejection.

## Phase 9 -- population-solver integration

GS maturity gates (`GET_FG_GS_MATURITY`, `GET_MOM_GS_MATURITY`) gate
SOLVE_FOR_POPS: warm-up forces LAMBDA iterations while the lagged
field is immature; the checkpoint metrics (MOM_GS_CHG_FILE/CAP_FILE)
reject a "complete but still-moving" restored field.  JFNK admission
is gated on USE_LAM_ES + maturity; trace-demand floors and
trace-polish sweeps are documented in the AGENTS/vault notes.  BA
uses the solver-native diagonal ALO (GET_FG_ALO_DIAG /
GET_MOM_ALO_DIAG) so the linearization matches the forward solve;
quarantined cells contribute zero response (their J is frozen) --
verified in do_var_cont.
Energy-balance/temperature updates are upstream stock and unchanged.

## Phase 10 -- grid revision and persistent state

Grid rebuild: expanded-grid rebuild remaps J_NEXT/H_NEXT by
MON_INTERP (clamped query points), preserving GS state; V/R/SIGMA
fingerprints trigger rebuild on drift >1e-8 rel or any V change;
BA_MATRIX_STORED cleared (stale depth map).  Checkpoint: full R_NM +
R_PNT_NM fingerprint + R/V endpoints + tag -> stale or foreign files
rejected and zeroed (runtime-verified for both FG T7 and the new MOM
ckpt case).  RESTORE_MOM_JNEXT_LATEST handles BA rollback with a
fingerprint check.  Restart-equivalence: verified bit-exact on the
manufactured problems; production model runs pending as always.

## Phase 11 -- regression suite

`run_all_tests.py`: 19/19 files pass after the drv_sob_prof fix,
~500 checks total:

- test_mom_solve.py (new, 85 checks): discrete-coefficient exactness,
  residual, admissibility, fixed-point contraction, zero-GAM,
  sigx/jump/strong/dec, outbc=1/2, inner BCs, G_ONLY, checkpoint.
- test_knop_ray.py (8): ray solve, ord/alt equivalence, freq match,
  stall, checkpoint, refresh, OpenMP determinism.
- test_sigma_true.py (16), test_rv_hydro_rdin.py (9),
  test_sob_tau_floor.py (14), test_vadat_keys.py (70),
  test_mom_sweep_table.py (8), test_knop_lagged_field.py,
  test_sob_rh_couple.py (9), test_edep_shock.py (29),
  test_interp_guards.py (57), test_jfnk_solver.py,
  test_jfnk_regressions.py, test_ng_abstain.py, test_rd_rv_file.py,
  test_rv_snhyd.py, test_rv_sniin.py, test_sob_prof.py (5, after fix).

Runtime defects found and fixed this session:

- drv_sob_prof stale call signature (missing ARAD_VEC/C4_VEC/ED_IN
  added by the DO_SOB_DEP_WINGS change) -> OOB crash at
  cmf_flux_sub_v5.f:2919.  Test-maintenance defect, not a production
  bug, but it shows the fixed-form call convention carries no
  interface checking.

## Phase 12 -- findings classification

Critical physics error: none found.  The Knop Paper-I discretization
is verified term-by-term in both solvers; the moment operator is
confirmed numerically bit-exact against an independent recomputation.

Critical numerical error: none found in the audited paths.

Incorrect/incomplete boundary treatment: none blocking; the red-edge
absorbing BC is a documented modelling choice (Phase 6); the
HBC/NBC counted floor is a defensible closure for unmeasurable
Eddington factors (logged, not silent).

Convergence weakness (open):
- Moment lagged map on strongly compressing production models can
  stall (~30% floor, documented in convergence_issues); the
  retained-field rollback contains it, so no wrong answer is
  committed, but the +refresh on the stiffest models is the known
  weak point.
- strong/jump manufactured cases show non-monotone (but decreasing)
  residuals -- consistent with the shock stiffness, not divergence.

Robustness issues:
- Zero-GAM floor: O(1e-2) systematic bias on |T1|<1e-2 cells +
  measure-zero exact-zero exception; field continuity verified, fix
  documented-but-unadopted (needs model-run regression to justify).
- KN_CHI_EXT/KN_ETA_EXT single-precision capture: ~1e-7 rel error in
  refresh coefficients (documented memory tradeoff).
- FG_JUMP_MAX torn reads across OMP siblings (benign, finite).
- drv_sob_prof stale signature: fixed (this session).  Root cause --
  fixed-form F77-style calls have no interface checking; new CALL
  sites should prefer module procedures where practical.

Performance issues: none new beyond documented capture-store costs.

Documentation/test gaps:
- Grid-refinement sweep for GAM/J convergence (Phase 5) -- property
  verified on power laws; a refinement scan is not yet a test.
- NCF=2/3 minimal-frequency runs (Phase 6 edge cases) -- index-audited
  but not runtime-driven.
- MOM refresh end-to-end on a manufactured compressing model
  (MOM_KNOP_REFRESH_V1 re-drive via the captured KN_M store) -- the
  store is populated by the driver; wiring a refresh case is the next
  test addition.
- Population-coupled test (Phase 9 end-to-end) remains model-run
  only, per repo policy.

## Acceptance criteria status

- Formal/moment consistency: VERIFIED (coefficient-exact moment
  operator; same-equation sweep/refresh; ord/alt agreement 2.6e-8).
- Zero-GAM: derived + runtime-validated; floor is a defensible
  regularization with a documented small bias -- PASS WITH NOTE.
- Frequency boundaries: audited + guarded; PASS.
- Velocity gradient: exact on power laws, refinement gap noted;
  PASS WITH GAP.
- Unphysical fields rejected, not silently accepted: PASS (all
  repairs counted + excluded from convergence).
- Lagged-field convergence independently demonstrated: PASS
  (manufactured models + refresh machinery; production stiff-model
  stall is documented, not hidden).
- Population updates gated on field maturity: PASS (gates verified;
  end-to-end is model-run only).
- Grid/restart consistency: PASS (bit-exact checkpoints, remap).
- Monotonic CMFGEN unchanged: PASS (flag-gated; stock tests
  bit-identical where applicable).
- Strongly non-monotonic models converge reproducibly: PARTIAL --
  manufactured strong compression converges; the production m15
  moment path's intrinsic stiffness is the documented open item.

Verdict: PASS WITH LIMITATIONS.  The implementation solves the
correct non-monotonic RT + moment problem; the residual risk is
convergence rate on the stiffest production models, not correctness.
