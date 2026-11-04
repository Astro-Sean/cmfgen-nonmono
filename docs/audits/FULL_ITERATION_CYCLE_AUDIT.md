# Full Iteration-Cycle Audit - CMFGEN Non-Monotonic Velocity Implementation

Audit baseline: `cmfgen_jfnk` git commit `30b596bd744efd7f54dbaf4b5b749796bde0dfb4`
(`main`, 30b596bd "readme: note the population audit diagnostic").
Scope: one complete great iteration of `new_main/cmfgen_sub.f` and every solver,
controller, and accelerator it dispatches to. Audited configurations:
SOB/KNOP x JFNK/non-JFNK x `USE_LAM_ES` T/F.

This document is an audit. No source changes were made.

---

## 1. Executive summary

The iteration architecture is a nested fixed-point scheme: an outer "great
iteration" that alternates a full frequency sweep of the radiation solve with
one population correction, and an inner Gauss-Seidel (GS) sweep inside each
Knop radiation call that converges a lagged sign-dependent frequency-coupling
field (`I_NEXT` on the formal path, `J_NEXT_FRQ`/`H_NEXT_FRQ` on the moment
path). On top of that sits a large, recently grown control stack in
`new_main/mod_subs/solve_for_pops.f`: a field-maturity gate, a population
freeze, marginal-release trials, a luminosity runaway watchdog, a
LAMBDA/BA mode controller, JFNK admission gates, BA residual acceptance,
BA matrix reuse, and NG acceleration.

Main conclusions:

1. The cycle *can* converge - the `USE_LAM_ES=T` SOB recipe and warm KNOP
   sequences have converged in production. The gating stack does connect
   radiation maturity to the population update, which upstream CMFGEN never
   needed.
2. However, the architecture has grown by accretion of failure-driven
   patches (23 "edited" revision dates in Sep-Nov 2026 alone). Many gates
   measure proxies for the quantity of interest rather than the quantity
   itself, and several of the most important consistency questions can only
   be answered at runtime, not by inspection.
3. The most serious structural finding: the final convergence gate
   `CONV_EXIT_MET` contains **no radiation-maturity term**. Exit is gated on
   population-side demand, exempted demand, JFNK-commit status, held pops,
   pending verdicts, and NG bookkeeping only. Field drift enters indirectly
   via `GS_POPS_HELD` and the LAMBDA demotion - both of which have
   configuration-dependent coverage (findings F-01, F-04).
4. JFNK is permanently unavailable on the `USE_LAM_ES=F` moment pathway
   (admission requires `USE_LAM_ES`, `solve_for_pops.f:2182`). That is a
   deliberate choice driven by the missing eval freeze on the MOM lagged
   store, but it means the pathway the moment solver was built for - the
   non-monotonic default under SOB - can never use Newton acceleration
   (finding F-02).
5. The frozen-J JFNK tier is provably impotent under SOB and in any regime
   where the SE residual is radiation-pinned: 186 admissions, 0 committed
   Newton steps on the tested model. The coupled tier is the only tier that
   can move that residual, and it costs one full radiation sweep per
   Jacobian-vector product (findings F-03, F-22).
6. Radiation-population coupling is internally consistent *by construction*
   only at the fixed point: the SE residual and the BA matrix are assembled
   on the in-sweep `RJ` before `FG/MOM_KNOP_REFRESH_V1` converges the lagged
   store (findings F-04, F-05). This is the correct order for a lagged
   scheme, but it means a population update can be committed against a
   coupling field that the refresh then materially changes. The resacc
   net is the only protection.
7. The LAMBDA/BA mode controller is stable in the tested regimes but has
   documented trap states (ba_demand_hold probe cadence, cnt_lam escape,
   saturated-demand demotion, release/refreeze bound counters). It is
   fragile by design: each transition must reset a dozen counters and the
   audit found the resets are in place but only because each was added after
   a production failure.
8. The moment pathway under SOB (`USE_LAM_ES=F`) is documented FAIL at
   runtime: bounded non-contraction, ~200K quarantined cells per sweep,
   refresh stalls at residual ~1.1. The GS freeze/maturity machinery is
   KNOP-scoped, so under SOB the immature moment field has no hold or
   demote gate beyond `MAX_CHNG`/resacc (finding F-16).

The complete picture: the fixed-point being solved is coherent, but the
evidence that each iteration advances toward it is mostly indirect
(demand contraction, churn censuses, drift share thresholds). Several of the
rejection/hold paths are correct by construction; several gates are tuned to
one model family's failure modes and are untested elsewhere.

---

## 2. Complete execution-flow diagram

```
CMFGEN_SUB great iteration (label 20000, cmfgen_sub.f:1463)
│
├─ iteration setup (1463-1497)
│    NUM_ITS_TO_DO--, MAIN_COUNTER++, LST_ITERATION, DBG-T guard (KNOP)
│
├─ population/level data: POPION, EVAL_LTE_V5, SET_TWO_PHOT_V3 (~1500)
│
├─ dT/dR section: mini continuum pass at last depth (~1530-1640)
│
├─ zero SE/BA accumulator arrays: SE%BA, STEQ_T, STEQ_T_EHB, BA_T, BA_ED (~1650)
│
├─ STEQ_MULTI_V10/STEQ_AUTO_V2 + STEQNE_V4: J-independent SE terms (~1700)
│
├─ dielectronic loop: COMP_OPAC + COMP_J_BLANK + SOBEW per recombination freq
│
├─ MAIN CONTINUUM FREQUENCY LOOP  DO 10000 ML=1,NCF  (~2050-3270)
│    ├─ QUAD_MULTI_V10: TDF opacity per frequency
│    ├─ SET_LINE_OPAC / COMP_OPAC: chi/eta assembly
│    ├─ NEG_OPAC floors: ESEC_CHK / SRCE_CHK / KN_NO_GAIN (2581-2597)
│    ├─ COMP_J_BLANK (2056/2639/4028):
│    │     USE_LAM_ES=T:
│    │       KNOP or SOB+DO_SOB_FG_KNOP -> FG_J_CMF_NONMONO_V1 (SIGMA_TRUE)
│    │       else                        -> FG_J_CMF_V13
│    │       RJ = TC  (or ES-ALI Jacobi update, DO_ES_ALI)
│    │     USE_LAM_ES=F (moment equations):
│    │       KNOP or SOB  -> MOM_J_CMF_NONMONO_V1
│    │       stock        -> MOM_J_CMF_V11 / MOM_J_DDT / MOM_JREL / PP
│    │     ES-ALI optional on both branches (comp_j_blank.f:789-1010)
│    ├─ STEQ_T += FQW*(chi*RJ - eta)   <- RJ is THIS frequency's new J
│    ├─ DO_VAR_CONT: BA fill dJ/dchi, dJ/deta (banded or diagonal ALO)
│    ├─ BA_EHB_BF_UPDATE_V1, UPDATE_BA_FOR_LINE, TWO_PHOT_RATE
│    ├─ NORDFLUX/COMP_OBS_V2, RLUMST luminosity accumulators
│    └─ COMP_FREE_FREE_V2, BA_EHB_FF_UPDATE_V1
│
├─ GS REFRESH  (3281-3293)
│    FG_KNOP_REFRESH_V1  (formal; converges lagged I_NEXT store)
│    MOM_KNOP_REFRESH_V1 (moment; converges lagged J_NEXT/H_NEXT store)
│    NOTE: refresh runs AFTER STEQ_T/BA were accumulated on pre-refresh RJ
│
├─ residual assembly: STEQ_BA_TWO_PHOT_RATE_V3, charge exchange, Penning,
│    STEQ_CO_MOV_DERIV_V3/STEQ_ADVEC_V4, EVAL_TEMP_DDT/EVAL_ADIABATIC,
│    shock power, radioactive decay, PREVENT_LOW_T, DEP_RAD_EQ=STEQ_T
│
├─ LINE LOOP  (~3772-, skipped for BLANK)
│    COMP_OPAC, COMP_J_BLANK, SUB_CMF_LINE or SUB_SOB_LINE_V3 (SIGMA_TRUE)
│    -> ZNET net rates, LLUMST
│
├─ luminosity diagnostics, EDDFACTOR -> RJ_ES (non-coherent ES)
│
├─ SOLVE_FOR_POPS (4556-4558)  --> see Section 3
│
├─ POPTOION, COMP_LEV_DIS_BLK, SUP_TO_FULL_V4
├─ DO_WIND_VEL_V1 velocity update -> RECOMPUTE_SIGMA_TRUE
├─ AUTO_CLUMP_REV; INIT_BA_DATA_PNT_V3; COMP_J_CONV_V2 -> RJ_ES
│
├─ label 9999 (~4954+): exit/transition logic
│    BA_TRACE_EVAL, BA_DEM_EXIT_OK,
│    CONV_EXIT_MET = (RD_LAMBDA&&!DO_LAMBDA_AUTO | LAST_LAMBDA!=MC)
│                    && !JFNK_UNCONVERGED_COMMIT && BA_EX_MAXCH<EPS
│                    && BA_DEM_EXIT_OK && !DONE_HYDRO_REVISION
│                    && !GS_POPS_HELD && !BAACC_PENDING && MC!=LAST_NG
│    LST_ITERATION -> write RVTJ/POPs/DC, STOP
│    else ramps, CONV_EXIT_MET->NUM_ITS_TO_DO=1, BA_STAGEB_ARM,
│    DO_T_AUTO, DO_LAM_AUTO handoff, SPECIFY_IT_CYCLE_V3, IN_ITS re-read,
│    RD_LAMBDA re-impose, BA_RDLAM_PROBE, BA_TRACE_PHASE -> GOTO 20000
│
└─ internal JFNK_RAD_SWEEP_V1 (5766-): replayed sweep for coupled-J evals,
     saves/restores ACCESS_F, UP/DOWN_RATE_TWO, STEQ_T_*, FOLD, MOM
     quarantine/relax state; FREEZE_I_NEXT arms the I_NEXT eval freeze.

SOLVE_FOR_POPS (solve_for_pops.f)
│
├─ entry resets: FORCE_LAMBDA/GS_FREEZE/TRIAL_FIRED/JFNK_UNCONVERGED_COMMIT/
│    JFNK_COMMITTED; RDLAM_PROBE_THIS consumed; NEWMOD resets counters
├─ luminosity runaway watchdog: LUM_DIV_FAC tiers, LUM_BLOCKED latch
├─ KNOP warm-up: RD_LAMBDA && <= N_LAM_RWARM -> force LAMBDA+FIXED_T
│    (skipped when a mature FGINEXT_GS checkpoint was restored)
├─ field maturity metrics:
│    GS_MAX_ALL=max(GS_MAXCHG,GS_SW_MAXCHG,GS_RF_DCHG)
│    GS_FG_IMMATURE (raw>0.05 or bulk share + unsafe backstop; churn term)
│    MOM_UNSAFE_B/MOM_BROAD/MOM_FREEZE_B (moment-side equivalents)
├─ demotion under RD_LAMBDA: immature field -> forced LAMBDA+FIXED_T
├─ GS_RAW_IMMATURE -> N_MATURE counter -> GS_FREEZE (KNOP only)
├─ GS_PRECOND: demand > GS_PRECOND_MAX postpones the freeze
├─ frozen call: forced LAMBDA (GS_LAM_CAP or MAX_LAM_COR by churn);
│    MOM_TRIAL/FG_TRIAL marginal releases (resacc + snapshot gated)
├─ stock cycle (!RD_LAMBDA): BA step kept at MAX_LIN_COR, or capped LAMBDA
├─ GS release damping (GS_REL_N window), GS_PRECOND_LIM / GS_FREEZE_LIM /
│    GS_REL_NOPROG_LIM structural-failure stops
├─ GS_HOLD_POPS + GS_FREEZE (or LUM_POPS_HELD): snapshot, GS_POPS_HELD
├─ JFNK admission (see Section 6) -> JFNK_SOLVE_CORR_V1 or SOLVEBA_V13
├─ post-solve: held-pops rollback, BA rejection/partial rejection handling,
│    NG_CLEAN bookkeeping, trial refund, GS_STATE publish
├─ mode controller: CNT_STABLE_LAM hysteresis, BA_ENTER_RAW gate + probes,
│    BA_SAT demotion, BA_DAMP_CAP damped-full, T_HELD_PHASE release
├─ SCR_RITE_V2 checkpoint; REVISE_R_GRID; UNDO_IT/AVE_FLIPS; NG acceleration
└─ next-call flags: LAMBDA_ITERATION, COMPUTE_BA, FIXED_T, FIX_IMPURITY
```

---

## 3. The four solver configurations

### A. SOB + no JFNK (`DO_NONMONO_SOB=T`, `USE_JFNK_SOLVER=F`)

Two radiation sub-paths:

- `USE_LAM_ES=T`: formal J via `FG_J_CMF_V13` (expansion-biased at sigma<-1)
  unless `DO_SOB_FG_KNOP` promotes it to `FG_J_CMF_NONMONO_V1`. Ray moments
  feed `GET_MOMS_NON_REL` for the `N_ON_J` Eddington factors. Sobolev line
  rates in `SUB_SOB_LINE_V3` use `SIGMA_TRUE`. This is the converged
  production recipe on m15 (it718).
- `USE_LAM_ES=F` (SOB default since 21-Oct-2026): continuum J via
  `MOM_J_CMF_NONMONO_V1` with `DO_KNOP_REFRESH` default on, paired with the
  fixed-J `GET_MOM_ALO_DIAG`/`GET_MOM_CELL_QUAR` recipe in `do_var_cont.f`.
  `FG_J_CMF_V13` still runs for ray moments. **Runtime status FAIL**:
  bounded non-contraction, ~200K quarantined cells/sweep, refresh stalls at
  residual ~1.1. No GS freeze machinery under SOB - `GS_FREEZE` is forced F
  (solve_for_pops.f:1518-1522) - so the only rollbacks are `MAX_CHNG` and
  resacc.

### B. SOB + JFNK (`DO_NONMONO_SOB=T`, `USE_JFNK_SOLVER=T`)

- Admission requires `USE_LAM_ES` (frozen tier only way in). Under SOB the
  coupled tier is never demand-admitted (`solve_for_pops.f:2314-2318`);
  `JFNK_COUPLED_SOB=T` arms a reach-fail fallback inside
  `JFNK_SOLVE_CORR_V1` that fires only after the frozen tier certifies
  unreachability.
- Frozen-J under SOB cannot move a radiation-pinned residual - the residual
  is the same fixed-J system SOLVEBA solves directly (verified: 186 admits,
  0 GMRES steps). The stable-LAMBDA relief (`JFNK_NREQ_EFF->1`) exists to
  cut admissions under alternating LAMBDA regimes.
- Under `USE_LAM_ES=F`, JFNK cannot be admitted at all (F-02).

### C. KNOP + no JFNK (`DO_NONMONO_KNOP=T`, `USE_JFNK_SOLVER=F`)

- `USE_LAM_ES=T`: `FG_J_CMF_NONMONO_V1` per frequency + `FG_KNOP_REFRESH_V1`
  after the continuum loop. Full GS maturity gate, warm-up (`N_LAM_RWARM`
  forced LAMBDA unless a mature `FGINEXT_GS` was restored), freeze,
  marginal-release trials, release damping, hold-pops option.
- `USE_LAM_ES=F`: `MOM_J_CMF_NONMONO_V1` + `MOM_KNOP_REFRESH_V1`. Same gate
  structure on moment-side metrics (`MOM_SW_*`, `MOM_RF_*`). Boundary-flag,
  quarantine, negative-J, capped-write and huge-drift censuses.
- Warm-up demotes population calls to damped LAMBDA while the lagged field
  re-converges from a cold or reallocated store.

### D. KNOP + JFNK (`DO_NONMONO_KNOP=T`, `USE_JFNK_SOLVER=T`)

- Two-tier admission: frozen tier at `JFNK_FRZ_GATE` (fallback `JFNK_COR_GATE`),
  coupled tier at `JFNK_COR_GATE`; `JFNK_GATE_NREQ` consecutive sub-gate hits;
  demand contraction via `BA_LAST_RAW`/`JFNK_LAST_DEM` (1.5 rise factor);
  blocked by `BA_RAW_REJECTED`/`BA_RAW_TRIAL_PENDING`.
- Newton is not admitted while `GS_FREEZE` (field immature or inside the
  `N_MATURE < GS_REL_STABLE` window); `JFNK_CPL_RELWAIT` additionally holds
  the coupled tier through the post-release damp window `N_REL_SINCE <=
  GS_REL_N`.
- `JFNK_COUPLED=T` admits the coupled tier once `JFNK_CPL_HITS >= NREQ`;
  each eval replays `JFNK_RAD_SWEEP_V1` at trial populations with the
  `FREEZE_I_NEXT` store freeze armed.
- `DO_BA_DIR_CHK` is diagnostic-only admission priced against the stored BA
  (`info=-2` exit); `JFNK_FRZ_FAIL_LIM`/`JFNK_FRZ_FAIL_HOLD` hold future
  admissions after no-commit streaks.

---

## 4. Detailed findings

Severity: Critical / High / Medium / Low. Class: confirmed defect (C),
potential defect (P), architectural limitation (A), performance (Perf),
correct behaviour to preserve (OK).

### Radiation field and refresh

**F-01 [High, P] Convergence exit has no field-maturity term.**
`cmfgen_sub.f` label 9999 (~5024, ~5452). `CONV_EXIT_MET` requires
`BA_EX_MAXCH<EPS`, acceptable exempted demand, no pending verdict, no held
pops, not a post-NG state - but no `GS_FG_IMMATURE`/`MOM_*` term. In the
free phase a frozen call demotes to LAMBDA, and `LAST_LAMBDA=MAIN_COUNTER`
blocks exit. Under `RD_LAMBDA && !DO_LAMBDA_AUTO` (a user-forced pure-LAMBDA
KNOP run) the first disjunct is true on every call, so a quiet LAMBDA demand
can exit while the lagged field is still immature. Add an explicit
field-maturity conjunct to the exit, or document that pure-LAMBDA under KNOP
is not a convergence target. Test: a forced-LAMBDA KNOP run with an
artificially persistent churn cell must not exit converged (TEST 14).

**F-04 [Medium, OK-by-design] SE residual is assembled on pre-refresh J.**
`cmfgen_sub.f:2658-2660` accumulate `STEQ_T` and BA inside the frequency loop;
`FG/MOM_KNOP_REFRESH_V1` runs at 3281-3293 after the loop. This is the
correct order for a lagged GS scheme (the refresh only updates the coupling
store for the next call). But it means a committed population step was
linearized about a coupling field the refresh can then materially change;
`GS_RF_DCHG` measures exactly that change and feeds `GS_MAX_ALL`. Residual
acceptance (`DO_BA_RES_ACC`) is the only net catching a step that was safe
on the old field and bad on the new. Preserve the order; do not add a
post-refresh residual re-eval unless measurements justify it (F-22).

**F-05 [Medium, P] Stock-cycle BA step kept on an immature field.**
`solve_for_pops.f:1911-1928`: under `!RD_LAMBDA && GS_FREEZE`, a scheduled
BA call proceeds at `MAX_LIN_COR` ("BA step kept; immature field does not
replace it"). Protection is the resacc net only; with `DO_BA_RES_ACC=F` or
a missing snapshot the immature-field BA commit is unchecked. The
`BA_SNAPSHOT_MISSING` demote path exists, but only routes to LAMBDA when
resacc was armed. Potential defect for resacc-off configs.

**F-16 [High, A - resolved 04-Nov-2026] Field maturity has no gate under
SOB + moment.**
`GS_FREEZE` is set only inside `IF(DO_NONMONO_KNOP)` (solve_for_pops.f:1518).
Under `DO_NONMONO_SOB && !USE_LAM_ES`, the moment lagged field can be
immature with no freeze, no demote, no hold - only `MAX_CHNG` and resacc.
This is exactly the configuration documented FAIL at runtime (200K
quarantined cells/sweep). Whether the moment store should share the KNOP
freeze machinery is an open architectural question; the current state is a
coverage gap.

*Resolution (T-02/T-14, 04-Nov-2026): the census block is now gated on
`HAS_LAGGED_FIELD`, which covers SOB+LAM_ES=F, so `MOM_RAN`/`MOM_UNSAFE_B`
demote a scheduled BA step and `LAM_FIELD_NOT_READY` blocks the
DO_LAM_AUTO handoff until `N_MATURE >= GS_REL_STABLE` clean sweeps.
`GS_FREEZE` itself stays KNOP-scoped; the freeze machinery was not
ported. See FULL_ITERATION_CYCLE_TODO T-14.*

**F-17 [Medium, P] Refresh failure vs maturity reporting.**
`FG_KNOP_REFRESH_V1`/`MOM_KNOP_REFRESH_V1` exit on stall/cap/plateau;
`GS_NBAD_BLK` best-pass census mitigates stale counts, and `MOM_RF_CONV`
distinguishes converged vs not. A capped or incomplete refresh still writes
the store partially; bulk-pass maturity can then release on a partially
refreshed field (TEST 11). Verified mitigations exist but the
incomplete-refresh path needs a dedicated runtime test - potential defect.

### Population solve and controller

**F-06 [Medium, OK] Held-population bookkeeping is consistent.**
`GS_POPS_HELD` calls restore `POPS_GS_SAVE` (unless the resacc rejection
already restored), report `MAXCH_APP=0`, keep raw `MAXCH` at the attempted
value for diagnostics, and reset `NG_CLEAN` so the duplicate SCRTEMP record
cannot enter the NG window. Verified at solve_for_pops.f:2698-2719, 2898-2903.

**F-07 [Medium, P] MAXCH carries three different semantics.**
Raw demanded correction (SOLVEBA), applied correction (JFNK commit), or
held-call attempted demand (held). The controller compensates case-by-case:
`BA_LAST_RAW` excludes JFNK calls, `MAXCH_CTRL=MAXCH_APP` under KNOP,
`JFNK_UNCONVERGED_COMMIT` blocks exit, `JFNK_COMMITTED` resets NG_CLEAN.
This is correct today but fragile - any new gate keyed on `MAXCH` must
re-derive which semantic it sees. Architectural fragility, not a defect.

**F-09 [Medium, P] JFNK demand gates read the wrong operator under LAMBDA
dominance.** `JFNK_DEM_IN = BA_LAST_RAW` only refreshes on non-LAMBDA calls;
before the first full call it falls back to `MAXCH`, which in alternating
regimes is a LAMBDA-attractor demand (documented admission on a 95% LAMBDA
demand while the last full demand was 4.8E4%). Mitigated by the contraction
term and `JFNK_GATE_NREQ` hysteresis; the SOB stable-LAMBDA relief needs the
`LAM_STABLE_PROVEN` latch, whose reset coverage the audit verified at every
demote branch (3345-3348, 3528, 3561, 3581, 3606, 3632, 3844-3847).

**F-11 [Medium, OK-reviewed] KN_NO_GAIN opacity floor.**
`cmfgen_sub.f:2581-2597`: under KNOP, `CHI` is floored at `CHI_SCAT`
(forbids negative net opacity = maser amplification) so the lagged map
stays non-expanding. This changes local physics on transient iterates;
opt-in and documented. Keep flagged as a physics-altering guard per policy.

**F-19 [Medium, P] Transient guards change the equations on bad iterates.**
`PREVENT_LOW_T`, `CHECK_LINE_OPAC`/`CHIL=1.0` resets, `VB_SIM`/`VC_SIM`/`SRAT`
zeroing in `SUB_SOB_LINE_V3`, `RLUMST` floors at 1E-10/1E-20, `ES_ALI`
`DEN=max(.,0.01)` (comp_j_blank.f ~995). Each is a documented transient-only
repair, but they make the residual a different function on bad iterates than
on good ones - JFNK FD JVPs taken across the boundary see a discontinuous
operator. Acceptable while flagged; each needs an activation-rate diagnostic
(they have one) and a bounded-influence argument.

### JFNK integration

**F-02 [High, A - resolved 04-Nov-2026] JFNK permanently gated to
`USE_LAM_ES`.**
`solve_for_pops.f:2182, 2496-2506`: the moment lagged store (`J_NEXT_FRQ`/
`H_NEXT_FRQ`) has no equivalent of `FREEZE_I_NEXT`, so each residual/JVP eval
would overwrite the coupling field the committed iterate relies on. The
moment pathway - the non-monotonic default under SOB and one of two KNOP
paths - can never use Newton. Closing this requires a MOM-side eval freeze
(or coupled evals that fully re-converge the store before committing).

*Resolution (T-03, 04-Nov-2026): the combination is closed, not ported.
A MOM eval freeze alone is small (~4 guard sites), but useful admission
under LAM_ES=F additionally needs a coupled-only admission rewrite (the
frozen tier is the same fixed-J system SOLVEBA solves - F-03), a new
direct coupled trigger for the SOB-mom route (the current trigger is
frozen-reach certification), and a determinism audit over the MOM
module state. That is a substantial architecture change for unproven
gain, so `CHECK_PARAM_CONSISTENCY` now stops on
`USE_JFNK_SOLVER .AND. .NOT. USE_LAM_ES` with a clear message instead
of the armed-but-always-deferred state. The defer reason in
`solve_for_pops` remains as unreachable backstop. See
FULL_ITERATION_CYCLE_TODO T-03.*

**F-03 [Medium, A] Frozen tier is useless under SOB.**
Frozen-J Newton solves the same fixed-J system SOLVEBA does; on a
radiation-pinned residual the reach gate certifies zero useful descent
(186 admits, 0 GMRES). `JFNK_COUPLED_SOB=T` is the only escape and it fires
only post-certification. Recommendation: either auto-suppress the frozen
tier under SOB once `JFNK_REACH_SKIP_STREAK` exceeds `JFNK_REL_REACH_MAX`
(half-implemented - the relief stands down but admissions continue), or make
`JFNK_COUPLED_SOB` the documented SOB default.

**F-22 [Medium, A/Perf] Coupled-J costs one full sweep per JVP.**
`JFNK_RAD_SWEEP_V1` replays the entire radiation pathway per residual eval;
FD Jacobian-vector products need one eval per Krylov direction. The two-tier
admission exists precisely because of this cost. No cheaper Jacobian option
exists on the tree (no matrix-free physics-based preconditioner reuse
beyond BA). Measure before proposing any further acceleration.

**F-20 [Medium, OK] Eval-state save/restore is comprehensive by inspection.**
`JFNK_RAD_SWEEP_V1` saves/restores ACCESS_F cursors, `UP/DOWN_RATE_TWO`,
`STEQ_T_SCL/NO_SCL`, `FOLD`, MOM quarantine age/counters, FG/MOM relax state
and residual history, and `FREEZE_I_NEXT` freezes the I_NEXT store. Eval
failures route through `IOS_EVAL` to a broken-eval abort rather than a
silent commit. Repeated-eval determinism (F(x) reproducible) is asserted in
comments and unit-tested at `test_jfnk_opint`/`test_jfnk_regressions` level,
not at full-iteration level.

### Freeze/release and NG

**F-08 [Medium, A] Population hold is opt-in and immature by default.**
`GS_HOLD_POPS` default F: frozen iterations apply a damped LAMBDA update to
populations, not a true transfer-only sweep. The true nested-iteration mode
(`GS_HOLD_POPS=T`) is available but unproven as default. Freeze is
field-maturity-driven (good) but the release path has grown four separate
machineries (raw release, MOM_TRIAL, FG_TRIAL, churn arm) each with its own
refund/budget/stop - the cycle-bound counters (`GS_FREEZE_LIM`,
`GS_PRECOND_LIM`, `GS_REL_NOPROG_LIM`, `GS_RELEASE_LIM`) terminate runs
rather than guarantee progress.

**F-10 [Medium, P] Release/refreeze cycle is bounded, not prevented.**
`GS_JUST_REL` audits release progress, `N_REL_NOPROG` bounds no-progress
releases, `GS_REL_N` damps the first live calls, accepted trials refund to
`TRIAL_MIN_USED` not zero. A persistent marginal field can still alternate
release->detonate->refreeze within those budgets until a STOP. This is by
design (stop loudly), but the budgets are model-tuned.

**NG [Low-Medium, OK] NG window hygiene is handled.**
`NG_CLEAN>=4` requires four clean committed iterates after any reject,
masked rollback, soft reject, armed-pending step, JFNK commit, or held-pops
record; `NG_DEFERRED` handles the stepped-over +1 target. `NG_ABSTAIN`
validity check on extrapolation size. `RECOMPUTE_SIGMA_TRUE` fires after NG
so `SIGMA_TRUE` tracks the extrapolated V. The remaining risk is conceptual:
NG extrapolates the LAMBDA/BA outer map, while JFNK commits belong to a
different map - the code excludes JFNK records correctly.

### Global convergence

**F-12 [High, A] "Converged" is population-demand converged plus bookkeeping.**
The exit gate measures `BA_EX_MAXCH` (applied-scale major demand) plus
exempted-demand channels and bookkeeping. It does not directly measure:
the SE residual (`BA_SE_RES` is diagnostic only), field drift, or luminosity
consistency. The defense is architectural - quiet demand under a mature
field implies self-consistency - but several documented cases had quiet
demand coexisting with a churning field (churn-arm trials exist for exactly
this). `JFNK_UNCONVERGED_COMMIT` and `GS_POPS_HELD` blocks cover their
specific cases; a general "no exit while field immature" term is absent
(F-01 is the concrete hole).

**F-13 [Medium, P] Restart coverage is partial.**
`FGINEXT_GS` restores a mature I_NEXT store (with `GS_FOREIGN_SEED` opt-in
adoption); warm-up skips only on validated maturity. A missing/stale
checkpoint restarts warm-up - safe but slow; a corrupted checkpoint is
rejected by the maturity getter. MOM-side restore exists (`RESTORE_MOM_*`)
but the pair is not exercised by an integration test (TEST 9).

**F-14 [Low, OK] Temperature coupling is staged, not inconsistent.**
`FIXED_T`/`RD_FIX_T`/`MOD_FIXED_T`/`T_HELD_PHASE`/`DO_T_AUTO` tier the T
release; `N_ACC_FULL` counts accepted full steps; `CON_SCL_T` forces BA
rebuilds when T was partly held. The electron-energy equation and
`STEQ_T_EHB` are assembled in the same sweep as the populations - no stale
energy residual was found.

### Tests

**F-23 [High, A] No full-iteration integration tests exist.**
`unit_tests/python/` has 26 tests; they cover solver internals (knop_ray,
mom_solve, jfnk_solver, jfnk_opint, sob_tau_floor, ng_abstain, pop_audit,
gs_drift_diag, sigma_true, vadat_keys) but nothing drives a complete great
iteration end-to-end, compares formal-vs-moment fields on the same model,
or exercises freeze/release cycling, restart restore, or coupled-J
determinism in-situ. TESTS 1-15 in the TODO are all gaps.

### Performance

**F-24 [Medium, Perf] Work is duplicated by construction.**
Per great iteration: NCF formal solves (or moment solves), one refresh pass
over the full ray store, up to a full extra sweep per JFNK residual eval,
plus a dielectronic mini-loop. `COMPUTE_BA` reuse and `DO_BA_RES_REUSE`
bound the BA cost. Nothing in the tree measures whether refresh passes or
JFNK evals dominate wall time - the instrumentation plan in the TODO is a
prerequisite for any further acceleration work.

---

## 5. Assessment vs the audit questions

- **Radiation and populations advance consistently?** Yes at the fixed
  point, conditionally en route. The gating stack ties population updates to
  field maturity, but the SE residual is always one refresh behind the
  committed coupling store (F-04) and the exit gate never directly checks
  field state (F-01, F-12).
- **Radiation maturity measured correctly?** The metrics are the right
  quantities (max drift, bad-cell share, churn share, refresh delta, capped
  writes) but thresholds are model-tuned and the freeze gate is KNOP-scoped
  only (F-16).
- **JFNK admitted at the correct stage?** Mostly yes on the tested paths
  (warm-up, freeze, gate hysteresis, contraction, release-window hold). The
  known wrong-operator admission under LAMBDA dominance is mitigated but
  not eliminated (F-09); USE_LAM_ES exclusion is deliberate (F-02).
- **Frozen-J and coupled-J solve the intended equations?** Frozen-J solves
  the same fixed-J SE system as BA - verified. Coupled-J replays the sweep
  at trial pops with a frozen store - verified by inspection; determinism is
  unit-tested, not integration-tested.
- **LAMBDA/BA controller robust?** Stable but complex; every trap state has
  a bound counter and the audit found the resets in place. Fragile to new
  gates (F-07).
- **Freeze/release effective?** Yes with the caveats in F-08/F-10.
- **Acceleration preserves the fixed point?** NG hygiene is correct
  (NG_CLEAN); JFNK commits are a different map and are correctly excluded
  from the NG window.
- **Final convergence = physical self-consistency?** Not directly proven by
  the gate; it is inferred. F-01 and F-12 are the actionable holes.

---

## 6. Priority recommendations (maps to TODO task IDs)

1. (T-01) Add an explicit field-maturity conjunct to `CONV_EXIT_MET`, or
   gate it behind the same `GS_*`/refresh state the freeze uses. [F-01/F-12]
2. (T-02) Decide whether SOB+LAM_ES=F shares the KNOP freeze machinery or
   gets its own gate; today it runs bare. [F-16]
3. (T-03) Close or formalize the JFNK/moment exclusion: either implement a
   MOM-side eval freeze or document `USE_LAM_ES=F` + JFNK as unsupported
   with a hard consistency-check rejection. [F-02]
4. (T-04) SOB + JFNK: make the frozen tier self-suppressing after certified
   unreachability (already half-done via the relief de-rate), and evaluate
   `JFNK_COUPLED_SOB=T` as the SOB default. [F-03]
5. (T-05) Build the 15-test integration suite, prioritizing TEST 1
   (monotonic control), TEST 5 (formal vs moment), TEST 9 (restart), TEST
   10 (failed Newton rollback), TEST 11 (failed refresh maturity), TEST 14
   (small-delta/large-residual no-exit). [F-23]
6. (T-06) Add the runtime instrumentation (sweeps, refreshes, BA
   builds/reuses, JFNK evals, rejections, freeze events, wall time, memory)
   before proposing any new acceleration. [F-22/F-24]
7. (T-07) Audit transient guards for bounded influence on committed
   iterates (activation-rate already logged; add a committed-iterate
   census). [F-19]
