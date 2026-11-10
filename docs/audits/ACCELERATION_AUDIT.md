# ACCELERATION AUDIT - DO_NONMONO_SOB great-iteration cost

Audit date: 10-Nov-2026. Source revision audited: `280c4a44`
(main, post fail-safe allocation fixes). Runtime evidence: job
12689655 `m15_2p82d_VAL_sob_jfnk` (the production recipe,
`dessart_m15_2p82d_SOB_l14_lte_np.yaml`), ITER_HISTORY 670 rows,
OUTGEN 143670 lines, wall 16:15:51 on 16 CPUs; plus
convergence_issues.md records for the LAM_ES=F partner
(job 12316397, cancelled it742) and the SOB-esali run (it1185,
EPS_TERM 1e-3).

## 1. Actual runtime solver pathway (Phase 1)

VADAT parsing (`rd_control_variables.f`):

- `DO_NONMONO_SOB=T`, `USE_LAM_ES=T` (explicit in VADAT; also the
  absent-key auto-set since 10-Nov-2026, lines ~3689-3696).
- `DO_KNOP_REFRESH` default expression (lines 3305-3307):
  `DO_NONMONO_KNOP .OR. (SOB .AND. .NOT.USE_LAM_ES) .OR.
  (SOB .AND. DO_SOB_FG_KNOP)` -> all three FALSE on this recipe
  -> `DO_KNOP_REFRESH=F`. The refresh machinery is armed but inert.
- `DO_JFNK=T`: admitted (requires SOB|KNOP + USE_LAM_ES=T;
  USE_LAM_ES=F is a hard stop in CHECK_PARAM_CONSISTENCY).
- `JFNK_COUPLED=F`, `JFNK_COUPLED_SOB=F`: absent keys. Under SOB a
  `JFNK_COUPLED=T` without `JFNK_COUPLED_SOB=T` is reset to F with a
  warning (lines 1586-1595).
- `DO_NG=T`, `BEG_NG=10`, `IBEG_NG=15`, `BW_NG=40`, `ITS/NG=15`;
  `NG_ABSTAIN=T` (SOB absent-key auto-set).
- `DO_AV=T`, `ITS/AV=8`, `NOSC_AV=4` (oscillation averaging, upstream).
- `DO_BA_RES_ACC=F` - explicit VADAT override; the SOB auto-set would
  have made it T.
- `DO_ES_ALI=F`, `DO_ES_ALI_BA=F`, `DO_FG13_ALO_TRI=F`: the three
  ES-scattering accelerators all OFF in this config.
- `DO_FG13_ALO=T`: auto-set under SOB+LAM_ES=T - the V13 diagonal
  OAB operator feeds the continuum dJ/dCHI,dJ/dETA diagonals.
- `SOL_METH=DIAG`, `NUM_LAM=2`, `LAM_ENTER=400`, `LAM_EXIT=5`.

Execution pathway per great iteration:

1. COMP_J_BLANK sweeps the frequency grid (8331 freqs): for each
   frequency `FG_J_CMF_V13` formal solve on clamped SIGMA, then the
   `RJ=TC` Lambda update (NUM_LAM passes). This is the ONLY
   radiation-field fixed-point map on this path - there is no
   moment-equation solve and no lagged J_NEXT store.
2. Sobolev line rates via `SUB_SOB_LINE_V3`/`SOBEW` on `SIGMA_TRUE`.
3. `SOLVE_FOR_POPS`: BA matrix assemble-or-recall (N_FIX_BA=2,
   FIX_BA=10 -> 420 of 669 non-LAMBDA rows ran with BA HELD),
   correction applied via trust-region gates; JFNK admission wrapper
   around the same frozen-J residual.
4. Scheduler: LAMBDA/FULL alternation, NG/AV/R-grid auxiliary writes.

**The two-level AA machinery in `MOM_J_CMF_NONMONO_V1` is not on the
production path.** It only runs under `USE_LAM_ES=F` (opt-in
experiment) or `DO_NONMONO_KNOP`. The audit questions about
MOM_RF_AA/MOM_AA_OUT therefore apply to the experiment arm, not to
the 722-iteration production run.

## 2. Where the time and iterations go (runtime evidence)

- 722 iterations, 16.27 h -> ~81 s/iteration.
- ITER_LEDGER on sibling runs (FRESH_sob, SN_A SOB archive):
  `t_rad ~ 78-80 s/iter` ~ 97% of t_iter. The formal solve over
  8331 frequencies IS the iteration cost. Population-solve time is
  a small remainder; refresh time is ~0 on this path.
- Mode histogram: 247 LAMBDA_TFIX + 422 FULL_TFIX; only 2 REJECTs.
- Convergence history: maxraw decays ~60 -> 5 over its 1-220
  (LAMBDA, ~0.97-0.99/it contraction), mid-run rebounds (maxraw
  pinned at 1E7 -> capped 100) around its 224-334, then a
  FULL_TFIX tail 0.5 -> 0.05 with repeated oscillations
  (0.15->0.62->0.10->0.81->0.09...). EPS_TERM=5e-2 reached at it722.

Accelerator census for the whole run:

| accelerator | attempts | committed | evidence |
|---|---|---|---|
| NG (DO_NG) | ~10 | 0 | every attempt: NUM_BAD_NG>=4 -> "Too many bad NG accelerations - cancelled" |
| JFNK (frozen tier) | 322 admissions | 0 | 322/322 "admission gate skipped Newton" (reach/census/FD-descent) |
| JFNK-CPL (coupled) | 0 | 0 | never armed (JFNK_COUPLED=F) |
| AV averaging | several | yes | AVE_FLIPS 6253 early, decaying; damping not acceleration |
| BA_RES_ACC | 0 | 0 | DO_BA_RES_ACC=F |
| MOM two-level AA | n/a | n/a | solver not on path |
| ES_ALI / ES_ALI_BA / FG13_ALO_TRI | 0 | 0 | all F in VADAT |

**Every accelerator armed on this recipe either never fired or was
structurally unable to help.** The 722 iterations are the bare
operator-split contraction.

## 3. Why the outer map contracts slowly (Phases 4 + 6)

Regime classification: F (outer population/SE iteration) dominated,
with a specific structure - the residual is *radiation-pinned*:

- The frozen-J STEQ residual |F| sits in rows whose radiative rates
  cannot move without re-solving J. Documented in the code itself
  (jfnk_solve_corr_v1.f:1040-1053): "|BArow|/|F| ~ 1e8-1e9 on the
  top rows, ~39k pop-floor regularized rows, cos(F,Jz) ~ 1e-4".
- Runtime confirms: REACH frac negative-huge, |Jz|/|F| ~ 1e9-1e13,
  cos(F,Jz) sign-flipping; the FD descent check also fails
  (JFNK_COS_MIN/FdCosMin). A Newton step on populations alone can
  only burn evals producing rejected ~0 steps - hence all 322
  skips. This is correct gating, not a JFNK defect.
- JFNK admission cost is still real: ~33 s/admission overhead on
  the audited era's build (~4.5 h over the run for zero committed
  steps, per the in-code audit comment; somewhat cheaper now via
  the census gate). Config-side, DO_JFNK=T on this recipe buys
  overhead only.

Why the residual is radiation-pinned: the continuum J is
ES-dominated (m15 2.82d dense shell). The `RJ=TC` Lambda update
converges an ES-scattering field at ~(1-albedo) per pass; NUM_LAM=2
passes/iteration leave J substantially lagged each call, so the
populations always see a stale field - the classic ALI slow mode
(Rybicki & Hummer 1991). The documented fixed-ratio decay is
~0.97/it on this model (convergence_issues.md, adp2 series), which
matches the observed maxraw decay rate.

NG failure mode: the banded secant fit needs 4 consecutive clean
linearly-contractive SCRTEMP records (NG_CLEAN). Under this map the
records carry oscillation/rebound structure; at >=4 depths the
extrapolation trips the LOCINC>10.1 / LOCDEC<0.09 / T<0 vetoes even
with NG_ABSTAIN holding per-component bad cases. NG_CLEAN itself
was not the blocker here (JFNK never committed, so JFNK_COMMITTED
never reset it); the NG fit quality was.

GS_FREEZE: forced F outside DO_NONMONO_KNOP - not holding anything
back on this path.

## 4. Two-level Anderson audit (Phase 2) - MOM/LAM_ES=F arm

Source state (mom_j_cmf_nonmono_v1.f):

- `MOM_RF_AA=T` (default), `MOM_AA_M=3` inner depth;
  `MOM_AA_OUT=T`, `MOM_AA_OUT_M=8` outer depth. Flags wired via
  SET_MOM_RF_AA / SET_MOM_KNOP_REFRESH bounds (5002-5004).
- Inner AA: Anderson on the alternating-direction two-pass
  composite map, applied at even pass boundaries (5946ff);
  weighted-residual Gram, non-affine cell exclusion
  (MOM_AA_GCAP), Tikhonov ridge, affine-subset-empty skip,
  ||corr||>||f|| rejection, per-cell 10*|f| cap, H dragged by the
  J ratio. History restarts each refresh call (the map changes
  with the source rebuild - correct).
- Outer AA: Anderson on the committed lagged J_NEXT map across
  great iterations (6544ff), same exclusion/rejection family,
  contractive-retry on leading m-1 windows.
- Memory: outer history = ND_NM*NCF*(2*MOM_AA_OUT_M+4)*8B.
  At ND=70, NCF~4-8k, M=8: ~22-45 MB - realistic; not an
  allocation risk on this grid. MOM_AA_OUT_FAIL records the triple
  and prints once on failure.

Runtime state on the failed LAM_ES=F run (job 12316397, recorded in
convergence_issues.md 22-Oct (d)):

- Both AA levels were active and evaluating every refresh.
- Inner AA: mostly "correction exceeds ||f||" - the safeguard
  correctly rejects extrapolating a non-affine map.
- Outer AA: rejected on every depth.
- Root cause is not accelerator deficiency: the lagged map is a
  bounded period-2 limit cycle pinned at the `1+sigma*K ~ 0`
  resonance band (~15K cells/pass on the GAM floor, ~200K
  quarantined/sweep ~11% of the field, freq-1 blue-edge J at the
  noise floor). AA on a non-contractive map can at best do nothing;
  the safeguards held, which is the correct behavior.

AA visibility: per-refresh counters print ('inner Aitken
acceleration: applied X, skipped Y, solve failures Z' + outer line
incl. allocation-failure flag). What OUTGEN/ITER_HISTORY does NOT
currently expose per iteration: AA-active/alloc-fail/reset state,
NG attempted/deferred/fired, JFNK committed-vs-admitted. The ledger
extension in ACCELERATION_TODO (D1) closes that gap.

## 5. Frequency-space coupling (Phase 3)

Measured on the actual model data (SN_HYDRO_DATA 70-pt R,V,Sigma;
CONT_FREQ 8331 freqs, median dLOG_NU ~1.0e-3 ~ 300 km/s):

- Characteristic drift per cell: |GAM|*dr = (V/(cR))*|1+sigma*K|*dr.
  Median ~ 0.001-0.003 bins; maximum ~1.3 bins at K=1/3 (2 cells),
  ~3.9 bins at the K=1 bound (9 cells, depths ~60-68).
- **Conclusion: the one-bin upwind approximation is NOT the
  convergence bottleneck on m15 2.82d.** Displacement >1 bin occurs
  in <13% of cells under the pessimistic bound, and those cells are
  not the documented residual argmax (which sits at the GAM-floor
  resonance band, where |1+sigma*K|~0 makes the drift SMALL).
- Caveat: CONT_FREQ contains zero-spacing bins (dLOG_NU min = 0,
  edge doublets); matched-span code uses real log(nu) once the
  cache fills, and MOM_DLN_MIN floors the denominator.
- MOM_FREQ_MATCH (auto-set T under LAM_ES=F) is the right
  discretization improvement where the drift DOES exceed a bin
  (Ceulemans et al. 2025 regime) - keep it, but it is not the m15
  bottleneck.
- MOM_RF_GS_DESC: default F; the frozen-read one-bin-per-pass
  propagation limit applies to the MOM refresh arm only. Since
  this path already fails on non-contraction (not slow
  propagation), a controlled GS_DESC comparison belongs to the
  LAM_ES=F revival experiment, not the production recipe.

## 6. JFNK coupled-SOB residual (Phase 5)

- Frozen tier: solves the same fixed-J system as SOLVEBA_V13 -
  consistent by construction, and certified unreachable by the
  radiation-pinned structure (322/322 skips).
- Coupled tier (`JFNK_COUPLED + JFNK_COUPLED_SOB`): each residual
  eval replays the formal solve AND `SUB_SOB_LINE_V3` against the
  perturbed state - the full dJ/dn coupling (Arramy et al. 2024).
  The A3 consistency gate enforces |F_eval-F_prod|_w <=
  JFNK_CPL_SOB_TOL=1e-5 before Newton may differentiate it, with a
  top-8 |diff| row dump (a systematic line-rate omission would show
  as a concentrated signature; re-solve truncation reads diffuse).
  Fires only as the reach-fail fallback, eval-budgeted via
  JFNK_CPL_SOB_MAX_EVALS. A failed enforced check detaches the
  callback permanently.
- The consistency machinery is adequate and never armed on this
  recipe (flags F). It is the only Newton tier that can touch the
  radiation-lag residual. Cost per eval ~ one full frequency sweep
  (~78 s), so an admission can cost 10-30 sweeps - expensive but
  bounded by the eval cap.

## 7. Findings summary

F1. The production SOB recipe is a pure operator-split Lambda
    iteration for the radiation field plus a held-BA/occasional-JFNK
    population correction. No accelerator committed a step in 722
    iterations.
F2. Wall time is ~97% formal frequency solve. Iteration COUNT is
    the only meaningful lever; per-iteration cost is fixed.
F3. The dominant slow mode is the ES-scattering J lag (documented
    ~0.97/it fixed ratio), which leaves the STEQ residual
    radiation-pinned and defeats every population-only accelerator.
F4. JFNK frozen tier provably cannot help (reach/census correct);
    it costs ~33 s/admission. On this recipe DO_JFNK=T is pure
    overhead unless the coupled tier is enabled.
F5. NG never gets a usable window; the map is not locally
    contractive. Not a scheduling defect.
F6. MOM two-level AA is implemented correctly with the right
    safeguards; it ran and correctly rejected on the LAM_ES=F arm,
    whose lagged map is a limit cycle (physics/discretization
    defect at the resonance band + blue-edge noise floor, not an
    accelerator problem).
F7. Frequency-drift >1 bin is rare on m15 (<=~9 cells); not the
    bottleneck. FREQ_MATCH remains the correct guard where it does
    occur.
F8. The accelerators designed for exactly this slow mode
    (DO_ES_ALI, DO_ES_ALI_BA, DO_FG13_ALO_TRI, DO_BA_RES_ACC) are
    all OFF in the production config.

## 8. Recommended experiments (Phase 7 matrix, single-change arms)

All arms: same m15 2.82d input, warm from the converged
`adp2_tauto` checkpoint where possible (a cold LTE start spends
~300 LAMBDA iterations re-deriving what a checkpoint already
knows), same EPS_TERM=5e-2, 16 CPUs. Record the full ITER_LEDGER +
ITER_HISTORY + wall time set.

- E1 (baseline warm): np recipe as-is, warm start.
- E2: E1 + DO_ES_ALI=T (ES_ALI_TOL=1e-4, ES_ALI_MAX=20). Verified
  mechanism for the ES J-mode; measures its effect on the outer
  count.
- E3: E2 + DO_ES_ALI_BA=T. Targets the documented ~0.97/it BA
  under-correction - the direct outer-map rate.
- E4: E3 + DO_BA_RES_ACC=T (the auto-set value; np explicitly
  disabled it).
- E5: E1 + JFNK_COUPLED=T + JFNK_COUPLED_SOB=T +
  JFNK_CPL_SOB_MAX_EVALS=20 (bounded probe). Tests whether the
  coupled Newton can buy outer iterations at ~1 sweep/eval.
- E6 (mom arm, separate question): LAM_ES=F + FREQ_MATCH + the
  AA census columns - only to characterize the limit cycle, not as
  a speed proposal.
- E7 (housekeeping): DO_JFNK=F on the np recipe recovers the
  ~3-4.5 h of inert admissions IF E5 shows the coupled tier cannot
  commit either.

Success criterion: same converged solution (populations within
EPS_TERM envelope, unchanged spectrum) at materially fewer great
iterations AND lower wall time. An arm that trades iterations for
more wall time is not a win (E5's sweep-cost evals are the risk).

## 9. Limitations

- The LAM_ES=F failure data is from the 22-Oct run record, not a
  re-run on this revision; the AA verdict for that arm relies on
  the documented counters.
- JFNK-CPL has zero runtime evidence anywhere; E5 is genuinely
  exploratory.
- The ~33 s/admission JFNK overhead figure predates the census
  gate; current cost is lower but the 322-skip outcome is unchanged
  in kind.
- K=1/3 vs K=1 bounds the drift diagnostic; K_ON_J is Eddington
  and depth-dependent. The conclusion (<1 bin almost everywhere)
  is insensitive to that bound.

## 10. D1 diagnostics landed (10-Nov-2026)

The audit's first deliverable is in place so every experiment arm
(E1-E7) reports accelerator state per call, not just per run:

- `ITER_LEDGER` column contract v2 appends after `rss_mb`:
  `ng_try,ng_ok,ng_def,ng_cln,jfnk_cmt,aa_app,aa_skip,aa_sng,aa_rst,
  aao_app,aao_rej,aao_sng,aao_fail,freqm`. New files get the v2
  header; `LEDG_ACTIVE` still gates the whole file off on stock
  paths. (`iter_ledger_mod.f90` LEDG_WRITE_ITER.)
- `ITER_SUM` gains an `accel` row: `ng:-|d|r|A` (due/committed/
  deferred/rejected states), `jfnk:-|a|A` (admitted vs committed
  Newton), `aa=off` when the MOM refresh did not run, else
  `aa=app/skip/sng/rst o=app/rej/sng` plus `AOFAIL` (outer-AA
  history allocation latch) and `freqm` (MOM_FREQ_MATCH armed).
  (`solve_for_pops.f`; same GET_MOM_AA_STATE feed as the ledger.)
- `MOM_AA_RST` counts damp-triggered inner-AA history restarts
  inside a refresh - the missing Phase-2 item 10 evidence
  ("how often is the history reset"), previously only visible in
  LUER text. New file-scope getter `GET_MOM_AA_STATE` exposes the
  armed flags plus last-refresh counters without module access.
  (`mom_j_cmf_nonmono_v1.f`.)

Tests: `test_iter_ledger.py` 31/31 (v2 header, per-column values,
committed-vs-admitted JFNK distinction, aao_fail/freqm flags,
stock-path silence). `solve_for_pops.o`, `iter_ledger_mod.o`,
`mom_j_cmf_nonmono_v1.o` compile clean into the archives.

What this makes visible in the next run: inner AA applied vs
skipped vs singular per refresh, outer AA applied vs rejected per
iteration, history-restart frequency, outer-history allocation
failure (previously only a one-line LUER message at init), NG due
vs deferred vs committed vs rejected per call, JFNK admitted vs
committed per call, and FREQ_MATCH arming - i.e. Phase-2 items
1-12 and the Phase-4 NG questions are answerable from the ledger
alone on the next np/LAM_ES=F run.

## 11. D1a diagnostic corrections (10-Nov-2026)

(Committed as `1a21f41b`..`c77114b5`; labelled "D2" in the first
draft of this section and in the pre-reconciliation AUDIT_LOG -
renamed to D1a so "D2" consistently means the outstanding OUTGEN
summary task in ACCELERATION_TODO.md.)

ChatGPT review of the D1 instrumentation found two reporting
defects; both corrected. No solver physics or algorithm change.

**Stale AA counters.** The `MOM_AA_*`/`MOM_AA_OUT_*` tallies reset
inside `MOM_KNOP_REFRESH_V1` *after* its self-guard early returns
(line ~5325 vs guards at ~5258). A call where the refresh
guarded out (or the formal-J path where it never runs) therefore
kept reporting the previous refresh's tallies. Fix: the getter
gates all seven activity counters on `MOM_RF_RAN` - the flag
reset at mom-solver entry and latched only after every guard
passes - so `aa_*=0,aa_ran=0` now reads 'no refresh this call'
rather than 'AA applied nothing'. Enabled flags and the
`MOM_AA_OUT_FAIL` latch are configuration state and are reported
regardless; the latch is never cleared by a non-running refresh.

**Enabled vs applied state.** `GET_MOM_AA_STATE` now returns
`RF_RAN` plus the live `AA_ON`/`AA_OUT_ON` enables. Ledger
contract v3 appends `aa_ran,aa_en,aao_en` after `freqm` (all v1/v2
columns keep their positions). The `ITER_SUM accel` row now reads
`rf=y|n aa:off|on[ app/skip/sng/rst | -] oa:off|on[ app/rej/sng |
-][+AOF] fm=on|off`, so 'enabled but idle', 'disabled' and
'alloc-failed' are three distinguishable states.

**Outer-AA gate decoupled (Task 3).** The outer AA block was
gated `MOM_RF_AA .AND. MOM_AA_OUT`. Source inspection shows zero
inner-AA state inside the block - it runs entirely on
XHIST_OUT/FHIST_OUT/RETJ_PREV/RETH_PREV against the committed
lagged-store map, a fixed point distinct from the inner two-pass
map - so the coupling was convenience, not a numerical
requirement. The gate is now `MOM_AA_OUT` alone. Defaults are
unchanged (both T), so default-path behaviour is identical; only
a `MOM_RF_AA=F, MOM_AA_OUT=T` config changes (it now receives the
outer acceleration it requested). VADAT help updated.

**Tests.** New `test_mom_aa_state.py` + `drv_mom_aa_state.f90`
exercise the production getter against real module state (links
the same production object set as test_knop_ray): 62/62 covering
cases A-G - ran+enabled, ran+inner-off, ran+outer-off,
guarded-out with stale counters (zeros), never-ran formal-J,
alloc-fail-without-refresh (AFAIL=T, counters 0, enables T), and
enabled-but-zero-applied. `test_iter_ledger` 31/31 with the v3
columns; `test_knop_ray` 41/41 unchanged.

**Refresh-state lifecycle.** `MOM_RF_RAN`: reset `.FALSE.` at
`MOM_J_CMF_NONMONO_V1` entry (~1755), set `.TRUE.` at ~6815 after
all refresh guards pass. Under `USE_LAM_ES=F` the solver runs
every call, so the flag is per-call accurate. Under `USE_LAM_ES=T`
the solver never runs, the flag stays F, and the getter can never
report moment-AA activity on the formal-J path - acceptance
criterion 7 by construction.

**Remaining limitation.** The counters are gated on read, not
cleared on the guarded-out path itself (chosen to leave the
refresh's internal bookkeeping untouched); any future consumer
must use `GET_MOM_AA_STATE`, not the module counters directly.
