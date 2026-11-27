# JFNK acceptance-path audit

Working the checklist `checklists/jfnk.md` (WP-0 .. WP-7) on branch
`jfnk-acceptance-audit`, baseline `a9b0b022` (+19 commits at work
start, HEAD `e2f9eb40`). Scope question: JFNK admissions rarely
commit Newton corrections - mostly info=3 exits, or best iterates
that never beat the 1e-4 commit floor in jfnk_solve_corr_v1.

Per-WP entries follow the checklist report format:
Question / Evidence plan / Measurement / Verdict / Tests /
Residual risk.

## WP-0: instrumentation (27-Nov-2026)

**Question.** Can every printed solver disposition be tied to the
branch actually taken, and can the eval budget be attributed by
call site? Before WP-0, `step_status` was computed and printed
before the trust-region decision block ran, so a row could claim
`ACCEPTED (enlarge TR)` for a step the chain then rejected or
re-routed through the line search / Cauchy fallback; ared/pred/rho
described the unscaled candidate, not the alpha-scaled step.

**What landed.**

- `step_status` removed; `jfnk_led_branch` is now set inside each
  TR branch (`nm_armijo_decrease`, `ls_accept_noTR`,
  `nm_increase_rejected`, `full_newton_noTR`, `accept_previous`,
  `reject_pred_negligible`, `reject_no_decrease`,
  `retry_smaller_tr`, `accept_full_newton`, `retry_larger_tr`,
  `accept_at_maxdel` - the previously unlabelled implicit accept
  when the radius is already at `mxdl`), at each reject/cycle
  site (`reject_runaway_step`, `reject_nan_inf`) and at every
  terminal/abort site. The `step disposition` print moved behind
  the decision and prints the taken label.
- On Armijo accept, `J_s`, `dir_deriv` (F-basis and the new
  PTC-folded G-basis `jfnk_led_dirg`), `norm_Js` are rescaled by
  `alpha_ls` and `F_pred`, `pred`, `rho`, `cos_angle` recomputed -
  the printed ratios now describe the step actually taken, with
  `alpha` printed alongside. Cauchy accepts are marked
  `cauchy_used` and print `non-Newton direction -- rho not
  meaningful` instead of a Newton-model ratio.
- `JFNK_LEDGER`: CSV, one row per decided Newton pass plus one
  terminal row per activation exit. Columns: `nits, fw, eta, gres,
  ref_norm, gmres_its, gmres_restarts, gmres_info, hess_res,
  true_res, true_res_ratio, dir_deriv_F, dir_deriv_G, cos_angle,
  alpha_ls, ls_steps, ls_accepted, nm_phi_ref, nm_phi_alpha,
  cauchy_used, del_before, del_after, n_reject, taken_branch,
  exit_reason, n_res_evals, wall_s`.
- Eval accounting (`jfnk_ev_*`/`jfnk_tm_*`): per-activation counts
  of entry/trial/LS/Cauchy `f()` calls, wall seconds per category,
  count of budget-sentinel (`all >= 1e19`) returns, LS evals
  wasted on TR-discarded accepts, reconciled against
  `JFNK_N_RESIDUAL_EVALS` (solver-side sum vs interface total; the
  difference is evals spent inside JVPs and `sub()`). Printed once
  per activation by `jfnk_eval_account_v1`.
- Wall clock is `jfnk_wtime_v1` (`system_clock`), not
  `omp_get_wtime`: the solver object stays linkable without
  `-qopenmp` (tests/jfnk_core links it bare).

**Tests.** `tests/jfnk_core/test_step_status.f90` parses the
ledger after a solve on the standard stiff family and asserts:
header contract; every `taken_branch` in the known label set;
nits nondecreasing with accept rows carrying the pre-increment
count and total matching the `sub()`-counted accepts; terminal
labels only in the last row; `0 < alpha <= 1` on `ls_accepted`
rows; the `n_res_evals` column reconciling exactly with both the
driver-side f-call count and the solver-side `ev_*` sum. The
driver's `f` increments the stub `JFNK_N_RESIDUAL_EVALS`, making
the reconciliation a real check rather than a tautology.

**Residual risk / owed.** The production admission-log harvest
(rule of engagement 3) is still owed - one live `DO_JFNK`
admission on m15 through `run_model_grid.py`, harvesting
`JFNK_LEDGER` + the accounting block, before WP-1+ verdicts cite
production behaviour. Unit-level coverage is green; no solver
behaviour changed (instrumentation only - same decision tree,
same eval counts).

## WP-1: non-monotone line-search acceptance (27-Nov-2026)

**Question.** The GLL line search accepts an alpha with |F| still
rising, then the trust-region chain discards it whenever del /= 0 -
and del == 0 is unreachable (JFNK_DEL=0 maps to -1 = auto TR at
newton_hook_v1:1009). Can a GLL-passed increase be committed
safely?

**What landed.**

- `JFNK_NM_ACCEPT` (default F): the `nm_increase_rejected` branch
  is split. With the key set and `jfnk_nm_credit_left > 0`, a
  GLL-passed increase commits (`nm_armijo_increase`): iterate
  advances, |F| reaches the nm ring via the normal commit push,
  del is left unchanged (never grown on an increase), one credit
  is spent. The credit refills to `JFNK_NM_CREDIT` (default 2,
  clamped 0..16) on any ared>0 commit; exhaustion falls back to
  the old reject path. Standard GLL safeguard.
- Bug fix regardless of the key: the discard branch now also
  counts `n_reject`, so it is bounded at n_reject_max like every
  other reject instead of only by `jfnk_del < mndl` (~70 halvings).
  It also gained the shared too_many_rejected_steps abort.
- The stale-history cap phi_ref <= 1.21*phi_0 (+10% in norm) stays
  where it was (solver:1843); the new branch documents it.
- Li & Fukushima citation corrected to Optim. Methods Softw. 13,
  181-201 (was Appl. Math. Comput.) at the nm comment block and
  rd_control_variables.f.
- Key plumbing: CVM/RCV/jfnk_solve_corr_v1 sync/drv_vadat echo/
  test_vadat_keys checks/NONMONO_KEYS row.

**Alternatives evaluated (WP-1.4).** A `JFNK_GLOB=LS|TR|HYBRID`
pure-LS globalisation was not taken: the LS ladder itself is not
the inconsistency - the discard was. The nm-accept repair keeps
the TR for all other branches and is the smaller diff; HYBRID =
current behaviour remains the only mode.

**Tests.** tests/jfnk_core/test_nm_accept.f90: homogeneous problem
F = x*r^-0.9 (exact J and J^-1), Newton overshoots 10x so the
alpha ladder probes +25%/+15% fails and a +4.1% GLL pass at
alpha=0.25 once the window holds an older larger norm. Run A
(default F): 18 nm_increase_rejected rows, 0 increase commits,
info=3 - the before picture. Run B (NM_ACCEPT=T, credit 2): 17
nm_armijo_increase commits, all LS trials, max consecutive run
2 (bound holds), tr_rejected delta equals genuine reject rows
only, iterate advanced 40 times. W3-equivalent (a GLL-accepted
increase survives to sub()) is B5 here. Full suite green incl.
test_step_status with the new label.

**Residual risk / owed.** Flag defaults F - production behaviour
is unchanged until a live DO_JFNK A/B (validation job) shows LS
evals converting into accepted iterates; promote the default to T
only after that harvest. Credit semantics bound consecutive
increases, not total fraction; a model that alternates
increase/decrease could still ride the cap indefinitely, by GLL
design.

## WP-2: graded Krylov-stagnation response (28-Nov-2026)

**Claim.** A `ginfo == 1` GMRES exit aborts Newton unless BOTH
pred_res and the measured true_res beat 0.9*ref_norm
(solver:1349-1371). Band [0.9,1.0)*ref directions - still a
residual decrease - are discarded; a stale-preconditioner stall
(true_res >= ref or unverifiable) aborts with no re-solve; and
the preconditioner-refresh flag is only consumed inside sub() on
accepted steps, so refresh requests raised during rejected passes
starve.

**Actions.**

- WP-2.1 (new key `JFNK_STAG_POLICY`, 'ABORT'|'GRADED', default
  'ABORT'): two flag-gated arms between the 0.9-tier accept and
  the abort. (a) measured true_res < ref_norm -> keep the inexact
  direction (ginfo=0); the dir_deriv > 0 half of the checklist
  gate is enforced by the existing descent gate (:1603-1638),
  which retries/hooksteps a non-descent direction before any
  trial eval. (b) true_res >= ref_norm or unverifiable -> consume
  `jfnk_need_precond_refresh` on the rejection path: call
  JFNK_FACTOR_PRECOND_V1(.TRUE.) directly, tighten gres 10x,
  re-arm the iteration budget, clear the stale jvp/pc flags
  (they would otherwise kill the retry's first matvec), and
  `cycle descent_retry`. After stag_retries = 2 the abort stands.
- WP-2.2 (budget decoupling): `JFNK_GMRES_RST` absent-key default
  is now ceil(JFNK_MAXGMRES/JFNK_M)-1, so the effective cap
  min(MAXGMRES,(RST+1)*M) equals the advertised budget (was a
  flat 1 -> 20 matvecs under 50/10, silently below the advertised
  50). Explicit key values win unchanged. This is a deliberate
  default change the checklist explicitly sanctions
  ("either default RST ... or print a warning"); documented in
  NONMONO_KEYS and the RCV comment.
- WP-2.3: `JFNK_GMRES_STAG` already exposed jfnk_gmres_stag_tol
  as a VADAT key (RCV:2019, sync solve_corr:408) - no new key
  needed. The "stall persists over a full restart cycle"
  refinement is documented as a follow-up consideration, not
  implemented.
- WP-2.4: covered by arm (b) - the refresh flag now has a
  consumption point on the rejection path via the direct
  JFNK_FACTOR_PRECOND_V1 call, not only inside sub().

**Tests.** tests/jfnk_core/test_stag_graded.f90: linear problem
F = (I+20S)(x-1), S the cyclic permutation - spectrum is a
circle, so a starved solve (m=3, rst=0, budget 2m) exits on the
restart cap at ~0.99*ref, squarely in the [0.9,1) band. Run A
(ABORT): aborts on the first in-band stall, nits=1, NFACTOR=0.
Run B (GRADED): arm (a) keeps the direction - 8 accepted iterates
vs 1, no refresh spent. Run C (GRADED, matvec fails on every
exit-measure call): arm (b) fires - NFACTOR=2 exactly, then the
abort stands (info=3). A stale jfnk_jvp_failed flag poisoned the
retry's first matvec before the flags were reset - arm (b) now
clears them.

**Regression note.** tests/jfnk_core/test_guards.f90 G14 fails
under -O3 builds ('trust_region_too_small' vs expected
'too_many_rejected_steps') on BOTH the HEAD and the WP-2 solver -
a pre-existing optimization-sensitivity boundary between the
pred_negligible and ared<=0 branches, not a WP-2 regression.
-O2 builds of either pass. Owed: pin the G14 surrogate to a
numerically robust arrangement.

**Residual risk / owed.** GRADED is default-off pending the
production A/B. Arm (b) spends one BA refactorization per retry -
expensive on a real model, which is why retries cap at 2. The
true_res < ref band relies on the exit measurement being
trustworthy under the FD matvec (WP-4 territory).

## WP-3 budget-exhaustion graceful stop

Land (unit-verified, all four sub-items):

- 3.1/3.2 `JFNK_EVAL_BUDGET_HIT` (new `JFNK_CMFGEN_INT_MOD` flag)
  is set in `JFNK_RESIDUAL_V1` at both sentinel sites (eval cap and
  wall cap), alongside a cached `JFNK_EVAL_CAP` refreshed at
  `JFNK_INIT_V1` and per eval. The solver tests the flag after
  every eval that can return the sentinel - entry eval, post-GMRES,
  trial eval, LS ladder, Cauchy sweep - and exits `info=2`
  `budget_stop`, keeping the last accepted (or entry) iterate.
  The driver maps `info=2` to the existing best-iterate commit
  (`jfnk_solve_corr_v1.f:2368`), so the old route - sentinel
  |Y|~1e22 landing in the NaN/contraction guards as a mismeasured
  trial - is gone. LS/Cauchy sites exit the loop keeping the best
  already-accepted trial.
- 3.3 Affordability gate at admission
  (`jfnk_solve_corr_v1.f` pre-factorization): when
  `JFNK_EVAL_CAP > 0` and the remaining budget cannot cover
  `m+2` evals (one restart cycle + entry + trial), the admission
  takes the existing `info=-2` skip instead of burning the LU.
- 3.4 `[JFNK_LS_MAX_EVALS_FRAC]` (real, default 0.25) caps the
  LS eval spend per Newton pass at a fraction of the remaining
  budget; 1.0 restores the unbounded ladder.

Test: `tests/jfnk_core/test_budget_stop.f90` - the residual
callback trips the flag + 1D20 sentinel after N calls. Run A
(mid-GMRES hit) and run B (entry-eval sentinel) both exit
`info=2 budget_stop` with a real iterate, not the sentinel.
jfnk_core 6/6 PASS in the new test, zero regressions in the
suite; test_vadat_keys 143/143.

**Residual risk / owed.** `JFNK_EVAL_BUDGET_HIT` stays set once
tripped (the budget does not grow), so the flag is one-shot per
activation - matches the sentinel semantics. `JFNK_SEC_HIT` now
also feeds the flag, so wall-clock exhaustion gets the same
graceful exit instead of sentinel-trials. Production A/B owed.

## WP-4: FD operator fixed linear map (28-Nov-2026)

The measured hess_res/true_res 150-450x split is the per-direction
FD step cap making the matvec a different operator every Arnoldi
step. Two flags, both default-off (stock h formula and cap path
unchanged):

- `JFNK_H_FREEZE` (F): the clip cap is computed from the first
  Arnoldi direction of each GMRES solve and frozen for the rest;
  gmres_m_v1 re-arms at solve entry (`jfnk_h_frozen_active`,
  `jfnk_h_freeze_nv`). Later matvecs whose fresh cap would have
  shrunk h are counted and logged once per solve. Frozen h makes
  the FD map exactly homogeneous within the solve -- the WP-4.2
  "fixed linear map" arm.
- `JFNK_H_FORM` ('INFNORM'): 'TWO_NORM' is the Knoll & Keyes 2004
  (Sec 3.2) h = eps*(1+|x|_2)/|v|_2; 'COMPONENTWISE' takes the
  scalar min over active components of eps*max(|x_i|,1)/|v_i|
  (a true vector h cannot act in a single directional
  difference, so the min is the conservative scalar bound).

WP-4.1 lintest probes (existing JFNK_LINTEST/JFNK_HSCAN machinery)
and the WP-4.4 BA-JVP A/B are production-harvest items: the A/B
verdict is hess_res vs true_res convergence at restart, read off
the ledger columns on the validation run.

## WP-5: preconditioner lag under PTC (28-Nov-2026)

`JFNK_PC_TAU_FAC` (real, default 4, 0=off). JFNK_FACTOR_PRECOND_V1
snapshots tau and x at factorization; per Newton pass, when
max(tau/tau_f, tau_f/tau) or the relative inf-norm x drift exceeds
the factor, `jfnk_need_precond_refresh` is raised and consumed by
the commit-path refold (diagnostic v1). Unit-verified both arms in
test_ptc_sign T1/T2 (white-box: set tau_at_factor / ptc_dtau past
the factor with the GMRES-its trigger disabled to isolate the
drift path).

## WP-6: merit/operator consistency under PTC (28-Nov-2026)

`JFNK_DESCENT_BASIS` (F|G|BOTH, default F). The gate tests the
folded <F,J_G s> (measured pre-strip at the SER anchor) when 'G',
either positive when 'BOTH'; the Armijo slope and Cauchy merit
follow the selected derivative. Telemetry (WP-6.1) is the
ledger's dir_deriv_F vs dir_deriv_G columns, live since WP-0.
WP-6.3 alpha-sweep diagnostic: on the first two Newton passes the
true ||F(x - a s)|| is evaluated along an alpha ladder and logged
against the model-space sweep in jfnk_newton_diagnostic_v1.

## Production harvest: the before picture (28-Nov-2026)

The D5 before/after table needs a baseline measured on production
admissions. Four finished m15 2.82d runs carry pre-WP binaries and
print per-admission `JFNK final diagnostics` blocks plus (in the
coupled arm) per-attempt solver exits; the pre-WP ITER_LEDGER also
carries per-iteration jfnk counters. Harvest script:
`docs/audits/jfnk_acceptance/jfnk_harvest.py`.

| run | admissions | Newton ran | skip (info=-2) | FD-guard defer | gmres_stagnation_abort | commits |
|---|---|---|---|---|---|---|
| `VAL_sob_jfnk_eps1em3` | 237 | 0 | 237 | 153 | 0 | 0 |
| `VAL_sob_jfnk_cplsob_eps1em3` | 237 | 133 | 104 | 130 | 133 | 0 |
| `VAL_sob_jfnk_ngoff` | 297 | 0 | 297 | 180 | 0 | 0 |
| `VAL_sob_jfnk_ngshadow` | 259 | 0 | 259 | 163 | 0 | 0 |
| total | 1030 | 133 | 897 | 626 | 133 | 0 |

The symptom is worse than "rarely accepted": in 1030 production
admissions not a single Newton correction committed. Two failure
classes split the population:

- **Frozen-tier reach/guard skip (897 admissions, info=-2).** The
  FD-image descent guard (`JFNK_FD_COS_MIN=0.05` + `JFNK_FD_AUTODEF`)
  vetoes the admission after the probe JVPs; the trial direction is
  certified a non-descent step before GMRES runs. Across the four
  runs cos(F,A_FDz) is < 0.05 on 43-54% of probed admissions and
  < 0 on 35-47%. Note cos(BAz,A_FDz) = 1.000 on every probe in every
  run: the FD image is exactly parallel to the stored-BA image, so
  the guard is not catching an FD-vs-BA disagreement - both
  operators agree the preconditioned direction does not descend F.
  This is the WP-6 regime (residual orthogonal to the reachable
  subspace), not a bug in the guard.
- **Coupled-attempt stagnation (133 admissions, all in cplsob).**
  The coupled tier runs real Newton: 96363 JVPs, 649 PC calls,
  32597 s of Newton time. Every attempt exits
  `gmres_stagnation_abort` - the restart-cap ginfo=1 path of WP-2,
  with the per-direction FD step cap inflating the FD operator
  ~300x over stored BA on the pinned rows
  (|diff|/|anal| = 303 on JVP dir 1; the h-sweep plateaus at
  ~18x only for h_eff > 1e-17, i.e. production h is deep in the
  cancellation regime). This is the WP-4 territory.

**A/B arms in flight.** `m15_2p82d_VAL_sob_jfnk` (job 13774528,
git_head a02daa7a, all flags at default): the WP-0/WP-7 telemetry
arm - JFNK_LEDGER + eval accounting + GATE_REC once it passes
LAMBDA. `m15_2p82d_VAL_sob_jfnk_cplsob_ab` (job 13778399, same
build): warm restart of the cplsob arm's own it666 checkpoint,
identical deck plus `JFNK_NM_ACCEPT=T`, `JFNK_STAG_POLICY=GRADED`,
`JFNK_H_FREEZE=T`, `JFNK_LINTEST=T` - the WP-1/2/4 after-arm on the
only configuration that reaches real Newton attempts.

## WP-7: admission-gate policy instrumentation (28-Nov-2026)

Report-only per the checklist; the data harvest rides the
validation run. New `JFNK_GATE_REC` line emitted once per
solve_for_pops call whenever USE_JFNK_SOLVER: great-iteration,
gated demand (BA_LAST_RAW basis), both gates, FRZ/CPL hit
counters, effective NREQ, contraction flag, BA_LAST_RAW age
(new CVM var BA_LAST_RAW_IT, stamped where BA_LAST_RAW refreshes
at the non-LAMBDA SOLVEBA exit), LAM_STABLE_PROVEN, resolved
tier (frz/cpl/-) and a one-word defer class (warmup, lam_es,
lambda, gs_freeze, rdlam_probe, dem_over, dem_rise, hits_lt_nreq,
frz_hold, trace_phase). The histogram and the 7.4 relative-vs-
absolute recommendation are owed in the validation report; no
defaults changed.

## Verdict

- WP-0: PARTIAL->VERIFIED pending production harvest (unit-verified;
  the checklist's harvest-before-logic-changes rule applies to
  WP-1+, not to this instrumentation WP).
- WP-1: PARTIAL->VERIFIED pending production harvest (unit-verified
  both arms; flag default F).
- WP-2: PARTIAL->VERIFIED pending production harvest (unit-verified
  all three arms incl. rejection-path refresh consumption;
  JFNK_STAG_POLICY default ABORT). JFNK_GMRES_RST default is now
  derived - a sanctioned default change; the advertised budget is
  now reachable. Known marginal: test_stag_graded B2 (and the
  earlier G14) flip under -O3 rebuilds of the same semantics -
  HEAD sources fail identically, i.e. pre-existing flag
  sensitivity, not a WP regression.
- WP-3: PARTIAL->VERIFIED pending production harvest (unit-verified:
  sentinel sites flag JFNK_EVAL_BUDGET_HIT, all five eval sites
  exit info=2 budget_stop, admission affordability gate, LS frac
  cap at default 0.25).
- WP-4: PARTIAL->VERIFIED pending production harvest (machinery
  committed; the hess_res/true_res A/B needs the live admission).
- WP-5: PARTIAL->VERIFIED pending production harvest (unit-verified
  both drift arms).
- WP-6: PARTIAL->VERIFIED pending production harvest (telemetry
  columns live since WP-0; basis gate + alpha sweep committed).
- WP-7: PARTIAL (instrumentation committed; histogram +
  recommendation owed to the validation-run report).

## Could not verify (D6, 28-Nov-2026)

- Per-Newton-step "before" detail: the pre-WP builds have no
  per-pass ledger (that is what WP-0 added). The before column of
  the D5 table is therefore admission-level (attempts, evals, JVPs,
  exit reason, commits, wall time), not LS-eval/Cauchy-eval level.
- `JFNK_H_FORM` additivity comparison: WP-4.3 asks for the defect
  under all three forms on one production admission. The cplsob_ab
  arm measures the frozen INFNORM map via JFNK_LINTEST;
  TWO_NORM/COMPONENTWISE need dedicated arms (one flag per run) or
  offline replay. Owed.
- Whether the FD-guard defers (897 of 1030 admissions) were the
  right call physically: cos(BAz,A_FDz)=1.000 on every probe means
  the guard never observed an FD-vs-BA disagreement - it vetoed on
  the residual itself being orthogonal to the reachable subspace.
  Whether a different (non-preconditioned) Newton direction could
  descend there is a WP-6 question the dir_deriv_F/dir_deriv_G
  ledger columns will answer on the telemetry arm.
- test_stag_graded B2 / test_guards G14 -O3 fragility: pre-existing
  flag sensitivity on HEAD sources, not a WP regression, but the
  surrogate needs a numerically robust arrangement before the
  suite is -O3-clean. Owed.
