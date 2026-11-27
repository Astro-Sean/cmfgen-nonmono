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

## WP-2 .. WP-7

Not yet worked this session. Owed in checklist order after the
WP-0 production harvest: WP-1 non-monotone accept (`JFNK_NM_ACCEPT`,
`nm_credit`, `n_reject` for discarded LS accepts),
WP-2 graded stagnation (`JFNK_STAG_POLICY`, GMRES budget
decoupling), WP-3 budget-exhaustion signal (`JFNK_EVAL_BUDGET_HIT`),
WP-4 FD-operator fixed linear map (`JFNK_H_FREEZE`/`JFNK_H_FORM`,
lintest A/B), WP-5 tau-drift preconditioner refresh
(`JFNK_PC_TAU_FAC`), WP-6 merit/operator consistency
(`dir_deriv_F` vs `dir_deriv_G` telemetry now lands via the
WP-0 ledger columns), WP-7 admission-gate report (histogram +
`BA_LAST_RAW` age, report-only).

## Verdict

- WP-0: PARTIAL->VERIFIED pending production harvest (unit-verified;
  the checklist's harvest-before-logic-changes rule applies to
  WP-1+, not to this instrumentation WP).
- WP-1..WP-7: OWED.
