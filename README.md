# CMFGEN with non-monotonic velocity fields and a Newton-Krylov population solve

Development fork of [CMFGEN](https://sites.pitt.edu/~hillier/web/CMFGEN.html)
(baseline `cur_cmf_18jun25`, `pitt-phyast/cmfsource`) for supernova
ejecta where the velocity law is not monotonic. The case we care about
is a Type IIn: a fast ejecta shell hits slow, dense circumstellar
material, and a reverse shock stalls the shell into a cold dense
shell. There `sigma = d ln V / d ln r` drops below -1, which the stock
code assumes cannot happen, so the transfer inside the compression
zone returns negative or NaN intensities.

The non-monotonic and Newton-Krylov solver paths are opt-in. With
their flags unset, they do not replace the stock solver. This tree
also has general correctness, build, and observation-frame fixes, so
it is not a byte-for-byte copy of upstream.


## The problem

CMFGEN alternates between a formal solution of the radiation field
and an update of the level populations, electron density, and
temperature. Two linear systems carry the iteration:

- The **transfer matrix**. In the comoving frame the intensity along a
  ray obeys an equation containing a frequency-shift term
  `gamma * dI/d ln nu`, with `gamma = (V/r)(1 + sigma*mu^2)`.
  Discretised on the frequency grid, this gives a matrix that is
  tridiagonal in frequency at each depth point. In homologous
  expansion `gamma > 0` everywhere, so each frequency couples only to
  the previous (bluer) one and the system sweeps outward in a fixed
  order. Lines in the expansion opacity can also be treated in the
  Sobolev approximation: each line meets the radiation field at
  exactly one resonance point per ray, where the optical depth is
  proportional to `1/(1 + sigma*mu^2)`.

- The **rate matrix** (the BA matrix in the code). At each depth, the
  statistical-equilibrium equations for the level populations plus the
  charge and energy equations form a dense `NT x NT` block. The stock
  update alternates fixed-point Lambda iterations with a Newton step
  on this block (`SOLVEBA`), whose largest demanded correction feeds
  the convergence metric.

Both break in an interacting supernova. When the flow decelerates into
a dense shell, `sigma` falls below -1 in the compression zone. The
comoving-frame frequency then folds back along a ray: one line
resonates at several points, `gamma` changes sign, and the upwind
direction of the frequency coupling reverses. A scheme hard-wired for
expansion pulls intensity from the wrong (uncomputed) neighbour, and
the Sobolev optical depth goes through a singularity at `gamma = 0`.
The second failure is in the population update: on a stalled model
the first full Newton step asks for order 1e7 per cent corrections,
mostly on physically empty trace levels that no step can satisfy.

## What the solver does

![Sobolev vs non-monotonic transfer](docs/knop_solver.png)


The formal solver follows Knop, Hauschildt & Baron (2009): the sign of
`gamma` picks which side counts as upwind. Expanding points couple to
the previous frequency, compressing points to the next one. The next
frequency has not been computed yet on the downward sweep, so its
intensity is stored and lagged by one iteration, kept only where
`sigma < -1`. This keeps the modified opacity positive and removes the
negative intensities.

A second, lighter path applies the Dessart, Audit & Hillier (2015)
Sobolev treatment directly in the compression zone: the unclamped
`sigma` drives the line opacity, and the optical-depth denominator is
floored so it cannot cross zero. Both paths are off by default.

## Iteration cycle

![Iteration cycle](docs/iteration_cycle.png)


Each pass is one model iteration: solve the radiation field, update
the populations, write a checkpoint, repeat until the populations stop
changing. Grey boxes are stock CMFGEN, blue are new code, amber are
existing routines with new switches. Everything new is off by default.

The new parts:

- The non-monotonic formal solver. In compression zones a point's
  intensity couples to last iteration's stored value instead of the
  next frequency.
- A refresh loop that re-solves the radiation field until the lagged
  values stop moving, so the populations are never updated against a
  stale field.
- A mode controller that keeps the model on safe fixed-point steps
  until the lagged field settles, then switches to the faster Newton
  or Newton-Krylov step.
- Checkpoints of the lagged fields, so a restart picks up where it
  stopped.

## Updated approach: demand floors and the trace polish

![Trace-polish flow](docs/stage_b_flow.png)


On interacting-supernova models many levels are physically empty at
some depths (a few 1e-20 or less of the local gas). They keep
demanding large relative corrections that no solver can satisfy -
cascade coupling refills them between solves - so the convergence
metric never drops below the exit threshold and the model stalls at
the finish.

This tree treats those rows as trace species and handles them in two
stages, both off by default:

1. **Demand floors.** An absolute population floor and a relative
   floor (fraction of the largest population at the depth) exempt
   trace rows from the convergence metric. Their corrections still
   commit - the floors only gate the metric. The largest exempted
   demand and the largest committed trace change are reported every
   call, on both the stock and Newton-Krylov solve paths, so hidden
   demand is never silent.

2. **Trace polish sweeps.** When the major species converge but
   exempted rows are still moving, the run arms a polish phase and
   loops for a bounded number of extra iterations instead of exiting.
   Each sweep takes a direct matrix step, commits only the exempted
   rows, and holds the majors, electron density, and temperature
   fixed. Trace species are one-way coupled - they absorb and emit but
   do not feed back on the field or the majors - so solving them at a
   fixed major state is the lower triangle of a block-triangular
   system and cannot oscillate against the majors.

The run exits when the trace motion falls below its threshold or the
sweep budget runs out; a spent budget never holds the run. What
remains is the unsatisfiable dead-stage class the floors were built
for.

## Code layout

| Path | Contents |
|---|---|
| `new_main/` | Main program and the VADAT input readers (`mod_subs/rd_control_variables.f` registers every new key; `mod_subs/control_variable_mod.f` declares them). |
| `new_main/subs/` | New-solver support routines (`recompute_sigma_true.f`, SOLVEBA variants, stage-B polish drivers). |
| `plane/` | The Knop formal solver `fg_j_cmf_nonmono_v1.f`: sign-aware frequency coupling, the lagged `I_NEXT`/`J_NEXT` maps, and their checkpoint/restore. |
| `subs/`, `newsubs/`, `com/`, `disp/` | Upstream CMFGEN source tree, carried with surgical edits only. |
| `obs/` | Observed-frame tools: `cmf_flux` / `obs_frame` sources and the `rd_cmf_flux_controls.f` key reader. |
| `unit_tests/` | Parser and solver-adjacent test drivers; `run_all_tests.py` runs the suite. |
| `docs/` | README figures and their generators plus the audit notes under `docs/audits/`. |
| `knop_solver_files/`, `jfnk_solver_files/` | Working copies of the files each feature branch touched, kept as a per-feature index; the build uses the tree sources, not these. |
| `exe/`, `mod/`, `lib/` | Build outputs (executables, `.mod` files, libraries). `exe_dbg/`, `mod_dbg/`, `lib_dbg/` are the debug-build equivalents. |

The exit test itself gates on the demanded correction, not the applied
step, so a capped or held update cannot release the exit while demand
remains. The test is evaluated before the final-iteration branch, so a
run whose iteration budget ends exactly on a converged state still
gets the convergence check and a report of which exit it took - clean
convergence, a saturated exempted demand, or a spent polish budget -
rather than writing out untested.

## Newton-Krylov population solve

The stock update alternates fixed-point iterations with a Newton step
on the depth-block rate matrix. On a stalled non-monotonic model that
step overshoots. The opt-in Krylov path treats statistical equilibrium
as `F(n) = 0` and solves it with GMRES - finite-difference
Jacobian-vector products, the rate matrix as preconditioner only
(Knoll & Keyes 2004; Arramy et al. 2024). Options allow a frozen or
recomputed radiation field during the residual evaluations, and the
lagged non-monotonic intensity store is held fixed for the duration of
each evaluation.

Each activation is bounded by iteration and restart caps plus two
budgets: `JFNK_MAX_EVALS` (residual evaluations) and `JFNK_MAX_SEC`
(wall-clock seconds per activation). Both default to 0, meaning off;
when `DO_JFNK` is set with no wall-clock budget the input reader
prints a warning that a stalled activation could consume the job.
Exhausting a budget forces the rejection path, so the iteration falls
back to the stock update rather than aborting.

## VADAT keys for the trace-exemption and polish machinery

All of these are opt-in and only take effect under
`DO_NONMONO_SOB`/`DO_NONMONO_KNOP` (or `DO_JFNK` where stated). With
the flags off none of this code runs and the defaults are inert.

| Key | Default | Meaning |
|---|---|---|
| `BA_DEM_PFLOOR` | -1 (auto) | Absolute population floor exempting trace rows from the demand metric. `-1` resolves to `1.0E-30` under the non-monotonic flags; `0` disables; `>0` is used as set. |
| `BA_DEM_RELFRAC` | -1 (auto) | Relative demand floor, as a fraction of the largest population row at the depth. `-1` resolves to `0.01` under the non-monotonic flags; `0` disables. |
| `BA_TR_MAXSWP` | 3 | Bound on Stage-B trace-polish sweeps. `0` disables the polish entirely; the reader no longer rewrites `0` to another value. |
| `BA_TR_ABSMAX` | 1.0E-08 | Trace convergence threshold on `BA_TR_ABS` (max committed `|dn|/POP_ATOM` over exempted rows). Under the non-monotonic flags an explicit `0` is auto-set to `1.0E-08`. |
| `KN_ADAPT_RELAX` | T | Adaptive write-back weight for the lagged-field refresh maps (`INEXT_RELAX`/`JNEXT_RELAX` adapt to the measured contraction). `F` keeps the fixed `0.5` weight. |
| `INEXT_RELAX` | 0.5 | Seed write-back weight for the lagged-field refresh; the adaptive omega moves it from this seed. Lower values (e.g. 0.3) damp oscillating restarts more strongly. |
| `KN_LUM_TRACE` | F | First-detection scan of `STEQ_T` after each energy-equation accumulation site; reports the site, frequency, depth, and operands at the first non-finite entry. Diagnostic only. |
| `KN_REPAIR_T` | F | Rebuild non-finite/non-positive `T` restored from a checkpoint at iteration entry, from the nearest finite depth. Lets a crash-poisoned `SCRTEMP` self-heal instead of stopping mid-sweep. |
| `GS_RAW_MATURE` | F | `T` requires the raw FG sweep max change to pass `GS_MATURE_TOL` for a sweep to count toward the `GS_REL_STABLE` release streak -- a bulk-pass sweep with a high max norm resets the count instead. Recommended for KNOP runs. |
| `JFNK_CPL_RELWAIT` | F | `T` defers the coupled JFNK tier while a call sits in the post-release damp window (`N_REL_SINCE <= GS_REL_N`); the call falls to the frozen tier or SOLVEBA. Note: a very large `GS_REL_N` defers the coupled tier indefinitely. |
| `JFNK_REL_REACH_MAX` | 5 | Consecutive reach-gate skips after which the SOB stable-LAMBDA relief stands down (stops relaxing `JFNK_GATE_NREQ`). `0` disables the stand-down. |
| `JFNK_RUNG_MAX` | 4 | Budget-ladder top rung (replaces the retired 3-strike pause/re-arm): counted failures raise the rung, a commit resets it; Newton is attempted on every admissible call at every rung. |
| `JFNK_RUNG_MIN_NITS` | 4 | Floor on Newton iterations per activation under the rung trim. |
| `JFNK_RUNG_MIN_GMRES` | 15 | Floor on GMRES iterations per Newton step under the rung trim. |
| `JFNK_RUNG_MIN_EVALS` | 12 | Floor on residual evals per activation under the rung trim. |
| `GS_RELEASE_LIM` | 0 | Rejected GS trial releases (`N_TRIALS_REJ`) before the run stops (KNOP); `0` disables. |
| `JFNK_FD_COS_MIN` | 0.05 | Descent floor on the finite-difference image of the clipped step, `cos(F,A_FD z)` -- the stored-BA REACH can pass a step the true Jacobian maps anti-parallel to F (KNOP_jfnkrel it35: REACH 0.735 vs FD cos -0.23). `0` disables. |
| `JFNK_FD_AUTODEF` | T | Defer any JFNK admission whose FD operator image is anti-descent (`cos(F,A_FD z) < 0`), even with all thresholds at 0. Arms the precheck block (one extra residual eval per admission). `F` recovers pre-guard behaviour. |
| `LUM_DIV_FAC` | 1.0E2 | KNOP only: per-call outer-luminosity growth factor counted as runaway. |
| `LUM_DIV_N` | 3 | KNOP only: consecutive runaway calls before `LUM_BLOCKED` fires (pops held, damped LAMBDA only). `0` disables the guard. |
| `LUM_DIV_FAC2` | 5.0 | KNOP only: lower per-call growth factor for the slow-runaway streak. |
| `LUM_DIV_N2` | 4 | KNOP only: consecutive slow-growth calls before `LUM_BLOCKED`. |
| `LUM_DIV_REC` | 1.0 | KNOP only: once latched, `LUM_BLOCKED` holds until L_out falls below `REC` times the pre-runaway level (`0` = never release). |
| `LUM_DIV_FAC3` | 1.005 | KNOP only: per-call growth factor for the drift tier (sustained small-factor growth). Off by default. |
| `LUM_DIV_N3` | 0 | KNOP only: consecutive drift-growth calls before `LUM_BLOCKED`; `0` disables the tier. |
| `LUM_DIV_CAP` | 1.0E15 | KNOP only: outer luminosity (L_sun) past which the run stops as non-convergent. The growth/cap checks are skipped during the `N_LAM_RWARM` warm-up window; the non-finite check stays armed. |
| `N_LAM_RWARM` | 4 | KNOP only: forced LAMBDA iterations at process start so no BA step runs against a stale (zeroed) `I_NEXT`; `0` disables. |
| `GS_FOREIGN_SEED` | F | KNOP only: on restart, accept a fingerprint-valid but tag-mismatched `FGINEXT_GS` (a checkpoint from a different run of the same model) as the lagged-field seed instead of cold-rebuilding `I_NEXT=0`. The seed stays immature (`I_NEXT_GS_OK=F`) so `N_LAM_RWARM` re-converges the map before any BA step. |

SOB-path acceleration probes (`DO_NONMONO_SOB`; all default OFF, so
the stock Lambda update is unchanged):

| Key | Default | Meaning |
|---|---|---|
| `DO_SOB_FG_KNOP` | F | Use the Knop sign-aware formal solver for the SOB+`USE_LAM_ES` J instead of `FG_J_CMF_V13`. |
| `DO_SOB_CLAMPED_LINES` | F | Clamp Sobolev line rates to the infrastructure `SIGMA` (consistency test only; biased physics). |
| `DO_FG13_ALO` | F (auto T under SOB+`USE_LAM_ES=T`) | `FG_J_CMF_V13` stores its diagonal approximate-Lambda operator for the continuum Jacobian. |
| `DO_ES_ALI` | F | Olson-Auer-Buchler Jacobi update of the coherent-ES Lambda iterate inside the `USE_LAM_ES` J loop. |
| `ES_ALI_TOL` | 1.0E-4 | `DO_ES_ALI` relative dJ/J convergence target per frequency. |
| `ES_ALI_MAX` | 20 | Cap on extra J passes per frequency under `DO_ES_ALI`. |
| `DO_ES_ALI_BA` | F | Amplify the stored ALO by the converged-ES Jacobi factor so the BA linearization sees the scattering-coupled J response. |
| `ES_ALI_BA_MIN` | 0.1 | Floor on the `DO_ES_ALI_BA` amplification denominator. |
| `DO_FG13_ALO_TRI` | F | Fill `TX(K,K+/-1)` with the stored off-diagonal Lambda* entries (effective only with `NUM_BNDS>=3`). |

Preconditioner-audit keys (`jfnk-precond-fixes`; all default OFF so
the default factorization path is unchanged; details in
`docs/audits/JFNK_PRECOND_AUDIT.md`):

| Key | Default | Meaning |
|---|---|---|
| `JFNK_PC_SIGN_DIAG` | F | Logging-only report at each factorization: fraction of rows with `d<0` and `sigma>0`, the `sigma/|d|` min/median/p99, and the near-cancellation count `|d+sigma| < 0.1|d|`. |
| `JFNK_PTC_SIGNED` | F | Sign-matched PTC shift: residual term, preconditioner fold, apply rescale, and ridge all use `sign(d_i)*sigma_i` so `|d_i+shift| = |d_i|+sigma_i` on loss-dominated rows. |
| `JFNK_REG_ZERO_COL` | F | Regularized (near-null/pop-floor) rows also get their column zeroed, so the factored block is the clean projected operator (per-depth and banded paths). |
| `JFNK_ROW_TINY_EQUIL` | F | Scale-free near-null row test `share_i = max_J |a_iJ|/colmax_J < JFNK_ROW_TINY` replacing the raw rowmax/blockmax test; logs reclassification counts per factorization. |
| `JFNK_PC_RCOND_GATE` | F | Replace the `x=1` forward-error micro-test by an explicit `rcond >= JFNK_PC_RCOND_MIN` acceptance threshold and cap the ridge ladder at `JFNK_PC_ITRY_CAP` attempts per depth. |
| `JFNK_PC_RCOND_MIN` | 1.0D-13 | rcond acceptance threshold under `JFNK_PC_RCOND_GATE`. |
| `JFNK_PC_ITRY_CAP` | 3 | Maximum factorization attempts per depth under `JFNK_PC_RCOND_GATE` (historical ladder allows up to 9). |
| `JFNK_PC_PAR_FACTOR` | F | Declared for the depth-parallel factorization option; pending the Phase-4 profiling outcome. |

Observation-frame (`cmf_flux`/`obs_frame`) keys:

| Key | Default | Meaning |
|---|---|---|
| `ALLOW_SIGMA_UNDEF` | F | `T` warns instead of aborting on cells with an undefined `SIGMA_TRUE`. |
| `SOB_SIGMA_MAX` | -1 (off) | If `>0`, clamp `|SIGMA_TRUE|` to this bound in the OBS estimator. Off keeps the same unclamped `SIGMA_TRUE` as the population solve. |
| `DO_SOB_DEP_WINGS` | F | Continue the deposited Sobolev line profile past 3 Doppler widths as Voigt damping wings out to `SOB_DEP_WING_EXT`; matches the full CMF solve, which carries wings the pure-Gaussian deposit drops. |
| `SOB_DEP_WING_EXT` | 2000 km/s | Wing extent for `DO_SOB_DEP_WINGS`. |
| `DO_SOB_DEP_FLOOR` | F | Keep the emissivity of masing/negative-opacity deposit cells with their opacity floored at `SOB_DEP_CHI_FLOOR`, instead of dropping them (they carry the P-Cygni wing flux). |
| `SOB_DEP_CHI_FLOOR` | 1.0 | Line opacity floor for masing deposit cells under `DO_SOB_DEP_FLOOR`. |

## Validation status

The Newton-Krylov path is off by default. Admission waits for a warmup
window, a mature and stable radiation field, and contracting matrix
demand; it is not triggered by a single low correction. A clean build
or a generated spectrum alone does not establish that the new paths
converge or produce a valid interacting-supernova spectrum.

## Numerical additions to CMFGEN

A short list of what this tree adds numerically; all are off by
default unless stated:

- Sign-aware comoving-frame frequency coupling for arbitrary velocity
  laws (Knop et al. 2009), so `sigma < -1` compression zones no
  longer produce negative/NaN intensities.
- Floored Sobolev optical-depth denominator in the compression zone
  (Dessart et al. 2015), removing the `1 + sigma*mu^2 -> 0`
  singularity.
- Lagged-field refresh sweeps and Gauss-Seidel checkpointing, so the
  stored `I_NEXT`/`J_NEXT` maps survive requeues and warm restarts.
- Jacobian-free Newton-Krylov population update (GMRES with the BA
  matrix as preconditioner only) with residual-based acceptance,
  admission gating, and bounded budgets.
- Trace-species demand floors and Stage-B polish sweeps, so
  physically empty levels stop stalling the convergence exit.
- ES-scattering ALI acceleration: Jacobi updates on the stored formal
  J diagonal plus BA-side amplification and optional tridiagonal
  Lambda* couplings.
- Luminosity-runaway guards (fast, slow, and drift growth tiers plus
  a hard cap) that hold the populations instead of letting a
  diverging field commit.
- Residual-acceptance rollback of failed Newton steps, per-equation
  quarantine, and per-depth-column rollback on partial failure.
- Observed-frame improvements: non-monotonic integration bounds,
  deposited-line Voigt wings, and emissivity retention at masing
  cells.

## Changes vs upstream

Solver behaviour is gated where stated; the tree also contains general
correctness and build fixes that are not controlled by solver flags:

- Non-monotonic Sobolev line treatment in compression zones, after
  Dessart, Audit & Hillier (2015). The unclamped velocity gradient is
  carried in the restart files and recomputed on restart.
- The Knop et al. (2009) formal solver with the lagged-intensity
  machinery, and a moment-equation variant for the same fields.
- A Jacobian-free Newton-Krylov population solve that replaces the
  stock Newton step on admitted iterations.
- Residual-based acceptance for Newton steps, including rollback of
  only the failed depth columns when the bulk residual remains
  acceptable, and a bound on repeated rejected steps.
- Per-equation quarantine for equation/depth pairs that fail
  repeatedly; fixed-point iterations still update those populations.
- Watchdog, checkpoint, and restart handling for the lagged-field
  state across job requeues.
- Rejection of inconsistent line-treatment settings by the input
  readers instead of running with mixed transfer assumptions.
- General fixes including the non-monotonic observed-frame
  integration bounds and headless build support.
- Robustness hardening on the Knop formal-solver path (degenerate
  ray/frequency spacings, short-ray end slopes, deallocation
  status checks, OpenMP private hygiene); itemized status in
  `docs/audits/KNOP_FORMAL_SOLVER_AUDIT.md`.

The `unit_tests/` suite (`unit_tests/run_all_tests.py`) compiles the
production parsers and solver-adjacent sources against small drivers;
`test_vadat_keys.py` covers the non-monotonic flag contract and the
JFNK VADAT keys.

## Building and running

```bash
module load intel mkl
./build.sh          # production build: main code + spectrum tools
./build_debug.sh    # bounds-checked debug build of the main code
```

The production build produces the model executable plus the flux and
observed-frame tools. All new behaviour is switched off unless the
corresponding input-file options are set; see `AGENTS.md` for the
build constraints and verification checklist.

## References

- Knop, Hauschildt & Baron 2009, *A&A* 496, 295,
  "A new formal solution of the radiative transfer in arbitrary
  velocity fields" - the `chi_hat` formulation.
  [DOI](https://doi.org/10.1051/0004-6361:200810290) |
  [ADS](https://ui.adsabs.harvard.edu/abs/2009A%26A...496..295K) |
  [arXiv:0902.1876](https://arxiv.org/abs/0902.1876)
- Knop, Hauschildt & Baron 2009, *A&A* 501, 813,
  "Comoving-frame radiative transfer in arbitrary velocity fields -
  II. Large scale applications" - the Gauss-Seidel solvers.
  [DOI](https://doi.org/10.1051/0004-6361/200911793) |
  [ADS](https://ui.adsabs.harvard.edu/abs/2009A%26A...501..813K)
- Dessart, Audit & Hillier 2015, *MNRAS* 449, 4304,
  "Numerical simulations of superluminous supernovae of type IIn" -
  the interaction model the Sobolev branch follows.
  [DOI](https://doi.org/10.1093/mnras/stv609) |
  [ADS](https://ui.adsabs.harvard.edu/abs/2015MNRAS.449.4304D) |
  [arXiv:1503.05463](https://arxiv.org/abs/1503.05463)
- Feldmeier & Nikutta 2006, *A&A* 446, 661,
  "Nonlocal radiative coupling in non monotonic stellar winds" -
  three-point coupling across velocity kinks.
  [DOI](https://doi.org/10.1051/0004-6361:20053380) |
  [ADS](https://ui.adsabs.harvard.edu/abs/2006A%26A...446..661F)
- Knoll & Keyes 2004, *J. Comput. Phys.* 193, 357,
  "Jacobian-free Newton-Krylov methods: a survey of approaches and
  applications" - the JFNK method and pseudo-transient continuation.
  [DOI](https://doi.org/10.1016/j.jcp.2003.08.010)
- Arramy, de la Cruz Rodriguez & Leenaarts 2024, *A&A* 690, A12,
  "Jacobian-Free Newton-Krylov method for multilevel NLTE radiative
  transfer problems" - JFNK on NLTE populations; the convergence
  requirements the solver here follows.
  [DOI](https://doi.org/10.1051/0004-6361/202449963) |
  [arXiv:2406.17234](https://arxiv.org/abs/2406.17234)
