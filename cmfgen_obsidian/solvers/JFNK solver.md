# JFNK solver

Jacobian-free Newton-Krylov solve for the population update, active
when `DO_JFNK=T` (`USE_JFNK` deprecated alias) and admitted by the
gates in [[JFNK admission]]. Works under `DO_NONMONO_KNOP` or
`DO_NONMONO_SOB`; under SOB only the frozen-J tier can be admitted -
`JFNK_COUPLED` is rejected at VADAT read because the eval sweep
(`JFNK_RAD_SWEEP_V1`) replays the KNOP continuum chain and does not
recompute `SUB_SOB_LINE_V3` line rates, so a coupled SOB residual
would be inconsistent with the outer equations.
Method follows Knoll & Keyes (2004, JCP 193, 357,
[DOI](https://doi.org/10.1016/j.jcp.2003.08.010)) and Arramy, de la
Cruz Rodriguez & Leenaarts (2024, A&A 690, A12,
[arXiv:2406.17234](https://arxiv.org/abs/2406.17234)).
Summary doc: [[JFNK_SOLVER]].

> [!quote] Arramy et al. 2024 (abstract)
> "the JFNK method can converge faster than our reference case by up
> to a factor 2" - measured in formal-solution evaluations vs their
> Rybicki & Hummer operator-splitting reference. A factor-2 ceiling on
> the *cheap* part of the problem is why this tree treats Newton as a
> gated endgame accelerator, not a loop replacement.

Intuition for "Jacobian-free": instead of assembling the Jacobian
matrix (huge, and under a non-monotonic sweep partly unknowable), the
code only ever asks for its *action* on a vector,
$Jv \approx [F(x+hv)-F(x)]/h$ - like probing a hillside's slope by
nudging a marble and watching where it rolls, rather than surveying
the whole terrain.

## Files

```
new_main/mod_subs/jfnk_cmfgen_int_mod.f      module: state, VADAT keys,
                                           JFNK_SWEEP_CB callback type
new_main/mod_subs/jfnk_cmfgen_interface_v1.f residual/JVP/preconditioner:
                                           JFNK_INIT_V1, JFNK_RESIDUAL_V1,
                                           JFNK_JVP_V1, JFNK_BUILD_BA_RAD_V1,
                                           JFNK_BA_RAD_MATVEC_V1,
                                           JFNK_FACTOR_PRECOND_V1,
                                           JFNK_PRECOND_V1, audits
new_main/mod_subs/jfnk_solve_corr_v1.f       driver JFNK_SOLVE_CORR_V1:
                                           F(0) test, Newton loop, commit,
                                           fallback chain
new_main/mod_subs/jfnk_solver_v1.f90         generic Newton/GMRES core:
                                           newton_hook_v1 (hookstep TR),
                                           gmres_m_v1, hook_step_v1
new_main/cmfgen_sub.f 5047-5945              JFNK_RAD_SWEEP_V1 - the
                                           coupled-eval sweep, registered
                                           as JFNK_SWEEP_CB before the loop
```

## What is being solved

Residual F(n) = statistical-equilibrium equations at trial
populations. Two modes with ~100x per-eval cost asymmetry:

- **Frozen** (`JFNK_COUPLED_THIS=F`):
  `F(x) = STEQ_coll(n(x)) + STEQ_RAD + CMAT_RAD*g(x)` where
  `CMAT_RAD = BA_full - BA_coll` (radiative-only Jacobian, no double
  counting) and `g(x) = exp(x)-1` under `JFNK_LOGX`. Exact at the
  baseline by construction; an eval costs seconds (SE assembly only).
- **Coupled** (`JFNK_COUPLED_THIS=T`): each residual/JVP eval freezes
  `I_NEXT` (`FREEZE_I_NEXT`), re-runs the full Knop continuum sweep at
  the perturbed populations (`JFNK_RAD_SWEEP_V1` via `JFNK_SWEEP_CB`),
  and restores all state - so the finite-difference JVP probes
  `dJ/dn`. ~250-400 s per eval on m15. The lagged `I_NEXT` store
  itself stays a parameter (lagged, not differentiated through -
  differentiating it would make every matvec a second coupling).

The depth-block BA matrix is the **preconditioner** in both modes -
exact-ish for the frozen residual, only the local part of the
coupled operator.

## SOB-seeded production stage

The m15 Type IIn workflow converges the fixed-hydro-temperature SOB
model first, then starts Knop/frozen-JFNK from its population checkpoint.
`cmfgen.warm_restart_from` names the SOB run directory and
`warm_restart_population_only: true` copies only `POINT1`, `POINT2`,
`SCRTEMP`, and `MODEL`. The source and target `SN_HYDRO_DATA` files must
be byte-identical and the MODEL NT/ND fingerprint must match. Radiation
state is intentionally excluded: `RVTJ`, `EDDFACTOR`, `FGINEXT_GS`, and
`MOMJNXT_GS` encode the source transfer operator and must be rebuilt by
Knop. The matched configurations keep `FIX_T=T`, `FIX_T_AUTO=F`,
`JFNK_FIXT=T`, diffusion inner boundary, thin outer boundary, and common
microphysics. This stages the nonlinear solve without treating an SOB
radiation field as a valid Knop fixed point. For the STELLA-seeded
path the zlaw14 yamls now set `JFNK_TAUTO=T` (see below): T stays
live but damped inside JFNK while the caller's `FIX_T`/`FIX_T_AUTO`
mode is unchanged outside JFNK admissions.

## Newton formulation (all on by default for a JFNK solve)

### `JFNK_LOGX` (default T) - logarithmic unknowns

Unknown for every level, ED, T: $x = \ln(n/n_0)$ so $n = n_0 e^x$.
Removes the `1e-100` floor that made F non-differentiable and broke
forward-difference JVPs crossing it. At `x=0` the derivative matches
the old fractional unknown, so the BA preconditioner is still the
right local Jacobian; away from x=0 the matvec carries `exp(x)` and
the PC apply multiplies by `exp(-x)` - kept consistent. Committed
step still limited to `MAX_JFNK` in ratio (`|x| <= ln(MAX_JFNK)`);
inside an eval x is clipped at +/-40 against exp overflow (clip
counted, difference retried).

### `JFNK_PTC` (default T) - pseudo-transient continuation

Knoll & Keyes (2004, JCP 193, 357 Sec 4.2); switched evolution
relaxation, Kelley & Keyes (1998, SIAM J. Numer. Anal. 35, 508).
Newton sees

$$ G(x) = F(x) + \frac{(x - x_k)\, s}{\Delta\tau} $$

anchored at the current iterate $x_k$ (re-anchored each accepted
step; baseline $x=0$ is the first anchor), with
$s_i = \max(|F_i|_{\rm baseline}, \text{floor})$ matching units.
The term vanishes at $x_k$, and every residual the Newton loop reads
is stripped back to F, so norms, ared/pred, EW forcing and
convergence all act on $|F|$ (06-Oct-2026: previously anchored at
$x=0$ with $|G|$ norms, which could declare convergence at
$F = -x s/\Delta\tau \neq 0$).
$\Delta\tau$ starts at 1, is multiplied by $|F_{old}|/|F_{new}|$ on
progress (capped 1e8 => pure Newton), halved on failure (floor 0.01).

Analogy: PTC is Newton with a flywheel. Far from the root, small dtau
makes the added term dominate and the iteration behaves like a damped
march - bounded, always moves downhill-ish. As the residual falls the
flywheel (dtau) spins up until the extra term vanishes and the method
*is* Newton. You get global behaviour early and quadratic endgame
convergence for free, without a handoff between two solvers.
The preconditioner rescales the BA solve by
`d'/(d'+sigma)`, where `d'` is the signed diagonal of the applied
preconditioner, including ridge or identity-row repair, and `sigma=s/dtau`.
The scalar limit is exact for either sign: large
dtau gives the BA solve; small dtau gives the PTC diagonal solve. A
zero `d'` uses the input residual divided by `sigma`; near a singular
`d'+sigma`, the unscaled BA result is retained.

### GMRES forcing floor

Eisenstat-Walker eta floored at `3e-4` - the bottom of Arramy et al.'s
useful band (3e-4..3e-2). Eta still loosens when the residual falls
fast so a distant iterate does not pay for a tight inner solve.
EW forcing computed in preconditioned units under left
preconditioning; right preconditioning forcibly disabled (hookstep
norm invalid there).

### JVP step and the noise floor

Forward difference `J*v ~ (F(x+h v) - F(x))/h` in `JFNK_JVP_V1`.
`JFNK_EPS` is a **lower bound on the perturbation norm
`h*||v||inf` - coupled mode only**: the sweep is iterative and the
residual carries evaluation noise, so a machine-eps step can put the
finite-difference quotient below the noise floor. Clip retries do not
shrink this norm below the configured floor. The check uses `h*||v||`
rather than scalar h, which varies inversely with `||v||`. Pernice &
Walker (1998, NITSOL) prescribe scaling h to the square root of
function noise. Frozen residuals use deterministic stored-array
arithmetic and may continue shrinking on clip retries.

### Preconditioner

`JFNK_FACTOR_PRECOND_V1`: per-depth equilibrated LU of
`BA_full + ridge`; null rows -> identity; Jacobi fallback on factor
failure. Rebuilt only on baseline changes (PC refreshes deliberately
rare - Arramy's Table 2 shows a rebuild ~doubles per-call cost).

Analogy: the depth-block BA preconditioner is a *street map*: exact
for local traffic (collisions, rates within one depth) but blind to
the highways (the nonlocal dJ/dn radiation coupling in coupled mode).
Perfect for the frozen residual where there are no highways; in
coupled mode it preconditions only the local part - which is why
Janett et al. 2024 (A&A 682, A68,
[DOI](https://doi.org/10.1051/0004-6361/202348048)) use a physics-based
"lightweight model" preconditioner and get 10-20x fewer iterations.
`POPS_WORK` is synchronized with `POPS_BASE` before the rebuild.
`JFNK_USE_BANDED` + `NUM_BNDS>=3` exist but untested (the
cross-depth radiative coupling is dropped self-consistently from both
the frozen matvec and the PC).

### `JFNK_TAUTO` (default = `VARFIXT` at VADAT read) - adaptive damped T

13-Oct-2026. Replaces the hard T mask inside the JFNK solve with a
damped-live-T mode, so a STELLA T guess (not in NLTE radiative
equilibrium) can seed the run without the first free-T Newton step
driving T to zero/overflow. When `JFNK_TAUTO` and the caller leaves
`FIXED_T=F`, `JFNK_TAUTO_ACTIVE` disengages the VARFIXT depth mask at
all six mask sites - F rows, JVP in/out, PC in/out, trial vectors -
so the energy equation enters the residual and the JVP keeps its
dF/dT columns. Step acceptance does the guarding, not the residual.

Per admission (in `JFNK_SOLVE_CORR_V1`):

- `jfnk_tauto_damp_v1` post-processes the Newton step and the Cauchy
  fallback direction: T components (flat index `(k-1)*NT+NT`, domain
  `JFNK_TAU_DOM_ST..END` = VARFIXT range, or 1..ND when that range is
  empty) are scaled by `JFNK_TAU_ALPHA` and capped at
  +/-`JFNK_TAU_DMAX` in log space; the commit path's existing
  exp-transform `T = T_base*exp(x_T)` keeps T positive.
- `JFNK_TAUTO_SCALE_V1` additionally scales PC-output T components by
  `JFNK_TAU_ALPHA` when `JFNK_TAU_PCW=T` (default T), so the left-PC
  Krylov basis does not propose oversized T directions that the line
  search then rejects.
- The energy residual `|RE|` over the live T rows is measured after
  each activation and `JFNK_TAU_RE0` latches the entry value.

Ladder state machine (runs once per activation, commit-gated):

- rung 1 alpha=`JFNK_TAU_A0` (0.05); 2 -> 0.2; 3 -> 0.5; 4 -> 1.0.
- Release rung->rung+1 only when `JFNK_TAU_NLOCK` committed
  activations have accumulated AND `|RE| <= JFNK_TAU_REL*RE0`
  (default 0.1) AND rung < 4. Prints `RELEASE`.
- Revert: if the committed log-T step `eps_T > JFNK_TAU_EPS`
  (default 0.1), rung drops one and the commit streak clears.
  Prints `REVERT`.
- `info=-2` reach-gate skips do not touch the ladder.

OUTGEN diagnostic per activation:
`JFNK TAUTO: |RE|=.. |RE0|=.. rung=N alpha=.. commits=..`.

Surrogate coverage (`cmfgen_jfnk_unit_test`, `shock` spec): T12
STELLA-shock mock (RE0~3.8e7, x*_T=3; damped first step 0.496 vs 3.0
undamped; bounded march to target), T13 ladder transitions
(lock/release/revert), T14 JVP T-column stays live under AUTO.

### Budgets and fallback

- `JFNK_MAX_EVALS` (40): per-activation eval budget; exhaustion marks
  the eval bad -> GMRES abort -> rejection path -> best-iterate or
  Lambda fallback. One coupled activation at nits=10, m=10 is ~110
  evals ~ 9 h worst case - the budget is load-bearing.
  08-Oct-2026 sizing note: the frozen tier is cheap per eval but a
  single GMRES cycle can spend ~60 JVPs (m=10, up to 20 iterations,
  2 restarts) plus per-Newton residual evals and line-search retries,
  so a budget of 40 starves mid-solve and the abort lands in the
  stagnation budget. The SOB zlaw14 config runs `JFNK_MAX_EVALS=200`.
- `JFNK_MAX_RESTART`: baseline/precond rebuild attempts per call.
- Fallback chain: bad F(0) / GMRES abort / trust-region rejection /
  non-finite step -> restore best iterate or LAMBDA-mode SOLVEBA.
  `SUCCESS=F` never commits a worse state.
- `JFNK_MAX_STAGNATE=3`: repeated failures self-disable JFNK
  (03-Oct-2026 fix: `JFNK_DISABLED` now clears `DO_JFNK_NOW` - it
  previously only set the LAMBDA flags and Newton ran anyway under a
  mutated operator). 04-Oct-2026: catastrophic failures (NaN, info=4,
  F(0) eval failure) now also count toward the disable budget - the
  old reset let a deterministic catastrophic mode retry forever.
  08-Oct-2026 re-arm (commit e56471f): a disable used to persist for
  the whole process while a requeue silently re-armed it - the same
  logical state behaved differently depending on SLURM timing.
  `JFNK_REARM_EVERY=10` admitted calls now re-arms JFNK
  deterministically (`JFNK re-armed ... retrying Newton`) and clears
  the stagnation count, so the retry policy is identical within a
  process and across restarts.

### CMAT_RAD activation caching (08-Oct-2026, commit e56471f)

`JFNK_BUILD_BA_RAD_V1` costs two full BA assemblies (BA_full +
BA_coll) per activation. The fixed point is invariant to `CMAT_RAD`:
it only multiplies the correction `x`, which vanishes at the root, so
a stale `CMAT_RAD` degrades the Newton path but never the answer.
`JFNK_CMAT_AGE`/`JFNK_CMAT_MAX_AGE=4` reuse it across up to 4
activations (`JFNK: reusing CMAT_RAD (age N)` in OUTGEN); a rebuild
is forced on dims change, unallocated matrix, or age limit. The
preconditioner is unaffected - `JFNK_FACTOR_PRECOND_V1` re-derives
`JFNK_BA_DIAG`/ridge from the live BA every activation.

### Failure semantics (04-Oct-2026 numerical review)

The parallel review of the opt run's it73 stall found the failure
handling deterministic but the retries bit-identical:

- `info=3` (Krylov/JVP/preconditioner abort) now exits Newton
  immediately instead of cycling. The preconditioner clip is a pure
  function of the unchanged LU factors and the input vector - the
  trust radius `del` is only consulted after a GMRES completes, and
  the `jfnk_need_precond_refresh` flag was only consumed via `sub()`
  on accepted steps, so each "shrink del, fresh GMRES" pass replayed
  the identical abort (observed: `res=3.1503E+01` frozen, 22 wasted
  coupled sweeps on one iteration).
- The driver restart loop now exits on `info=3` for the same reason
  (re-factoring the same BA at the same point can only replay).
- `jfnk_del <= 0` now trips the `mndl` abort: a zero `snrm_tr`
  (GMRES returned no step) used to produce `jfnk_del=0`, which
  slipped the old `>0` guard and entered the `del==0`
  no-trust-region branches mid-solve.
- The `info==2` hookstep recalc now requires `jj>0`: early GMRES
  returns (DGELSY error, near-zero beta, JVP/clip abort) used to
  reach it with a stale or uninitialized `jj` and index a Krylov
  space from a different point.
- `JFNK_PRECOND_V1` zeroes fixed-unknown components (T under
  FIX_T/FIX_T range, Ne under FIX_NE) **before** the `|z|>1e8` clip
  test: DGETRS mixes all NT rows of a depth block, so a masked row
  is a genuine solve output of arbitrary magnitude that the pack-time
  mask discards anyway - it must not count as a clip event.
- The stagnation-fallback step now uses `MAX_LIN_COR` (factor ~10);
  it previously carried `MAX_LAM_COR` (factor ~1e5) while the rarer
  catastrophic path had the safe cap - inverted.

### Preconditioner conditioning (04-Oct-2026 audit)

- `JFNK_ROW_TINY` (default 1e-30, opt config now 1e-10): near-null
  per-depth rows are identity-replaced before equilibration. 1e-30
  only caught essentially-zero rows; rows with rowmax ~1e-23*blockmax
  survived, equilibrated to scales ~1e22, and produced |z|~1e8-1e22
  clips on any input. `JFNK_BA_RIDGE_FAC` is NOT the lever - it is
  only consulted on the banded retry path, not the depth-LU ladder.
- The clip test still only guards the LU path; the Jacobi-fallback
  and banded paths and the post-solve PTC rescale are unclipped
  (reported, not yet changed - masked components in those paths
  are ~0 by construction).

### Null-row projection replaces identity-row solve (08-Oct-2026)

Live zlaw14 run, it31/33/35: every GMRES entry aborted on the
`|z|>1e8` clip because ~6300 regularized rows solved `z_i=v_i` with
residual components ~1e12. Identity replacement made the clip a
certainty on rough states, so JFNK never reached a second residual
eval and stagnated out after 3 activations.

- `JFNK_ROW_REG(J,K)` now persists the factor-time regularization
  mask (LU near-null rows, banded zero/near-null rows and zero cols).
  At apply time masked components are set to zero - a linear
  projection, so `M^{-1}` stays consistent and no clip event fires.
  Those directions get no Newton update; LAMBDA/globalization still
  owns them.
- Applied in every branch: LU solve, LU-fail Jacobi fallback,
  DGETRS-fail fallback, pure-Jacobi, banded DGBTRS output, and the
  log/PTC post-processing loop (where the `|D|<tiny` branch would
  otherwise return `RHS/sigma` and resurrect a masked component).
- Rowmax detection reads the raw `BA_A` row, not the ridged
  `LU_FAC`, so a null row cannot escape regularization on ridge
  retries (ITRY>0) and invert to an oversized z.
- Diagnostic: `Reg-row projections` counter (`jfnk_n_reg_proj`)
  printed in the JFNK diagnostics block.
- Residual note: after the projection fix landed, GMRES still clipped
  because rows with rowmax between `1e-30` and `1e-10`*blockmax
  escaped regularization and inverted to `|z|>1e8` anyway.  The
  documented lever is `JFNK_ROW_TINY=1.0D-10` (opt config, now set in
  the zlaw14_jfnk yaml + live VADAT): those rows join `JFNK_ROW_REG`
  and are projected instead of inverted.

### |z|>1e8 clip removed -- finite pass-through (08-Oct-2026)

Even with `ROW_TINY=1e-10` (64729 rows projected at it64, ~half the
matrix), the forcing-norm apply still clipped: a handful of
*non-regularized* rows produced `|z|>1e8` on the ~1e12 residual.
The audit's `||z||_inf=5.8e6` is measured post-PTC rescale, so it
masked the pre-rescale clip.  On reflection the clip itself is the
bug: with null rows projected, a huge finite component is faithful
`M^{-1}` output on a badly conditioned block, and clipping it is
itself the nonlinearity the guard exists to detect.  Fix (commit
ff3933d): finite components of any magnitude pass through unclipped
in both the LU and Jacobi apply paths; only nonfinite (NaN/Inf)
output is zeroed and flagged as `jfnk_pc_clip_event`.  `|z|>ZMAX`
remains a diagnostic counter.  Downstream bounds are unaffected:
`|s|>10` step scaling, hookstep TR, line search, and the +/-40
log-clip on x.  Worst case reverts to stall->LAMBDA fallback, the
same end state as the abort but with a chance of progress.

### Stored-BAMAT recall crash on warm restart (08-Oct-2026)

zlaw14_jfnk requeued at it46 (READY=1); the new process restored the
checkpoint-stored BA matrix (`READ_BA_DATA_V3` on `SCRATCH/BAMATPNT`,
SUCCESS=T) so the it47 FULL call ran `COMPUTE_BA=F`.  The recall path
in `CMF_BLKBAND_V3` then read `SCRATCH/CSCRATCH0001` and STOPped --
job failed exit 6.

Invariant broken: "stored BAMAT implies stored LU-factor cache".
`STORE_BA_DATA_V3` writes BAMAT on any `COMPUTE_BA .AND. WRBAMAT`
non-LAMBDA call, but `WRITE_BCD_MAT` only runs inside
`CMF_BLKBAND_V3` solves with `WR_BA_INV=T`.  All 46 iterations of this
run went through LAMBDA (`WR_BA_INV=F`) or the JFNK internal
factorization, so no `*SCRATCH` factor file was ever written; the
first post-restart recall had nothing to read.  Stock runs hit a
WR_BA_INV solve by it1 so the cache always exists.

Fix in `cmfgen_sub.f` at the `.NOT. NEWMOD` read site: after
`READ_BA_DATA_V3` succeeds, `INQUIRE` the depth-1 and depth-ND files
of the type `METH_SOL` needs (`CSCRATCH` for DIAG,
`BCSCRATCH`+`DSCRATCH` for TRI); if absent, `SUCCESS=.FALSE.` so the
first call recomputes BA and rewrites the cache.  Same pattern as the
02-Oct SUCCESS-initialization fix one block up.

Recovery for the live run: `warm_restart: true` in the yaml preserves
RVTJ/POP*/POINT1/2/SCRTEMP but wipes SCRATCH/ entirely, which removes
BAMATPNT+BA_ASCI+CSCRATCH and forces `COMPUTE_BA=T` on resume --
equivalent outcome to the code fix, no checkpoint loss.

### Cross-activation state disarm (09-Oct-2026)

Two module-state leaks found by a pass over everything that survives
between `JFNK_SOLVE_CORR_V1` calls.  Both break the assumption that
each activation starts from a clean slate.

**PTC scale/anchor.** `jfnk_ptc_scale_ready`, `jfnk_ptc_anchor_ready`,
`jfnk_ptc_scale` and `jfnk_ptc_x0` persist across activations, but the
F(0) identity-test eval runs before `JFNK_SET_PTC_SCALE_V1` re-arms
them.  A stale nonzero anchor adds `(0 - x0_old)*s/dtau` to F(0):
frozen mode fails the identity test on every second and later
activation (JFNK permanently disabled), coupled mode adopts the
contaminated baseline so every FD JVP inherits `stale/h`.  A stale
`jfnk_ptc_scale` sized for a previous N (R-grid revision) is also an
out-of-bounds index.  Fix in `jfnk_cmfgen_interface_v1.f`: INIT clears
both readiness flags before the F(0) eval, and
`JFNK_REBUILD_BASELINE_V1` clears `scale_ready` before the rebuild
eval; `JFNK_SET_PTC_SCALE_V1` then re-arms with the current baseline.

**GMRES hookstep basis.** `jj` and `beta` were SAVE locals in
`gmres_m_v1` - the only state that must survive a call, for the
info==2 hookstep recompute that recombines the caller's h/v.  But
`JFNK_DIAGNOSTIC=T` runs `jfnk_gmres_dim_experiment_v1`, which calls
the same routine with m_test up to 40 (> production m=10).  The last
test solve left `jj~40` behind; a later production recompute then
indexed `y(1:jj)`, `h(1:jj+1,1:jj)`, `v(:,1:jj)` past m-sized arrays -
silent heap corruption in the trust-region step.  Fix in
`jfnk_solver_v1.f90`: the pair promoted to module state
`jfnk_gmres_jj`/`jfnk_gmres_beta`, checkpointed at each completed
Arnoldi cycle; the info==2 branch guards `jj <= m` and falls through
to a fresh solve otherwise; the dim experiment saves and restores the
pair so its test solves cannot poison production recompute.

`./build.sh`: solver tests 21/21, interface tests 9/9, harness 3/3.

## State save/restore and eval failures (order-of-operations audit, 03-Oct-2026)

Every path that calls `JFNK_COMPUTE_STEQ_V1` or an eval sweep is now
bracketed by `JFNK_SAVE_MOD_STEQ_V1` / `JFNK_RESTORE_MOD_STEQ_V1`
(shared helpers): the residual eval, `JFNK_BUILD_BA_RAD_V1`, the
frozen-mode radiative split in the driver setup, and the frozen
baseline rebuild.

05-Oct-2026 (state leaks and evaluation failure handling, fixed):

- `FREEZE_I_NEXT` was defeated inside evals - the FG solver cleared it
  on every INIT call, so each eval re-advanced the lagged I_NEXT store
  at perturbed populations. `F(x)` was eval-order-dependent; the FD
  JVP inherited a drift term that does not shrink with h (the ~1e18
  JVP vs ~1e13 BA.v signature). INIT-clear now gated on `JFNK_IN_EVAL`
  (`fg_j_cmf_nonmono_v1.f`).
- `UP_RATE_TWO`/`DOWN_RATE_TWO` accumulate (`+=`) inside
  `TWO_PHOT_RATE` and are only re-zeroed by `SET_TWO_PHOT_V3` at the
  outer-iteration head - eval k saw production leftovers plus k eval
  contributions. Now saved/zeroed/restored around the eval sweep in
  `JFNK_RAD_SWEEP_V1` (`cmfgen_sub.f`).
- `BAACC_PENDING` is disarmed at JFNK entry (`jfnk_solve_corr_v1.f`):
  adjudication of a prior armed BA step against post-JFNK residuals
  could silently roll back a JFNK commit or falsely accept a bad step.
  Its raw-demand trial is marked unresolved so it cannot reopen the
  Lambda-exit gate as if a BA residual test had accepted it.
- An applied best iterate, rebuilt baseline, or clipped Newton commit is
  not a converged JFNK solve. `JFNK_UNCONVERGED_COMMIT` prevents
  CMFGEN_SUB's MAXCH-based temperature release, staged source/solver
  ramps, and final-iteration exit from treating it as one. Only an
  unclipped `JFNK_INFO=0` commit may satisfy those outer gates.
- If JFNK is disabled after stagnation, the SOLVEBA fallback keeps the
  caller's correction cap, bounded by `MAX_LIN_COR`; it does not reset
  the cap to `MAX_LAM_COR`.
- A successful JFNK call returns applied correction in MAXCH; it no
  longer writes into `BA_LAST_RAW` (a raw-demand slot) in
  `SOLVE_FOR_POPS` - the last BA-measured demand is kept until a real
  SOLVEBA re-measures.
- `JFNK_RAD_SWEEP_V1` saves/restores the EDDINGTON selector and
  `ACCESS_F` even when no continuum cursor was captured. `FOLD` is
  cold-started inside each trial sweep, then restored so trial iteration
  state cannot seed the production sweep. The two two-photon accumulators
  are independently guarded. `COMP_J_BLANK` converts its
  EDDFACTOR/ES_J_CONV record mismatch and mean-intensity blowup STOPs
  into `JFNK_EVAL_BROKEN` during trial evals. The frequency loop exits
  immediately; the sweep epilogue marks the species residual invalid,
  restores the production state, and lets the existing JFNK failure
  path fall back to Lambda/BA. Outside `JFNK_IN_EVAL`, the original
  STOP behavior remains.

Coverage includes `SE%STEQ_ADV` and `SE%T_EHB` in
addition to `SE%STEQ` / `STEQ_ED` / `STEQ_T` / `STEQ_T_EHB` - the
advective and EHB arrays previously leaked perturbed values into
`GENERATE_FULL_MATRIX_V3` (the SOLVEBA right-hand side) and the
resacc acceptance metric.  In frozen mode the module STEQ now keeps
the full outer-iteration value for the whole activation, so the
LAMBDA-mode `SOLVEBA_V13` fallback solves the same residual the
outer iteration measured.  A coupled eval admitted without a
registered `JFNK_SWEEP_CB` returns the budget-exhaustion sentinel
(forces F(0) abort -> LAMBDA fallback) instead of silently
evaluating a collisional-only residual.

## The F(0) identity test

Before Newton runs, one eval at x=0 validates the residual plumbing:

- Frozen mode must reproduce the outer `STEQ` (stored-term split must
  satisfy `F(0) = STEQ_full`); mismatch -> abort.
- Coupled mode **adopts the eval residual as the Newton baseline** -
  the lagged-store offset is real physics (the store moved during the
  sweep), not a bug. Guarded by a >N/2 spread check.
- Aborts on NaN or `jfnk_resid_nbad > 0` (components capped at 1e20,
  not zeroed - the earlier guard overwrote the whole vector with 1e20
  and masked divergence).

## Measured failure signature (the reason for the gates)

First live coupled fire, m15 it11 (job 12066716):

```text
|F|_w = 2.7e13;  linearity probe F(3v) -> 5e29
FD-JVP vs eps 1e-4..1e-8: ||Jv|| flat ~1e32 (never converges)
|FD| vs |BA v| rel diff 7e15; PC ||Mz-F||/||F|| ~ 3e2
=> admitted far outside the basin; 8 PC rebuilds x 8 clip retries,
   GMRES never past its first matvec, ~3.2 h wasted
```

Newton is only worth firing near the root. Hence the two-tier demand
gates - [[JFNK admission]].

Second live coupled attempt, m15 it73 (job 12073580): admitted on a
95% LAMBDA demand while the last full demand was 4.8e4% and
|F|_w=2.2e12 vs target 2.2e8 (the gate measured the wrong quantity -
fixed 04-Oct-2026, now `BA_LAST_RAW`). 13,826/128,730 near-null BA
rows regularized; surviving degenerate rows clipped on every retry.
~1.5 h burned before the cancel.

## Deliberately not implemented

- Full Chen & Shen (2006) preconditioner update from Krylov products.
- Larger `JFNK_M` by default (m=10; hybrid config uses 30 for coupled).
- Differentiating through `I_NEXT`; Anderson acceleration of the outer
  store-update map (belongs to the outer fixed point; needs a
  contractive inner map first).
- Line search on the committed Newton step - the admission gates +
  PTC + trust region + best-iterate fallback substitute for it so far;
  a backtracking line search on the commit is the candidate upgrade if
  activations fail rather than defer.

## Related

- [[JFNK admission]] - when this is allowed to fire
- [[Knop formal solver]] - the sweep `JFNK_RAD_SWEEP_V1` replays
- [[BA matrix and SOLVEBA]] - the fallback and the preconditioner
- [[JFNK_SOLVER_AUDIT]] - the full audit

## Reference audit vs Arramy, de la Cruz Rodriguez & Leenaarts 2024 (A&A 690, A12; arXiv:2406.17234)

Audit 08-Oct-2026 of the CMFGEN JFNK against the reference JFNK-NLTE
implementation.

### Points of agreement

- GMRES as the Krylov solver (paper finds GMRES strictly superior to
  BiCGSTAB for these problems).
- Jacobian-free matvec by finite-difference residual calls, one extra
  F evaluation per JVP - same cost model.
- Left preconditioning: (P^-1 J) ds = -P^-1 F.
- Preconditioner is the per-depth block of the (approximate)
  analytical Jacobian, refreshed per outer iteration, cost ~ one
  residual call. Our BA block is the same object; we factor it with
  equilibrated LU rather than inverting a point/block Jacobi, and the
  optional banded path adds the nonlocal cross-depth terms the paper
  leaves as future work.
- Krylov tolerance: paper recommends rtol < 1e-2 for stability
  (optimal 3e-4..3e-2). We use adaptive Eisenstat-Walker eta with
  floor 3e-4 - inside their recommended regime.
- Convergence tracked on the nonlinear residual norm, not the inner
  GMRES residual alone.

### Deliberate deviations

- Complex-step JVP (Im F(n+i eps v)/eps, eps~1e-200, O(eps^2)) vs
  our real forward FD with h = eps*max(1,|x|inf)/|v|inf. The paper
  needs complex-step because raw n spans 1e4..1e17; we answer the
  same problem with the ln(n/n0) unknown map plus scale-aware eps.
  Complexifying CMFGEN's RT code is not feasible; real FD is the
  standard Knoll & Keyes formulation the paper also cites. Harness
  T1 confirms min error ~4e-9 at eps~5.6e-9 (~sqrt(eps_mach)).
- Globalization: paper uses a plain update n += dn. We add hookstep
  trust region + line search + pseudo-transient continuation - the
  paper itself warns coarse steps may not converge; SN ejecta states
  start far outside the Newton basin.
- Singular/preconditioner rows: not discussed in the paper. BA has
  structurally null rows; we identity-replace at factor and project
  to zero at apply (JFNK_ROW_REG), which is a linear map and keeps
  M^-1 consistent.
