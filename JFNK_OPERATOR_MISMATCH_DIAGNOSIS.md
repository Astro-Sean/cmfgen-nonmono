# JFNK Operator-Mismatch Diagnosis Report

**Purpose**: self-contained technical report on why the Jacobian-Free
Newton-Krylov (JFNK) solver in a modified CMFGEN cannot make progress on
the frozen-J NLTE statistical-equilibrium residual, despite an LU
preconditioner that inverts the stored Jacobian nearly exactly.

**Status**: diagnosis complete; awaiting a decision on the fix path.
The model continues to converge via ordinary Lambda iterations
(demand 2.3% and falling at great iteration 368).

---

## 1. System under study

- Code: CMFGEN (NLTE statistical equilibrium + radiative transfer),
  locally extended with a JFNK solver tier in `cmfgen_jfnk/`.
- Model: `m15_2p82d_SOB_jfnk_zlaw14_fresh`, STELLA snapshot at 2.82
  days, 70 depth shells.
- Unknown vector `x`: ln-space population corrections, `n_i = n0_i * e^{x_i}`,
  per depth block of NT = 1814 entries (levels of all ions + ED + T rows).
  Total N = 1814 * 70 = 126 980. Temperature is held fixed this run.
- Residual (frozen-J tier):
  `F(x) = STEQ_MAT(n0*e^x) + STEQ_rad_stored + CMAT_RAD*(e^x - 1) + PTC(x)`
  - `STEQ_MAT`: collisional/rate terms, evaluated live at perturbed pops.
  - `STEQ_rad_stored`: radiative terms frozen at the outer-iteration state.
  - `CMAT_RAD*(e^x-1)`: stored block-diagonal linearization of the
    radiative terms vs population (exact derivative at x=0).
  - `PTC(x)`: pseudo-transient continuation term.
- Inner solver: right-preconditioned restarted GMRES (m=30) in a
  weighted inner product; globalization by hookstep trust region +
  line search; forcing tolerance eta ~ 1e-3.
- Preconditioner `M`: per-depth dense block of the stored `SE%BA`
  matrix assembled by `GENERATE_FULL_MATRIX_V3`, equilibrated
  (DGEEQU) and LU-factored (DGBTRF/DGETRF) with ridge fallback.

## 2. Original symptom

At JFNK admission, GMRES reduced the residual only ~1-2% before
stalling, then "lucky breakdown" (Arnoldi `h(j+1,j) ~ 0` at
iterations 2-4), Newton abort, Lambda fallback. Repeatable at every
admission.

|F| at admission ~ 3.03e11, dominated by HI rows at depths 64-70
where coll and rad terms each ~1e12 and nearly cancel (F_i ~ 1e11).

## 3. Campaign summary — preconditioner quality

| PC | Result |
| --- | --- |
| LSQ (normal eqns) | ~0% reachability — conditioning kills it |
| JACOBI | ~7% residual reduction, then stalls |
| Equilibrated LU of stored BA | solves `BA z = F` to 99.98% (`REACH frac=0.9999`, `cos(F,BAz)=0.9994`) — the stored system is fully reachable |

Conclusion: the stored BA system is solvable. The failure is that the
**actual Jacobian** (as seen by the finite-difference JVP) is not the
stored BA.

## 4. Bugs found and fixed during diagnosis

### 4.1 Residual save/restore leaks (fixed earlier)

`JFNK_RESIDUAL_V1` brackets module state around each eval. Audit found
missing coverage: `SE(ID)%STEQ_ADV`, `SE(ID)%T_EHB`, and the two-photon
accumulators `UP_RATE_TWO`/`DOWN_RATE_TWO` (coupled path). All added.
Post-fix evals are deterministic.

### 4.2 Baseline-offset bug (fixed; the main artifact source)

`jfnk_fx` retained the *outer-iteration* residual while `JFNK_JVP_V1`
computed `JV = (F_eval(x+h v) - jfnk_fx)/h`. Since
`F_eval(0) != F_outer` at the 1.1e-4 absolute level, every matvec
carried a constant offset `(F_eval(0) - F_outer)/h ~ 1.1e-4/1.5e-8 ~
8e3` on a handful of rows.

This single bug explained all earlier confusing observations:

- apparent **cross-depth** FD column entries (~8e3 scale),
- inflated soft-column norms (2e4 vs true ~5e2),
- `|A_FD z|` up to 59x|F|,
- repeat==pristine determinism (the offset is eval-invariant),
- the same response rows across admissions.

Fix: adopt `F_eval(0)` as `jfnk_fx` whenever the eval is finite, in
both frozen and coupled modes. Verified: pristine vs repeated column
probe now differs by exactly 0.0.

### 4.3 Two-photon call signature (fixed earlier)

`STEQ_BA_TWO_PHOT_RATE_V3` was called with 6 args instead of 7 (a
logical landed in an integer output-unit slot). Corrected.

## 5. Clean measurements (job 12153906, admission at it 367)

### 5.1 Stored-system reachability

```text
REACH: frac=9.9878E-01 |Jz|/|F|=1.0007 cos(F,Jz)=9.9939E-01
```

The LU-factored stored BA solves `BA z = F` nearly exactly.

### 5.2 FD-JVP action on that same z

```text
REACH-FD: |A_FD z|/|F| = 2.0824  cos(F,A_FDz) = 0.176
```
The true Jacobian acting on the stored-Newton direction produces a
response ~2x|F| but 82% orthogonal to F. GMRES has nothing useful to
build on: Arnoldi terminates at dimension ~2 (the preconditioned
operator maps the Krylov sequence into a near-invariant 2D subspace)
and stalls at 1.6% residual reduction.

### 5.3 Column-by-column probe (`COLPROBE`)

`JVP(e_i)` returns true Jacobian column i; compared against stored
`BA(:,i)`:

- **Stiff column** (top residual row, HI lev22 @ K70):
  `|FDcol|=1.1884e19 |BAcol|=1.2650e19 |diff|/|BA|=6.4e-2`.
  Stiff directions are captured well.
- **Soft column** (top |z| component, J=316 of NIV @ K60, z_i=-46.6):
  `|FDcol|=5.49e2 |BAcol|=1.23e3 |diff|/|BA|=0.90`.
  - All significant diffs at K=60 (same depth as the column). **Zero
    cross-depth entries post-fix.**
  - Off-diagonal entries agree to ~2% (e.g. -6.50 vs -6.48, -5.56 vs -5.55).
  - The dominant diff is the **diagonal itself**: `FD=3.98e-2` vs
    `BA=-1.10e3`. The true residual row barely responds to that
    level's own population; the stored BA carries an explicit rate
    diagonal ~3e4x larger.

### 5.4 Structural checks on that row (J=316, NIV, K60)

From the run's own WARNINGS equation-selection report:

- NIV g.s. equation (`EQXzV`, J=272) **is** replaced by the ionization
  conservation equation at all 70 depths — identical decision in the
  residual and the stored BA (same code, same baseline state).
- FIXPOP fixes zero levels at that ion/depth.
- So row 316 is an ordinary rate equation in **both** operators.
  The 0.04-vs-(-1104) gap is not a structural/flag mismatch.

### 5.5 Ruled out with direct evidence

| Hypothesis | Evidence against |
| --- | --- |
| Eval contamination / stale state | repeat-vs-pristine `diff=0.0` exactly |
| Finite-difference nonlinearity | `JVP(0.01 z)/0.01` bit-identical to `JVP(z)`; also `h` auto-normalizes, so identical output is guaranteed — the residual was separately shown linear at these scales by the exact-column probes |
| Cross-depth physics (advection, nonlocal rad) | `INC_ADV=F`; K69-only probe response at K70 top row = `0.0000`; all colprobe diffs same-depth |
| REPLACE flag flip on this row | row 316 is not an `EQXzV`/`EQ_SPECIES` slot; NIV g.s. replaced in both operators |
| FIXPOP identity rows | no levels fixed (WARNINGS) |
| OpenMP race / nondeterminism | no `!$OMP` in STEQ path; bit-identical repeats |
| Stored-system singularity | LU solves `BAz=F` at 99.98%; canonical solves `v=e1` also fail only in the null-adjacent directions |

### 5.6 The named missing-derivative channel

`JFNK_SYNC_ATM_FROM_POPS` (`jfnk_cmfgen_interface_v1.f:571-598`)
re-derives, from the perturbed `POPS`, on **every** residual eval:

1. `POPTOION` — repacks populations into `ATM(ID)%XzV`, `DXzV`
   (also feeds the population-dependent `REPLACE` test).
2. `POPION` — per-depth total ion population.
3. `COMP_LEV_DIS_BLK` — level-dissolution (occupation-probability)
   factors.
4. `SUP_TO_FULL_V4` — LTE reference and full-atom populations from
   super-levels.

`GENERATE_FULL_MATRIX_V3` builds `BA` treating all of those as
constants: it writes only the explicit rate-derivative terms
(d(rate)/dn through the stored rate coefficients). The true Jacobian
of `F` therefore contains additional `dF/dn` couplings through
`POPION`, dissolution factors, LTE reference populations, and the
super-level expansion — concentrated within each depth block.

For a level whose inflow and outflow both scale with `n_i` through
these factors (near-proportional balance), the **true self-derivative
is ~0** while BA's explicit diagonal is the full outflow rate sum
(-1104 here). These near-null true directions are exactly where the
LU solve parks large `z` components (|z_i| up to ~147 measured),
because BA prices them as real stiffness.

### 5.7 Net effect on the Newton step

`z = M^{-1} F` solves the stored system almost exactly, but its large
components live on columns whose true Jacobian is near-null or
mis-modeled. `A_FD z` lands 82% orthogonal to `F`; the Krylov space
`{r, AM^{-1}r, ...}` exhausts at dim ~2. An optimal damped step along
`z` alone yields only ~3% residual reduction
(`lambda* = (F.Az)/|Az|^2 ~ 0.085`), matching GMRES's ~1.6%.

## 6. Why this is intrinsic to the formulation

CMFGEN's `BA` was designed as a **damped-iteration matrix** for the
production solver (its `SOLVEBA` path does fixed-point/lambda steps
with clamping), not as a Newton Jacobian. Dropping the implicit-factor
derivatives is harmless for a damped iteration but fatal when `M^{-1}`
is used inside Krylov-Newton: the preconditioner manufactures "free"
directions the true residual does not honor, and conversely the true
operator's null-adjacent structure is invisible to `M`.

Secondary structural note: `REPLACE` and `FIXPOP/ZERO_STEQ` decisions
are recomputed inside every eval (at the *perturbed* state), while the
stored BA froze them at the baseline. Perturbations near the `FAC`
threshold can discontinuously switch a row between rate and
conservation equations — a residual non-smoothness the FD JVP sees as
an O(1/h) artifact. This is a real but secondary gap; it did not
produce the measured dominant diff (row 316 isn't replaceable), but it
is worth removing regardless.

## 7. Current runtime state

- Job 12153906 (interactive partition): great iteration 369.
- MAXCH demand 5.2% -> 2.3% across iterations 367-368; Lambda
  iterations converging the model unaided.
- JFNK machinery itself works correctly: admission gating, strike
  counter, Lambda fallback, warm restart, diagnostics all healthy.
- JFNK simply cannot accelerate: 3 strikes -> cooldown, Lambda keeps
  making progress.

## 8. Options on the table

**A. Accept frozen-J limits.** Let Lambda converge (in progress, ~2.3%
demand). JFNK stays as bounded, harmless attempts. Cost: nothing.
Benefit of JFNK on this residual: none proven.

**B. Structural-consistency fix (cheap, ~1-2 h).** Save the baseline
`REPLACE`/`ZERO_STEQ` arrays at preconditioner build and pass them to
the eval via `USE_PASSED_REP=T`. Removes the discontinuous
equation-structure flips under perturbation — a genuine JVP-quality
improvement. Will **not** close the dominant implicit-factor gap.

**C. Augment BA with implicit-factor derivatives (research-level).**
Add `d(POPION)/dn`, `d(dissolution)/dn`, `d(LTE-ref)/dn`,
`d(superlevel-map)/dn` terms — either analytically into
`GENERATE_FULL_MATRIX_V3`'s BA writes, or numerically as a low-rank
update to the stored blocks. Would make `M` a true Jacobian. Effort:
significant; physics-code risk: high; payoff: uncertain.

**D. Alternative solver strategies.**

- Damped application of `z = M^{-1}F` directly (max ~3% per step —
  poor).
- Different Krylov method (BiCGStab/FGMRES) — same operator pair,
  same fundamental limit.
- Restrict the Newton solve to the stiff-row subspace where BA is
  trustworthy (~6% column error) — but the `|z|<1` restricted probe
  already shows `cos(F, A_FD z_restricted) = -0.29`: even the trusted
  components' aggregate response is anti-correlated, because per-column
  6% errors accumulate across a direction built from cancellation-
  dominated rows.
- Reformulate the residual/unknowns so BA *is* a faithful Jacobian
  (e.g., solve in the exact equation structure BA models, with
  factor-dependent terms moved into the frozen part) — a
  formulation-level redesign, related to option C.

## 9. Questions for the advising AI

1. Is there a standard treatment in NLTE solvers for the
   implicit-factor derivative gap (POPION/dissolution/LTE-reference)
   short of full analytic Jacobian assembly?
2. Is the deeper problem that the **frozen-J residual itself** is
   near-unsolvable at the dominant rows? The top rows are
   `coll ~ -rad ~ 1e12` cancellations with `F ~ 1e11`; a pure
   population Newton cannot change the frozen radiative terms, so the
   feasible correction manifold may be nearly empty there. If so, is
   any population-space Newton step doomed regardless of
   preconditioner quality?
3. Does the diagnostic evidence support a *partial* Newton (solve
   only the depth blocks/rows where BA column error < threshold), or
   does the anti-correlated restricted-probe result (`cos=-0.29`)
   refute that too?
4. Given Lambda iterations are converging the model at ~0.5-2% demand
   reduction per iteration, is the rational choice simply option A+B
   and declaring JFNK-inapplicable for frozen-J residuals of this
   structure?
5. Any standard low-rank/deflation technique to repair `A M^{-1}`'s
   effective rank-2 behavior when the missing Jacobian terms are
   within-block dense couplings?

## 10. Key file references

- Residual/JVP/preconditioner driver:
  `cmfgen_jfnk/new_main/mod_subs/jfnk_cmfgen_interface_v1.f`
  - eval state re-derivation: `JFNK_SYNC_ATM_FROM_POPS` (lines ~571-598)
  - residual assembly: `JFNK_ASSEMBLE_STEQ_V1` (~525-559),
    `JFNK_RESIDUAL_V1` (~1600-1660)
  - JVP: `JFNK_JVP_V1` (forward difference, h auto-normalized)
  - weighted inner product: `JFNK_DOT_V1` (~4265)
- Newton/GMRES/hookstep:
  `cmfgen_jfnk/new_main/mod_subs/jfnk_solver_v1.f90` (`gmres_m_v1` ~1651)
- Solver control / diagnostics / admission:
  `cmfgen_jfnk/new_main/mod_subs/jfnk_solve_corr_v1.f`
- BA construction:
  `cmfgen_jfnk/new_main/subs/generate_full_matrix_v3.f`
  (`REPLACE` decision ~325-339, `FIXPOP` call ~366-374,
  column scaling by `POPS(J,L)` ~384-389)
- Fixed-level handling:
  `cmfgen_jfnk/new_main/subs/fixpop_in_ba_v3.f`
- Doc: `cmfgen_jfnk/cmfgen_obsidian/solvers/JFNK solver.md`

## 11. Reproduction

Validation run:
`cmfgen_validation/m15_2p82d_SOB_jfnk_zlaw14_fresh/cmfgen/2.82days/OUTGEN`,
SLURM job 12153906, admission block near line 1274-1600. All numbers
above are quoted from that block. The executable is
`cmfgen_jfnk/exe/cmfgen_dev.exe` (built ~16:22, unit tests T1-T3 pass).

## 12. Follow-up on the fresh OUTGEN (job 12153985, 11-Oct-2026)

This follow-up corrects the earlier runtime summary: Lambda fallback did
not permanently prevent FULL-BA calls. After the cooldown window, FULL-BA
resumed at it399 and reduced demanded MAXCH from 32.5% to 1.35% by it402.
It then hovered near 1.4-1.5% through it404; `EPS_TERM=0.05%` was still
not reached. The outer BA path works but has a repeatable endgame plateau.

The frozen-J operator mismatch persisted on the attempts at it405-407.
At it384, the soft column had `|diff|/|BA|=0.831` and the top-F column
had `0.064`; `cos(F,A_FD z)=0.0405`, followed by GMRES stagnation. The
frozen tier is not a dependable finisher in the current formulation.
The configured `JFNK_FRZ_GATE=100%` admitted it at demands of 4.03% as
well as near 1.5%; the 5% coupled gate does not help this SOB run because
`JFNK_COUPLED` is unavailable under `DO_NONMONO_SOB`.

### Fallback correction and live regression

The old non-catastrophic failure and cooldown branches forced Lambda mode,
set `COMPUTE_BA=T`, and caused repeated full BA rebuilds. Cooldown calls
cost about 366 seconds, compared with about 70 seconds for ordinary FULL
calls in the same OUTGEN. The run later returned to FULL mode, so the
correct conclusion is avoidable delay and mode churn, not a permanent
lockout.

A pass-through fallback was built and picked up by the requeued process.
It preserved the caller's FULL mode, but its first cooldown fallback
reached `SOLVEBA_V13` with an uninitialized local change cap and stopped
at it408 (`Maximum change for normal iteration must be > 1`, job exit 6).
The cap is now initialized from the saved caller input before dispatch;
the corrected build has not yet been validated in a model run. Job
12153986 was canceled by its failed dependency.

### Implicit-factor prototype and limits

Source review narrows the mismatch claim. `STEQ_MULTI_V10` has explicit
BA columns for direct electron-density and ion-density effects. The
candidate missing derivative is through the aggregate `POPION` and the
level-dissolution coefficients in `MOD_LEV_DIS_BLK`; full-level
populations are still rebuilt from the trial superlevels. A new opt-in
`JFNK_LAG_FACTORS` mode snapshots `POPION` and the dissolution arrays for
one frozen-J activation while leaving direct ED, T, and ion-density
columns live. This is a modified residual, not an exact Newton Jacobian
for the physical outer problem. The outer iteration refreshes factors
between activations, so only a controlled CMFGEN run can establish
whether it reduces the measured production column gap and improves
true outer convergence.

The surrogate T18 analogue passes: live FD matches its exact Jacobian,
but the live Jacobian differs from the baseline BA by `2.11e4` in the
chosen direction; with the aggregate factor held at baseline, FD matches
BA to `5.0e-7`. This validates the test model only, not CMFGEN.

### VADAT and temperature policy

The run used `FIX_T=T`, `FIX_T_AUTO=F`, no `DO_T_AUTO` override,
`TAU_SCL_T=0`, `JFNK_TAUTO=T`, and `N_FULL_TREL=0`. `FIX_T_AUTO` is the
`VARFIXT` depth-mask behavior; by itself it does not release global
`FIX_T`. `DO_T_AUTO=T` releases global fixed T when `MAXCH_CTRL<50%`;
with the current ~1.5% population demand, that is an early release, not
a gradual thermal safeguard. `JFNK_TAUTO` only damps JFNK's T components
when T is already free; normal FULL-BA temperature steps do not use that
ladder. The prepared YAML now sets `FIX_T_AUTO=T`, `DO_T_AUTO=T`,
`JFNK_FRZ_GATE=3%`, two eligible contracting calls, and
`JFNK_LAG_FACTORS=T`. These are proposed next-run settings, not values
used by the failed job.

Status: production build succeeded after the cap fix; solver tests passed
21/21, interface tests 9/9, and surrogate T15-T19 passed. The full
surrogate suite retains four pre-existing T3/T4 failures. Config
resolution passed through `run_model_grid.py --print-resolved-config`.
No new production run has been submitted; fallback behavior and
`JFNK_LAG_FACTORS` remain UNVERIFIED in CMFGEN, and no spectrum was
produced by the failed job.
