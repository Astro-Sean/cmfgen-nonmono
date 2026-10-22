# JFNK_CONVERGENCE_BAKEOFF.md

Branch: `convergence-bakeoff` (off `development`). Date: 22-Oct-2026.

Evidence base: the synthetic test bed in `tests/convergence_bakeoff/`
(toy coupled NLTE populations + mean intensity + lagged store), which
links the production `jfnk_solver_v1.f90` unchanged. Nothing here was
run on a real CMFGEN model; the production-side statements are
flag-gated code + code reading, and the toy statements are measured.
Section 4 lists exactly what needs a real run.

## 1. Review findings, verified

| # | Claim | Verdict | Evidence |
|---|-------|---------|----------|
| 1 | Merit is the raw Euclidean norm; SQRT_STE weights freeze at the baseline | CONFIRMED | `JFNK_COMPUTE_RSCALE_V1`, `jfnk_cmfgen_interface_v1.f:478` -- `SQRT_STE` weights are `1/sqrt(max(|F0_i|, floor))` from the admission-time `jfnk_fx`; later calls replay them. |
| 2 | `DIAG_RDG` branch unreachable (`(1:4)=='DIAG'` shadows it) | CONFIRMED + FIXED | `jfnk_cmfgen_interface_v1.f:490-512`; reorder committed as `bc7fb139`. |
| 3 | Lagged intensity store held fixed inside each eval; outer Picard refresh caps the rate | CONFIRMED | `FREEZE_I_NEXT=.TRUE.` set at `jfnk_solve_corr_v1.f:473`, cleared at `:2104`; `JFNK_RAD_SWEEP_V1` (`cmfgen_sub.f:5537`) brackets quarantine/adaptive-refresh state per eval. Toy reproduces the cap: `rho_lag` instances converge at the outer rate regardless of inner method. |
| 4 | Frozen-J re-solves SOLVEBA's problem; coupled mode is the potential win; no sweeps-per-residual measurement existed | CONFIRMED | Toy: `s2f` costs the same sweep-equivalents as `s0` when it converges (med 34.4 vs 28.7) and converges less often (16/21). `s2c` costs ~1 transfer solve per JVP: 163 JVPs median, budget out on 18/21. On the toy the coupled tier is never competitive; it pays for information the frozen tier gets free. |
| 5 | FD JVP uses sqrt(eps_mach)-scale h on a noisy residual; noise source unidentified | PARTIALLY CONFIRMED | `EPS0 = max(sqrt(eps_mach), JFNK_EPS_JVP_LOC)` at `jfnk_cmfgen_interface_v1.f:2740-2758`. On the toy, deterministic evaluation noise (iterative J-solve count changes, population clips, cancellation at gross-flux stiffness to 1e18) produces a clean FD window to h=1e-12 -- the 36x |FD|/|BA v| discrepancy seen in production is NOT reproducible deterministically. Remaining suspects: OMP reduction order and unbracketed sweep state. `JFNK_NOISE_DIAG` now measures this on real runs (see runbook). |
| 6 | Trace species are one-way coupled, sit in the Krylov vector, and drive the "unsatisfiable demand" machinery | CONFIRMED | `JFNK_N=NT*ND` (`jfnk_solve_corr_v1.f:499`): every level row is a Krylov unknown. Toy: `s2f` stagnations all land on trace rows (worst row = first trace row at every stalled activation; majors already at their fixed point -- the `|s0-s2f|=1.4e-11` sanity diff on i10 proves majors converged while the trace residual stayed ~1). Eliminating them analytically takes every method from mixed results to 21/21. |

Additional Phase-0 facts:

- Ng interaction: `do_ng_accel_band_v4` needs 4 consecutive clean
  SCRTEMP records; `JFNK_COMMITTED` (`solve_for_pops.f:2720`) resets
  that window on every committed Newton step, so Ng and JFNK are
  effectively exclusive in a run.
- Coupled-eval cost: one `JFNK_RAD_SWEEP_V1` per eval = one full
  continuum formal solve (not an inner loop); confirmed at
  `cmfgen_sub.f:5537`.
- Trace rows are linear in trace populations at fixed majors and J on
  the toy by construction (`atr` block + linear source `bsrc`); the
  production one-way coupling is approximate but standard.

## 2. Held-out results (21 instances, 8 regimes)

Success = gross-flux-normalised residual <= 1e-6 inside 200
sweep-equivalents (enforced at c_ba=1). Rule committed before the
held-out run: `tests/convergence_bakeoff/DECISION_RULE.md`
(`7c07dbe8`).

| method | ok | med SE @0.3 | @0.1 | @1.0 | IQR @0.3 | worst @0.3 | med cheap | med JVP | failures |
|--------|-----|------------|------|------|----------|------------|-----------|---------|----------|
| s0 (stock) | 21/21 | 28.7 | 24.9 | 42.0 | 26.1-33.9 | 45.6 | 39 | 0 | -- |
| and2+a | 21/21 | 12.8 | 11.6 | 17.0 | 12.8-12.8 | 14.1 | 19 | 0 | -- |
| s2f+b+c1+a | 21/21 | 12.0 | 10.0 | 19.0 | 12.0-12.8 | 13.6 | 169 | 43 | -- |
| s2f+a | 21/21 | 14.3 | 12.1 | 22.0 | 12.0-14.3 | 15.6 | 208 | 103 | -- |
| s2f+b | 20/21 | 34.4 | 26.8 | 61.0 | 28.0-36.6 | 52.6 | 626 | 321 | stag:1 |
| and2 | 19/21 | 36.5 | 31.5 | 54.0 | 30.0-41.7 | 51.8 | 76 | 0 | stag:2 |
| s2f+c1 | 18/21 | 34.4 | 26.8 | 61.0 | 25.9-36.0 | 56.8 | 1018 | 309 | stag:3 |
| s2f | 16/21 | 34.4 | 26.8 | 61.0 | 19.1-47.6 | 55.2 | 5054 | 2719 | stag:5 |
| and4 | 18/21 | 40.9 | 35.3 | 60.5 | 24.6-48.9 | 66.1 | 79 | 0 | stag:3 |
| and8 | 18/21 | 34.5 | 29.9 | 51.0 | 23.9-48.5 | 56.4 | 73 | 0 | stag:3 |
| and2+b | 18/21 | 53.1 | 45.7 | 79.0 | 41.2-55.7 | 63.5 | 103 | 0 | stag:3 |
| s2c | 3/21 | 108.0 | 106.0 | 115.0 | 90-199 | 199 | 7 | 163 | budget:18 |

Qualified (>=95%): `s0`, `and2+a`, `s2f+a`, `s2f+b+c1+a` (`s2f+b` sits
exactly at 20/21 but gains nothing over stock).

Per the rule, beating stock by >=1.5x median: `s2f+b+c1+a` 2.39x,
`and2+a` 2.24x, `s2f+a` 2.0x -- all clear the bar at every c_ba
weight, and all three carry `+a`. No qualified method without `+a`
beats stock at all: `s2f` converges less often for the same cost.

Failure-mode classification: all non-`+a` failures are `stagnation`
on trace rows (worst-row diagnostic); `s2c` failures are `budget`
(eval cap inside the JVP). No divergence, timeout or crash in the
held-out set.

Sensitivity: the ranking is unchanged at c_ba in {0.1, 0.3, 1.0}
(cheap-eval-heavy methods shift slightly but never cross `+a`).
Anderson depth is fragile in the wrong direction (4,8 worse than 2);
the c1 theta knob changed nothing alone. `+a` has no knobs.

## 3. Recommendation (toy evidence only)

**Recommended: `and2+a` -- Anderson(2) on majors with exact trace
elimination.** 12.8 median SE vs stock 28.7 (2.24x), tightest IQR of
any method (12.8-12.8), zero JVPs, one knob (depth=2), no Newton
machinery at all: it is `s0` with two extra ingredients (a 2-term
secant combination over the stock map + traces resolved exactly per
eval). By the tie-break it beats `s2f+b+c1+a` (12.0, but three knobs +
the whole JFNK stack) -- the 0.8-SE median gap does not justify the
complexity difference.

Ranked fallbacks:
1. `s2f+b+c1+a` (12.0 med) if a Newton solve is wanted anyway for
   other reasons (e.g. coupled tier later).
2. `s2f+a` (14.3 med) -- plain frozen JFNK + trace elimination.
3. `s0` stock -- always qualified; the honest default until real-model
   validation.

What is NOT supported: bare JFNK (`s2f` loses to stock on
robustness at equal cost), bare Anderson (19/21 or worse),
coupled-JVP (`s2c` pays a formal solve per JVP and budgets out).

## 4. Needs human validation on a real model

The production port committed here is `JFNK_TRACE_ELIM` (VADAT, default
F): the Gauss-Seidel form of `+a` -- exempt rows (same
`BA_DEM_PFLOOR` predicate) masked out of the Newton residual and
correction per activation via `JFNK_MASK_FIXED_UNKNOWNS`; traces are
left to the existing Stage-B `BA_TRACE_PHASE` polish. This is the
minimal, reviewable port. The Anderson variant (`and2+a`) is NOT yet
ported; it would wrap `SOLVEBA` calls in an Anderson-m=2 combination
with the residual-decrease safeguard -- a driver-level addition, not
attempted because trace-elimination semantics needed verification
first.

Runbook (a stalled non-monotonic KNOP checkpoint, `DO_JFNK=T`):

1. Baseline stock: `DO_JFNK=F`. Log sweeps and `MAXCH` per iteration.
2. Recommended candidate: `DO_JFNK=T`, `JFNK_TRACE_ELIM=T`. Same
   gates otherwise. Compare: sweeps per decade of `MAXCH`, final
   `MAXCH`, count of JFNK fallbacks to SOLVEBA, Stage-B polish entry
   (`BA_TRACE_PHASE` iterations) and whether trace rows settle.
3. Noise diagnostic: `DO_JFNK=T`, `JFNK_COUPLED=T` (KNOP only),
   `JFNK_NOISE_DIAG=T`. Grep `JFNK_NOISE eval=` lines: `|dY|` at
   `|dX|=0` is the eval non-determinism floor; repeat with
   `OMP_NUM_THREADS=1` vs production thread count to separate OMP
   reduction-order noise from stateful sweep decisions. This
   decides whether the smooth-F candidate (D) is worth building.
4. Acceptance criteria for `JFNK_TRACE_ELIM`: (a) runs with the flag
   reach the same `MAXCH` floor as stock on the same checkpoint;
   (b) sweep count to that floor does not increase vs stock;
   (c) `BA_TR_ABS` converges in the Stage-B polish after the masked
   activation commits (i.e. traces actually get resolved);
   (d) no new `jfnk_resid_nbad` spikes vs baseline.

Overturn criteria: if masked activations converge the majors but the
trace polish then fails to contract (BA_TR_ABS stalls), the GS
variant is insufficient and the in-eval trace solve (variable
projection) is needed instead -- a larger port deferred deliberately.

## 5. What the toy cannot represent

- The toy's `+a` is variable projection on a *linear* trace block;
  production trace coupling is only approximately one-way, and the
  GS port relies on Stage-B polish converging rather than an exact
  solve per eval. If real trace rows are nonlinear in themselves,
  the toy overstates the win.
- The toy BA is a genuinely imperfect model of its Jacobian (AL*),
  but production BA has far more structure (per-ion blocks, banded
  coupling); the ~2x margin may not transfer.
- The toy residual is deterministic and non-smooth only through
  Jacobi-count switching; production's FD noise floor (finding 5)
  was NOT reproduced -- the toy cannot say whether JFNK JVPs are
  trustable on the real model. `JFNK_NOISE_DIAG` exists to answer
  exactly this.
- Lagged-store semantics are a scalar blend in the toy
  (`rho_lag`); production I_NEXT quarantine/repair logic is
  stateful in ways the toy does not model.
- S1 (stock Ng) is represented by the Anderson family, not by the
  actual `do_ng_accel_band_v4` code path.
- Ensemble size: 32 instances, 21 held-out -- the failure-rate
  resolution is ~5%.

## 6. Reproduction

```bash
cd tests/convergence_bakeoff
make                       # links production jfnk_solver_v1.f90 unchanged
python3 gen_ensemble.py    # regenerates instances/ + split.txt
python3 run_bakeoff.py --set HELDOUT --methods-file methods.txt \
    --out results_heldout.csv
python3 analyze.py results_heldout.csv
./bakeoff sanity instances/i10_thin.in   # same-fixed-point check
DBG=1 ./bakeoff s2f instances/i10_thin.in  # worst-row diagnostics
```

Committed artifacts: `results_heldout.csv`, `results_tune.csv`,
`heldout_summary.txt`, `DECISION_RULE.md` (pre-registered at
`7c07dbe8`), `analyze.py`, `gen_ensemble.py`, `run_bakeoff.py`.
