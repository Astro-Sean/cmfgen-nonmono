# DECISION_RULE.md -- pre-registered comparison rule

Committed before any HELDOUT-set run. The tuning set (11 instances,
`instances/split.txt` tag `TUNE`) may be used to pick hyperparameters;
held-out numbers are the only evidence reported in the final comparison.
Amendments, if any, land as separate commits with a stated reason, and
both versions' results are reported.

## Success definition

A run succeeds iff the gross-flux-normalised residual
`r_i = |F_i| / G_i`, where `G_i` is the exact row gross flux
(inflow + outflow terms accumulated in `eval_F` alongside `F_i`, the
analog of accumulating `sum_j |STEQ_ij| n_j` during assembly),
satisfies `max_i r_i <= max(nml_tol, 10 * floor_i)` within 200
sweep-equivalents. `floor_i` is the deterministic row noise floor
`k * eps_mach` with `k = 100`. On the toy the deterministic floor is
far below `nml_tol` everywhere it was measured (`run_noise`, h-scan
to 1e-12), so the criterion reduces to `rmax <= nml_tol` plus the
budget -- but the floor clause is retained because it is the
production-relevant part of the rule.

Status `converged` in the RESULT record = met the residual criterion
inside the outer-iteration cap; `budget`, `stagnation`, `diverged`,
`timeout`, `crash` are failures.

## Sweep-equivalent cost

`sweep_eq = n_sweep + c_ba * n_ba`, reported at
`c_ba in {0.1, 0.3, 1.0}`. The harness enforces the 200-SE budget at
the strictest weight (`n_sweep + n_ba > 200` => `budget`), so a run
reported as converged is within budget at every `c_ba`. Cheap evals
(frozen-mode residuals, exact trace solves, Anderson least-squares
evals) and JVPs are reported raw but do not add sweeps: a frozen-JVP
shares the eval's lagged transfer operator, matching the production
frozen tier where the JVP is a cheap stored-BA response. A
coupled-JVP costs one full transfer eval in the toy and is counted
through `n_sweep` inside `eval_F` mode 0.

## Disqualification

Success rate < 95% on the held-out set (i.e. fewer than 20 of 21
instances) disqualifies the method outright.

## Ranking

Among qualified methods, rank by median sweep-equivalents to success.
A method must beat the best of `s0` by >= 1.5x (median at `c_ba=0.3`)
to justify its added complexity. Report IQR and worst case alongside
the median.

Tie-break, in order: fewer new tuning knobs; less code; no dependence
on per-instance tuning.

## Method matrix (ablations)

`methods.txt`: `s0`, `and2`, `and4`, `and8`, `s2f`, `s2f+b`, `s2f+c1`,
`s2f+a`, `s2f+b+c1+a`, `and2+b`, `and2+a`, `s2c`.

Mapping to the spec: `s0` = S0 stock; `s2f`/`s2c` = S2 frozen/coupled;
`+b` = B weighted residual; `+c1` = C1 secant-h latch
(`theta=0.1` namelist); `and*` = C2 Anderson depth 2/4/8; `+a` = A
reduced-space trace elimination. S1 (stock Ng) is represented by the
`and*` family: the production Ng band accelerator is an Anderson-m=2
equivalent and cannot be called from the toy harness -- stated as a
limitation in the audit doc.

## Sensitivity

- Report all statistics at `c_ba in {0.1, 0.3, 1.0}`.
- Hyperparameter fragility: a method that only qualifies at one tuned
  knob value is called fragile. On the tuning set the relevant knobs
  are `theta` (c1 latch, {0.01, 0.1, 0.5}) and Anderson depth
  ({2, 4, 8}); held-out runs use the single committed method matrix
  above.

## Failure-mode taxonomy

`diverged` = rmax > 1e3 or NaN; `budget` = outer-iteration cap
(sweep-equivalent proxy for the 200-SE budget) exhausted with
rmax > tol; `stagnation` = no 5% residual improvement over the
stagnation window; `timeout`/`crash` = harness failures (counted as
method failures in the success rate).
