# Convergence strategy

The techniques this tree uses to make a non-monotonic CMFGEN model
converge, organized by method rather than by code module. The
component-level documentation lives in `solvers/` and `control/`;
this folder is the "how do we get it to converge" view.

## The one-paragraph version

Convergence here is not one algorithm - it is a sequence of
increasingly aggressive methods, each gated on evidence that the
cheaper, safer stage has done its job. Lambda iterations (bounded
fixed-point steps) run while the radiation field and populations are
still transient. Once the field is a demonstrated fixed point, the
quasi-Newton BA step is released under a residual-acceptance rollback
net. Only when the demanded correction has contracted does JFNK fire,
and only inside its trust region. Acceleration (Ng, Jacobian reuse)
fills the gaps; limiters, caps, and rollback make every step
recoverable.

```text
cheap, safe, slow                    expensive, fast, narrow basin

Lambda iteration -> BA (SOLVEBA) -> JFNK (frozen -> coupled)
        |               |                  |
     cap-limited    resacc net          trust region
     always on      rollback           + PC + PTC

        ----> acceleration layered on top: Ng, AV averaging,
              FIX_BA Jacobian reuse, warmup/restart machinery <----
```

## Why staged

- The field solve and the population update are only consistent at a
  fixed point; an aggressive update on a transient field re-poisons
  the map (the ordering rule in [[Field maturity and GS freeze]]).
- Newton is quadratically fast only inside its basin. Measured waste
  outside it: ~3.2 h and ~1.5 h on the two live coupled fires
  ([[JFNK solver]]). Hence sustained-demand gates, not a schedule.
- Every layer keeps the previous layer as fallback: JFNK -> BA ->
  Lambda -> hard stop.

## Initial LAMBDA handoff

For non-monotonic SOB and Knop runs, `DO_LAM_AUTO` no longer releases
full iterations on one `MAXCH < 50` call alone. `SOLVE_FOR_POPS` counts
consecutive calls below `LAM_EXIT`, using applied `MAXCH_CTRL`; forced
LAMBDA calls reset the count. Luminosity ratio and drift do not enter
this decision because shocks permit physical deposition between the
boundaries. The handoff occurs after `LAM_N_STABLE` qualifying calls. The stock
50% handoff remains unchanged when both non-monotonic flags are off.
This prevents the initial handoff from bypassing the same stability
controls used by later LAMBDA exits. After handoff, a raw BA demand
cannot reopen the LAMBDA or JFNK gate while its residual check is
pending or after the trial was rejected, masked, or soft-rejected; a
periodic probe must pass residual acceptance first. These changes do
not prove the BA direction converges. Runtime/model validation is
pending.

The gate alone has no escape when the LAMBDA attractor floors above
`LAM_EXIT` (m15 SOB: demand parked ~45-49% with `stable=0/3`). The
`BA_PROBE_LAM` rd_lam probe closes that: every `BA_PROBE_LAM` Lambda
calls with demand in `[LAM_EXIT, LAM_ENTER)` fire one `BA_DAMP_CAP`-
bounded full call inside the `RD_LAMBDA` phase; residual acceptance
judges it, a clean verdict pins `CNT_STABLE_LAM=N_STABLE_LAM` so the
ordinary `DO_LAM_AUTO` switch hands off, and a rejection rolls back
and waits out the next cadence. See [[Lambda iteration]].

## Map of techniques

| Technique | Note | Where it lives |
|---|---|---|
| Bounded fixed-point update, mode control | [[Lambda iteration]] | `SOLVEBA_V13`, mode controller in `solve_for_pops.f` |
| Ng / AV acceleration, Jacobian reuse | [[Acceleration methods]] | `ngaccel.f`, `do_ng_accel_band_v4.f`, `ave_flips` |
| Caps, limiters, trust regions, rollback | [[Step limiting and rollback]] | `solve_for_pops.f`, resacc, `INEXT_RELAX` |
| Warmup, T release, restart | [[Warmup and staged release]] | `N_LAM_RWARM`, `DO_T_AUTO`, checkpoints |
| Exit tests, metrics, stall detection | [[Convergence criteria]] | `MAXCH<EPS`, `MOM_SW_*`, watchdogs |

## Component notes these build on

- [[The great iteration]] - the call sequence everything hangs off
- [[BA matrix and SOLVEBA]] - the linearized update and resacc
- [[Field maturity and GS freeze]] - the ordering enforcement
- [[JFNK admission]] - the Newton gating chain
- [[Failure modes and guards]] - the observed blow-ups each guard answers

## The standing open problem

As of Oct-2026 neither the T path nor the F path has demonstrated a
monotonically contracting endgame on the Type IIn model: BA kicks are
bounded and absorbed but not shrinking ([[convergence_issues]]
Sect. 1.1, 1.9). The techniques in this folder bound the damage and
order the stages; whether they converge faster than pure Lambda is
the unresolved question.

## Recommended hardening steps (06-Oct-2026 review)

Mirror of [[convergence_issues]] item 22 (architecture-level facts only).

- Relative (rate-normalised) SE residual for the acceptance test in
  [[BA matrix and SOLVEBA]] so the max channel is not saturated by the
  inner-boundary equations.
- Implemented: a degraded `BA_PROBE_LAM` trial restores its full
  pre-probe state on the LAMBDA follow-up and skips the stale solve;
  shared by SOB and Knop/JFNK (see [[Step limiting and rollback]]).
- Residual-judged NG, then safeguarded Anderson acceleration as the
  NG replacement ([[Acceleration methods]]).
- OBS J for SOB-converged populations from the Knop formal solver
  rather than the clamped-sigma stock FG solver ([[Sobolev non-mono
  solver]], [[Knop formal solver]]); COH_ES=F for the final model.
- Runtime validation of [[JFNK solver]]: no run has yet reached the
  Newton loop; the SOB-seeded frozen-J stage is the first target.
