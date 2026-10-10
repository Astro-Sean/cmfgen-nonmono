# Step limiting and rollback

The damping layer: every update path has a cap on what it may commit,
and every aggressive step has a net to catch it. The design invariant
([[JFNK admission]]): **every unchecked step needs a net** - resacc
for BA, best-iterate + commit checks for JFNK, caps for Lambda.

## Committed-step caps

| Cap | Applies to | Scale |
|---|---|---|
| `MAX_LAM` (`MAX_LAM_COR`) | Lambda update | production profile `1e10`, effectively unclamped |
| `MAX_LIN` (`MAX_LIN_COR`) | full BA update | production profile `1e1` |
| `MAX_JFNK` | committed JFNK step | ratio bound, `|x| <= ln(MAX_JFNK)` |
| `GS_LAM_CAP` | Lambda during near-mature freeze | factor 3 |
| `BA_DAMP_CAP` | damped-full leash while demand saturated | SOB path |
| `MAX_dT_COR` | temperature correction | 0.2 |
| `T_MIN` | temperature floor | 0.3 |

- The 1e7% sentinel: demanded corrections saturate there; a demand
  still at the sentinel after a committed step means marching, not
  contracting - `BA_SAT_RAW`/`BA_SAT_LIM` demote on it and
  `BA_ENTER_RAW` holds the Lambda->full exit.
- `SCALE_OPT=MAJOR` (the MAJOR limiter): the step scale is set only by
  levels above `1e-10*Ne`, so a trace level screaming for a huge
  correction cannot throttle the whole model.
- `LAM_SCALE_OPT`: `NONE` (unlimited Lambda corrections) is what let
  the LTE-at-T_gas cold start explode on m15 ([[convergence_issues]]
  Sect. 6). The current production recipe nevertheless keeps
  `NONE` with `MAX_LAM=1e10`: small Lambda caps made the LTE-to-NLTE
  cascade impractically slow, so safety is supplied by resacc, the
  demand gates and the damped-full leash instead of by throttling
  every Lambda step.
- Multiplicative update `POPS *= (1 - STEQ*SCALE)` with `BIG_LIM < 1`
  keeps populations positive by construction.

## Under-relaxation on the lagged field

- `INEXT_RELAX` (0.5): damped writes to the `I_NEXT` store - the
  lagged field is a fixed-point map, and undamped writes make it
  oscillate.
- `JNEXT_RELAX`: same damping applied to `H_NEXT_FRQ`/`J_NEXT_FRQ` in
  the MOM path.

## Residual acceptance ("resacc")

`DO_BA_RES_ACC` - mandatory on KNOP runs. Snapshot `POPS`, take the BA
step, re-measure the SE residual at the next call, roll back what got
worse. Per-depth masked rollback under `BA_REJ_FRAC` (default 0.25):
failing columns are restored and their corrections zeroed while the
rest of the step stands; bulk failure or too many columns rolls back
the whole state and demotes to Lambda. Since 13-Oct-2026 the
`BA_PROBE_LAM` probe path uses the same tiers - a few-column failure
is masked (kept progress at the good columns), only a bulk failure
restores the whole pre-probe state and skips the Lambda solve. Trust
radius `BA_TRUST_CUR` halves on rejection; `BA_NREJ_LIM` hard-stops a
run of rejected steps. Since 07-Oct-2026 the counter is progress- and
hold-aware: calls under `GS_POPS_HELD` defer judgement to the first
population-moving call, and a masked/soft step that improved either
residual channel resets `BA_NREJ_CONSEC` instead of incrementing it.
`DO_BA_QUAR` (default F) can additionally park individual failing
(equation,depth) pairs - up to `BA_QUAR_MAX` - so one stiff equation
no longer vetoes its whole depth column. Full mechanics in
[[BA matrix and SOLVEBA]].

## Field-side rollback

- `MOM_KNOP_REFRESH` exits on the stall test (two consecutive
  same-direction residual ratios >0.90 after pass 8) and **restores
  the best pre-refresh field** - an unconverged refresh keeps the
  converging sweep field rather than its own last iterate.
- `MOM_RF_KEEPNEG` quarantine: repaired negative-J cells keep their
  prior stored value instead of importing a neighbour's (the
  bluer-frequency fallback ratcheted it1 L 2.6e10 -> 2.7e20; removed
  `cd45cf2`).
- The `I_NEXT` store write cap `NONMONO_INEXT_CAP` (1e20) drops
  non-finite/exploding writes and counts them (`GS_CAPPED`).

## Newton-side damping

- PTC (`JFNK_PTC`): pseudo-transient continuation adds `x*s/dtau` to
  the residual - small `dtau` makes early iterates a bounded damped
  march, and `dtau` grows toward pure Newton as the residual falls
  (Knoll & Keyes 2004). See [[JFNK solver]].
- Hookstep trust region on the GMRES step; best-iterate restore and
  Lambda-mode SOLVEBA fallback on rejection; `exp(x)` clipped at
  +/-40 inside evals.
- `JFNK_MAX_EVALS` budget and `JFNK_MAX_STAGNATE` self-disable bound
  the cost of a bad activation.

## Related

- [[BA matrix and SOLVEBA]] - resacc internals, sentinel handling
- [[JFNK solver]] - PTC, trust region, clip semantics
- [[Field maturity and GS freeze]] - `GS_LAM_CAP`, `GS_LAM_TRANSIENT`
- [[Failure modes and guards]] - the blow-ups each cap answers
