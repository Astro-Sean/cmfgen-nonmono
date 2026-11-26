# KNOP remediation WP-4: store line search

26-Nov-2026. Fourth work package of the Knop non-monotonic solver
remediation (spec: `checklists/knop_solver.md`). WP-3 made the
maturity gate read the committed state; WP-4 replaces the fixed
`FG_COMMIT_NC` store bound on non-certified commits with a
residual-descent line search that uses the WP-3 committed verify
as its merit function.

## The gap WP-4 closes

On a non-certified exit the legacy path wrote a provisional step
bounded by the fixed `FG_COMMIT_NC` (5%): the store ratchets
toward the fixed point at <= NC per refresh. On the fixture the
fresh field after a stalled solve already reads ~3e3 better than
the retained store; the cap throws that descent away and needs
O(100) refreshes to close it. The other direction is worse: an
expansive step whose fresh field is *worse* than the retained
store still commits 5% of it -- the cap is a magnitude bound,
not a descent test, so it ratchets bad states too.

## Changes

### `KN_STORE_LS` -- residual-descent backtracking (VADAT, default T)

Inside `KNOP_RAY_SOLVE_V1`, the non-certified commit branch now
runs a line search on the committed-state residual:

1. `KN_LS_VER` (the WP-3 verify, refactored to take the store as
   an argument) measures the baseline residual `R_LS0` at the
   retained store.
2. For each rung alpha in {1, 0.5, 0.25, 0.125, NC-floor}:
   `KN_LS_BLEND` builds a scratch candidate
   `stored + alpha*(fresh - stored)` for both march directions,
   and `KN_LS_VER` measures its committed residual `R_LS1`.
3. The first rung with `R_LS1 < R_LS0` commits that blend and
   records `FG_RAY_LSA(LS)=alpha`. `FG_COMMIT_NC` survives only
   as the rung floor (a rung below it is skipped, and
   0 < NC < 0.125 appends itself as rung 5).
4. If every rung measures and none descends, the previous store
   is retained: `FG_RAY_LSA=-2`, `FG_RAY_VC=2`, `FG_LS_REJ` +1.

Certified commits (`WRITE_STORE`), the `FG_COMMIT_JUMP` veto, and
the non-finite/cap guards are unchanged; `KN_STORE_LS=F` restores
the fixed-cap `WRITE_STORE_CAP` path (including the NC=0 retain)
exactly.

Census: `FG_LS_HIST(5)` bins the accepted rung, `FG_LS_REJ`
counts rejections, `FG_LS_R0`/`FG_LS_R1` sum residual
before/after -- printed in the refresh epilogue when any LS ran.
Per-ray verdicts land in `FG_RAY_LSA` (accepted alpha, -1 never
ran, -2 rejected).

### Helpers added in the solve's CONTAINS

- `KN_LS_VER(STM,STP,QM,QP,NB,NT)`: the WP-3 committed verify
  with the store passed in (was `KN_COMMIT_VER`, same body).
- `KN_LS_BLEND(FLD,STORE,CAND,ALPHA,DCH)`: candidate store +
  change-norm for one march direction.
- `SEED_FLD2(FLD,STORE)`: `SEED_FLD` variant that reads an
  arbitrary store array instead of the module store (used by
  `KN_LS_VER` to rebuild each candidate's committed point).

## Fixture measurements (test_knop_commit_gate)

| case | commit path | line search | verdict |
|---|---|---|---|
| lsok | stalled ordered solve, LS on | baseline 1.4e5 -> alpha=1 accepts, CRES=3.95e-2 | full step |
| lsrj | solve1 accepts alpha=1; solve2 over-relaxed + capped | all rungs worse than 3.95e-2, LSA=-2 | rejected, fresh1 kept |

`lsok` is the spec's forward direction: the same stall that
`nccm` parks at CRES=1.4e5 under the 5% cap now commits the fresh
field in one refresh -- committed residual 3.95e-2, matching the
solve's own iterate. `lsrj` manufactures the spec's expansive
step: solve 1's LS parks fresh1 in the store; solve 2 runs
`KN_RAY_OMEGA=1.36` pass-capped after one march, so its iterate
exits non-certified kicked off the fixed point (RESID=1.0). Every
blend toward it reads worse than fresh1's 3.95e-2, the search
rejects, and the committed metrics still report fresh1 -- the
expansive map cannot ratchet.

## Design notes

- On this fixture's contracting map the iterate seeds from the
  retained store, so the fresh field is always nearer the fixed
  point and alpha=1 always descends -- rejection needs an
  expansive step, which `KN_RAY_OMEGA>1` manufactures (the raw
  over-relaxation grows the orbit; the pass cap catches it before
  the divergence verdict, so the ray still reaches the LS).
- An unmeasurable baseline (`R_LS0 < 0`: verify failed, NOPROBE,
  or flag off inside the verify) falls through to the legacy NC
  path untouched -- the search never blocks a write it cannot
  judge.
- `drv_knop_ray` case 80 (`cer0`, `FG_COMMIT_NC=0`) now disarms
  the LS explicitly: under the new default the NC=0 case takes
  the search (floor 0 -> alpha=1 accepts) rather than the legacy
  retain the case exists to cover.
- `nccm`/`nclg` in `drv_knop_ordering` also disarm `FG_STORE_LS`
  -- their signature is the *legacy* provisional write the WP-3
  gate flag test needs.

## Verification

- `test_knop_commit_gate` 18/0 (was 10/0): WP-3 checks unchanged
  plus the `lsok`/`lsrj` accept/reject signatures and census.
- `test_vadat_keys` 124/0 (was 122/0): default T, explicit F
  lands.
- `test_knop_ray` 50/1 (T18 pre-existing), `test_knop_ordering`
  50/0, `test_knop_operator` 32/0, `test_knop_oneway` 8/0,
  `test_knop_p2_state` 19/0, `test_knop_replay_p2` 29/0,
  `test_knop_spv` 93/0, `test_iter_conv_gate` 18/0,
  `test_conv_exit` 53/0, `test_fg_gs` 36/0,
  `test_mom_aa_state`/`test_mom_augres`/`test_mom_gate` green.

## Verdict

VERIFIED. Non-certified commits now descend the committed
residual: the stalled fixture store reaches ~4e-2 in one refresh
where the 5% cap needed O(100), and a manufactured expansive step
is rejected outright (LSA=-2) with the prior store retained.
Next: WP-5 (signed adaptive relaxation).
