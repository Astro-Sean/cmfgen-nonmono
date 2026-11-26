# KNOP remediation WP-3: committed-state gate

26-Nov-2026. Third work package of the Knop non-monotonic solver
remediation (spec: `checklists/knop_solver.md`). WP-2 made the
donor freshness contract measurable; WP-3 makes the same contract
apply to what the store actually keeps, so the maturity gate can
never certify a stale or bounded write.

## The gap WP-3 closes

Until now the maturity gate read `GS_MAXCHG` -- the solve's own
residual -- and `GS_NBAD_EFF` -- the best-pass bad-cell census.
Both measure the solve's iterate, not what `WRITE_STORE` or
`WRITE_STORE_CAP` wrote. A non-certified exit commits a
provisional `FG_COMMIT_NC`-bounded step from the previous store;
on a cold or stale store that write parks the committed field far
off the fixed point while the best-pass census reads clean. The
gate then certifies on the solve's report while the store holds
uncommitted error -- exactly the period-2/divergence scenario the
spec calls "best pass < tol but committed pass > tol".

## Changes

### `KN_GATE_COMMIT_RES` -- committed-state verify (VADAT, default T)

After each ray's commit decision writes a store (`WRITE_STORE` on
a certified commit, `WRITE_STORE_CAP` on a provisional one, or a
rejected write), `KN_COMMIT_VER` inside `KNOP_RAY_SOLVE_V1`
rebuilds the committed point -- the production field at unstored
cells, the committed store at stored cells, via `SEED_FLD` --
and runs one L_QP verify pass over it:

- `FG_RF_COMMIT_RES`: the committed max-norm residual.
- `FG_NBAD_CM_RAY`/`FG_NTOT_CM_RAY`: per-ray committed bad-cell
  census (module SAVE arrays on the `FG_RAY_VC` pattern, written
  under the ray's own LS slot), summed to `FG_CM_NBAD`/
  `FG_CM_NTOT` after the refresh loop.

`GET_FG_GS_MATURITY` exports both. The `GS_MAX_ALL`/`GS_NBAD_EFF`
feed in `solve_for_pops` reads the committed metrics when the
census exists, and falls back to the pre-commit residual and
best-pass census only when the verify never ran (flag off, or no
store written that refresh -- nothing committed, nothing to
verify). Sentinels are -1 on both sides of the flag so a missing
verify can never read as a clean census.

## Fixture measurements (test_knop_commit_gate)

| case | commit path | committed verify | verdict |
|---|---|---|---|
| capok | certified full write | CRES=4.2e-16, CM_NB=0/55800 | certifies |
| nccm | NC-bounded provisional write | CRES=1.4e5, CM_NB=55552/55800 | holds |
| nclg | same write, flag off | CRES=-1, slots -1 | legacy only |

`nccm` is the spec's wedge: the ordered march stalls (RESID=4e-2,
NXIT=6), the exit is not certified, and the provisional commit
bounds the write to 5% of the 1e-6 store seed -- the committed
store sits ~1e-6 while the map's fixed point is ~4. The committed
verify reads that gap directly: CRES=1.4e5 and 55552 committed
bad cells, against a legacy census of zero. `capok` proves the
honest path unaffected: a certified commit stores the fixed point
and the verify certifies it. `nclg` proves the flag disarms
cleanly.

## Design notes

- One L_QP pass per committed ray per refresh; skipped on
  NOPROBE, on rejects that wrote nothing, and under the flag.
- The verify reuses the pre-commit save/restore machinery (the
  `Q_M`/`Q_P` block at the sweep tail); `KN_COMMIT_VER` is
  contained in the solve so the scratch arrays are host-visible.
- `NONMONO_INEXT_CAP` applies inside the sweep itself, so a
  cap-clipped committed state is a fixed point of the capped
  operator -- the verify cannot flag it. The honest wedge for a
  committed-off-fixed-point store is the NC-bounded write; the
  fixture uses that path (the cap construction was dropped after
  it measured a consistent 0).

## Verification

- `test_knop_commit_gate` 10/0 (new): the table above.
- `test_vadat_keys` 122/0 (was 120/0): default T, explicit F
  lands.
- `test_knop_ordering` 50/0, `test_knop_ray` 50/1 (T18
  adaptive-omega diag remains the only failure, pre-existing),
  `test_iter_conv_gate` 18/0, `test_conv_exit` 53/0,
  `test_fg_gs` 36/0.

## Verdict

VERIFIED. The maturity gate now reads the committed state, not
the solve's report: a provisional write that parks the store off
the fixed point holds the gate at CRES=1.4e5 with a full
committed census, where the legacy best-pass census reported
clean. Next: WP-4 (the spec's remaining invariants).
