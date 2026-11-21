# Live M-state capture design

Status: LARGELY ALREADY ARMED. Audit conclusion: the existing capture
points already sit at the requested boundaries; the gaps are specific.

## Required capture points vs what exists

| Point | Requirement | Existing | Gap |
|---|---|---|---|
| A: pre-M | initial field, stores | REFRESH_START_SEED_VIEW + RAY_SOLVE_START stage rows at 10694-99 | none |
| B: post-M | M output, RES_M, IXIT_M, IM_NI, M boundary data | IM_NI in deck phase-1 (10790); KN_PA col 1; KNOMG IDIR=2 omega history; KNAPPLY per-pass | full FLD field and M end-stage row absent |
| C: pre-P | IM_NI, P seed, march inputs | deck phase-1 + NEXT_REFRESH_SEED/RAY_SOLVE_START (10773-76) | none |
| D: post-P | FLDP_END, exit state | deck phase-2 + FLDP_END_DUMP (10819) | none |

## Implemented addition (diagnostic-only)

`M_SOLVE_END` stage row added after the IM_NI assignment (~10739):
records `fld(ni:ni,1:ncf)` - bit-identical to IM_NI - with RES_M and
terminal OMEGA_RAY in the scalar slots and NP_M as the pass field.
Syntax-checked clean; build_v5 carries it. Production stage files now
contain the live boundary value AND the deck contains it, so replay
comparison has two independent production references.

Still missing (deferred): full M field dump (M_FLD_FULL) if IM_NI
comparison proves insufficient, and KN_HOT_SPENT-at-ray-entry (reveals
whether production's shared budget was already consumed).

## Comparison metrics

Same RELERR denominator as two_channel_replay_metrics.json
(REL_FLOOR=1e-300); integers/logicals exact-equality; counts and
exclusion reasons reported; M and P channels separate. IM_NI
comparison: replay-computed vs deck-recorded, per frequency index,
exact bit comparison plus distributional metrics - the cheapest and
highest-value single check when runtime opens.

## Blockers

All runtime. No production capture exists with these records; no
allocation for replay execution. Held jobs unchanged.
