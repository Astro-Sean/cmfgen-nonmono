# M-state comparison plan and results

Status: PLANNED ONLY - no runtime data. Order of operations when an
approved allocation and a captured run exist:

1. Replay-IM_NI vs deck-IM_NI (bitwise + distributional). Cheap; needs
   no new capture.
2. Replay M trajectory vs production KNOMG IDIR=2 rows (per-pass RES,
   RATIO, NOSC, applied omega) once captured.
3. Replay NP_M/IXIT_M vs production census row (KN_PA col 1) once
   captured.
4. KN_HOT_SPENT at ray entry: if production entered the ray with spend
   > 0 and replay hot-pass count differs, classify as budget-channel
   divergence (not input error, not FP).

Divergence taxonomy for any M mismatch: initial field / retained store
/ coefficients / donor-topology / boundary / ordering / FP order /
missing module state / hot-budget state.

Results table (all null until runtime):

| ray | stage | max_abs | max_rel | med_rel | p90 | p99 | n | excl | verdict |
|81|IM_NI| - | - | - | - | - | - | - | NOT_RUN |
|82|IM_NI| - | - | - | - | - | - | - | NOT_RUN |
|84|IM_NI| - | - | - | - | - | - | - | NOT_RUN |
