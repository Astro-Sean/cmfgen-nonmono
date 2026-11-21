# First-divergence report - production vs replay

Status: NO RUNTIME DATA. This is the ranked hypothesis list and the
evidence each candidate requires. The comparison machinery (KNSTAGE
20261127 records at REFRESH_START, seeds, PASS_START, raw march,
relaxed result, march result, FLDP_END, Q scratch states,
classification, commit, next-refresh seed) exists in source; the
production side has never been captured.

## Ranked candidate causes for the ~2.4e-4 ray-81/82 P maxdiff

R1. M input-state gap (leading hypothesis). The v2 deck captures
post-M state; replay recomputes M from retained store columns.
If replay-M differs from production-M, the P boundary data (IM_NI)
and inner boundary differ from pass 1 onward. Evidence required:
production-M field and IM_NI dump, or a deck written pre-M. Detectable
today: compare replay-M residual/npass/ixit (KN_PA_* col 2) against
the deck's recorded production M diagnostics when captured.

R2. Relaxation schedule. Refuted as sole cause by the verify defect
(~1.0 on 81/82 at the production end state under omega=1), but a
schedule mismatch could still contribute along the trajectory. Mode B
with a captured KNOMG2 schedule tests it; seed-sweep does not.

R3. Commit-state / next-refresh seed. Cannot explain an FLDP_END
difference (dump precedes cap) but can explain a stored-state or
later-refresh difference. Stage join on COMMIT_END / REFRESH_START.

R4. FP/ordering noise. Bounded by ray-84 control (3.5e-8) and O0/O3
driver identity; cannot be ruled out for grazing rays until the
solver-level O0/O3 field comparison runs.

R5. Controller-state side channel. Divergent NOSC/RATIO/exit decisions
(see controller_state_attribution.md) alter pass count, not the
applied omega chain. Detectable via KNAPPLY proposed-vs-applied and
IXIT mismatch.

## First-divergence procedure when capture exists

Join production vs replay stage rows on (refresh serial, ray, channel,
pass). Scan stages in order: seed -> M formal -> M boundary (IM_NI) ->
PASS_START -> per-pass KNAPPLY/KNOMG -> march result -> FLDP_END ->
verify -> classification -> commit -> next seed. First stage past the
documented FP tolerance (see two_channel_replay_metrics.json schema;
RELERR = |a-b|/max(|a|,|b|,REL_FLOOR)) names the divergence class.
A divergence at M formal or IM_NI confirms R1; only appearing at
PASS_START with matching M formal implicates seed/store state;
appearing first at KNAPPLY applied column would indicate an injection
defect (must not happen - records are bit-checked).
