# Hot-budget telemetry design (Phase 2.27F)

Status: designed, not implemented. Diagnostic-only; must not change
budget values, decisions, or scheduling.

## What must be captured

Per refresh: KN_HOT_SPENT initial (0) and final (may exceed 480 by
overshoot), plus for every budget event:

- refresh serial (KN_DIAG_REFRESH), ray LS, OMP thread no., channel
  sweep (M=2/P=1), pass NPO, NHOT after increment
- KN_HOT_SPENT before/after the consume
- L_HOT_OVER verdict; IXIT=4 when fired
- TRY_HOT admissions: LS, sweep, NPO, W_LO/W_HI window, and whether
  the admission consumed on an exhausted budget (stale-read entry)
- sweep end: executed hot count (NHOT at exit), exit reason

## Mechanism

Per-thread lock-free buffers, not writes in the hot loop:

- The consume already sits inside CRITICAL (KN_HOT_SPD). Append the
  event record to a thread-private buffer INSIDE the same critical:
  the buffer append is then ordered with the counter, giving the
  true synchronised decision order for free - no extra atomics, no
  added serialisation outside the existing critical extent (~ns).
- TRY_HOT admission logging is outside the critical (like the read);
  mark those records `entry` vs `consume`.
- Buffers flush after the PARALLEL DO region (per-thread files
  KNBT_T<tid>.txt or one merged keyed by a seq assigned at flush
  is NOT acceptable for ordering - keep the critical-order seq: an
  INTEGER*8 seq incremented inside the same critical is already
  ordered by construction; it adds one shared increment inside an
  existing critical, scheduling impact negligible).
- Flush target: `./data/` per candidate layout; fields for rays
  81/82/84 still come from the existing KN_RAY_DUMP/stage writers.

## Scheduling-perturbation assessment

- Added work inside KN_HOT_SPD: one buffer store + one seq
  increment; same order as the existing increment. Critical hold
  time grows <1 us/event; hot events are <= ~600/refresh. Verdict:
  does not meaningfully alter DYNAMIC scheduling.
- No new global critical around the ray solve; no serialisation.
- Alternative rejected: OMP ATOMIC seq for ordering - seq order
  would then be decoupled from spend order and could misattribute
  the boundary allocation.

## Record format (text sidecar, mirrors KNOMG style)

```
RF LS TID IDIR NPO NHOT SPENT_BEFORE SPENT_AFTER OVER IXIT
```
plus `ENTRY LS TID IDIR NPO W_LO W_HI` and
`END LS IDIR NHOT EXECUTED IXIT` rows.
