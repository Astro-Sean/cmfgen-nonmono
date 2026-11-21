# KN_HOT_SPENT OpenMP audit (Phase 2.27F)

Status: source-reviewed.

## Parallel context

- Ray loop: `!$OMP PARALLEL DO NUM_THREADS(N_TH) SCHEDULE(DYNAMIC)`
  at 9129 (FG_KNOP_REFRESH_V1). DYNAMIC: rays are handed to threads
  as they finish, so which rays run concurrently - and hence which
  rays are in HOTMODE at the same wall-clock window - is
  scheduling-dependent run to run.
- KN_HOT_SPENT is a module SAVE variable: one shared instance.
- W_LO/W_HI/NHOT/HOTMODE/L_HOT_OVER live in KNOP_RAY_SOLVE_V1 or its
  sweeps: per-call, thread-private. No race on these.

## Synchronisation inventory

| Access | Site | Sync | Class |
|---|---|---|---|
| spend + over-cap test | 11676-79, 12400-03 | CRITICAL (KN_HOT_SPD) | atomic test-and-consume |
| eligibility read | 15024 (TRY_HOT) | none | unsynchronised shared read |
| KN_HOT_CNT update | 10754-58, 10844-49 | CRITICAL (KN_HOT_UPD) | atomic |
| reset | 8905 | before parallel region | safe (implicit barrier at PARALLEL DO entry) |
| JFNK push/pop | 2621, 2674 | serial driver context | safe |
| summary read | 9526 | after parallel region | safe (implicit barrier) |

## Classification: B

Atomic counter update but non-atomic test-and-enter.

- The authoritative consume is correct: inside CRITICAL, increment
  and L_HOT_OVER are one atomic unit. Two threads cannot both claim
  the last unit; whoever increments past 480 sees L_HOT_OVER=T and
  exits IXIT=4 before marching.
- The entry test `IF(KN_HOT_SPENT .GE. KN_HOT_TOTAL_PASSES)RETURN`
  at 15024 races with concurrent spends. A ray can read a stale
  below-cap value and enter HOTMODE when the budget is actually
  exhausted. Bounded consequence: the next pass-top consume charges
  one unit and exits it before any restricted march - at most one
  wasted budget unit per stale admission, no extra solve passes.
- INTEGER*8 read is practically atomic on x86-64; even a torn or
  stale value only affects the advisory entry, never the spend.

## Determinism verdict

NOT deterministic. With SCHEDULE(DYNAMIC) and N_TH > 1, identical
input + executable can allocate different per-sweep hot allowances
across runs whenever >= 2 sweeps compete for the budget tail. For
rays whose hot run ends on the per-sweep 120 cap without global
contention, allocation is deterministic. Whether rays 81/82's
production allowance was contended is knowable only from telemetry
(per-ray executed hot count vs 480 total consumed).

Replay note: production itself is an ensemble over thread schedules
for contended hot rays - 'replay fidelity' can at best target the
captured realisation, not a unique production trajectory.
