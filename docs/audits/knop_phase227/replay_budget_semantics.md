# Replay budget semantics (Phase 2.27F)

Status: documented; none implemented beyond existing Modes A/B.

## Mode A - fresh-budget ray replay (exists)

drv_knop_replay pokes deck state and calls KNOP_RAY_SOLVE_V1 in a
serial driver. KN_HOT_SPENT resets per refresh call site... NOTE: in
the driver, the refresh prologue at 8905 does not run per ray the
same way; the driver replays a single ray with its own fresh
counter (0), giving the ray up to min(120/sweep, 480) hot passes.
Production may have allocated it fewer (contention) or zero.
Mode A is a ray-local diagnostic; NOT a production-refresh
reproduction. For rays 81/82 (stall-class, hot-mode candidates) its
pass count can exceed production's.

## Mode B - recorded-schedule replay (exists)

KNOMG schedule injected at pass top; strict validation fails on any
pass beyond the recorded count. The recorded schedule length
implicitly encodes production's pass count INCLUDING hot passes
(armed-ray KNOMG rows cover every pass). Mode B therefore reproduces
production's pass count for armed rays without needing the budget -
the schedule itself is the binding constraint. Replay's own
KN_HOT_SPENT still increments but cannot truncate earlier than the
schedule (fresh 480 >= any recorded length <= 240/ray max).

## Mode C - production-allocation replay (not implemented)

Minimum sufficient state: per (ray, sweep) executed hot-pass count
from production telemetry (the END rows above), enforced as a
per-sweep NHOT ceiling equal to production's executed count. The
global interleave order is NOT required: atomic monotone consume
means each sweep only ever observed the cap through its own consume
verdicts; replaying the same per-sweep executed count reproduces
every decision the budget influenced.

Required telemetry first: the Part-5 buffers. Do NOT claim Mode C
exists until per-sweep executed hot counts are captured and the
replay driver enforces them.

## Interaction with the three-ray deck set

A 3-ray replay cannot reproduce budget competition; that is fine
for B (schedule-bound) and C (allowance-bound). It is a known
limitation of A only.
