# Replay experiments - corrected evidence status

Phase 2.27B: new isolated normal/debug solver-object builds and synthetic
Fortran IO fixtures are prepared, not executed. omf has been changed to one
uninterrupted ray call; earlier one-pass mode limitations describe the
superseded implementation, not a new test result. See runtime_test_results.md
for executable hashes and blocked controls. The new builds vary solver
optimization as well as driver optimization; historical driver-only
agreement does not predict their bytewise result.

Provenance and source positions: ../CMFGEN_PHASE227_ADAPTIVE_OMEGA_AND_RAY_CENSUS.md.
Historical inputs are job 13446721 v2 decks, block 0, with SHA256s in
../knop_phase226/deck_sha256.txt. No new captured omega history is available.

Historical baseline P maxdiff: 81=2.656368e-4, 82=2.442209e-4,
84=3.487919e-8. O0/O3 printed metrics agree; full-field bytewise comparison
and new regressions have not been executed in this continuation.

The seed-omega sweep tests initial values only. It does not settle whether
actual production omega histories change M median/p90/p99 or P trajectories.
No claim of schedule irrelevance or measured zero floating-point error is
supported. The production-flags-plus-history arm is not complete.

Concrete omf defects: the path was only 32 characters; filename text was
also parsed as a numeric loop count; incomplete histories silently fell
back to base omega; duplicates overwrote passes; a non-short-circuit AND
could evaluate SIZE on an unallocated array in other loop modes. These are
hardened in source; syntax check only, installed binaries not replaced.
Repeated one-pass full solves still reset solve-local histories and repeat
M. That algorithmic limitation must be addressed before calling omf faithful.

Deck FLDP_END is pre-Q/pre-certification/pre-commit P, not a retained store.
KR_DBG_FLDP is also a local solve field. The deck additionally carries the
post-M field, whereas baseline replay computes M again. Locate the first
divergent stage before testing alternative relaxation or cap operations.
Historical verify mode is not independently validated as a one-pass map on
an identical full state; do not call its scalar result a fixed-point error.

Required output contract for each matched refresh/pass/stage/channel:
max abs; max/median/p90/p99 relative; cell count and exclusions; residual;
class and reason; field/input/executable hashes. Use a declared common
normalization, e.g. abs(test-reference)/max(abs(reference), declared floor),
with a separately reported bright subset. Choose and record the floor
before comparing arms; do not tune it to reduce discrepancy. Nonfinite
values are failures, not silently excluded cells. M and P stay separate.
No unavailable percentile values are filled with zeros.

Regression command pending allocated execution: run each O0/O3 executable
on each immutable deck in separate new working directories, hash and compare
KN_REPLAY_END.bin bytewise and numerically, then run malformed-sidecar cases
with nonzero-exit assertions. Require ray 84's historical 3.487919e-8 result
to agree to the original printed precision and document any deviation.
