# Phase 2.27B lifecycle comparison

Status: SOURCE-REVIEWED, isolated normal/debug compile/link PASS WITH
LIMITATIONS. No new production trajectory has been captured. No runtime
fidelity claim is made.

## Previous versus revised replay

The previous omf set KN_RAY_MAXPASS=1 and repeatedly called the complete ray
solver, writing P output directly into sparse store slots between calls.
That restarted controller/stopping/acceleration state, recomputed M and its
P boundary, invoked certification/commit repeatedly, and bypassed P commit
capping through caller overwrite. It was a different operator.

Revised omf loads one checked schedule and calls KNOP_RAY_SOLVE_V1 ONCE with
the original multi-pass limits and termination policy. A diagnostic-only
hook at ordered pass entry substitutes the recorded applied omega AFTER
the hot-budget precheck and BEFORE the first field update. M and P each
have their normal uninterrupted loop. Q verification does not consume
production schedule records. Residual, oscillation, AA/P2 and stopping
histories remain inside the same production call. The controller still
calculates the next proposed omega; KNAPPLY records both proposed and
applied values. No controller formula, bound or acceptance test changed.

Commands after allocated testing is separately authorized:

    <isolated>/knray_replay <deck> stage
    <isolated>/knray_replay <deck> omf <canonical-history> <refresh>

Run each command in a NEW isolated directory. The driver maps the deck ray
to solver ray 1 but labels diagnostic output with the original ray ID.
The v1/v2 parser is unchanged. Output now calls the deck omega
`deck_omega`, not an inferred replay terminal value.

## Canonical schedule

Header: `KNOMG2 original_ray refresh np_P np_M`.
Rows: `refresh ray channel direction pass proposed applied applied_bits`.
Channel 1=P, 2=M; direction 0=ordered spatial march with frequency order
selected locally by GAM. There is no single expanding/compressing sign for
an entire mixed-sign ray. The old unordered ASC/DSC telemetry is rejected.
`applied_bits` is the exact 16-digit hexadecimal IEEE binary64 encoding.
All required rows must be finite, in range, correctly ordered within each
channel and present exactly once. Extra records/tokens, identity mismatch,
wrong direction and a decimal-to-bits mismatch fail explicitly. Runtime
consumption must equal both declared channel counts. No base fallback.

KNAPPLY_Lxxx.txt records each actually entered ordered production march at
the application boundary, including passes that later terminate by tolerance
or divergence. It uses ES26.17E3 plus bit encoding. Legacy KNOMG remains for
controller-end context but is not sufficient as the authoritative schedule.
Use phase227_records.normalize_omega with expected per-channel counts from
matching authoritative capture records. A bare observed maximum pass does
not prove that a missing terminal row never existed.

## Stage files and meanings

KNSTAGE_Lxxx.bin uses paired native Fortran unformatted records. Header:
magic 20261127, character(64) stage, KN_DIAG_REFRESH (integer kind 8),
ray/channel/pass and two default-integer dimensions, residual and context
scalar (LDP). The next record is a Fortran-order LDP field array.
KN_DIAG_REFRESH increments once per eligible armed outer refresh and is
process-local, never reset by the great-iteration ledger. LEDG_N_RF resets
in LEDG_ITER_RESET and therefore cannot identify whole-run appended blocks.
KNRDUMP and census filenames now use the full diagnostic serial with I0,
without modulo wrap, and KNOMG/KNAPPLY/KNCELL/stage rows share that serial.
The v2 binary deck record schema is unchanged. omf uses the canonical
history's explicit serial; baseline stage mode labels its one local ray
invocation as serial 1 and must be aligned by pinned input identity, not
silently treated as an outer-refresh reproduction.
KNCELL_Lxxx.txt holds selected-cell events at row min(2,NI), frequency
min(18497,NCF), with raw VNEW or the actually assigned relaxed value.

| Stage | Actual state represented |
|---|---|
| REFRESH_START_SEED_VIEW | per-ray M seed view, NOT a global refresh-entry checkpoint |
| NEXT_REFRESH_SEED | actual SEED_FLD result read on this refresh, before solve |
| RAY_SOLVE_START | seeded M or P field at channel entry |
| PASS_START | live pass-entry field; Q also has surrounding VERIFY labels |
| FORMAL_RESULT | selected interior cell's raw VNEW before sanitize/cap |
| RELAXATION_APPLIED | selected cell's actual assigned value after ordinary/P2 branch |
| FIELD_UPDATE | full field after march, before Aitken/AA and termination handling |
| FLDP_END_DUMP | P solve field before Q, classification and commit |
| VERIFY_START / VERIFY_END | scratch field before/after independent verification |
| CERTIFICATION | unchanged solve field; context scalar carries authoritative class |
| COMMIT_CAP | OLD sparse retained store on entry to capped branch; no claim of an applied cap yet |
| COMMIT_END | actual post-decision sparse retained store, including rejected/no-op decisions |

Sparse store records have depth-slot x (NCF-1) shape and frequency slot l-1;
full ray fields have NI x NCF shape. They must not be compared as identical
arrays or stages. CERTIFICATION records Q as residual context, not the
production pass residual. COMMIT_END is emitted on the ordinary successful
finalization path; divergence/heal and allocation-failure paths still need
explicit stage coverage. FORMAL_RESULT/RELAXATION_APPLIED cover selected
interior cells, not boundary rows or a hypothetical simultaneous raw map.
A true global REFRESH_START now records selected rays' retained sparse
stores after refresh eligibility guards and before coefficient rebuild.
The ray-only replay does not execute this outer refresh setup; it cannot
claim global-refresh equivalence from channel-entry agreement alone.
These limitations must be closed or explicitly excluded before all stage
requirements pass.

## Acceptance still outstanding

Run baseline and recorded-history modes with both isolated optimization
variants. Reproduce ray 84's approximately 3.487919e-8 P max difference to
its original printed precision; compare O0/O3 fields bytewise. For 81/82,
compare post-M first, then identical P pass/stage fields. Relative metrics
must use a fixed documented denominator and floor and report counts and
nonfinite failures, separately for both channels. Median/p90/p99 results,
production per-pass agreement and next-refresh seed comparisons are NOT
available yet. Same-refresh commit capping cannot explain a difference
between two pre-cap FLDP_END fields.

## Equivalence checklist (2.27C review)

| Requirement | Production site | Replay | Source verdict | Runtime gate |
|---|---|---|---|---|
| M solved once, before P | SWEEP_FLD_ORD INWARD, RES_M/IXIT_M before P results | same single call, same order | YES | pending |
| M boundary data retained through P | IM_NI + stores into march (~11413) | deck-poked identical inputs | PARTIAL: live-vs-deck IM_NI unverified | pending |
| P consumes same M info | march reads IM_NI/stores | same modules/data | YES in source | pending |
| Controller history across passes | sweep-local RES/RATIO/NOSC | alive across single call | YES | pending |
| Stopping at same point | pass-end checks | unchanged | YES in source | pending |
| omega substitution only | pass boundary OMEGA_RAY | KN_RAY_OMEGA override, adaptation still runs | YES in source; KNAPPLY provenance not yet produced | pending |
| Replay-only branch affects numerics | none intended | verify at runtime | UNVERIFIED | pending |
| Persistent state coverage | replay_state_inventory.md | module init verified; per-pass solve-local equivalence assumed from shared routine | PARTIAL | pending |

