# Adaptive omega source audit

Phase 2.27B update: omf now uses a single multi-pass ray call and a strict
both-channel KNOMG2 schedule. Applied values are taken from KNAPPLY at the
ordered application boundary, with proposed/applied binary64 text and
applied bit encoding. Q does not consume production schedule rows. This
source repair compiles in isolated normal/debug builds but has no runtime
fidelity result. Earlier one-pass omf descriptions below are historical;
see replay_state_inventory.md and replay_lifecycle_comparison.md. Neither
legacy KNOMG nor seed-omega sensitivity establishes schedule independence.

Status: PARTIAL. Provenance, commands and source-stage table are in
../CMFGEN_PHASE227_ADAPTIVE_OMEGA_AND_RAY_CENSUS.md (inspection
2026-11-21T05:55:59Z, HEAD e28833fe plus dirty diagnostic edits).

OMEGA_RAY is scalar within a channel solve. Its pass-top snapshot precedes
field updates; march sites apply VOLD+OMEGA_RAY*(VNEW-VOLD). Pass-end
residual-ratio and oscillation logic calculate the next value, bounded by
controller branches at 0.45 and 1.0. Solve-local histories reset on entry.
The ordered call writes MERGE(2,1,INWARD): 1=P, 2=M. The unordered writer
instead uses ascending/descending frequency direction; the formats are
not interchangeable. Actual capture records have not yet validated either.

KNOMG records ref, ray, pass, direction/channel, used omega, next omega,
residual, ratio and oscillation counter. Used omega is sampled at pass top,
written later. Early-exit placement must be audited for complete coverage.
The armed binary's ES14.7 retains only eight significant digits, not exact
binary64 state. Source now uses ES26.17E3 (syntax checked, no binary rebuild).
Capture also needs an explicit solve/channel identifier,
and a record for every applied pass, including terminal passes. Retain base
omega and controller state separately. No raw history is available yet.

Q verification branches exit before the later ordinary telemetry call in
the inspected path; the previous assertion that Q rows necessarily append
was not established. A parser must reject ambiguous resets rather than
silently overwrite a pass.

The v2 scalar record includes BOTH configured KN_RAY_OMEGA and INEXT_RELAX;
terminal OMEGA_RAY is in the outcome record. No complete pass schedule is
stored in v2. INEXT_RELAX is a separate generation-level mechanism, not
an omitted scalar deck input. WRITE_STORE does not apply it a second time.

omf is exploratory, NOT a validated trajectory replay. It resets the full
ray solve on every pass and re-solves M. This pass adds explicit failure on
missing/malformed/nonfinite/out-of-range/wrong-ray/noncontiguous P history,
uses production NPO from the deck, supports long paths, and avoids unsafe
SIZE evaluation on an unallocated array. These source changes passed an
ifort syntax check only; existing executables remain unchanged. M history
injection and uninterrupted P-loop state preservation remain outstanding.
