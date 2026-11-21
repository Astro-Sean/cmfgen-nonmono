# Part G - free-matter post-flag validation design (NOT SUBMITTED)

Purpose: the A1/A2/B1/B2 arms all inherit KN_MAT_FREEZE=T, which holds
populations by design and makes the release/freeze transcript
bookkeeping-only. A free-matter post-flag arm is needed to test the
normal release controller and the final population-acceptance
sequence. This is a design record only - submit only after the
pre-checks below are re-verified at submission time.

## Design

Clone the B2 deck `dessart_m15_2p82d_KNOP_B2_POSTFLAG` exactly, with a
single VADAT change:

    KN_MAT_FREEZE: F        # was T

Everything else identical: same model, same initial iteration, same
probe/diagnostic arming. Rationale: B2 is the post-flag frozen-matter
arm, so a B2-clone with free matter isolates the matter-freeze factor
against an already-staged, hash-verified sibling.

## Verified inputs (as of 2026-11-20)

| item | value | status |
|---|---|---|
| cmfgen_dev.exe | exe_frozen/58a699d0f37864e6_abe4725807474143 (sha256 58a699d0 f37864e6..., source commit 595f9ae7) | pinned path verified; hash re-verify at submit |
| VADAT_INIT | sha256 f238af97 9b9613d3... | identical to all arms |
| SN_HYDRO_DATA | sha256 347eda19 ce4f645c... | identical |
| SCRTEMP (pops) | sha256 5704dc48 6d6f5a16... | iso_M it157 donor - identical initial populations |
| FGINEXT_GS (field) | sha256 48077bb3 65c20760... | identical initial radiation field |
| MODEL | sha256 c187873d 6db028fc... | identical donor copy |
| warm_restart_from | cmfgen_validation/m15_2p82d_KNOP_iso_M/cmfgen/2.82days | same donor as B2 |
| KN_MAT_FREEZE | F | the only intended diff |

Note on the executable choice: the live exe path currently holds
d50dcacd (Phase 2.22 dev build, contains KN_RDUMP_*). The pinned
58a699d0 is the correct binary for this experiment because (a) it is
B2's exact binary - the comparison isolates KN_MAT_FREEZE alone;
(b) it is a frozen path, immune to live-tree rebuilds. If a separate
replay-deck experiment is wanted later, that is a different job on
the d50dcacd-class binary with KN_RDUMP_LS armed.

## Required measurements

- population release history (GS_STATE release events, counts)
- freeze/hold history (same counters, now meaningful)
- FG_RF_ALL_CONV per refresh (post-flag certification)
- GS_RF_DCHG per refresh
- population MAXCH per great iteration (demand AND applied)
- final convergence status + termination reason
- ray-exit census and commit verdicts (same extraction as A1)

## Comparison logic

- vs B2 (frozen, post-flag): isolates the release controller + matter
  coupling under identical field physics.
- vs A1/A2 (frozen, pre-flag): NOT a clean comparison for population
  behaviour (different binary + freeze state); still usable for
  field-trajectory context.
- The ten-release question: needs a matched frozen-matter control
  (B2) reaching the same iteration range - B2 is pending, so the
  comparison arm exists but has no data yet.

## Preserve

All pending jobs untouched (A1/A2 running, B1/B2/diagcap/obs queued).
Existing staged dirs not modified. No new job submitted by this
document.
