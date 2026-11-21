# Phase 2.23 experiment arm design (not submitted - QOSGrpCpuLimit)

All arms clone the diagdump_i deck (interactive partition where
possible) so results share the A2-identical initial state and the
frozen pre-flag solve semantics. Pending-job policy: no submissions
while B1/B2/diagcap/obs are queued on QOSGrpCpuLimit except small
interactive diagnostic captures, which use a different QOS.

## Experiment A - baseline reference (already exists)
- A1/A2 frozen-matter baseline + diagcap_i/diagdump_i PSD captures.
- No new job needed; serves as the control for B/C.

## Experiment B - KN_ORD_P2DAMP=T only
- Config: diagdump_i deck + KN_ORD_P2DAMP=T (P2N=2, P2REL=3 defaults).
- Everything else identical: same checkpoint, same armed rays.
- Record: res/ixit/verdict per ray per refresh; NP2_CNT damped-commit
  counter (printed ~9833 under LUER); KN_CS census; capped count;
  GS_RF_DCHG; q residual census.
- Success: rays 81P/82P reach ixit 1-2 AND q_resid drops below 1E-2
  without verdict-3 commits persisting. Failure modes to watch:
  damped commits that lower RES while q stays high (metric artifact),
  damped cells releasing and re-entering a cycle.

## Experiment C - FG_ORD_AA=T only
- Config: same + KN_ORD_AA=T (AA_M=2 default).
- Record: KN_AA_APP/VETO/REJ/SNG counters per refresh, same
  certification outputs as B.
- Success: same as B plus rollback rate low (REJ << APP); a high REJ
  rate means AA mixes are harmful on this state.

## Experiment D - combined (only if B or C helps)
- Same + both keys. DO_AIT suppressed by AA anyway.

## Purity note
- B/C/D need only the pinned pre-flag binary 93364cac IF it carries
  the option code: strings check required before submission. P2DAMP
  and AA exist in current source; the pinned trees' source also
  carries them (verified by grep on exe_frozen trees) - but the
  compiled binary must be checked for the symbols.
- If a v2 deck is wanted for replay, the armed run must use a build
  carrying the 20261123 writer (current source = uncommitted Phase
  2.22 work) - a different executable than the A/B baseline; flag it
  as diagnostic-capture, not a pure config A/B.

## Free-matter arm (recap from Phase 2.22 design)
- Clone B2 deck: pinned 58a699d0, KN_MAT_FREEZE=F, identical donor
  hashes (VADAT_INIT f238af97, SN_HYDRO_DATA 347eda19, SCRTEMP
  5704dc48, FGINEXT_GS 48077bb3). Record release/freeze history,
  FG_RF_ALL_CONV, GS_RF_DCHG, population MAXCH, final status.
  NOT submitted: B2 remains its frozen control; the same QOS limit
  blocks it.
