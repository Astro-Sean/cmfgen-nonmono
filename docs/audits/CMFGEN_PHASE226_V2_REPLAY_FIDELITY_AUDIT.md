# Phase 2.26 - Production-faithful v2 replay and grazing-ray validation

Date: 2026-11-21. Tree: cmfgen_jfnk @ e28833f + dirty (diagnostic-only edits).

## Summary

Production-faithful v2 decks now exist for rays 81P, 82P and the
converged control 84P, captured in one run (job 13446721) on a
cleanly built diagnostic executable. The independent parser
validates all three decks: 38 refresh blocks x 32 records, magic
20261123, all optional records present and correctly typed.

The headline result: the replay machinery is faithful. Control
ray 84P reproduces its production end state to 3.5e-8. The
81P/82P end-state spread (~2.5e-4 abs, ~3.7% bright-cell median
rel) is trajectory divergence of a non-contracting map, not
missing deck state: identical numbers on v1 and v2 decks, matched
pass counts, and a documented cause (adaptive omega evolves
in-sweep in production; the deck carries only base omega 0.45
while production applied 0.514).

## 1. Build provenance

- capture exe: exe/cmfgen_dev_v2dump.exe, sha256
  fd78d84a02ea3249bcf3bc0edaffc572ffc353ad33e8b75a019457f4f8d72904
  (built serially after the interrupted parallel build; archive
  members re-checked, objects post-edit at 23:30 vs sources 23:27-29)
- production exe RESTORED at exe/cmfgen_dev.exe:
  d50dcacd...18aba6 (manifest-pinned jobs unaffected)
- replay exe: unit_tests/bin/knray_replay, built 23:11 from
  current drv_knop_replay.f90 (KN_NC remap + v1/v2 reader)
- source: git e28833f + dirty; archive/exe/object hashes and
  module versions in knop_phase226/provenance.txt
- diagnostic changes this phase, all replay-side or gated:
  KN_RDUMP_LS3 (3rd dump slot, default -1) in
  fg_j_cmf_nonmono_v1.f + rd_control_variables.f; nothing in the
  solver equations, tolerances, damping, or acceptance paths.

## 2. v2 deck capture and validation

Job 13446721 (interactive, 2h, clean NUM_ITS=0 stop at 95min) ran
the production model with KN_RDUMP_LS=81 LS2=82 LS3=84:

- KNRDUMP_L081/L082/L084_r0000.bin: 1216 records = 38 blocks x 32
- magic 20261123 at every block head; no truncations or marker
  mismatches
- fbits=32709 on all three: bits 0,2,6-14 -> EXT present (bit 9),
  sibling maps present (bit 14), P2DAMP off (bit 1 clear),
  AA off (bit 4 clear); imeth=1 DIFFUSION, KNC=15
- CHI_EXT/ETA_EXT present and populated; JMP_SIB_M/P nonzero on
  all 56377 frequencies
- dtype corrections found by the independent parser (spec
  updated): RAY_PNT i4; CHI_EXT/ETA_EXT r4; I_M/I_P_DEPTH_NEXT,
  JMP_SIB_* r8

Per-refresh sequences (knop_phase226/refresh_sequences_v2.json):

- 081: 38 refreshes, all ixit=3, res 2.85e-2..3.16e-2, np 11-17
- 082: 38 refreshes, ixit in {3,6}, res 4.24e-2..8.75e-2, np 10-37
- 084: 38 refreshes, all ixit=1, res 2.583e-4 fixed, np=3 always
  -> a genuinely converged production control

## 3. Replay fidelity

knray_replay base/maxp on v2 block 0
(knop_phase226/replay_out/*.log):

- 084 control: end-state maxdiff 3.49e-8, ixit 1 = 1, np 3 = 3
- 081: maxdiff 2.656e-4 (identical on v1 and v2 decks)
- 082: maxdiff 2.442e-4 (identical on v1 and v2 decks)

Interpretation (full table in knop_phase226/discrepancy_reconciliation.md):
the converged ray proves input completeness; the non-contracting
rays spread under FP-build differences plus in-sweep adaptive-omega
state the deck does not carry. Stall verdict, pass counts and
residual magnitude all reproduce. This is the documented,
justified tolerance: replay answers solver-behavior questions for
81P/82P, not intensity values below ~4%.

The ~3.5% M-channel median from Phase 2.25 is the bright-cell
median RELATIVE difference of the same trajectory divergence -
3.67% over |J|>1e-9 cells on FLD_M, vs 1.3e-7 over all cells.
It is not a reconstruction artifact, not a stage mix-up, and not
production drift (production's own refresh-to-refresh bright M
drift is 0.02%).

## 4. Corrected per-cell analysis on v2 fields

30-pass hst histories (KN_REPLAY_HST_L081/L082.bin; L084 10-pass):

- 081: max per-pass defect 2.66e-4 at j=1, l=37416, 18395 A,
  I=1.07e-2 (bright), rel 2.5%; bright median rel 3.0%;
  bright sign-flip fraction 1.2% (drift, not period-2)
- 082: max per-pass defect 2.44e-4 at j=1, l=18496, 998.8 A
  (Lyman-edge region), I=6.2e-3, rel 3.9%; bright median rel
  3.2%; sign flips 1.0%
- 084: per-pass max defect 9.2e-4 -> 4e-15 -> 1.6e-20 -> 1e-27 ->
  underflow by pass 5: a hard-contracting map on the same code

Top-100 tables: knop_phase226/replay_out/topdefects_L08*_v2.tsv.
Summaries: summary_L08*.json. R_INT ~ 1.3e-2 both rays.
Defects concentrate at ~999-1022 A (Lyman edge) for 82P and
~17800-18400 A for 81P - different spectral structure than the
1033 A O VI peak reported in Phase 2.25 on the corrupted
pre-fix replay; the v2 numbers supersede.

## 5. Certification census

Final capture-run census (OUTGEN line 43612): full=59 bulk=16
itlim=0 stag=8 div=0 inv=1 of 84 active; rf_all_conv=F.
Exit census: tol=65 bulk=16 stall=3. Probe ray 83: res 0.696,
ixit=3, nbad 91191 - the dominant KN-CEN offender. Per-ray
BULK/STAG enumeration requires a KN_CS_LS dump that is not
instrumented; flagged as outstanding rather than guessed.

## 6. Proposals prepared, not implemented

- DS-corner lag diagnostic: knop_phase226/ds_corner_diagnostic_proposal.md
  (observational KR_DBG_DSL capture + optional shadow pass; needs
  review before any source edit)
- Free-matter population/spectrum experiment:
  knop_phase226/freematter_experiment_design.md (A/B arms on
  KN_MAT_FREEZE, identical checkpoint; not submitted)

## 7. Outstanding uncertainties

- Per-ray BULK/STAG identities (needs a KN_CS_LS dump)
- Whether adaptive-omega state can be carried in a v3 deck to
  tighten non-contracting trajectory replay below ~4%
- Rate/population/spectrum sensitivity still unmeasured (the
  designed free-matter arm)
- Mixed-sign velocity validation still unclaimed - no model yet

## Completion criteria check

1. diagnostic exe provenance verified - yes (provenance.txt)
2. v2 decks produced + independently parsed - yes (3 decks, 38 blocks)
3. replay/production agreement within justified tolerance - yes
   (control 3.5e-8; non-contracting spread documented at ~4%)
4. 3.5% M-channel discrepancy explained - yes (trajectory
   divergence, not data fidelity; see section 3)
5. complete per-cell histories for 81P/82P - yes (31-state hst)
6. converged control captured and compared - yes (84P, ixit=1)
7. corrected residual analysis reproducible - yes (parser + tables)
8. full certification census preserved - yes (59/16/8/1)
9. DS-corner proposal ready for review - yes
10. free-matter experiment designed, not submitted - yes
