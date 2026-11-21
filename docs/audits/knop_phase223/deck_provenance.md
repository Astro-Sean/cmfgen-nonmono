# Phase 2.23 replay-deck provenance and replay status

## Decks (diagdump_i job 13433532, interactive, TIMEOUT after 2h)
- KNRDUMP_L081_r0000.bin (14.0 GB, ~34 refreshes appended)
- KNRDUMP_L082_r0000.bin (13.9 GB)
- Armed via VADAT KN_RDUMP_LS=81, KN_RDUMP_LS2=82 in
  m15_2p82d_KNOP_diagdump_i/config_resolved.yaml; repo dev exe (d50dcacd).
- Format v1 (magic 20261122). FBITS=0x1C5 both decks:
  FREQ_MATCH + ORD_SWEEP + BULK_CONTR + ADAPT_OMEGA_V2 + PLATEAU_DET set;
  P2DAMP, AITKEN, AA, probe, CHI_EXT-present, JMP_SIB-present NOT set.
- ray 81: NI=5 (tangent depth 5); ray 82: NI=4 (tangent depth 4).
  KN_INNER_METH=DIFFUSION.

## Version gap (measured)
- Live exe d50dcacd (built Nov 20 ~17:26): KNRDUMP=2 strings -> writes v1.
- Pinned pre-flag exe 93364cac / post-flag 58a699d0: KNRDUMP=0 strings.
- Current working source writes v2 (magic 20261123, +KN_CHI_EXT/KN_ETA_EXT
  +JMP_SIB records). Replay driver drv_knop_replay.f90 requires 20261123.
- v1 deck lacks KN_CHI_EXT/KN_ETA_EXT which BUILD_FREQ's on-node branch
  reads (fg_j_cmf_nonmono_v1.f ~13146-13151): replay of v1 decks
  segfaults inside KNOP_RAY_SOLVE_V1 for rays with on-node cells
  (job 13444839). v1 decks are not replayable end-to-end.
- knray_replay relink blocked by incomplete lib/libsubs.a (117 members;
  unix objects incl. term_io.o missing after the Nov 20 18:15 rebuild).
  prod_obj cache holds the needed objects for a manual link.

## Path to a replayable deck
1. Rebuild cmfgen_dev.exe from current source (v2 writer; also restores
   the archive member set). Safe for pending jobs: B1/B2/diagcap use
   frozen trees, observers use cmf_flux.exe, A1/A2 already launched.
2. Armed interactive capture (clone of diagdump_i, 1-2 refreshes is
   enough; deck ~0.4 GB/refresh) for rays 81/82 (+second arm for a
   converged grazing control such as ray 78).
3. Relink knray_replay (v2) and run base/verify/loop/om modes.
