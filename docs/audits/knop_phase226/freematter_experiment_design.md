# Free-matter population/spectrum experiment design (NOT submitted)

## Question

Do the 81P/82P non-contracting ray residuals (bright-cell per-pass
defects ~3% local, integrated ~1e-2 R_INT) measurably change the
population solution or the observer spectrum?

## Design

Two production runs from the SAME checkpoint, identical inputs:

- arm A: KN_MAT_FREEZE=T (baseline - matter frozen, as in the
  diagdump config). Radiation field evolves, populations cannot
  respond. Serves as the transfer-only reference.
- arm B: KN_MAT_FREEZE=F (free matter). The ray fields feed
  J at depths 1-5 -> radiative rates -> population equations.
  All certification (KN_CS_*, FG_RF_ALL_CONV, LAM_FIELD_NOT_READY)
  stays on; no ray is artificially accepted.

Identical: executable (v2 build or later validated), VADAT,
SN_HYDRO_DATA, atomic data, restart checkpoint, OMP/thread count,
job resources. Both arms on the same v2-writer build so decks are
available if a follow-up replay is needed.

## Comparisons

- population: max |db| per level per depth vs iteration;
  final CONV_EXIT_MET, FG_RF_ALL_CONV status, census counts
- rates: photoionization/bb radiative rates at depths 1-5 in the
  ~999-18395A defect regions vs baseline
- spectrum: identical observer formal-solve inputs, identical
  pops-or-frozen-pops; max fractional flux difference, UV/EUV
  bands around 999A and 1033A (O VI), 18395A region; state the
  denominator wherever reference flux ~ 0.

## Provenance required

- exe sha256 for each arm, source commit + dirty patch,
  VADAT hash, checkpoint sha256, job ids, node/compiler/env.
- The experiment may only run after the replay-fidelity gate has
  passed (done: control ray 3.5e-8) and CPU/QOS headroom exists.
  Existing pending jobs must not be displaced.

## Constraints

- No acceptance-policy change, no physical-significance mask,
  no forced ray acceptance. If arm B never reaches
  FG_RF_ALL_CONV, that is a result, not a failure.
- Observer solver unchanged.
