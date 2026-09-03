# Unresolved Assumptions Requiring Confirmation

**Status legend:** RESOLVED (verified against `cmfgen_source/` this session) /
OPEN (needs a real CMFGEN run or further source trace).

## RESOLVED this session (citations added to the report)

### #1 — Does LTEPOP_WLD_V2 use the species `_IN` T column or T_IN?  → RESOLVED

**Answer:** It uses the **model T from T_IN**, not the species `_IN` T column.

**Evidence:** `set_new_model_estimates.f` sets T from T_IN first
(`REGRID_T_ED_V3(...,'T_IN')` line 163; `REGRID_T_ED_V2(...,'T_IN')` line 206),
then calls `LTEPOP_WLD_V2(..., ED, T, ND)` at line 367-370 with that model T.
`REGRID_LOG_DC_V1` reads the species `_IN` T only as `OLD_T`, the interpolation
reference coordinate (lines 149-166), never as the physical LTE temperature.

**Consequence for the proxy:** writing T_local to the species `_IN` T column is
safe and conservative (it matches the physical T, so DC interpolation is
self-consistent). Writing T_ion would change only the interpolation coordinate,
not the physical LTE T — so B5's caution is upheld by defaulting to T_local.

### #2 — Does the species `_IN` ED column override GAMMAS_IN?  → RESOLVED

**Answer:** No. ED is set from GAMMAS_IN / T_IN, not the species `_IN` ED.

**Evidence:** `GETELEC_V2(...,'GAMMAS_IN')` (lines 189, 236) and
`REGRID_T_ED(...,'T_IN')` (line 163/206) set the model ED. `REGRID_LOG_DC_V1`
reads the species `_IN` ED only as `OLD_ED` (interpolation reference, line 150).
`LTEPOP_WLD_V2` at line 370 uses the model ED.

**Consequence:** v2 writing the GAMMAS_IN ED into the species `_IN` ED column is
consistent and safe.

### #3 — F_TO_S: does the proxy write superlevel or full-level DCs?  → RESOLVED (NF safe)

**Answer:** v2 writes **NF (full-level)** DCs by default and pre-applies F_TO_S
replication (lowest full level's DC → all full levels in the superlevel),
mirroring `regrid_log_dc_v1.f:609-619`. This is self-consistent regardless of
DC_METH. `--write-levels NS` is available for the legacy RTX path but not
recommended unless NS==NF.

**Evidence:** `regrid_log_dc_v1.f:609-619` (F_TO_S pass); round-trip test
`test_f_to_s_replication` confirms the replication logic.

## OPEN — require a real CMFGEN run or further source trace

### #4 — What happens when NOLD < N (proxy writes fewer levels than the model)?

When F_OSCDAT has fewer levels than MODEL_SPEC's NF, v2 writes NOLD = min(NF,
n_levels_F_OSCDAT) < N. `REGRID_LOG_DC_V1` lines 557-602 extrapolate the missing
levels using an excitation temperature. This could introduce error in the
extrapolated levels' DCs.

**Needed:** a real CMFGEN run on a species where NOLD < N, comparing the
first-iteration population corrections against the LTE baseline. Alternatively,
trace `regrid_log_dc_v1.f:557-602` to confirm the extrapolation is stable for
the proxy's DC=1 (LTE) case (where extrapolation of DC=1 is trivially exact).

**Mitigation already in place:** v2 reports `nlev` per species in
`proxy_audit.json`; any species with nlev < NF is visible. For LTE mode (DC=1)
the extrapolation is exact regardless of NOLD, so this only matters for GREY/
CONTINUATION modes.

### #5 — Demonstrated improvement (lower first-iteration CMFGEN corrections)

The acceptance criterion "any claimed improvement is demonstrated by lower
first-iteration CMFGEN population corrections or faster convergence on a real
test model" has NOT been demonstrated yet. The round-trip test confirms
**internal consistency** (reconstruction == target), not that the proxy reduces
CMFGEN's first-iteration corrections relative to the old proxy.

**Needed:** run CMFGEN (via `run_model_grid.py`) on the v2 LTE `_IN` files and
on the v1 `_IN` files for the same model, and compare:
- maximum / 95th / 99th percentile population correction factors at iteration 1;
- number of iterations to convergence;
- charge and abundance residuals after iteration 1.

The CORRECTION_SUM / CORRECTION_LINK files from a CMFGEN run provide the
per-iteration correction diagnostics. This is the only remaining acceptance
criterion that requires compute time rather than source inspection.

### #6 — Occupation probability W (level dissolution)

`LTEPOP_WLD_V2` calls `OCCUPATION_PROB` (line 59) which computes W_C2 (level
dissolution / occupation probabilities) from ED, edges, and Z. v2 defaults
W_occ = 1.0 (documented as audit assumption 5). This means the proxy's LTE
populations omit level dissolution, so CMFGEN's first-iteration LTE (which
includes OCCUPATION_PROB) will differ from the proxy's by the W factor in
high-density layers.

**Needed:** confirm whether OCCUPATION_PROB significantly changes the LTE
populations at the densities of the target epoch (0.5 d has high density). If
it does, v2 could emulate OCCUPATION_PROB; if not, W=1 is acceptable. This is a
**physical-accuracy** gap, not a correctness bug — the round-trip is still exact
because both sides use W=1 in the emulation.

### #7 — HDKT literal vs computed

`cmfgen.f:117-122` computes HDKT from physical constants and notes the old
literal `4.7994145`. v2 uses the literal `4.7994145` to match the routines that
read our files. If the CMFGEN build in use computes HDKT differently (e.g.,
updated constants), there would be a tiny systematic offset in the LTE
populations.

**Needed:** confirm which HDKT path the production CMFGEN binary uses
(`cmfgen_base/batch.sh` build). Low impact (<1e-4 relative) but listed for
completeness.

## Summary

| # | Assumption | Status |
|---|---|---|
| 1 | Species `_IN` T column vs T_IN for LTEPOP | RESOLVED (T_IN) |
| 2 | Species `_IN` ED column vs GAMMAS_IN | RESOLVED (GAMMAS_IN) |
| 3 | F_TO_S: NF vs NS DCs | RESOLVED (NF safe) |
| 4 | NOLD < N extrapolation | OPEN (low risk for LTE) |
| 5 | Demonstrated improvement vs old proxy | OPEN (needs CMFGEN run) |
| 6 | Occupation probability W=1 omission | OPEN (physical-accuracy gap) |
| 7 | HDKT literal vs computed | OPEN (low impact) |

The two BLOCKING-priority assumptions (#1, #2) that affected the B5 temperature
decision are now RESOLVED by source inspection. The remaining OPEN items are
either low-risk for the default LTE mode (#4, #7) or require compute time / a
physical-accuracy study (#5, #6) rather than further source reading.
