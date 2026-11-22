# Phase 2.27G - E1/E2 runtime results (21-Nov-2026)

First runtime experiments of the phase series. Executed with the
offline `knray_replay` driver (build_v6, ifort -O3 production flags,
OMP_NUM_THREADS=1) against the archived v2 decks
`KNRDUMP_L081/082/084_r0000.bin`. Not a production run; no jobs used.
Working dir: cmfgen_validation/phase227b_20261121_0604/e1e2/.

## E0 - baseline reproduction (base mode = Mode A)

RUNTIME CONFIRMED: replay reproduces the historical numbers exactly.

| Ray | replay_vs_prod maxdiff | NPO replay/prod | IXIT replay/prod |
|-----|----------------------|-----------------|------------------|
| 81  | 2.656368E-04         | 23 / 12         | 3 / 3            |
| 82  | 2.442209E-04         | 27 / 24         | 3 / 3            |
| 84  | 3.487919E-08         | 3 / 3           | 1 / 1            |

Discrete check (no new instrumentation): replay runs MORE passes than
production on 81/82 - a step-type trajectory divergence, not smooth
amplification of a tiny seed. Ray 84 matches pass-for-pass.

## E1 - replay M vs deck (production) M

RUNTIME CONFIRMED: the recomputed M field differs materially.
RELERR=|a-b|/MAX(|a|,|b|,1e-300).

| Ray | FLDM med | FLDM p99 | IM_NI med | IM_NI p99 | bitexact |
|-----|----------|----------|-----------|-----------|----------|
| 81  | 5.7e-3   | 3.7e-2   | 3.1e-2    | 9.8e-2    | 1 col    |
| 82  | 2.7e-3   | 4.3e-2   | 3.6e-2    | 1.0e-1    | 1 col    |
| 84  | 4.1e-6   | 7.7e-6   | 4.1e-6    | 8.4e-6    | 1 col    |

Exactly one column (NCF=56377 cells) is bitwise identical on every
ray - likely a pass-invariant column; the rest of the field diverges.
Even the control ray shows a real ~4e-6 M mismatch (not FP noise).

## E2 - imposed-M experiment (deck IM_NI forced into the P march)

New solver hook KN_IMPOSE_IMNI/KN_IMPOSE_ROW + driver mode 'impose'
(diagnostic-only, off in production). Engagement verified: e1-vs-e2
FLDP differ everywhere (0 bitexact, maxdiff 1.7e-4) and pass counts
shift (81: 23->20, 82: 27->26).

| Ray | base maxdiff | impose maxdiff | collapse? |
|-----|--------------|----------------|-----------|
| 81  | 2.656e-4     | 2.625e-4       | NO        |
| 82  | 2.442e-4     | 2.438e-4       | NO        |
| 84  | 3.488e-8     | 1.272e-8       | partial   |

## Attribution

RUNTIME CONFIRMED: the M reconstruction gap is real but NOT causal
for the P-end discrepancy. A bit-exact production boundary changes
the answer by ~1%. The 2.4e-4 originates downstream of IM_NI:
P-trajectory divergence (pass count, window/stall decisions) or
state not captured in the deck.

Supporting evidence that the stall-regime end state has ~1e-4-scale
sensitivity: imposing a 3e-2 boundary change moved the end state
only 1.7e-4 (effective gain ~6e-3); a 2.4e-4 gap is consistent with
decision flips (pass counts, hot-window timing) rather than
continuous amplification of a boundary error.

Ranking after E1/E2:
1. P-side decision/state divergence (pass count differs even with
   identical boundary - first candidate for per-pass telemetry).
2. Unrecorded input state the deck lacks (anything production had
   at P start that is not in the v2 record).
3. FP build differences (bounded by ray 84: cannot produce 1e4x).
4. M reconstruction - REAL but exonerated as the P-gap cause.

## Next capture requirements (unchanged need, now proven)

- Per-pass rows for armed rays: RES, applied omega, HOTMODE,
  W_LO/W_HI, NSTALL, IXIT, pass kind (add window bounds - currently
  not recorded anywhere; replay recomputes them from its own
  residuals, a discontinuous function).
- Production-vs-production repeat in the same capture: if two
  production refreshes differ ~1e-4 on rays 81/82, the 2.4e-4 lies
  inside the production ensemble width and no single-ray replay can
  remove it. Record IM_NI on both to measure run-to-run M
  variability directly.
- Hash-verify deck inputs at capture (coefficient packs, stores,
  seed, IM_NI).

Phase 2.27 remains PARTIAL.
