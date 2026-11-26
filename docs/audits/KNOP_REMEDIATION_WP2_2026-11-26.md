# KNOP remediation WP-2: freshness contract and order modes

26-Nov-2026. Second work package of the Knop non-monotonic solver
remediation (spec: `checklists/knop_solver.md`). WP-0 manufactured
the stale-donor fixture and flipped the ordered march off by
default; WP-1 landed the solver defaults and the tolerance-ladder
invariants. WP-2 makes the freshness claim measurable and adds the
per-ray order-mode selector the spec's fallback plan needs.

## What WP-0 measured and what the comment claimed

The `SWEEP_FLD_ORD` header claimed the only lagged reads were the
DS-corner source preview and "compressing interpolant donors JPR on
the unmarched side of the front" -- and then stated "the main |Q|
wavelength link never waits a pass". The second clause is false:
the off-row JPR share of the |Q| donor IS a full pass stale, and on
the WINX=0.95 fixture that share is 0.943 of the compressing donor
weight on the I_M pass. The comment is corrected at source.

## Changes

### `KN_RAY_ORDER_MODE` -- order-mode selector (VADAT, default -1)

-1 is a sentinel resolved from `KN_RAY_ORDERED` (T -> 0, F -> 1),
so drivers that poke `FG_ORD_SWEEP` directly keep their behaviour.
Concrete modes:

- **0** ordered march (`SWEEP_FLD_ORD`, the legacy default pre-WP-1).
- **1** frequency-major alternating (`SWEEP_FLD`, the WP-1 default).
- **2** frequency-major, first pass down the dominant coupling
  sign. A uniform-compressing ray (some GAM<0, none >0) starts
  descending, so the L+NFM donor column is fresh on pass 1 where
  mode 1 spends a fully stale ascending pass. Every other ray
  keeps the mode-1 start, so mode 2 reduces to mode 1 exactly on
  mixed and expanding-only rays.
- **3** SCC/topological ordering. Reserved: the spec only lands it
  if mode 2 underperforms on mixed rays. Stops in
  `CHECK_PARAM_CONSISTENCY` until then; out-of-range values stop
  too. An explicit mode and an explicit `KN_RAY_ORDERED` that
  disagree resolve to the mode with a warning.

### `KN_FRESH_DIAG` -- stale-donor census (VADAT, default F, diagnostic)

Per pass, the share of |Q| donor weight the march read from a
not-yet-updated cell:

- `SWEEP_FLD` (alternating): a row's whole donor weight is stale
  when the pass direction leaves the coupled column pending --
  compressing rows on ascending passes, expanding rows on
  descending passes.
- `SWEEP_FLD_ORD`: the JPR share on the unmarched side of the
  front (`KN_FRESH_ACC` mirrors `KN_ORD_SPV_IMPL`'s donor branch;
  the ABI-frozen worker does not export NNF2, so it recomputes
  it). Expanding own-row donors and self-mapped JPR=own-row shares
  are always fresh.

Stored per (pass, dir, ray) in `KN_FRESH_FR` (dir 1=I_M, 2=I_P;
-1 marks a pass that never ran), printed one line per pass on
LUER. `KN_ORD_SPV_IMPL`'s ABI is untouched; the census adds no
numerics, only counters.

## Fixture measurements (test_knop_ordering)

| case | march | fresh census | verdict |
|---|---|---|---|
| ord | mode 0 | I_M stale 0.943/pass, I_P 0 | stalls, NXIT=6 |
| uno | mode 1 | passes: 1.0, 0.0, 1.0 | clean, 3 passes |
| m2u | mode 2 | passes: 0.0, 1.0 | clean, 2 passes |
| m2m | mode 2 | identical to unom | bitwise identical |

The ordered march's stale share is the WINX share on the I_M
direction exactly -- the lag lives entirely where the front
ascends into unmarched deeper rows. Mode 2 recovers the fresh
first pass mode 1 wastes, and on this uniform-compressing fixture
that halves the pass count (2 vs 3 at exact zero residual).

## Verification

- `test_knop_ordering` 50/0 (was 39/0): ord fresh-m stale ~0.943 /
  fresh-p ~0; uno [1.0, 0.0, 1.0]; m2u [0.0, 1.0] clean; m2m
  bitwise = unom at all three expanding fractions.
- `test_vadat_keys` 120/0 (was 114/0): default resolve to 1,
  mode=2 lands, mode=3 and out-of-range stop, conflict warns and
  overrides, KN_FRESH_DIAG lands.
- `test_knop_ray` 50/1 (T18 adaptive-omega diag remains the only
  failure, pre-existing), `test_knop_spv` 93/0 (ABI unchanged),
  `test_iter_conv_gate` 18/0, `test_conv_exit` 53/0, `test_fg_gs`
  36/0, `test_knop_oneway` 8/0.

## Verdict

VERIFIED. The freshness contract is now measured, not claimed: the
ordered march's |Q| lag is the WINX share on the ascending front,
the alternating sweep pays one stale direction per pair, and mode
2 eliminates the wasted opening pass on uniform-compressing rays.
Next: WP-3 (store freshness) measures whether the same contract
holds on the production route.
