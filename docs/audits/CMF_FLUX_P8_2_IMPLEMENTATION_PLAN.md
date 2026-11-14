# P8.2 ES redistribution convergence diagnostic - implementation plan

Companion to `CMF_FLUX_P8_2_ES_CONVERGENCE_AUDIT.md` (which holds the
source call-chain audit and the experiment results). This file records
what was changed, how the diagnostic is built, and how to run it.

Scope: diagnostic only. No production array, operator, iteration
order, observer-frame transformation, or default control value was
changed. `ES_RESID_DIAG=F` (default) leaves runs bit-identical.

## Files changed

| File | Change |
|---|---|
| `obs/es_resid_diag.f` | NEW. Module `ES_RESID_DIAG_MOD`: `ES_RESID_SEED/UPDATE/RENDER/FINAL` wrappers plus the pure `ES_RESID_COMPARE` core and an `ES_J_CONV` reader using the production direct-access record layout. |
| `obs/cmf_flux_cntrl_var_mod.f` | Three module variables: `DO_ES_RESID_DIAG=.FALSE.`, `ES_DIAG_AFLOOR=1.0D-12`, `ES_DIAG_DEPTHS=.TRUE.` |
| `obs/rd_cmf_flux_controls.f` | Optional-key reads via `RD_STORE_LOG`/`RD_STORE_DBLE` with `L_FALSE` (absent key keeps the default). |
| `obs/cmf_flux_sub_v5.f` | `USE ES_RESID_DIAG_MOD`; three gated call sites: `ES_RESID_SEED` (pre-sweep, `SOB_PHASE==1`, incoherent only), `ES_RESID_UPDATE` (immediately after the end-of-sweep `ES_J_CONV` write+close), `ES_RESID_RENDER`+`ES_RESID_FINAL` (after the final in-memory convolve, before the `ESEC*RJ_CMF_ST` fold). No production line moved or duplicated. |
| `obs/Makefile` | `es_resid_diag.o` added to `OBJ1` and to the `cmf_flux_sub_v5.o` dependency line. |
| `unit_tests/f90/drv_es_resid.f90` | NEW. Manufactured driver: `compare A..J` synthetic pairs, `seq`/`seqstale`/`seqcold` file-level chains, `omp` determinism mode. |
| `unit_tests/python/test_es_resid.py` | NEW. 20-check battery covering spec section 11 cases 1-14. |
| `scripts/p8_es_sweep_compare.py` | NEW. Sweep post-processor: parses `ES_RESID_DIAG` + `obs_fin_10` grids, emits `p82_residual_history.csv`, `p82_final_audit.csv`, `p82_spectrum_metrics.csv`, `p82_line_metrics.csv`, and three PNGs. |
| `docs/audits/CMF_FLUX_P8_2_ES_CONVERGENCE_AUDIT.md` | Source audit, sequences, design, results. |
| `docs/TODO.md`, `README.md` | Status entries. |

## Diagnostic design

- Generation counting: a pre-existing `ES_J_CONV` at seed time is
  generation 0; each completed end-of-sweep `COMP_J_CONV_V2` write is
  the next generation. The final render consumes the in-memory
  `RJ_CMF_ST` (verified against the latest stored file generation by
  `ES_RESID_RENDER`).
- Residuals: `ES_RESID_COMPARE` implements the spec's `S_K`, `D_K`,
  `R_K`, `R_GLOBAL` exactly, with explicit zero-field rules
  (`S=0,D=0 -> R=0`; `S=0,D!=0 -> R=1`) and the active-frequency mask
  `max(|NEW|,|OLD|) > AFLOOR * S_K` (`AFLOOR=1D-12` default). It is a
  pure routine - no file I/O, no module state - so the manufactured
  driver exercises production mathematics directly.
- IEEE safety: sequential `IEEE_IS_FINITE` checks on both fields
  before any ordered comparison or division (Fortran does not
  short-circuit `.AND.`). Invalid nodes mark the generation
  `INVALID` and are located (K,M,which); they are never reported as
  converged. Negatives are counted, not clamped.
- `ES_J_CONV` read-back: `READ_DIRECT_INFO_V3` supplies `ND`, `RECL`,
  `INIT_REC` from the sidecar `INFO` file, so the reader matches the
  production layout written by `OUT_RV_TO_EDDFACTOR`.
- Phase coverage: `SOB_PHASE==1` only; phase 2 is coherent and
  produces no file generations. Under `ACCURATE` the render grid can
  differ from the file grid - the elementwise render check is then
  skipped and reported as such.
- Thread safety: diagnostic code runs only in sequential post-sweep
  regions; no OMP directives. The `omp` test case is byte-identical
  under 1 vs 4 threads.

## Test design

`test_es_resid.py` builds `drv_es_resid` against the production
objects (`prod_obj` on `obs/es_resid_diag.o`; nothing is copied) and
checks:

- compare cases A-J: identical (R=0), 1% uniform (R~0.01), single-node
  localization, zero/zero (R=0), zero-to-nonzero (R=1), masked weak
  node, NaN-in-new, +Inf-in-old, negative counting, contractive series
  q=0.2/0.8/0.98 (residual ratios recover q), oscillatory series
  (R stays O(1), not flagged converged).
- `seq`: seed file -> two generations -> render check on the latest
  field (must match, `FINAL_RENDER_INPUT_GENERATION=2`).
- `seqstale`: render check against a deliberately stale field (must
  fail: `FIELD_MATCH=F`) - the check distinguishes generations.
- `seqcold`: no file at seed -> `SEED GEN=NONE`, sweep 1 reports
  `FIRST_GENERATION` with no invented comparison.
- `omp`: byte-identical output at 1 vs 4 threads.

Result at implementation time: 20/20 PASS.

## Build and run

```
cd obs && make -j4            # builds es_resid_diag.o + cmf_flux.exe
cd unit_tests && python3 python/test_es_resid.py
```

In `CMF_FLUX_PARAM_INIT` (optional keys; defaults shown):

```
F [ES_RESID_DIAG]        ! master switch
1.0D-12 [ES_DIAG_AFLOOR] ! active-node mask floor
T [ES_DIAG_DEPTHS]       ! per-depth rows
```

With `ES_RESID_DIAG=T` the run writes `ES_RESID_DIAG` in the obs
workdir; `batobs.sh` already removes `ES_J_CONV` at teardown but
leaves `ES_RESID_DIAG`.

Sweep experiment layout used here:

```
<model dir>/obs_p82_<sob|knop>_n<N>/   sibling clones of obs/
  CMF_FLUX_PARAM_INIT                  (NUM_ES=N, ES_RESID_DIAG=T,
                                        KNOP: DO_NONMONO_KNOP=T)
```

`scripts/p8_es_sweep_compare.py --root <model dir> --out <outdir>`
collects residuals, spectra, final-generation audit and runtimes.

## Deviations from the prompt

- Key names: implemented as `ES_DIAG_AFLOOR`/`ES_DIAG_DEPTHS`
  (module-variable style consistent with `ES_RESID_DIAG`), not the
  prompt's suggested `ES_DIAG_ACTIVE_FLOOR`/`ES_DIAG_DEPTH_ROWS`.
  Semantics identical.
- `TOTAL_MASKED` is not on the `GLOBAL` line; masked counts are
  carried per depth (`NMASK`) and summed in the analysis script.
- Render-consumption test: the strongest production-linked check is
  the in-module `ES_RESID_RENDER` elementwise comparison of
  `RJ_CMF_ST` against the stored latest file generation, exercised by
  `seq`/`seqstale`. A full observer-frame render inside the unit test
  is not practical (needs the whole model state); the limitation is
  recorded in the audit.
- KNOP sweep uses NUM_ES=2..12 (the production floor of 2 is enforced
  by `rd_cmf_flux_controls.f`), SOB uses 1..12.
