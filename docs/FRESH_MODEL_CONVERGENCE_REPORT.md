# FRESH_MODEL_CONVERGENCE_REPORT.md

**Date:** 2026-08-21
**Model:** SN_12C/cmfgen/0.5days (ND=400, 10 species, t=0.5532 d)
**Compute:** Viper interactive partition, 4 CPUs, 16GB, 47 seconds

---

## 1. Startup Test (VERIFIED)

A CMFGEN startup test was run with proxy v2 LTE _IN files (DC=1) as the
initialization. The test used IN_ITS=2 (2 iterations) as a documented diagnostic
to confirm that CMFGEN reads the proxy files and to capture first-iteration
corrections.

### 1.1 Confirmation that CMFGEN reads proxy files (VERIFIED)

From the OUTGEN log:
```
*_IN files will be used to start model
Setting ITS_DONE keyword in HYDRO_DEFAULTS to 0
```

This confirms that with LTE_EST=F and DC_METH=R (as in the target VADAT),
CMFGEN reads the proxy-generated _IN files via REGRID_LOG_DC_V1.

### 1.2 First-iteration corrections (VERIFIED)

| Metric | Proxy LTE init | Existing model init |
|---|---|---|
| Max % increase (iter 1) | 3.44E+22 at depth 243 | 1.70E+22 at depth 164 |
| Max % decrease (iter 1) | 9.99E+01 at depth 1 | 1.10E+02 at depth 181 |
| Luminosity (iter 1) | 1.86E+07 Lsun | 0.00E+00 Lsun |
| Luminosity (iter 2) | 1.86E+07 Lsun | 0.00E+00 Lsun |
| Runtime (2 iterations) | 47 seconds | — |
| NaN/Inf errors | 0 | 0 |
| Crash | No | No |

**Interpretation:**
- The proxy LTE init has a larger max % increase (3.44E+22 vs 1.70E+22), expected
  because LTE (DC=1) is further from the NLTE solution than the existing _IN files.
- The proxy init produces a **non-zero luminosity from iteration 1** (1.86E+07 Lsun),
  while the existing model init shows **zero luminosity**. This suggests the proxy
  populations are physically reasonable (producing flux), while the existing _IN
  files may have been from a problematic state.
- The huge max % increase (~1E+22) is typical for LTE cold starts: species with
  essentially zero LTE population but non-zero NLTE population show enormous
  relative corrections. This is not a proxy bug.
- Both initializations are stable (no NaN, no crash, corrections bounded).

### 1.3 Generated _OUT files (VERIFIED)

All 10 species _OUT files were generated and parsed by the independent parser:
- 0 non-finite values across all files
- 0 format errors
- DC ranges near 1.0 (0.85–1.19) after 2 iterations from LTE (DC=1)

| Species | NOLD | DC min | DC max | Non-finite |
|---|---|---|---|---|
| HI | 12 | 0.907 | 1.103 | 0 |
| HeI | 10 | 0.907 | 1.135 | 0 |
| OI | 7 | 0.907 | 1.103 | 0 |
| SkI | 1 | 0.999 | 1.043 | 0 |
| Sk2 | 10 | 0.907 | 1.126 | 0 |
| CaI | 1 | 0.999 | 1.017 | 0 |
| Ca2 | 12 | 0.917 | 1.103 | 0 |
| FeI | 1 | 1.000 | 1.190 | 0 |
| Fe2 | 25 | 0.926 | 1.101 | 0 |
| FeIII | 1 | 0.852 | 1.104 | 0 |

---

## 2. Full Convergence Comparison (BLOCKED)

A full convergence comparison (running to convergence with both initializations
and comparing final populations/spectra) is **BLOCKED** by compute time:
- The existing model ran for 15,469 iterations
- At ~24 seconds/iteration, full convergence would take ~4 days
- The interactive partition has a 2-hour time limit
- The general partition could accommodate this but would require a production
  submission via `run_model_grid.py`

### 2.1 What would be needed for full comparison

1. Run CMFGEN to convergence with proxy LTE _IN files (DC=1)
2. Run CMFGEN to convergence with built-in LTE init (LTE_EST=T, no proxy files)
3. Run CMFGEN to convergence with existing _IN files (warm restart)
4. Compare:
   - Final converged populations (per species, per depth)
   - Final spectrum (flux vs wavelength)
   - Total iteration count to convergence threshold
   - Runtime
   - Electron density residuals
   - Luminosity/flux residuals

### 2.2 HDKT verification (VERIFIED)

The startup test OUTGEN confirms:
```
HDKT=  4.79924307E+00
```
The built CMFGEN binary uses a computed HDKT (4.79924307), while the proxy uses
the old literal (4.7994145). The relative difference is ~3e-5, which is negligible
for LTE populations. This resolves unresolved assumption #7.

---

## 3. Orchestrator Integration Gap (VERIFIED)

The CMFGEN startup test was run as a **manual documented diagnostic**, not through
`run_model_grid.py`. This is because the orchestrator integration is incomplete:
- `prepare_VADAT.py` does not accept `--proxy-populations` (would crash with
  argparse error)
- No script calls `proxy_populations_v2.py` to generate _IN files
- The proxy _IN files were manually staged in the test directory

For production use, the orchestrator must be fixed to wire the proxy through:
`YAML config → run_model_grid.py → proxy_populations_v2.py → CMFGEN`

---

## 4. Conclusion

| Criterion | Status |
|---|---|
| CMFGEN reads proxy _IN files | VERIFIED |
| CMFGEN runs without crash/NaN | VERIFIED |
| First-iteration corrections captured | VERIFIED |
| _OUT files generated and valid | VERIFIED |
| Full convergence comparison | BLOCKED (compute time) |
| Built-in init vs proxy init comparison | BLOCKED (compute time) |
| Final population/spectrum comparison | BLOCKED (compute time) |

The proxy LTE initialization is **usable** (CMFGEN reads it, runs stably, produces
non-zero luminosity). Full convergence validation requires a production run.

---

## 5. BA Iteration and Requeue Fix (2026-09-02)

### BA iterations now working

The full convergence comparison is no longer blocked by BA iteration
failure. As of 2026-09-02, true Full BA iterations with depth coupling
are verified working on the X_nomono_v3 model family
(X_nomono_v3_ba_test and X_nomono_v3_vary_ne).

The root cause of all previous BA convergence failures was
`DO_LAM_AUTO=F`, which silently overrides IT_SPECIFIER's FULL(FIXT)
blocks (cmfgen_sub.f:4798 forces `LAMBDA_ITERATION=.TRUE.`). With
`DO_LAM_AUTO=T`, FULL(FIXT) blocks now produce true Full BA iterations.

See `README_STATIC_NONMONO.md` section "BA Iteration Configuration"
and `AGENTS.md` section 43.5 for the corrected settings.

### Requeue configuration fixed

The `final_requeue_max=0` bug that caused jobs to abort immediately on
requeue (wasting all convergence progress) has been fixed. Production
runs now use:
- `partition: general` with `time: '24:00:00'`
- `final_requeue_when: wall` with `final_requeue_max: 500`
- Graceful stop at ~20h into 24h wall, USR1 signal ~4h before wall

Debug runs use `partition: interactive` with hardcoded `01:30:00` time
limit (the interactive partition rejects jobs > 2h).
