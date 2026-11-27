# WP3 bit-identity + WP6 T1 census render audit

Dated 27-Nov-2026.  Rendered the converged `m15_2p82d_VAL_sob`
model (it414, 24-Nov) on the committed-state worktree exe
(`scratch/cmfgen_sobcheck_wt` @ `f69772a7`, ifort -O3,
`cmf_flux.exe` md5 `a9228c50b5aaad31b31b2d98877bdfc7`) in a
scratch copy (`scratch/render_valsob_*`).  Baselines are the
24-Nov OBS job products.

## WP3: unarmed render must be bit-identical (verdict: PASS)

| product | baseline md5 | new md5 | verdict |
|---|---|---|---|
| `obs/obs_fin_10` | a60eddc0a4d10419683c3a2cc5db26b3 | a60eddc0a4d10419683c3a2cc5db26b3 | IDENTICAL |
| `obs/OBSFLUX` | 531ebe61083bdf1660a3b500f9ccc1d6 | 531ebe61083bdf1660a3b500f9ccc1d6 | IDENTICAL |
| `obs/dFR_DATA` | - | - | IDENTICAL |
| `obs/EDDFACTOR_INFO` | - | - | IDENTICAL |
| `obs/OUT_FLUX` | 58d3f4de... | differs | 2 lines only: added `[WARN][SET_TWO_PHOT]` / `[INFO][TWO_PHOT]` log tags (print-tagging commits, not SOB WP code) |
| `obs/EDDFACTOR` | efec4214... | differs | last-ulp: maxabs 6.8e-16, median rel 4.4e-16 - OpenMP reduction-order noise |
| `obs/IP_DATA` | 6b12d153... | differs | last-ulp: maxabs 1.3e-17, median rel 4.8e-16 - same |
| `obs/EWDATA` | dangling symlink | dangling symlink | unchanged (never written under DO_SOB_LINES=F) |

The spectral products (`obs_fin_10`, `OBSFLUX`) are bit-identical:
the census instrumentation + EXPONX bound do not alter the
default render path.  EDDFACTOR/IP_DATA differences are
thread-order nondeterminism: a second render of the same exe on
the same inputs (job 13775543) reproduces obs_fin_10/OBSFLUX/
OUT_FLUX exactly while EDDFACTOR/IP_DATA differ at the same
scale (maxabs 6.7e-16 / 1.3e-17) - run-to-run FP noise, not a
code-path change.  VERDICT: PASS.

## WP6: armed census render (first attempt)

`CMF_SOB_T1_CENSUS=1` armed render on the same staged copy
produced bit-identical `obs_fin_10`/`OBSFLUX` (armed
instrumentation does not perturb the solve) but ZERO census
lines.  Root cause found: the staged `CMF_FLUX_PARAM_INIT`
carries `DO_SOB_LINES=F` and `DO_NONMONO_SOB=F`, so
`CMF_FLUX_SUB_V5` returns at the `DO_SOBOLEV_LINES` guard (line
2875) before the EW line loop, the `SOB_T1_CEN_ACC` hooks, and
the `SOB_T1_CENSUS_REPORT` call (line 3238).  The EW estimator
- and therefore the T1/GAM census - only runs under
`DO_SOB_LINES=T`; the `SIGMA_TRUE` branch needs
`DO_NONMONO_SOB=T` in the OBS param.

Fix: resubmitted (job 13775472) with `DO_SOB_LINES=T` +
`DO_NONMONO_SOB=T` in the staged INIT.  The env arming mechanism
itself was verified by a standalone driver against the worktree
`libsubs` (CEN_ARMED=T under `CMF_SOB_T1_CENSUS=1`, report
prints).

## WP6: armed census output (job 13775749, VERIFIED)

Staged copy rendered with `DO_SOB_LINES=T`, `DO_NONMONO_SOB=T`,
`GLOBAL_LINE=SOB`, `CMF_SOB_T1_CENSUS=1`.  The EW estimator ran
the unclamped `SIGMA_TRUE` branch of SOBEW_GRAD_V2.  Census block
(`obs/OUT_FLUX` lines 240-247 of the staged render):

```
SOB_T1_CENSUS: samples =    420926590
SOB_T1_CENSUS: min T1 =    7.393938E-03  at ray     36 depth      1
SOB_T1_CENSUS: min GAM =   -2.738882E+00  at ray     52 depth     15
SOB_T1_CENSUS: max GAM =    7.205230E+07  at ray     47 depth     25
SOB_T1_CENSUS: sigma range =   -1.556105E+02 ..   -1.036970E+00
SOB_T1_CENSUS: undefined sigma =        0
SOB_T1_CENSUS: GAM NaN =        0  Inf =        0
SOB_T1_CENSUS: EXPONX NaN =        0  Inf =        0
```

Read: 4.21e8 (ray,depth,freq,line) samples; the 1e-4 T1 floor is
never engaged (min T1 = 7.4e-3); sigma reaches -155.6 in the
compression shell (deeper than the -82 first estimate); negative
GAMs exist down to -2.74 but EXPONX stays finite (the X < -700
bound is never approached); the largest GAM (7.2e7) sits in the
1/X branch with no overflow.  No NaN/Inf anywhere.

The EXPONX-guard counter line is absent because the worktree exe
predates the `sobew_grad_v2.f` fold; from the zero overflow
counts the bound did not engage.  A worktree rebuild + rerun
will add the printed counters for completeness.

## WP6: guard-counter verification (job 13783512, VERIFIED)

Worktree rebuilt at `118f6779` (`cmf_flux.exe` md5
`9c403faa5fb7412fda09ae2e9d574764`) and the same staged render
rerun.  All census lines reproduce the 13775749 numbers exactly;
the new line reports

```
SOB_T1_CENSUS: EXPONX guard =        0  d_EXPONX guard =        0
```

The X < -700 bound never engaged on m15, matching the WP3
analysis (the dangerous sigma_true window (-1.0005, -1) is never
entered).  OUT_FLUX vs the 13775749 copy differs only in two
last-ulp T1 warning digits (OpenMP order noise, ~1e-14 rel) plus
the added counter line.  The `obs_fin_10` in this deck is the
DO_SOB_LINES=T reduced EW product (815 frequencies), not the
production spectrum - expected for this staged INIT.

Regression artifact: the staged render dir
`scratch/render_valsob_census/` (OUT_FLUX lines 240-248; pre-fold
copy saved as `OUT_FLUX_13775749`) is the armed-census
reference; reproduce with
`scratch/render_valsob_census/run_render.sh`.
