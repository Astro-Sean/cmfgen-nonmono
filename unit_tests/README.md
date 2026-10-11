# unit_tests -- production-file regression tests for the modified CMFGEN/JFNK tree

Every test is a Python script that compiles the **exact production
Fortran files** from `cmfgen_jfnk/` (and, where the signature allows,
the pristine reference in `cmfgen_reference/cur_cmf/`) together with a
small CLI driver in `f90/`, runs the binary, and asserts on its
output.  No Python reimplementation of production code.

Run everything:

```bash
python3 unit_tests/run_all_tests.py        # exit status = # failed files
python3 unit_tests/python/test_<name>.py   # single test
```

Binaries land in `bin/` (rebuilt only when sources change), fixtures
and scratch in `work/`, diagnostic plots in `plots/`.  Expected aborts
(`STOP`) are run with `check=False` and asserted via return code +
error text.

## Test inventory

| test file | production file(s) under test | what is asserted |
|---|---|---|
| `test_sigma_true.py` | `new_main/subs/recompute_sigma_true.f` | SIGMA_TRUE_CORE flag gating, unclamped `R/V dV/dR-1`, boundary cells, stagnation and `ND<3` edge cases, new-vs-reference identity when flag off |
| `test_sob_tau_floor.py` | `subs/sobjbar_sim.f`, `sobew.f`, `sobew_grad_v2.f`, `sobew_grad_v3.f` | sigma>=0 and sigma=-0.5 bit-identical to reference; sigma<-1 hits the Dessart+2015 tau floor; analytic BETA check against SOBJBAR_SIM; finite output |
| `test_interp_guards.py` | `subs/lin_interp.f`, `mon_interp.f`, `mon_interp_fast.f`, `mon_interp_fast_v2.f`, `mon_int_funs_v2.f` | valid-path bit-identity vs reference for all five routines; aborts on non-monotonic/repeated/reversed grids, out-of-range targets, `ND`/`NIN` too small |
| `test_rd_rv_file.py` | `newsubs/rd_rv_file_v2.f` | `NONMONO_V` VINF check against `MAXVAL(ABS(V))` (not `V(1)`); interior-peak profiles; 10% tolerance edges; `VINF<0.1` bypass; all-zero V, non-decreasing R and R-scaling rejections; old-behaviour adapter comparison |
| `test_sob_rh_couple.py` | `subs/sob_rh_couple.f` (new file) | zero pairs under homology; pair count/symmetry vs independent geometry; `mu^2` identity; DRC net-rate term vs recomputed value; suppression of `V<=0`/negative-opacity rows; `BB_COR=0`/`ETAL=0` safety |
| `test_ng_abstain.py` | `new_main/subs/do_ng_accel_band_v4.f` | `NG_ABSTAIN` per-component hold on zero-diff and non-positive extrapolation; siblings still accelerate; upstream whole-depth veto preserved when flag off; history restoration; vacuous guard |
| `test_rv_sniin.py` | `new_main/subs/rv_sn_model_sniin.f` | RDINR beta-law path bit-identical new-vs-ref; corrected SIGMA centred-difference denominator verified exactly against the analytic formula (reference matches the old formula, `|delta|~12`); uninitialized-V coordinate into LIN_INTERP aborts via the strict guards |
| `test_rv_snhyd.py` | `new_main/subs/rv_sn_model_snhyd.f` (new file, VEL_LAW=13) | all new guards abort (ND, missing/malformed SN_HYDRO_DATA, non-decreasing input R, `V<=0`, non-positive expansion factor, non-decreasing expanded R, TAU_OPT=1 with no positive Kappa); endpoints pinned to expanded hydro profile; monotone-cubic V on nodes; homologous V -> SIGMA=0; beta-law -> SIGMA~BETA-1; velocity reversal keeps `SIGMA<-1` under both TAU_OPT |
| `test_edep_shock.py` | `new_main/subs/get_edep_shock_power.f` | flags-off path bit-identical to reference including non-monotonic V (legacy tolerance); `DO_NONMONO_SOB`/`DO_NONMONO_KNOP`/`both` abort on non-positive or non-decreasing V and accept monotone V; Gaussian deposition profile shape/location; 4*pi normalization vs prescribed power; upstream `clump_fac=0` and `sum(wrk)=0` guards |
| `test_vadat_keys.py` | `new_main/mod_subs/rd_control_variables.f` (+ real `rd_store_log.f`, `control_variable_mod.f`, `jfnk_solver_v1.f90`, `check_param_consistency.f`, `fg_j_cmf_nonmono_v1.f`, `mom_j_cmf_nonmono_v1.f`) | a production VADAT parses; every new non-mono/LAM/BA/JFNK key round-trips into the module store; mutual exclusion; GLOBAL_LINE/USE_LAM_ES/LAM_SET prerequisites abort; DO_JFNK non-mono gate; JFNK_COUPLED KNOP-only; USE_JFNK deprecated alias; JFNK_LAG_FACTORS/JFNK_FRZ_GATE/JFNK_GATE_NREQ round-trip; DO_NONLOC_SOB auto-off; DO_CL/GAMRAY_TRANS legacy-coordinate rejections; IB_METH whitelist + case normalization; VEL_LAW=13 -> SN_HYDRO_MODEL |
| `test_jfnk_solver.py` | `new_main/mod_subs/jfnk_solver_v1.f90` (via `tests/jfnk_core` and `cmfgen_jfnk_unit_test`) | Newton/GMRES/hookstep core T1-T3 on an analytic stiff problem (PTC off, SER anchor, far-start trust region); surrogate-SE systems T0-T14: FD-JVP error curve, demanded-step SVD, PC variants, GMRES, end-to-end Newton, masked-T rows, tauto ladder/FSM.  Known author-flagged diagnostics are XFAIL |

## Supporting files

- `common.py` -- compiler discovery (ifort preferred, gfortran
  fallback), per-test object/link caching, process runner,
  denormal-tolerant comparison, `Report` pass/fail collector.
- `f90/drv_*.f90` -- one CLI driver per test (read fixtures, call the
  production routine, print machine-parseable output).
- `f90/adapt_*.f` -- signature adapters where new and reference
  routines differ (`rd_rv_file_v2` `NONMONO_V`, `do_ng_accel_band_v4`
  `NG_ABSTAIN`).
- `f90/stub_*.f` -- no-op stubs for externals that exist only on
  code paths the driver never reaches (`SCR_READ_V2`, KNOP solver
  internals).
- `f90/drv_edep_ref.f90`, `drv_edep.f90` -- paired drivers for the
  reference (no `DO_NONMONO_*` variables) and new trees.

## Generated plots (`plots/`)

`sigma_true_nonmono.png`, `sob_tau_floor.png`, `interp_guards.png`,
`rd_rv_vinf_profiles.png`, `sob_rh_couple.png`, `ng_abstain.png`,
`rv_sniin.png`, `rv_snhyd.png`, `edep_shock.png`

Plots are diagnostics only -- every assertion is numeric.

## Notes on upstream quirks the tests had to accommodate

- `rv_sn_model_sniin.f` reads `V_HYDRO` as a whole array inside its
  row loop and calls `LIN_INTERP` with the uninitialized output `V` as
  the interpolation coordinate.  The first quirk means a data row must
  carry `NHYDRO+6` values; the second is now caught loudly by the new
  strict `LIN_INTERP` monotonicity guards (verified by
  `test_rv_sniin.py` vfill=Z/I cases).
- The same file's R/V "reversal" loop swaps each pair twice and is
  therefore a no-op; `V_HYDRO` direction expectations in
  `test_rv_sniin.py` account for this.
- Several interpolation routines treat `NQ`/`NV`/`ND` arguments as
  leading-dimension strides, not logical sizes; drivers declare arrays
  with the exact dimensions the routines expect.
- `rd_control_variables.f` rejects duplicate VADAT keys at the
  RD-store level, so test mutations may replace an existing line or
  append a missing optional key, never duplicate one.
