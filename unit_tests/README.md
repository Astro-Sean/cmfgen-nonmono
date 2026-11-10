# unit_tests -- production-object regression tests for the modified CMFGEN/JFNK tree

Every test is a Python script that exercises the **exact production
object files** of this tree (the repo root, parent of `unit_tests/`;
`common.py:TREE`): production sources resolve to the `.o` the
production build produced - a sibling `.o` where the build keeps one
(e.g. `obs/`), else the member extracted verbatim from `lib/*.a`
(`common.py:prod_obj`, cached under `bin/prod_obj/`).  Only test-only
code compiles fresh: `f90/` drivers, stubs and signature adapters.
Reference-tree sources (`common.py:REF` =
`../cmfgen_reference/cur_cmf/`) compile fresh with the probed
production flag set (`prod_comp_flags`) so new-vs-reference numerical
comparisons are apples-to-apples; last-bit differences between a
shipped `-O3` object and a separately compiled reference are asserted
through `ref_close` (ULP-level), never silently.

Prerequisites: a completed production build (`make DEBUG=nodebug all`
or `./build.sh`) so `lib/*.a`, sibling `.o`s and `mod/*.mod` exist and
are fresh, plus `ifort` (the production compiler - gfortran cannot
consume the ifort objects/module files).  Stale objects fail loudly:
`prod_obj` compares member/source timestamps from the archive headers
and `require_fresh_archives` guards tests that link whole archives.

Run everything:

```bash
python3 unit_tests/run_all_tests.py        # exit status = # failed files
python3 unit_tests/python/test_<name>.py   # single test
bash unit_tests/build.sh                   # rebuild the se_driver harness
bash tests/f90/build_and_run.sh            # JFNK solver/interface suites
```

Binaries land in `bin/` (rebuilt only when sources change), fixtures
and scratch in `work/`, diagnostic plots in `plots/`.  Expected aborts
(`STOP`) are run with `check=False` and asserted via return code +
error text.

## Test inventory

| test file | production file(s) under test | what is asserted |
|---|---|---|
| `test_cmf_blkband_pin.py` | `new_main/subs/cmf_blkband_pin_null_v1.f` | Pins only exactly disconnected zero equations with zero RHS; preserves nonzero residual rows, zero-column equations, and near-zero live entries |
| `test_sigma_true.py` | `new_main/subs/recompute_sigma_true.f` | SIGMA_TRUE_CORE flag gating, unclamped `R/V dV/dR-1`, boundary cells, stagnation and `ND<3` edge cases, new-vs-reference agreement when flag off |
| `test_velocity_profiles.py` | `new_main/subs/recompute_sigma_true.f` (`SIGMA_TRUE_CORE`, `RECOMPUTE_SIGMA_TRUE`) | V001-V012 manufactured profile registry (`velocity_profiles.py`): analytic dV/dR vs complex-step; production sigma vs analytic sigma and vs an independent stencil mirror (bitwise); sign-change and compression census on every sweep (resolution, shock position/width/strength, close reversals, seeded perturbation, linear grid); h^2 refinement rate; degenerate inputs (V=0, collapsed spacing, reversed grid, infall pocket) flagged loudly via N_UNDEF/aborts |
| `test_sob_tau_floor.py` | `subs/sobjbar_sim.f`, `sobew.f`, `sobew_grad_v2.f`, `sobew_grad_v3.f` | sigma>=0 and sigma=-0.5 match reference within ULP; sigma<-1 hits the Dessart+2015 tau floor; analytic BETA check against SOBJBAR_SIM; finite output |
| `test_interp_guards.py` | `subs/lin_interp.f`, `mon_interp.f`, `mon_interp_fast.f`, `mon_interp_fast_v2.f`, `mon_int_funs_v2.f` | valid-path agreement with reference within ULP for all five routines; aborts on non-monotonic/repeated/reversed grids, out-of-range targets, `ND`/`NIN` too small |
| `test_rd_rv_file.py` | `newsubs/rd_rv_file_v2.f` | `NONMONO_V` VINF check against `MAXVAL(ABS(V))` (not `V(1)`); interior-peak profiles; 10% tolerance edges; `VINF<0.1` bypass; all-zero V, non-decreasing R and R-scaling rejections; old-behaviour adapter comparison |
| `test_sob_rh_couple.py` | `subs/sob_rh_couple.f` (new file) | zero pairs under homology; pair count/symmetry vs independent geometry; `mu^2` identity; DRC net-rate term vs recomputed value; suppression of `V<=0`/negative-opacity rows; `BB_COR=0`/`ETAL=0` safety |
| `test_sub_sob_line.py` | `new_main/mod_subs/sub_sob_line_v3.f`, `new_main/mod_subs/line_mod.f` (`DRC_SIM`) via `drv_sub_sob_line.f90` | full manufactured module-state fixture driving the production object: `DO_NONMONO_SOB` routes `SIGMA_TRUE` into `SOBJBAR_SIM` while `DO_SOB_CLAMPED_LINES` restores flag-off results bitwise; `VB_SIM`/`VC_SIM`/`SRAT` zeroed on `ETAL_MAT<0`, `NEG_OPACITY` and `BB_COR<=0` rows; valid-row `VB_SIM=-JBAR/ETAL_MAT`, `VC_SIM=JBAR*CHIL_MAT/ETAL_MAT^2`; zero-population and zero `U_STAR_RATIO` BA rows stay finite; `DO_NONLOC_SOB` fills `DRC_SIM` (and shifts the STEQ residual) while flag-off leaves it zero; `LST_ITERATION` `SOBEW` EW differs under `SIGMA_TRUE` vs clamped `SIGMA` |
| `test_do_wind_vel.py` | `new_main/subs/do_wind_vel_v1.f` via `drv_do_wind_vel.f90` | `DO_V_PARAMS` presence + `DO_NONMONO_SOB`/`DO_NONMONO_KNOP` STOPs with the documented incompatibility message; flag-on with no file returns silently (default path unchanged); flags off with a valid file parses the store and early-returns without a V revision |
| `test_ng_abstain.py` | `new_main/subs/do_ng_accel_band_v4.f` | `NG_ABSTAIN` per-component hold on zero-diff and non-positive extrapolation; siblings still accelerate; upstream whole-depth veto preserved when flag off; history restoration; vacuous guard |
| `test_rv_sniin.py` | `new_main/subs/rv_sn_model_sniin.f` | RDINR beta-law path matches reference within ULP; corrected SIGMA centred-difference denominator verified exactly against the analytic formula (reference matches the old formula, `|delta|~12`); uninitialized-V coordinate into LIN_INTERP aborts via the strict guards |
| `test_rv_snhyd.py` | `new_main/subs/rv_sn_model_snhyd.f` (new file, VEL_LAW=13) | all new guards abort (ND, missing/malformed SN_HYDRO_DATA, non-decreasing input R, `V<=0`, non-positive expansion factor, non-decreasing expanded R, TAU_OPT=1 with no positive Kappa); endpoints pinned to expanded hydro profile; monotone-cubic V on nodes; homologous V -> SIGMA=0; beta-law -> SIGMA~BETA-1; velocity reversal keeps `SIGMA<-1` under both TAU_OPT |
| `test_edep_shock.py` | `new_main/subs/get_edep_shock_power.f` | flags-off path matches reference within ULP including non-monotonic V (legacy tolerance); `DO_NONMONO_SOB`/`DO_NONMONO_KNOP`/`both` abort on non-positive or non-decreasing V and accept monotone V; Gaussian deposition profile shape/location; 4*pi normalization vs prescribed power; upstream `clump_fac=0` and `sum(wrk)=0` guards |
| `test_sob_prof.py` | `obs/cmf_flux_sub_v5.f` (`ADD_SOB_LINE_PROF`) | stored Doppler core integrates to CHIL*CLUMP; off-grid SOB line named and skipped; >10% masing-depth drop warns per line and counts N_DROP; non-monotonic NU aborts before deposition |
| `test_obs_nonmono_dispatch.py` | `obs/cmf_flux_sub_v5.f`, `obs/COMP_JCONT_V4.INC`, `obs/obs_frame_sub_v9.f`, `obs/rd_cmf_flux_controls.f`, `obs/cmf_flux_cntrl_var_mod.f` (source-text pins) | INC_GRID abort on genuinely non-monotonic V (field-detected, not flag-only); SOB+LAM_ES=F moment dispatch to MOM_J_CMF_NONMONO_V1 (SET_MOM_NU_NEXT, SIGMA_TRUE, FREQ_INDX/NCF, boundary controls) with stock V8 preserved; FIRST_FREQ one-sided dLOG_NU; DO_NONMONO_FREQ_MATCH absent-key auto-set / explicit-F override forwarded to SET_FG_FREQ_MATCH+SET_MOM_FREQ_MATCH; OBS_FRAME edge-clamp census at all three clamp sites + report; OBS_FREQ strictly-decreasing abort; CHK_OBS_FLUX at OBSFLUX/OBSFRAME phase-1/2 writes; two-phase SOB_FLUX_OBS invariants (phase-2 reload branch, coherent-ES/EDDFAC forcing, OBSFLUX phase-1 only, maser-drop counters) |
| `test_vadat_keys.py` | `new_main/mod_subs/rd_control_variables.f` (+ real `rd_store_log.f`, `control_variable_mod.f`, `jfnk_solver_v1.f90`, `check_param_consistency.f`, `fg_j_cmf_nonmono_v1.f`, `mom_j_cmf_nonmono_v1.f`) | a production VADAT parses; every new non-mono/LAM/BA/JFNK key round-trips into the module store; mutual exclusion; GLOBAL_LINE/USE_LAM_ES/LAM_SET prerequisites abort; DO_JFNK non-mono gate; JFNK_COUPLED KNOP-only; USE_JFNK deprecated alias; JFNK_LAG_FACTORS/JFNK_FRZ_GATE/JFNK_GATE_NREQ round-trip; DO_NONLOC_SOB auto-off; DO_CL/GAMRAY_TRANS legacy-coordinate rejections; IB_METH whitelist + case normalization; VEL_LAW=13 -> SN_HYDRO_MODEL |
| `test_jfnk_solver.py` | `new_main/mod_subs/jfnk_solver_v1.f90` (via `tests/jfnk_core` and the `unit_tests/se_driver` harness) | Newton/GMRES/hookstep core T1-T3 on an analytic stiff problem (PTC off, SER anchor, far-start trust region); surrogate-SE systems T0-T22: FD-JVP error curve, demanded-step SVD, PC machinery (T4) + per-depth LU residual census (T21) + JVP homogeneity/additivity (T22), GMRES, end-to-end Newton, masked-T rows, tauto ladder/FSM.  Known author-flagged diagnostics are XFAIL |
| `test_jfnk_opint.py` | `new_main/mod_subs/jfnk_op_probe_v1.f90` (+ adapters in `jfnk_cmfgen_interface_v1.f`, call site in `jfnk_solve_corr_v1.f`) | The shared operator-consistency battery (OPPROBE lines) run on the surrogate via `se_driver seprobe`: bitwise residual repeatability, 13 direction classes, FD-vs-exact Jacobian errors with per-direction h-sweeps, additivity/homogeneity at the truncation floor; module compiles standalone; production call site fenced under JFNK_HSCAN |
| `test_knop_ray.py` | `plane/fg_j_cmf_nonmono_v1.f` (`KNOP_RAY_SOLVE_V1`) | first runtime coverage of the refresh fixed-point solver on a manufactured core ray (NI=8, NCF=6, one compressing band): alternating and ordered sweeps converge (NXIT=1), committed stores agree to ~3e-8, `DO_NONMONO_FREQ_MATCH` fires multi-bin slots under strong drift, forced stall/cap exits visibly (NXIT=4), near-fixed-point seed converges in <=4 passes |
| `test_jfnk_regressions.py` | `new_main/cmfgen_sub.f`, `new_main/mod_subs/jfnk_*.f*` (source-text pins) | BUG-001 true-residual ordering pins and BUG-006 save-ordering pin (every mutating save site arms before the first mutation); guard ordering that needs CMFGEN state to run directly |
| `test_knop_lagged_field.py` | `plane/fg_j_cmf_nonmono_v1.f` (source-text pins) | fail-closed refresh commit defaults, NEW_FREQ sweep-write guards and the lagged-field drift fixes present in source |
| `test_knop_status_text.py` | `plane/*nonmono*.f`, KNOP_STATUS_MOD (source-text pins) | sentence-based status narration present; no stale format strings |
| `test_mat5pen.py` | `subs/mat5pen.f` | TRI block-Thomas kernel: all four ENOTZERO/DNOTZERO combinations match A+D(F.C-E)-B.C to ~2e-16 vs NumPy; production m[k]=C[k]-B[k].phi[k-1] update shape matches DGEMM |
| `test_mom_solve.py` | `plane/mom_j_cmf_nonmono_v1.f` via `drv_mom_solve.f` | discrete moment operator on manufactured atmospheres with controlled compression bands; doublet-spacing limit (compressing rows copy J_upwind), multi-reversal and strong-compression convergence |
| `test_mom_sweep_table.py` | `plane/mom_j_cmf_nonmono_v1.f` | GS-sweep report formats: MOM_MIN_DRA/MOM_MAX_COH argument counts match their descriptors (pre-fix formats dropped values) |
| `test_rv_hydro_rdin.py` | `new_main/subs/set_rv_hydro_model_v3.f` (VEL_LAW=11 RDINR branch) | R-grid read from RDINR, rescale to hydro inner edge, V rebuild; malformed-input rejections |
| `test_lam_handoff.py` | `new_main/mod_subs/lam_handoff_mod.f90` (via `drv_lam_handoff.f90`) + source pins in `solve_for_pops.f`, `cmfgen_sub.f`, `rd_control_variables.f`, `check_param_consistency.f` | LAM_FIELD_GATE held/frozen/raw-immature/maturity-window legs on and off the lagged-field routes; LAM_STABLE_STEP earn/reset/pin semantics incl. field-not-ready and forced/blocked pin vetoes; LAM_HANDOFF_OK stock-vs-nonmono predicates (stock keeps single `MAXCH_CTRL < 50`); HAS_LAGGED_FIELD census gate covers SOB+LAM_ES=F; latch clears, T-hold arm, blocked-handoff diagnostic, dead BA_PROBE_LAM consistency warning |
| `test_iter_ledger.py` | `new_main/mod_subs/iter_ledger_mod.f90` (via `drv_iter_ledger.f90`) | per-iteration work ledger CSV: header once, one row per write, counters/timers/state reconcile with caller values, `n_dem` counts `*demote` tokens, JFNK eval/JVP deltas vs reset snapshot, positive RSS on Linux, no file before `LEDG_ITER_RESET` (stock path silent) |
| `test_sup_to_full.py` | `subs/sup_to_full_v3b.f` (via `drv_sup_to_full.f90`) | super-level -> full-atom mapping invariants: exact population conservation, departure-coefficient sharing (`INT_SEQ=0`), LTE input returns LTE populations, zero-population super level yields zeros without NaN, extreme contrasts stay finite |
| `test_comp_obs.py` | `obs/comp_obs_v2.f` (via `drv_comp_obs.f90`) | COMP_OBS_V2 flux audit: zero-velocity limit equals CMF result, relativistic transform vs independently evaluated expectation, coverage/non-finite diagnostics STOP loudly, Jy normalization scales with RMAX |
| `test_comp_obs_retain.py` | `obs/comp_obs_v2.f` | NU_STORE retention fix: rays with FREQ_CONV_FAC>1 keep the CMF frequencies they need (old ~2-frequency store STOPped 'CMF_FREQ outside range'); fluxes still match independent evaluation |
| `test_gs_drift_diag.py` | `new_main/subs/solve_for_pops.f`, GS relax accessors (source-text pins) | GS_DRIFT_DIAG VADAT key declared/parsed; diagnostic separates raw sweep drift, refresh residual and retained-field delta with omega/residual-history/argmax/census/freeze state; FG_RF_DCHG reset at sweep entry; MOM_SW_* snapshot cannot go stale |
| `test_iter_conv_gate.py` | `new_main/subs/solve_for_pops.f`, `cmfgen_sub.f` (source-text pins) | convergence exit cannot fire while the lagged radiation field is immature: both CONV_EXIT_MET blocks conjunct `.NOT. LAM_FIELD_NOT_READY`; census block opens under HAS_LAGGED_FIELD covering SOB+LAM_ES=F |
| `test_pop_audit.py` | `new_main/mod_subs/control_variable_mod.f`, `rd_control_variables.f`, `solve_for_pops.f` (source-text pins) | POP_AUDIT VADAT keys declared/parsed; audit called under flag; routine iterates SPECIES_BEG_ID..SPECIES_END_ID, treats terminal XzV_PRES=F ion slot as one density row, sums every stage, compares against POP_SPECIES |
| `test_popoxy_conservation.py` | model `POP*`/POPOXY output files (no Fortran build) | per-species population file conservation: finite non-negative populations, per-stage sums consistent with totals, ion fractions in [0,1], closed species accounting at every depth; skips cleanly without a file |
| `test_comp_obs_v3_retain.py` | `obs/comp_obs_v3.f` (via `drv_comp_obs_v3_retain.f90`) | same NU_STORE retention fix as V2, exercised through COMP_OBS_V3: all observer bins produced, retained-window fluxes match independent evaluation |
| `test_sp_mon_interp.py` | `pgplt/do_vec_op.f` (`SP_MON_INTERP`), `pgplt/mon_interp_sp.f` (`MON_INTERP_SP`) via `drv_sp_mon.f90` | 05-Oct-2026 grid-validation guards on both single-precision monotone interpolators: valid increasing/decreasing grids complete, repeated endpoints, reversed/non-monotonic input, unordered output and out-of-range targets all STOP; linear column reproduced to ~2e-7 |
| `test_get_ibound.py` | `plane/get_ibound.f` via `drv_get_ibound.f90` | INCIDENT_INTENSITY file guards: accepts a valid descending grid + angles, rejects `NCF<2`, non-positive/non-finite and non-decreasing frequencies, out-of-range angles and out-of-order requests; linear frequency interpolation verified; out-of-range/NaN requests return zero |
| `test_eddfactor.py` | `new_main/subs/open_rw_eddfactor.f` via `drv_eddfactor.f90` | filename gate applies only to EDDFACTOR/ES_J_CONV (upstream compared FILENAME to itself); IOSTAT on the RV pointer read treats a truncated file as empty and rewrites the pointer; fresh creation writes EDDFACTOR + _INFO with correct access/compute flags |
| `test_rd_sn_data.py` | `new_main/subs/rd_sn_data.f` via `drv_rd_sn_data.f90` | NISO=-1 sentinel on the SN_HYDRO_DATA header gate: absent isotope header still STOPs, explicit `Number of isotopes: 0` now passes the gate (upstream rejected both identically), missing age header still caught |
| `test_bcd_mat.py` | `new_main/subs/read_bcd_mat.f` via `drv_bcd_mat.f90` | unrecognized TYPE now STOPs in READ/WRITE_BCD_MAT (upstream returned silently); WRITE error names the correct routine; BC/C/D scratch round-trips bitwise; missing scratch file STOPs on open |
| `test_gff.py` | `subs/gff.f` via `drv_gff.f90` | GFF_VEC bad-T diagnostic fires once at the corrupt depth for non-positive and NaN temperatures; clean runs produce finite positive Gaunt factors; scalar GFF agrees with the vector path |
| `test_jh_prev.py` | `plane/get_jh_at_prev_time_step.f` via `drv_jh_prev.f90` | IOSTAT on header/R-V/J-H records fails loudly on a truncated JH_AT_OLD_TIME; OLD_V and V must be positive and strictly decreasing before LOG(V) interpolation; manufactured 100-frequency file verifies exact-interpolation, frequency clamps, OGREY branch and boundary/time-step guards |
| `test_scr_read.py` | `subs/scr_read_v2.f` via `drv_scr.f90` | WRITFAIL flag on IOSTAT write errors (torn POINT/SCRTEMP pairs used to resume silently); multi-record bitwise round-trips across two iterations incl. explicit IREC_RD addressing; NEWMOD on absent files; WRITE_RVSIG=F keeps header R/V but skips per-iteration copies; GS checkpoint hooks linked in |
| `test_ins_line.py` | `subs/ins_line_v6.f`, `subs/ins_line_obs_v5.f` via `drv_ins_line.f90` | DO WHILE scan OOB fixes on the TRANS_TYPE blanketed-line search: all-SOB decks return cleanly (upstream read past N_LINES), BLANK lines merge into the grid with correct indices and monotonic frequencies, INCL_ALL_LINES and NU_MIN/NU_MAX handling verified; documents a still-latent OOB in the unpatched profile K-loop |
| `test_mom_gs.py` | `plane/mom_j_cmf_nonmono_v1.f` via `drv_mom_gs.f90` | MOMJNXT_GS checkpoint with real allocated GS state: atomic .tmp rename, bitwise J_NEXT/H_NEXT/MOM_QUAR_AGE round trip, tag/dims/R-V-endpoint/expanded-grid fingerprints, negative-field value rejection, truncated-file read rejection, cold-start and missing-file fallbacks; RESTORE_MOM_JNEXT_LATEST rollback restores tag-free and leaves the store untouched on fingerprint mismatch |
| `test_fg_gs.py` | `plane/fg_j_cmf_nonmono_v1.f` via `drv_fg_gs.f90` | FGINEXT_GS checkpoint with real allocated GS state: atomic rename, bitwise I_P/I_M_DEPTH_NEXT round trip, full R/V/SIGMA + COMP_DEPTH_IDX fingerprints, KN_FREQ grid sanity, stale-tag rejection vs GS_FOREIGN_SEED immature-accept (I_NEXT_GS_OK stays F), negative-intensity zeroing, truncation; RESTORE_FG_INEXT_LATEST rollback semantics |
| `test_get_pops_v5.py` | `new_main/subs/get_pops_at_prev_time_step_v5.f` via `drv_get_pops_v5.f90` | Positive/strictly-decreasing V guards on both current and OLD_MODEL_DATA velocity grids before LOG(V) interpolation; ND_OLD/=ND mixed-grid run verifies MON_INTERP remap + endpoints; inner-boundary mismatch STOPs; outer-boundary warning still dies downstream in MON_INTERP (upstream semantics) |
| `test_rd_cont_j.py` | `new_main/mod_subs/rd_cont_j.f` via `drv_rd_cont_j.f90` | IOSTAT on EDDFACTOR direct-access reads: torn RJ_HIGH/RJ_LOW records fall back to zero bounds instead of crashing the LAMBDA warm-up; mid-walk failure clamps LOW_FREQ=HIGH_FREQ keeping the last good field; linear interpolation, hi/lo clamps, ACCESS_F bookkeeping, RJ/K_MOM/SOURCE consistency; ACCURATE=T STOPs |
| `test_check_ions.py` | `new_main/subs/check_ions_present.f` via `drv_check_ions.f90` | DPTH_INDX=ND/ED_VAL=ED(ND) default when no T exceeds T_VAL: zero, positive and negative T decks all complete, fort.4 warning emitted, summary prints the inner-boundary depth instead of an undefined index |
| `test_vturb.py` | `new_main/mod_subs/set_abund_clump.f` via `drv_vturb.f90` | DO_NONMONO_SOB/KNOP VTURB laws use ABS(V)/MAXVAL(ABS(V)) for POW and MIN(ABS(V),VTURB_VEND)/VTURB_VEND for TRUNC_POW; zero-velocity range falls back to VTURB_MIN; flag-off keeps signed-V formulas incl. legacy NaN on negative V; unrecognised law STOPs |
| `test_sup_to_full_v4.py` | `new_main/mod_subs/sup_to_full_v4.f` via `drv_sup_to_full_v4.f90` | MOD_HOLD_LEV_DIS gate: flag-off recomputes POPION from Z_POP>0.01 POPS rows and refreshes B/A/X_LEV_DIS; flag-on keeps incoming POPION and dissolution constants while the super->full mapping still conserves populations |
| `test_add_dust.py` | `plane/add_dust_eta.f` via `drv_add_dust.f90` | L<RAY(IP)%NZ bound on the innermost-point scan: a tangent ray with all R_RAY above R(1) no longer reads R_RAY(NZ+1); outer points carry the R^-4 dilution fill, quadrature-weight accumulation verified at exact grid abscissae |
| `test_characteristics_v2.py` | `plane/characteristics_v2.f`, `plane/characteristics_v2_gam.f` via `drv_characteristics{,_gam}.f90` | prod-vs-ref comparison: tangent-impact clamp (0.9999*R_GRID) keeps Z real so the clamped interval refines (ref skips it); ABS(dV) insertion refines velocity-dip flanks the signed test missed (+12% NZ); well-signed dV reproduces ref NZ within 2% (clamp perturbs the tangent interval only); all rays finite, |MU|<=1 |
| `test_set_upper_axis.py` | `pgplt/set_upper_axis.f` via `drv_set_upper_axis.f90` | ED positive/strictly-increasing guard before the XV inversion: flat, decreasing, negative and zero ED all STOP with the labelled error; a valid deck returns NXED and the label |
| `test_par_opac.py` | `opac/par_opac_sub.f` via `drv_par_opac.f90` + `stub_par_opac_deps.f90` | non-monotonic V rejected at entry (stock CMF solvers only): flat, reversed and negative V all STOP with the error in OUT_FLUX; monotonic V passes the guard (dies later on missing model state, as expected) |
| `test_tau_regrid.py` | `new_main/subs/do_tau_regrid_v2.f` via `drv_tau_regrid.f90` | prod-vs-ref: FLUX_MEAN(1)=0 leaves TAU_OLD(1) at the floored TA value (ref puts LOG10(0)=-Inf in and dies); endpoints preserved, ordered finite grid out; identical result on a normal deck |
| `test_adiabatic.py` | `new_main/subs/eval_adiabatic_v3.f` via `drv_adiab.f90` | SIGMA_TRUE selection under DO_NONMONO_SOB/KNOP: the AD_CR_V ratio (3+sigma) proves the flag-on paths use SIGMA_TRUE, an unallocated SIGMA_TRUE falls back to SIGMA identical to flag-off, and the flag-off deck stays finite |
| `test_ave_dc.py` | `misc/ave_dc.f` via shipped `exe/ave_dc.exe` | V-coordinate guards on the shipped binary with real OSCDAT fixtures: nonpositive/flat V1 and V2 each die with their labelled error; identical V grids complete and write ED=sqrt(ED1*ED2) at every depth |
| `test_alteration_coverage.py` | every build input vs `../cmfgen_reference/cur_cmf` (no Fortran build) | self-checking manifest: each altered source is registered as exec/pin/gap/exempt; unregistered alterations, stale entries and missing covering tests all fail; known gaps printed loudly |

## Supporting files

- `common.py` -- compiler discovery (ifort required for production
  objects), production-object resolution (`prod_obj`/`prod_objs`:
  sibling `.o` or `lib/*.a` member with archive-header timestamp
  freshness), production flag/link probing from the tree Makefiles
  (`prod_comp_flags`, `prod_link_flags`, `link_prod`),
  `compile_and_link_prod` mixed-source builds, whole-archive
  freshness (`require_fresh_archives`), per-test object/link caching,
  process runner, `ref_close` ULP comparison, `Report` pass/fail
  collector.
- `velocity_profiles.py` -- common manufactured velocity-profile
  library (V001-V012) shared by every non-mono solver test:
  analytic R/V/dV/dR/sigma, gradient sign changes, compression
  regions, parameter sweeps, seeded perturbations, rvfile writer,
  and a Python mirror of the production hybrid sigma stencil.
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
`rv_sniin.png`, `rv_snhyd.png`, `edep_shock.png`,
`velocity_profiles.png`

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
