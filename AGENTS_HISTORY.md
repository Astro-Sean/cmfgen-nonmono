# AGENTS_HISTORY.md - cmfgen_jfnk status log and reference inventory

This file holds the sections moved out of AGENTS.md on 10-Oct-2026
so AGENTS.md stays under the ~16 kB agent rule-injection limit.
Update conventions are unchanged: dated status entries (newest first
within the Status section) and the diff-vs-pristine inventory live
here. AGENTS.md section 7 points at this file.


## 7. Full diff vs pristine (as of 29-Sep-2026)

Modified: `.gitignore`, `Makefile`, `Makefile_definitions`,
`new_main/cmfgen.f`, `new_main/cmfgen_sub.f`,
`new_main/mod_subs/Makefile`,
`new_main/mod_subs/comp_j_blank.f`,
`new_main/mod_subs/control_variable_mod.f`,
`new_main/mod_subs/do_var_cont.f`,
`new_main/mod_subs/line_mod.f`,
`new_main/mod_subs/mod_cmfgen.f`,
`new_main/mod_subs/rd_control_variables.f`,
`new_main/mod_subs/solve_for_pops.f`,
`new_main/mod_subs/sub_sob_line_v3.f`,
`new_main/subs/Makefile`,
`new_main/subs/auto_add_ion.f`,
`new_main/subs/check_ions_present.f`,
`new_main/subs/do_ng_accel_band_v4.f`,
`new_main/subs/rd_sn_data.f` (NISO=-1 sentinel so
`Number of isotopes: 0` is accepted - same 2-line fix as V6;
upstream conflated "header missing" with "zero isotopes"),
`new_main/subs/rv_sn_model_sniin.f` (SIGMA finite-difference
denominator fix, 29-Sep-2026, see status note below),
`new_main/subs/solveba_v13.f`,
`obs/cmf_flux_cntrl_var_mod.f`, `obs/cmf_flux_sub_v5.f`,
`obs/cmf_flux_v5.f`, `obs/COMP_JCONT_V4.INC`, `obs/mod_cmf_obs.f`,
`obs/obs_frame_sub_v9.f`, `obs/rd_cmf_flux_controls.f`,
`plane/get_moms_non_rel.f`, `plane/Makefile`,
`spec_plt/Makefile`, `subs/Makefile`, `subs/gff.f`,
`subs/scr_read_v2.f`,
`subs/sobew.f`, `subs/sobew_grad_v2.f`, `subs/sobew_grad_v3.f`,
`subs/sobjbar_sim.f`.

OBS / spectrum code (18-Sep-2026, DO_NONMONO_SOB plumbing):
`obs/cmf_flux_cntrl_var_mod.f` (flag, init F),
`obs/mod_cmf_obs.f` (SIGMA_TRUE array),
`obs/rd_cmf_flux_controls.f` (reads [DO_NONMONO_SOB], echoes mode +
warns DO_CMF_EW keeps clamped sigma),
`obs/cmf_flux_v5.f` (allocates SIGMA_TRUE, re-derives it from R,V after
RD_RVTJ_VEC - RVTJ stores clamped sigma only),
`obs/cmf_flux_sub_v5.f` (SOBEW_GRAD_V2 gets SIGMA_TRUE when flag set;
COMP_EW_ABS_EW keeps clamped SIGMA - it divides by GAM and is not
sigma < -1 safe),
`obs_template/CMF_FLUX_PARAM_INIT` (new F [DO_NONMONO_SOB] key).

OBS / Knop formal-solver plug-in (19-Sep-2026, DO_NONMONO_KNOP):
`obs/cmf_flux_cntrl_var_mod.f` (flag, init F),
`obs/rd_cmf_flux_controls.f` (reads [DO_NONMONO_KNOP], echoes mode),
`obs/COMP_JCONT_V4.INC` (standard-grid branch calls
FG_J_CMF_NONMONO_V1 with SIGMA_TRUE+NCF when flag set; extended-grid
site stops with an explicit error -- SIGMA_TRUE_EXT not implemented),
`obs/obs_frame_sub_v9.f` (NI_MAX bound uses sum|dV|+max|V| instead of
V(1): non-monotonic rays traverse a velocity range multiple times;
also per-thread OpenMP work arrays, bounded band hunts, sqrt-radicand
guards -- see 30-Sep-2026 robustness entry),
`obs_template/CMF_FLUX_PARAM_INIT` (new F [DO_NONMONO_KNOP] key).
Pipeline: prepare_CMF_PARAM.py knop_mode_active() stamps
DO_NONMONO_KNOP=T only when non_mono_mode=KNOP.

OBS I_NEXT restart (Oct-2026): cmf_flux builds its own frequency grid,
so the FGINEXT_GS checkpoint cannot be restored by position.
FGINEXT_GS format VER=3 adds the writer's KN_FREQ grid ahead of the
intensity records; CMF_FLUX_SUB_V5 calls SET_FG_GS_OBS_GRID(NU,NCF)
before its frequency loop, and the solver init dispatches to
RESTORE_FG_INEXT_OBS, which log-log remaps the stored field onto the
OBS grid (slot s = field of freq s+1 on either grid).  (KC,LS) indices
carry over: OBS reads NP from RVTJ and builds rays with the same IMPAR
on the same R.  VER=2 files still restore positionally in the FINAL
path but are rejected in OBS (no grid to remap).  The OBS workdir is
MODEL/obs/, so the restore tries FGINEXT_GS then ../FGINEXT_GS.

Added: `new_main/subs/recompute_sigma_true.f`,
`plane/fg_j_cmf_nonmono_v1.f`, `plane/mom_j_cmf_nonmono_v1.f`,
`subs/sob_rh_couple.f`, `pgplot_stubs/`,
`obs_template/` (pipeline OBS-job template, copied from V6 with
`NON_MONO_SOLVER` stripped - key does not exist in dev7),
`build.sh`, `build_debug.sh` (29-Sep-2026; strict-exit wrappers --
nodebug now builds `make all`, see Sect. 5),
`docs/`, `cmfgen_obsidian/`, `knop_solver_files/`,
`knop_solver_files_full/` (review snapshots), this file plus
`README.md`. Project documentation lives in `cmfgen_obsidian/`
(audits, history, solvers, control, reference); keep new `.md`
notes inside the vault.

Also built by `./build.sh`: `exe/cmf_flux.exe`, `exe/obs_frame.exe`
(the `all` target orders `d_obs` after `d_dust`/`d_disp`; the narrower
`cmfgen` target does not build obs).

## 8. Status

### JFNK frozen-tier endgame review (11-Oct-2026)

OUTGEN for `m15_2p82d_SOB_jfnk_zlaw14_fresh` (job 12153985) confirmed
repeated frozen-J GMRES stagnation and a soft-column gap (`|diff|/|BA|`
about 0.83 at it384; stiff top-F column about 0.064). The initial
interpretation that fallback permanently prevented FULL-BA was too
strong: FULL-BA resumed at it399 after the cooldown and reduced demanded
MAXCH from 32.5% to 1.35% by it402, then plateaued around 1.4-1.5%.

The prior fallback forced Lambda and a fresh BA build for finite
stagnation/cooldown; observed cooldown calls took about 366 s versus
about 70 s for most FULL calls. A caller-mode pass-through build reached
the it408 cooldown path but passed an uninitialized `T1` to SOLVEBA and
stopped (job exit 6). The cap is now initialized from `COR_LIM_IN` at
entry. Job 12153986 was canceled by its failed dependency. The corrected
binary has not been exercised in a model run.

Source review found that `STEQ_MULTI_V10` includes direct ED and ion
density BA columns. The opt-in `JFNK_LAG_FACTORS` prototype therefore
freezes `POPION` and `MOD_LEV_DIS_BLK` factors only; it does not freeze
ED/T or direct ion-density columns. It is an experimental residual
change, not verified as the physical Jacobian. Surrogate T18 passes with
live/BA gap 2.11e4 and lagged FD-vs-BA error 5.0e-7; this does not
establish CMFGEN production behavior.

The prepared YAML uses `JFNK_FRZ_GATE=3%`, two contracting solves, and
`JFNK_LAG_FACTORS=T`. It also sets `FIX_T_AUTO=T` and `DO_T_AUTO=T` per
user choice. `DO_T_AUTO` releases global FIX_T at `MAXCH_CTRL<50%`, so
this is an early temperature release at current population demand;
`JFNK_TAUTO` only damps T corrections during JFNK, not ordinary FULL-BA.
Source snapshot committed as `b063fb7`; residual-count telemetry committed
as `8380f24`; the production-file regression suite committed as `8e16f96`;
the prepared YAML committed in the parent repository as `a7696c68`. No
redeployment has occurred.
Re-run via `scripts/run_model_grid.py` only after explicit confirmation.

Iteration 407 exposed a telemetry defect: `jfnk_total_f_evals` counts
Newton-core calls but omits residual callbacks inside `JFNK_JVP_V1`, so
OUTGEN said F evaluations=0 with 12 JVPs. The reports now print both
core F calls and `JFNK_N_RESIDUAL_EVALS` (all residual callbacks); T20
guards the log fields.

Verification: production `./build.sh` succeeded after the cap and
telemetry fixes; standalone solver tests 21/21, interface tests 9/9,
harness T15-T20 passed, and the separate production-file suite passed
11/11. Full surrogate suite still reports the existing T3-mild demand
threshold failure and three T4 LU-accuracy failures; these do not touch
the new factor path. Config resolution via `--print-resolved-config`
passed. Runtime behavior of the new fallback and factor mode remains
UNVERIFIED.

- Builds clean; `DO_NONMONO_SOB` live in binary: VERIFIED.
- Optimized and debug builds both complete from a clean tree
  (29-Sep-2026; `exe/` 13:25, `exe_dbg/` 13:22). See the
  29-Sep-2026 status entry for the SIGMA denominator fix and the
  `d_lte` link failure that the OBJ_FILES restriction works around.
- Runtime behaviour (flag off = unchanged; flag on = stable Sobolev on
  a non-monotonic model): UNVERIFIED until a controlled model run.
- VEL_LAW=12 (USE_HYDRO_GRID) runtime validation: BLOCKED -- the
  interactive partition rejected submission (AssocMaxJobsLimit, job
  12023231 cancelled). No `skip_stella` + `vel_law: 12` config exists
  yet; see the 29-Sep-2026 entry.
- `DO_NONMONO_KNOP` in cmf_flux (19-Sep-2026): VERIFIED at runtime on
  the TypeIIn_10d model (OBS job 11842074). FG_J_CMF_NONMONO_V1 ran the
  full formal J (~4 calls/freq, JNU_CLAMP_DIAG 0 events); OBSFLUX is a
  clean 93969-bin spectrum. LIMITATION (fixed same day, see below):
  obs_fin_10 (OBSFRAME ray product) carried ~13% negative/NaN bins.

### obs_frame_sub_v9 bidirectional frequency hunt (19-Sep-2026)

`obs_fin_10` corruption root cause: the per-ray frequency interpolation
hunts the comoving-frequency table with a warm-started index
(`K=INDX(I-1)`) but only walks K DOWN (`DO WHILE(T1.GT.SM_FREQ_CMF(K-1))
K=K-1`).  That assumes T1 = nu_obs*(1-VMU_RAY)*GAM is monotonic along
the ray -- true for monotonic V (VMU swept by Z), false for a
non-monotonic V law where VMU_RAY reverses and T1 jumps back up in
frequency.  The stale bracket then extrapolates ETA/CHI far outside the
table -> negative/NaN intensities (~13% of bins; OBSFLUX unaffected
because it does not do the per-ray transform).

Fix: added the upward hunt (`DO WHILE(T1.LT.SM_FREQ_CMF(K)) K=K+1`)
immediately after the existing down-hunt, so the bracket always
satisfies SM_FREQ_CMF(K-1) >= T1 >= SM_FREQ_CMF(K).  For monotonic V
the added loop never executes (unchanged stock behaviour).  Table-edge
cases still hit the existing STOP guards; the band window already
carries a 1.2*MAX_VMU margin + 10-index pad.  Verified the mirrored-ray
extension (R,V same at same R, Z/VMU negated) is exact for any
spherical V(R), and that the point-1 search is safe (ray point 1 is
always the outermost near-side point, VMU(1)>0).  Same routine serves
both cmf_flux.exe and obs_frame.exe (Makefile OBJ1/OBJ2).

Rebuilt cmf_flux.exe + obs_frame.exe; OBS resubmitted through
run_model_grid.py on TypeIIn_KNOP/config_resolved.yaml (skip_final=T
-> OBS only) as job 11842754.

Follow-up (same day): job 11842754 STOPped at 'Bad K indx (loc 2b)' --
T1 landed 3.3 ppm below the lowest tabulated CMF frequency (edge
discretization undershoot at the receding ray endpoint, ML=93969
I=214).  The up-hunt correctly detected it but aborted.  User added
0.1% edge-tolerance clamps at all three boundary hits (loc 1, loc 2
top -> K=2, loc 2b bottom -> K=SM_NCF): within-tolerance T1 clamps to
the edge bracket (tiny extrapolation of a smooth continuum); genuine
off-grid T1 still aborts with diagnostics.  cmf_flux.exe/obs_frame.exe
rebuilt clean (forced recompile after stale-.o suspicion -- on
NFS/lustre always `rm the .o` before make to rule out clock skew).

INFRASTRUCTURE BUG FOUND + FIXED (run_model_grid.py): the resubmit
wipe's restart_files preserve list had no POINT1/POINT2 -- the
SCR_READ_V2 checkpoint files (R,V,SIGMA,POPS + nonmono GS radiation
state).  They are the ONLY warm-restart state; losing them forces a
cold start.  Worse, POP<SPECIES> files (required by cmf_flux -- no
POP* => zero species => N_LINE_FREQ=0 => segfault in INDEXX) are only
rewritten at CMFGEN *completion*.  After repeated OBS-only
resubmissions wiped the TypeIIn_KNOP run dir, OBS jobs failed with
the INDEXX segfault (11843842).  Fix: POINT1/POINT2 + POP* (already
globbed) + VADAT (skip_final only) now all preserved; verified live
on an M15 OBS-only submission.

Recovery: no POP*/POINT* backups exist for TypeIIn_KNOP, so the only
path is a FINAL rerun to completion.  Per user request: resubmitted
configs/TypeIIn_10d_DEV7_KNOP.yaml with init_num_its=100 (FINAL
completes at it100 unconverged -> writes POP* -> chained OBS tests
spectral generation on the velocity-reversal state).  FINAL=11844120,
OBS=11844121; IN_ITS NUM_ITS set to 100 directly (the slurm only
re-instates FINAL_NUM_ITS on requeue/zero, the static IN_ITS=60 would
have capped at 60).

Outcome (19-Sep-2026, second half):

- FINAL ran cleanly through 3 requeue cycles (it1 cold start, ~4-5
  min/it, zero capped/NaN GS writes throughout, cells>0.5 collapsed
  706M -> ~10M).  User stopped the population solver manually at it62
  (IN_ITS NUM_ITS=0); CMFGEN finalised, wrote all 8 POP* + RVTJ +
  POINT1/2 and exited COMPLETED.  Chained OBS 11844121 ran but
  STOPped at 'Bad K indx (loc 2b)' with a NEW signature: T1=3.4954e-3
  vs SM_FREQ_CMF(SM_NCF)=3.4990e-3 at I=308, i.e. 0.101% below the
  edge -- just over the 0.1% tolerance.

- Root cause: not roundoff this time.  FREQ_CMF(NCF)=3.49897e-3 is the
  absolute grid bottom (END_INDX(93)=NCF, pad exhausted): the CMF
  frequency grid built for the SE solve under-covers the
  observer-frame Doppler reach at deep far-side ray points on the
  non-monotonic field.  dev6's equivalent up-hunt
  (K.LT.SM_NCF guard, no error check) silently extrapolated below the
  table -- one source of its corrupted bins.

- Fix (obs_frame_sub_v9.f): replaced edge extrapolation with a
  boundary clamp at all three edge sites.  T1 within 2% of the edge
  is clamped to the edge frequency so T2=0 (or 1 at the top) selects
  the boundary opacity exactly -- bounded, cannot go negative, and
  the deep-ray-point contribution is negligible anyway.  Grossly
  off-grid T1 (>2%) still aborts with diagnostics.  Comments updated
  to state the true mechanism (grid undercoverage, not roundoff).

- OBS resubmitted OBS-only via SKIP_FINAL=1 (job 11850776): all
  preserved files verified restored (POP*x8, RVTJ, VADAT, OUTGEN,
  POINT1/2).  COMPLETED in 15:47.

- RESULT: obs_fin_10 = 93,906 bins, ZERO negative / NaN / Inf / zero
  bins (previous run: 5,986 negative + 6,640 NaN ~= 13% corrupted).
  Zero spike bins >5x the 21-bin median.  8 dips <0.2x cluster at
  4923.5 A (He I 4922) and 5018 A (Fe II 5018) -- physical narrow
  absorption, not artifacts.  Balmer series strong: Ha peak 2.74e7 Jy
  (9.4x local continuum, dominant line), Hb 9.4x, Hd 4.8x.
  Qualitatively Type IIn-like: strong Balmer emission on a
  line-rich continuum, dominant Ha.

- Caveats: spectrum is from the it62 unconverged state (diagnostic
  only); MEANOPAC 'zero or negative flux-mean opacity' warning is
  non-fatal at this convergence level; the edge clamp used boundary
  opacity for the small set of off-grid deep ray points; line
  attribution not yet performed.

### obs_fin vs OBSFLUX discrepancy: broad blue emission is REAL; CV_BOUND under-converged (19-Sep-2026)

User report: obs_fin_10 spectrum shows an additional broad blue
emission component on many lines (e.g. Halpha wing at v ~ -2000 to
-5000 km/s), suspected gridding bug.  Investigation with the per-ray
diagnostics (WR_IP -> IP_DATA, WR_CMF_IP -> IP_CMF_DATA, IP_FG_DATA
file-existence trigger -> solver CV_BOUND per ray per freq):

- On monotonic M15, obs_fin/OBSFLUX agree at ratio ~1 everywhere.
  On the non-monotonic TypeIIn they differed by 55-94x in the Ha
  blue wing -- a non-monotonic-path effect only.

- Per-ray comparison (cont 5200A/7000A): IP_DATA (obs_frame
  independent ray integration) vs CV_BOUND*2 (=IPLUS -> IP_CMF_DATA
  -> OBSFLUX).  Agreement ~1.0 for rays with p > ~3e5 (units 1e10
  cm); interior rays p < R_shell ~ 1.4e5 differ by 300-500x.
  The discrepancy boundary is exactly the dense interaction shell
  (rho jumps to ~5e-13, V jumps 160->5800 km/s at R~1.4e15 cm,
  tau_es~1 to the shell top).

- DECISIVE test -- bolometric boundary flux.  Model luminosity
  L=9.85e8 Lsun => int_0^1 I mu dmu = 2H = L/(8*pi*R1^2) ~ 5.4e8
  cgs.  IP_DATA integrates to ~1.3e8 (order-consistent; residual
  factor ~2-4 from unconverged it62 state + quadrature).  CV_BOUND
  integrates to ~6.4e6 -- ~43-84x BELOW the luminosity the
  atmosphere must emit.  Energy conservation is violated by the
  solver's boundary intensity; obs_frame does not over-produce.

- Cross-check: int CV_BOUND dmu = RJ(1) = 7.54e-9 at 5200A exactly
  -> the solver is internally consistent (CV_BOUND is the same
  intensity field that fed J for SE).  RJ depth profile shows
  radiation trapped below the compressing shell: RJ(75)~8.5e-6
  vs RJ(1)~7.5e-9 (~1100x contrast).  obs_frame's interior
  intensity ~1e-5 ~ shell-depth J -- physically correct emergent
  intensity for a tau~1 line of sight to the hot shell.

- Mechanism: in compressing zones (GAM<0, 70/150 depths incl. the
  shell DeltaV~5600 km/s jump) the effective opacity is
  chi_hat = chi + |GAM|/dLOG_NU -- enormous across the jump, so
  DTAU_hat walls off the interior.  The compensating source
  |Q|*I_NEXT (lagged Gauss-Seidel upwind intensity, INEXT_RELAX=0.5)
  re-injects photons advected in from the adjacent frequency.  At
  it62 the lagged fields were still far from converged -> the shell
  acts as a wall -> CV_BOUND shows only outer-wind emission
  (~2e-8 ~ surface S) and OBSFLUX shows spurious deep P-Cygni-like
  absorption where obs_fin correctly shows shell emission.

- VERDICT: the broad blue emission component in obs_fin_10 is REAL
  -- it is line emission from the approaching shocked shell + fast
  inner ejecta (V ~ 4500-8000 km/s toward the observer), not a
  gridding artifact.  obs_fin_10 (OBS_FRAME_SUB_V9) is the
  trustworthy product.  OBSFLUX (CV_BOUND from FG_J_CMF_NONMONO_V1)
  is the defective product until the GS coupling converges.

- Pipeline bug found+fixed: FGINEXT_GS / MOMJNXT_GS (the GS lagged
  I_NEXT/J_NEXT checkpoints) were NOT in the run_model_grid.py
  preserve list -- every resubmission cold-started the coupling,
  preventing convergence through the workflow.  Added both to
  restart_files (POINT1/2 do NOT contain this state despite the
  comment claiming so).

- Follow-up: FINAL resubmitted (11851691, NUM_ITS=200, cold start as
  single_final always is) to let the compressing coupling converge;
  OBS 11851692 chained.  When it runs, re-create IP_FG_DATA (touch)
  in obs/ before the formal solve to re-check CV_BOUND bolometric
  recovery vs IP_DATA.

### SN-model parameter audit (18-Sep-2026)

Parameters checked for SN-snapshot correctness in the generated
VADAT/IN_ITS (M15_HeOCa sanity run):

- `IT_ON_T`, `DO_GT_AUTO`: grey T overwrite of hydro T - OFF (see
  note below).
- `FIX_NE`: was T (pinned hydro Saha estimate). Now F - hydro Ne is a
  rough estimate; CMFGEN solves Ne from charge conservation. Changed
  in the single_final + final profile dicts and the sanity config.
- `INIT_TEMP_V2`/`T_INIT_TAU`: dead path for `SN_T_OPT=USE_HYDRO`
  (T/ED set by SET_ABUND_CLUMP, untouched). Would re-interpolate T on
  a tau scale if SN_T_OPT=USE_T_IN - avoid that option.
- `FIX_T=F` + `FIX_T_AUTO=T` + `DO_T_AUTO=T` (IN_ITS): staged
  release - T held fixed until MAXCH<50%, then freed. Reasonable for
  a hydro-sourced T; OUTGEN "Temperature held fixed" is this
  mechanism, not a bug.
- `LTE_EST=T`: cold-start LTE guess only; populations are iterated.
  Fine.
- `ST_INT=40`/`END_INT=55`, `ND_QUAD=0`: extended-grid refinement
  range - dead because `N_INS=0`. If N_INS is ever >0 these restrict
  refinement to depths 40-55; for ejecta the whole grid should be
  refined instead.
- `GREY_TAU=10`: only consumed by the IT_ON_T grey loop - inert.
- `D2_EQ_D1=T` (IN_ITS): pins outermost BA depths - standard outer
  BC helper. Fine.
- `T_EXC=0.7` = T_EXCITE_MIN: floor for LTE departure coefficients.
  Fine.
- `ADJUST_R_DEFAULTS`: live tau-regridding control file (required
  when REV_RGRID=T, ignored when F). Now excluded from run-dir copies
  with ADJUST_CORRECTIONS/IT_SPECIFIER/table*.dat - a REV_RGRID=T run
  without it fails loudly instead of silently regridding.
- `SN_DEBUG=T`, `SN_DEBUG_LUM=T`: verbose SN diagnostics only -
  cosmetic.
- `DIF=F` + `IB_METH=ZERO_FLUX` + `LUM_METH=ZERO` + `LSTAR=0`:
  consistent zero-flux inner boundary for ejecta with no photosphere.
- `DO_NOMONO_KNOP=T`/`DO_NOMONO_SOB=F` stamped by the pipeline: V6
  keys, ignored by dev7 readers (dev7 flag is DO_NONMONO_SOB).

- OBS spectrum transfer for SOB populations (06-Oct-2026):
  `prepare_CMF_PARAM.py` stamps both DO_NONMONO_KNOP and DO_NONMONO_SOB,
  with GLOBAL_LINE=SOB and all TRANS_*=SOB. The OBS reader dispatches
  Knop formal J for the continuum/electron-scattering source and retains
  SIGMA_TRUE Sobolev EW/frequency handling. Main VADAT population modes
  remain mutually exclusive. `COH_ES=F` is used in the matched SOB and
  frozen-J configs; prepare_CMF_PARAM propagates it and ensures NUM_ES>=2
  plus IIn wing extents. Build verified; runtime OBS validation pending.

### Sanity run note (18-Sep-2026, jobs 11810955/11810956)

The config `configs/M15_HeOCa_25d_DEV7_sanity.yaml` (Stockholm
M15-He-O-Ca-evo snapshot, t=25.10 d, ND=150, monotonic homologous
ejecta) gives a **good starting setup for Stockholm-style models**:
CMFGEN passes VADAT checks, reads `SN_HYDRO_DATA` (incl. NISO=0 after
the rd_sn_data fix), starts Lambda iterations and does not immediately
blow up. Key settings that get it this far: `STOP_IF_BP=F` (pipeline's
REL_OBS=F/USE_LAM_ES=T trips the upstream SN check otherwise),
`LSTAR=0`, `IB_METH=ZERO_FLUX`, `VTURB=50`, `IT_ON_T=F`,
`DO_GT_AUTO=F`, `REV_RGRID=F`, `N_OB_INS=0`, `DIF=F`, 8 CPUs on the
interactive partition. Convergence to a finished spectrum is still
UNVERIFIED.

Note IT_ON_T = F works much better

**Why IT_ON_T must stay off for LSTAR=0 ejecta models (investigated
18-Sep-2026):** `IT_ON_T=T` runs a grey-atmosphere T re-estimate in
`SET_NEW_MODEL_ESTIMATES` before iteration 1. `COMP_GREY_V4` builds
TGREY from `HFLUX propto LUM`, i.e. the VADAT `LSTAR`. With `LSTAR=0`
(pure ejecta, no photospheric luminosity) the grey solution has
`HFLUX=0` and TGREY collapses to the `0.95*T_MIN` floor (~2850 K).
The correction `T = w*TGREY + (1-w)*T_hydro` applies with weight
`w = 1-exp(-TA/GREY_TAU) ~= 1` at every optically thick depth, so the
correct hydro T is replaced by the collapsed floor. Electron density,
LTE departure coefficients and initial populations are then computed
on the wrong T, and the first BA corrections explode. IT_ON_T is only
valid for photospheric models with a constant luminosity diffusing
through the atmosphere; for ejecta snapshots leave it F.

### TRIDIAG audit (18-Sep-2026, job 11812352)

Question: is SOL_METH=TRIDIAG (NUM_BNDS=3) stable with the current
static-snapshot setup? Answer: under this configuration the BA matrix
has EXACTLY ZERO off-diagonal depth blocks, so the block-Thomas solve
in CMF_BLKBAND_V3 reduces to independent per-depth solves - TRIDIAG
is numerically identical to DIAG.

Why all off-diagonals vanish here:

- VJ = 0: do_var_cont.f line ~331 `ELSE IF (USE_LAM_ES)THEN` is an
  EMPTY branch (upstream design, same in V6). USE_LAM_ES=T means TX/TVX
  (dJ/dchi,deta matrices) are never computed, so dJ_LOC=0 and VJ=0.
  Every J-coupling term in BA_UPDATE_V7 is VJ-multiplied.
- dZ (CMF resonance-zone line coupling, update_ba_for_line.f): off-
  diagonal entries are built only from dJ_LOC -> zero.
- sub_sob_line_v3 BA off-diagonals: VJ-multiplied -> zero.
- steq_advec_v4: returns immediately (INC_ADV=F).
- Diffusion T-coupling: dead (DIF=F; would need NUM_BNDS=ND anyway).
- All steq_* writes go to band M=DIAG_INDX only.
- CONT_VEL=T -> the F2DA/FC ray-by-ray path (which could fill VJ) is
  not taken.

Verified empirically: TRIDIAG job 11812352 and DIAG job 11811962 on
the same M15_HeOCa_25d snapshot produce bitwise-identical correction
trajectories (it5: 21892.4377234943 both; it6: 135.113165529951 vs
...52, a 1-ulp storage-path difference).

Implications:

- The earlier "TRIDIAG is unstable" observation was confounded by the
  IT_ON_T/DO_GT_AUTO/FIX_NE bugs fixed since - not by the solver.
- TRIDIAG adds no stability benefit in this config either; it is a
  no-op. Its only value here is exercising the banded-solver path.
- For TRIDIAG to carry REAL coupling you would need USE_LAM_ES=F
  (moment-equation J via MOM_J_CMF with clamped sigma) or INC_ADV=T -
  a different, riskier J solution. Not advisable for the non-mono SOB
  path where the moment solver is intentionally bypassed.
- The stale comment in prepare_VADAT.py claiming "TRIDIAG forces
  USE_SELECTIVE_JVAR -> VAR_MOM_J_CMF_NONMONO_V1" describes V6-only
  machinery that does not exist in dev7. Dev7 TRIDIAG uses the
  standard upstream CMF_BLKBAND_V3.
- prepare_VADAT SOB forced-dict still stamps SOL_METH=DIAG
  unconditionally - harmless under this config (equivalent), kept as
  conservative policy until real coupling exists (Knop plug-in).

Standalone test submitted: configs/M15_HeOCa_25d_DEV7_tridiag.yaml
(depth_independent_solver: false), model dir
dev7_sanity/M15_HeOCa_25d_DEV7_TRIDIAG. Generated files verified:
VADAT SOL_METH=TRI, MODEL_SPEC NUM_BNDS=3, N_TYPE=N_ON_J,
IT_ON_T=F, FIX_NE=F, DO_GT_AUTO=F. Jobs 11812352 (FINAL, running) /
11812354 (OBS, pending). Early iterations match the DIAG baseline.

### SN_HYDRO_DATA optional headers (18-Sep-2026)

Stockholm/ejecta-only SN_HYDRO_DATA files carry no "Energy Lum [lsun]"
or inner-boundary headers. The pipeline no longer requires them:

- proxy_populations_v2.parse_sn_hydro_data: derives r_inner_cm and
  t_inner_k from the innermost grid point, and estimates the ejecta
  luminosity as the blackbody at the tau=2/3 surface built from the
  kappa column (_ejecta_lum_lsun; falls back to the inner-boundary
  blackbody when kappa is unusable). Matches the header's established
  photospheric-flux meaning; the STELLA converter wrote |flux| at
  tau=2/3. For M15_HeOCa_25d this gives 3.2e9 Lsun (vs CMFGEN
  converged ~6.8e8 - right order for a proxy T_rad anchor).
- convert_sn_hydro_to_superlite.py: the bogus 1e42 Lsun default
  replaced by the same photospheric estimate.
- smooth_sigma_hydro_data.py: fixed-position header/section parsing
  replaced with keyword scans (no blank-separator assumption);
  energy_lum=None when absent (writer already skips it).
- run_model_grid.py SEDONA prep already falls back to the
  inner-boundary blackbody when the header is missing - unchanged.

### USE_LAM_ES=F / dev6 A/B test (18-Sep-2026)

Audit conclusion on dev6: its STANDARD solver path is bit-for-bit
upstream (do_var_cont.f, cmf_blkband_v3.f identical). USE_LAM_ES=F works
in both trees via the moment-equation J path; dev6's instability traces
to the parallel nonmono stack (SOLVE_FOR_POPS_NONMONO_V1 et al., gated
by NON_MONO_SOLVER = DO_NOMONO_KNOP .OR. DO_NOMONO_SOB), which this tree
deliberately does not import.

USE_LAM_ES=F requirements (upstream check_param_consistency.f): for an
SN_MODEL, one of USE_J_REL / USE_DJDT_RTE / USE_LAM_ES must be T. The
Stockholm reference (tmp/dev6_obs_test/M15_model, converged ~40 its)
ran USE_J_REL=T + REL_OBS=T + REL_CMF=T + TRIDIAG + NUM_BNDS=3 +
IT_ON_T=F - real tridiagonal coupling via VAR_JREL_V5 filling TX/TVX.
IT_ON_T=F was the decisive stability fix: with LSTAR=0, IT_ON_T=T runs
a grey-T re-estimate that collapses the temperature before iteration 1.

Pipeline changes enabling the test:

- prepare_VADAT.py enforce_sn_iin_physics() now takes user_overrides;
  keys the user sets explicitly in vadat_overrides are preserved even
  when not in user_overridable (previously DO_NOMONO_KNOP=T and
  REL_OBS/REL_CMF=F were stamped unconditionally).
- run_model_grid.py enforce_cmfgen_sn_iin_settings(): REL_OBS/REL_CMF
  stamps use _force_unless_user; an explicit YAML vadat_overrides value
  wins.
- plot_sn_hydro_cmfgen_inputs.py: keyword-scan headers; missing Energy
  Lum falls back to the inner-boundary blackbody for the title only.

A/B runs submitted (identical inputs, only the binary differs):

- dev6: configs/M15_HeOCa_25d_DEV6_jrel.yaml, install_dir
  cmfgen_nonmono_dev_V6, model dir dev6_sanity/M15_HeOCa_25d_DEV6_JREL.
  VADAT verified: USE_J_REL=T, USE_LAM_ES=F, REL_OBS=T, REL_CMF=T,
  SOL_METH=TRI, NUM_BNDS=3, DO_NOMONO_*=F (standard stack),
  IT_ON_T=F, FIX_NE=F. Jobs 11813665 (FINAL) / 11813666 (OBS).
  OUTGEN confirms "Calling MOM_JREL_V9".
- dev7: configs/M15_HeOCa_25d_DEV7_jrel.yaml, same solver config, model
  dir dev7_sanity/M15_HeOCa_25d_DEV7_JREL. Jobs 11813746 / 11813750.

Caveat: USE_LAM_ES=F puts J on the moment solver, whose sigma clamp is
inconsistent with non-monotonic fields - for the DO_NONMONO_SOB path
keep USE_LAM_ES=T + DIAG. The JREL+TRIDIAG config is the correct
reference/sanity configuration only for monotonic models.

### Frozen-Jacobian stall diagnosis + fix (18-Sep-2026)

The DIAG sanity run plateaued at MAXCH ~= 0.357 for ~175 iterations and
the job exited at it427 without converging. Mechanism:

- `COMP_BA=F` + `STORE_BA_INV=T` + `N_FIX_BA=0` + `FIX_BA=15`: once
  MAXCH < 15, `COMPUTE_BA` is permanently disabled and every iteration
  solves with the stored inverse of the last Jacobian (a stale-Jacobian
  fixed point). OUTGEN shows `read BAMATPNT` + `BA matrix **NOT**
  computed`.
- `DO_NG=T` could not rescue it: Ng acceleration requires
  MAXCH < BEG_NG (=0.1) and MAXCH sat at 0.357, so NEXT_NG stayed at
  1500 forever.
- Dominant residual was at depth 1 (~0.65% decrease) and depth 42
  (~0.16% increase); the stale inverse never refreshes the coupling
  that drives them.

Fix (YAML only, no source change): `FIX_BA=1.0D0` + `N_FIX_BA=10` in
vadat_overrides. N_FIX_BA=10 selects the periodic-recompute branch
(fresh BA when MAXCH > VAL_FIX_BA on the 10-iteration check, plus the
MAXCH_SUM escape) and FIX_BA=1 delays any freeze until below 1%.

Resubmitted (fresh start - the orchestrator refresh does not preserve
population `_IN` files, so these runs re-converge from scratch, a
cleaner test of the fix):

- DIAG: configs/M15_HeOCa_25d_DEV7_sanity.yaml, jobs 11816623 (FINAL) /
  11816624 (OBS). VADAT verified: N_FIX_BA=10, FIX_BA=1.0D0.
- TRIDIAG: configs/M15_HeOCa_25d_DEV7_tridiag.yaml, jobs 11816688 /
  11816689. VADAT SOL_METH=TRI, MODEL_SPEC NUM_BNDS=3 verified.
- Old BA stores were stashed under dev7_sanity/ba_stash/, not deleted.

### JREL+TRI segfault: upstream bug, not a dev-tree regression (18-Sep-2026)

The dev6/dev7 JREL A/B runs (11813665, 11813746) crashed identically:
~27-33 clean Lambda iterations, first real TRIDIAG BA solve returns
exactly 0.0000 corrections, segfault next iteration
(`forrtl: severe (174): SIGSEGV` + `array index out of bounds`).

Bounds-checked debug build (`make DEBUG=debug` -> exe_dbg/, ifort -O0
-check bounds -traceback; job 11815551 on tmp/jrel_dbg):

- Root cause (VERIFIED): `Error in LINPAC (or BLAS) routine DGETRF_2
  in CMF_BLKBAND, depth=1, IFAIL=11` -> `Error in CMF_BLKBAND_V3 -
  shutting code down`. The LU factorization of the banded matrix hits
  a zero pivot at depth 1 ON THE FIRST TRI SOLVE (debug run started
  in BA mode at it1; production runs reached it only after ~33
  Lambda its). The block-Thomas banded matrix is structurally
  singular at the outer boundary when JREL off-diagonal coupling
  populates NUM_BNDS=3.
- Production fallout (consistent): the failed TRI solve returns
  zero/garbage corrections ("0.0000" before each crash), corrupting
  later state -> `gff.f:101` bounds abort with `IGB=-2147483648`
  because `LOG10(T(I))` evaluated on `T(ND) <= 0` (NaN -> INT
  overflow) via genopaeta_v10.f:146 -> comp_opac -> cmfgen_sub.f:1485
  (DTDR inner-boundary opacity block).
- All TRI-path files (cmf_blkband_v3, solveba_v13, ba_update_v7,
  do_var_cont, var_jrel_v5, generate_full_matrix_v3) are
  bit-for-bit upstream in dev7 - this is an upstream bug in
  USE_J_REL=T + TRIDIAG, consistent with the JREL path being
  unexercised on current upstream. Fix is open; JREL runs should use
  SOL_METH=DIAG until the depth-1 band structure is understood.
- Latent bug found en route (kept, minimal): `check_ions_present.f`
  used uninitialized `DPTH_INDX` to index `POP_SPECIES` when all
  T(K) <= T_VAL on the startup path. Now initialized before the loop.
- `Error opening BAMATPNT in READBA IOSTAT=29` is benign (missing
  saved-inverse file on first start; handled).
- Temporary diagnostic prints in new_main/cmfgen_sub.f (iteration
  start + DTDR entry T(ND)/minT) were REMOVED 18-Sep-2026 during the
  Knop integration (they printed unconditionally on the default path).
  The subs/gff.f bad-T report remains in the working tree - it is in
  the -check bounds build only (exe_dbg); the optimized exe/ build is
  unaffected.

### DO_NONMONO_KNOP integration (18-Sep-2026)

Knop et al. (2009) non-monotonic formal solver ported from dev6 as a
minimal, opt-in path. Only the formal solver was ported - the V6
population solver fork (SOLVE_FOR_POPS_NONMONO_V1, JFNK, runaway
guards, moment solver MOM_J_CMF_NONMONO_V1) was deliberately NOT
imported.

Activation: `T [DO_NONMONO_KNOP]` in VADAT (orchestrator: set
`DO_NONMONO_KNOP: T` in cmfgen.vadat_overrides plus
`non_mono_mode: KNOP` - non_mono_mode alone does not activate the
enforcement block). Main VADAT population mode is mutually exclusive
with DO_NONMONO_SOB (Fortran STOP if both T). Requires GLOBAL_LINE=BLANK
(or NONE with every present ion TRANS_*=BLANK), LAM_SET=F, INC_GRID=F,
FG_OPT=INT/INS, no FIXED_J/rel/plane-parallel modes - all enforced in
rd_control_variables.f; GLOBAL_LINE and LAM_SET gates were added
03-Oct-2026. OBS uses the same prerequisites for Knop-only spectra, but
allows both OBS flags T for Knop formal J with SOB line frequencies.
DO_BA_RES_ACC=F warns.

J-update modes (both valid since 19-Sep-2026):
- USE_LAM_ES=T : J from the formal solution (Lambda iteration).
- USE_LAM_ES=F : J from MOM_J_CMF_NONMONO_V1 - the non-monotonic
  variable-Eddington-factor moment solve (ported from dev6):
  chi_hat = chi+|GAM|/dlnu with GAM=(V/R)c^-1(1+SIGMA_TRUE*K_ON_J),
  PSI/W diagonal uses |1+sigma*f| (always positive), signed coupling
  in the RHS via J_upwind = RSQJNU_PREV (Q>=0) or J_NEXT_FRQ(:,k+1)
  (Q<0, lagged Gauss-Seidel over great iterations).  J_NEXT_FRQ is
  checkpointed to MOMJNXT_GS by DUMP_MOM_JNEXT_GS (SCR_RITE_V2) and
  validated on restore (tag/dims/expanded-grid fingerprint) by
  RESTORE_MOM_JNEXT_GS inside the solver's first call; a failed
  restore leaves J_NEXT=0 and the N_LAM_RWARM warm-up covers it.
  The SOLVE_FOR_POPS warm-up skip now requires BOTH the FG I_NEXT
  and (when the mom solver ran) the MOM J_NEXT checkpoints to be
  mature (GET_MOM_GS_MATURITY accessor, same build-order workaround
  as GET_FG_GS_MATURITY).

  NOT ported (deliberate): VAR_MOM_J_CMF_NONMONO_V1 /
  EDD_J_VAR_NONMONO_V1 / J_NEXT_MD / STORE_J_NEXT_MD - the V6
  linearization solver for the selective-JVAR ALO stack.  Under
  DO_NONMONO_KNOP/DO_NONMONO_SOB the stock VAR_MOM_J_CMF_V12 is
  ALSO bypassed (19-Sep-2026): do_var_cont.f gained a fixed-J +
  diagonal-ALO branch that zeroes TX/TVX and fills only
  TX(K,K,1)=ALO_DJDCHI(K), TX(K,K,2)=ALO_DJDETA(K) from the formal
  solver's own per-frequency dJ/dchi, dJ/deta (GET_FG_ALO_DIAG).
  V12 linearizes the monotonic (clamped-sigma) moment equation and
  produced an all-NaN TX at the first post-warm-up BA iteration
  (M15 it9) - same failure V6's ~5000-line replacement stack was
  built to avoid.  NUM_BNDS>1 under KNOP gets no cross-depth
  J-coupling; port the V6 selective-JVAR stack properly if TRI is
  ever needed rather than re-enabling V12.

Files touched (diff vs cmfgen_reference/cur_cmf):
- plane/fg_j_cmf_nonmono_v1.f    NEW (copy of dev6 file + added
  GET_FG_GS_MATURITY accessor; new_main compiles before plane so a
  module USE is not possible there)
- plane/mom_j_cmf_nonmono_v1.f   NEW (dev6 file minus J_NEXT_MD /
  STORE_J_NEXT_MD which only serve VAR_MOM_J_CMF_NONMONO_V1; added
  GET_MOM_GS_MATURITY accessor)
- plane/Makefile                 fg_j_cmf_nonmono_v1.o +
  mom_j_cmf_nonmono_v1.o + rules
- new_main/mod_subs/control_variable_mod.f  DO_NONMONO_KNOP (def F),
  N_LAM_RWARM (def 4)
- new_main/mod_subs/mod_cmfgen.f NU_NEXT_FREQ scalar
- new_main/mod_subs/rd_control_variables.f  reads + guards;
  USE_LAM_ES=F under KNOP prints an informational note (non-mono
  moment solver) instead of STOPping since 19-Sep-2026
- new_main/mod_subs/comp_j_blank.f  Knop-vs-V13 dispatch,
  GET_MOMS_NON_REL under USE_LAM_ES, FIRST_FREQ one-sided dLOG_NU,
  MOM_J_CMF_NONMONO_V1-vs-V11 dispatch under USE_LAM_ES=F -
  all gated on DO_NONMONO_KNOP
- new_main/cmfgen_sub.f  SET_FG_GS_TAG(NITSF) + SET_MOM_GS_TAG(NITSF)
  after SCR_READ, NU_NEXT_FREQ in both freq loops,
  RECOMPUTE_SIGMA_TRUE gates extended to KNOP
- new_main/mod_subs/solve_for_pops.f  ITS_SINCE_START counter +
  N_LAM_RWARM LAMBDA warm-up (skipped on mature FGINEXT_GS +
  MOMJNXT_GS restores), FORCE_LAMBDA guards the mode controller,
  SIGMA_TRUE refresh gates
- subs/scr_read_v2.f  DUMP_FG_INEXT_GS + DUMP_MOM_JNEXT_GS at end of
  SCR_RITE_V2 (self-guard on unallocated arrays; zero cost stock)

Untested/known limitations:
- N_COMP=0 path only so far (monotonic sanity model): the
  compression-zone I_NEXT machinery is compiled but not exercised.
- USE_LAM_ES=F is now supported under KNOP via MOM_J_CMF_NONMONO_V1
  (ported 19-Sep-2026).  Before the port this config hard-STOPped -
  verified at runtime (18-Sep-2026, job 11818525 exit 6, message
  `DO_NONMONO_KNOP=T requires USE_LAM_ES=T`).  The pipeline does not
  force-correct USE_LAM_ES in KNOP mode (it stays forced in SOB
  mode), so the VADAT carries the user's choice explicitly.
  Validation status: M15_HeOCa_25d_DEV7_knop_MOM (monotonic, Q>=0
  branch only) runs and the MOM GS sweep max change decays;
  TypeIIn_10d_DEV7_KNOP_MOM exercises the Q<0 / J_NEXT branch at
  ~85% of depths - convergence there is the open question.
- FG_J_CMF_NONMONO_V1 warns `SIGMA_TRUE_IN is all zeros` at most once
  per process (SIGMA_ZERO_WARNED flag) when the field is identically
  zero - expected on a homologous model; retained as a diagnostic for
  genuinely non-monotonic inputs where a zero profile indicates a
  missing RECOMPUTE_SIGMA_TRUE.  The warning caught a real bug
  (19-Sep-2026): RECOMPUTE_SIGMA_TRUE's internal gate was
  `DO_NONMONO_SOB`-only, so KNOP-only runs kept the zero-initialised
  SIGMA_TRUE (masked on the homologous sanity model where sigma=0 is
  the true answer).  Gate fixed to
  `.NOT. DO_NONMONO_SOB .AND. .NOT. DO_NONMONO_KNOP`.

Validation (18-Sep-2026): `configs/M15_HeOCa_25d_DEV7_knop.yaml`
(monotonic M15 sanity snapshot + DO_NONMONO_KNOP=T) submitted via
run_model_grid.py, job 11817508. Confirmed in OUTGEN:
`Calling FG_J_CMF_NONMONO_V1`, chi_hat path, USE_LAM_ES branch,
`no compressing points - GS storage skipped`, GS restart fallback to
LAMBDA rebuild, `SOLVE_FOR_POPS: KNOP GS warm-up` forcing LAMBDA for
the first 4 its. Convergence of a full non-monotonic model is still
to be demonstrated (e.g. the TypeIIn 10d snapshot that diverged on
the SOB path).

### TypeIIn 10d non-monotonic SOB test (18-Sep-2026)

Config `configs/TypeIIn_10d_DEV7_SOB.yaml`, model dir
`TypeIIn/cmfgen/10.00days`, snapshot ~9.84 d post-explosion from the
existing TypeIIn STELLA .swd. ND=169 (adaptive grid), sigma_min=-98.06,
38.7% of zones sigma<-1 - genuinely non-monotonic, the intended
DO_NONMONO_SOB test case. Jobs 11817076 (FINAL) / 11817077 (OBS) on
partition `small` (24h): the interactive partition has
MaxMemPerNode=32000 and the estimator asked 35G (ND=169), so the job
sat PENDING(MaxMemPerLimit) until moved. Debug submissions on
interactive must keep the memory request <= ~31G.

Pipeline lesson (non_mono_mode activation): `cmfgen.non_mono_mode: SOB`
alone does NOT activate the orchestrator's non-monotonic block in
run_model_grid.py - activation requires a `DO_NOMONO_*` or
`NON_MONO_SOLVER` key set T in `cmfgen.vadat_overrides`
(run_model_grid.py ~8126). With only non_mono_mode set, the block is
skipped, `GLOBAL_LINE` is stamped BLANK (SN_IIN default) into the
resolved overrides, `n_ob_ins` is not forced to 0, and the forced SOB
key table never lands. The config therefore carries
`DO_NOMONO_SOB: T` + `DO_NONMONO_SOB: T` explicitly. Resolved overrides
verified: GLOBAL_LINE=SOB, all TRANS_*=SOB, DO_SOB_LINES=T,
SOB_FREQ_IN_OBS=T, NON_MONO_SOLVER=T, DO_NOMONO_KNOP=F, USE_LAM_ES=T,
SOL_METH=DIAG, VEL_LAW=11, N_OB_INS=0, DO_LAM_IT=T, DO_LAM_AUTO=T,
FIX_BA=1.0, N_FIX_BA=10. In-job INIT writes VADAT; verify
`T [DO_NONMONO_SOB]` once the FINAL starts.

Other settings: `FIX_T=T` + `FIX_T_AUTO=F` (hold shocked hydro T),
`DIF=T` + `IB_METH=DIFFUSION` (inner-boundary tau ~28 at conversion;
prepare_VADAT auto-selects ZERO_FLUX if tau<50 in the generated file -
check VADAT when written), `VTURB=20`, `LSTAR=0`, `n_ob_ins: 0`
(VEL_LAW=11/RDINR grid match), `proxy_populations: true`.

First submission (jobs 11817076/11817077) aborted at it9: LAMBDA
corrections grew ~1e33 -> ~1e109 %/iteration from it1
(depths 1-4 dominant, `FG_J_CMF_V13 - SIGMA .LT. -1.0` warning at
depth 1, sigma=-1.005) -> "bad initial population guesses". A
`LAM_SCALE_OPT=LIMIT` resubmit was prepared but reverted per user
directive: the clip only masks the divergence - if the system diverges
the cause is elsewhere (J at the reversal depths, the boundary
sigma=-1.005 cell, or the initial populations) and must surface.
Diagnose the divergence, do not damp it.

Non-local Sobolev coupling (19-Sep-2026) - Rybicki & Hummer 1978:

New file subs/sob_rh_couple.f implements the RH78 (ApJ 219, 654)
non-local coupling between distant Sobolev resonance points, in the
geometry-free form of Feldmeier & Nikutta (2006, A&A 446, 661;
astro-ph/0510806).  Motivation: the Dessart floor is an ad-hoc patch;
the physically correct treatment gives photons trapped at sigma ~ -1
a real escape/reprocessing channel at remote resonance points.

Method: for each line l and depth r,

  dJ_l(r) = 0.5 * INT dr' r' v' GAMH_l(r') S_l(r') K(r,r')
  K(r,r') = EXF(t)*EXF(t') * 1/sqrt[(r^2-r'^2)(r^2v'^2-r'^2v^2)]
  t = GAMH_l(r)/|1+mu^2 sigma(r)|, mu^2 from RH Eq.12

applied as dZNET_l = -dJ_l*CHIL_l/(ETAL_l*BB_COR), matching the
SOBJBAR_SIM convention ZNET_l = 1 - Jbar/(S_l*BB_COR) (same
conversion as BETAC_SIM).  Two-point kernel only - 3-point
directions would need the extra exp(-t'') extinction (bounded
overestimate without it).  Kernel uses UNCLAMPED |1+mu^2 sigma|;
the Dessart floor is only for the local beta integral.

VADAT flag: DO_NONLOC_SOB (default F).  Requires DO_NONMONO_SOB=T
(uses SIGMA_TRUE for the resonance geometry); rd_control_variables
warns and disables it otherwise.  Lagged (lambda-iteration)
correction - no Jacobian derivatives added.

Call site: sub_sob_line_v3.f after SOBJBAR_SIM, gated on
DO_NONMONO_SOB .AND. DO_NONLOC_SOB.  Geometry (partner pairs, mu^2,
Upsilon1, weights) is line-independent and cached in SOB_RH_MOD;
rebuilt only when the R/V signature changes.

Verified standalone (test_sob_rh driver): v~r gives 0 pairs
(monotonic LOS velocity along every ray); v~1/r couples all pairs
with dJ>0/ZNET down (correct sign); plateau+shock-rise profile
couples 3376 pairs with the largest correction in the dense rising
region.  Bookkeeping dZNET = dJ*chi/eta/BB_COR verified exactly.

First run: TypeIIn_10d_DEV7_SOB.yaml + DO_NONLOC_SOB=T, job 11820632.
Watch for: "SOB_RH: N non-local coupling pairs" in OUTGEN, whether
the STEQ explosion at sigma~ -1.005 depths is damped, and TUNE cost.

Result (job 11820632): SOB_RH fired - 9964 coupling pairs over 148
depths - but the run still aborted at it11 (~1e107) with the SAME
signature as floor-only: var 193 (Fe-group population eqn) dominant
at all depths, corrections growing ~10 orders/iter.  OBS dep 11820633
auto-cancelled.  Conclusion: two-point non-local coupling is a modest
(few-%) correction to Jbar and cannot de-stiffen net rates that are
catastrophically wrong in the sigma ~ -1 region.  The STEQ runaway
is not rescued by RH78 coupling alone.  Next options: (a) quantify
dJ/Jbar in production (add diagnostic), (b) 3-point coupling +
kink handling per Feldmeier-Nikutta, (c) damped/limited correction,
(d) accept Knop as the non-mono path for this model class.

RH78 revision 2 (19-Sep-2026, job 11821347) - Jacobian-clean split:

Diagnostic run 11820830 revealed the first implementation's real
defect: folding dJ into ZNET_SIM injected the remote term into the
BA Jacobian via ZNET_SIM/POPS(NUP) and (under LAMBDA_ITERATION) via
VB_SIM=-JBAR/ETAL.  The lagged remote term has no local-population
derivative, so both channels were wrong.  Fixed: SOB_RH_COUPLE now
returns DRC_SIM = CHIL*dJ/BB_COR (rate units, eta*ZNET scale) in a
separate line_mod array; sub_sob_line_v3 applies it only in the STEQ
residual T2 = ETAL*ZNET - DRC_SIM.  ZNET_SIM and JBAR untouched.
Verified against FN06 Eq.10-16 (kernel, mu^2, tau, pair cuts all
match); standalone driver re-verified all three velocity cases.

Result: STILL DIVERGES - it11 abort ~1e112, same ~10 orders/iter
growth.  Diagnostics: max dJ/J grows 62 -> 1e13 across iterations
(tracks diverging remote S), max dRATE 4e6 -> 1.7e8.  Conclusion:
the Jacobian artifact was real but not the driver.  The Sobolev
line-rate map itself is anti-contractive here: remote absorption
DRC demands local populations that cannot balance it, so each
lambda iterate amplifies.  Adding physically-correct nonlocal
coupling strengthens, not damps, the channel because dJ follows
the exploding source functions.  Consistent with the Knop+BLANK
control converging on the same hydro - the escape-probability
line formulation is the diverger for this sigma structure.

Files: subs/sob_rh_couple.f (DRC interface), line_mod.f DRC_SIM,
sub_sob_line_v3.f zeroes DRC_SIM per call + T2 subtraction only.
DO_NONLOC_SOB remains default-F.  NOTE: sob_rh_couple.f was twice
found reverted on disk mid-session (stale editor buffer?); if
editing it, verify `grep "SUBROUTINE SOB_RH_COUPLE"` shows the
13-arg (DRC,ZNET_SIM,...) signature before building.

### KNOP_MOM stability fixes + TypeIIn rerun (19-Sep-2026)

Three fixes applied after auditing the TypeIIn_10d_DEV7_KNOP_MOM
divergence (job 11832268, L_outer ~10x/iter to 6e14 by it11):

1. get_moms_non_rel.f - degenerate-J guards on every J-denominator
   ratio.  K/J -> 1/3 (isotropic limit, same threshold and philosophy
   as the FEDD guard in fg_j_cmf_nonmono_v1.f), H/J and N/J -> 0, and
   NMID_ON_J -> 0 where RSQ_JNU is sub-tiny (1.0E-30_LDP).  Unguarded
   quotients at EUV frequencies (JNU_STORE~1e-40) produced garbage
   Eddington factors that QFROMF integrated into grid-wide degenerate
   QF (bounded at 723/724 depths) -> near-singular tridiagonal.
   Verified: rerun reports ZERO "QF_NM bounded" lines (was hundreds).
2. Under-relaxation on lagged GS fields (omega=0.5):
   fg_j_cmf_nonmono_v1.f INEXT_RELAX applied to I_P/I_M_DEPTH_NEXT
   writes; mom_j_cmf_nonmono_v1.f JNEXT_RELAX applied to the
   J_NEXT_FRQ update.  Form: OLD + w*(NEW-OLD); cells with OLD<=0
   fill directly (cold start unaffected).  Fixed point unchanged;
   damps the lagged-feedback gain.  Knop-path only - stock solvers
   untouched.
3. comp_j_blank.f - upstream string bug fixed at both post-solve
   inner-flux blocks: `'ZERO FLUX'` (space) -> `'ZERO_FLUX'`
   (underscore), matching VADAT/solver spelling.  Previously the
   zero-flux branch never matched and HFLUX_AT_IB always used the
   diffusion formula.  Verified: inner luminosity now exactly 0.
   (Bug also exists in cmfgen_reference - flag upstream.)

Infra fix (cmfgen_base/submit_cmfgen.sh, both FINAL+OBS templates):
check_outgen_for_fatal's exclusion list covered the FG GS diagnostic
line but not the MOM one - `MOM_J_CMF_NONMONO_V1 ... capped/NaN
writes= 0 NCF=` matched the `NaN` fatal pattern, so every MOM-mode
job FAILED instead of requeuing at graceful stop (job 11834123,
exit 3:0).  Added exclusion `MOM_J_CMF_NONMONO_V1.*capped/NaN
writes= *0 NCF=` (only the all-clear line; capped>0 still trips).

Result (job 11834123, resubmitted fresh - diverged it11 checkpoints
moved aside first): STILL DIVERGES.  L_outer 2.9e8 -> 1.8e8 (it1-3,
healthy decay) then 6.2e8 -> 3.9e9 -> 3.9e10 -> 5.5e11 -> 5.6e12 ->
7.5e13 -> 1.2e15 -> ... -> 1.5e18 at it15, tracking the unguarded
run almost exactly.  Negative-J resets localize at expanded depths
~192-196 (sigma_true ~ -37, chi ~ 1.7e-5 EUV).  Conclusion: QF
degeneracy was a symptom amplifier; the driver is the compressing-
branch lagged feedback itself - |PSIPREV|*J_NEXT(k+1) injected once
per great iteration with loop gain >1 at sigma<<-1 depths.  omega=0.5
halves the increment but cannot change the sign.
GS max-change stayed ~1e5-1e6 throughout - the lagged field never
converged and warm-up LAMBDA remained forced, so the growth is in
the transfer solve itself, not the BA step.

### PSI_NM compressing-weight fix - divergence cured (19-Sep-2026)

Root cause pinned down by comparing with the formal solver's
compressing source: FG uses SOURCE_PRIME = S + |Q_k|*I_NEXT where the
explicit upwind weight is the SAME local |Q_k| that sits on the
implicit side.  The moment solver instead used
|PSIPREV| = DTAUONQ*|GAM*(1+sigma*K_ON_J_PREV)| - the explicit J_{k+1}
coupling weighted by the k-1 Eddington factor.  The mismatched weight
can exceed the diagonal's |Q_k|, so the lagged fixed-point map
J_next = F(J_next) is EXPANSIVE at strong-compression cells and J
grows ~10x per great iteration.

Fix: at every compressing injection site use PSI_NM
(= DTAUONQ*|Q_k|, already computed for the diagonal) as the J_{k+1}
coefficient instead of ABS(PSIPREV).  Equal explicit/implicit weights
make the RHS a convex combination - J_k <= max(S, J_{k+1}, neighbours)
- so the lagged map is non-expansive and cannot amplify.  Sites:
interior G_ONLY + else branches, both outer-boundary branches, and
the negative-solution rebuild path (5 places; expanding branch still
uses signed PSIPREV, unchanged).  This also mirrors the V11
convention's spirit (velocity factors at k, Eddington at the coupled
frequency) while giving the FG-consistent local-|Q| choice.

Result (job 11837622, fresh start): STABLE + CONTRACTING.
L_outer 2.86e8 -> 1.5e8 over its 1-5 and pinned at ~1.5-1.6e8 through
it12+ (prior runs were at 3e9..1e18 by then).  MAXCH bounded ~300-1100
at the shock-ramp depths - the same limit-cycle regime as the
USE_LAM_ES=T formal run.  GS max-change decaying (1e5 -> 5e4);
warm-up LAMBDA still forced while the GS field matures, so the BA
phase with the diagonal-ALO linearization is yet to be exercised on
this model.  QF_NM bounded count remains 0; inner L = 0 (ZERO_FLUX).

Also fixed earlier today (same session): get_moms_non_rel degenerate-J
guards (K/J->1/3, H/J,N/J->0, NMID_ON_J->0 at JNU<1e-30);
INEXT_RELAX/JNEXT_RELAX=0.5 under-relaxation on the lagged writes;
comp_j_blank 'ZERO FLUX'->'ZERO_FLUX' upstream typo; and the
submit_cmfgen.sh check_outgen_for_fatal exclusion for the
MOM_J_CMF_NONMONO_V1 'capped/NaN writes= 0 NCF=' line which had been
failing every MOM-mode requeue.

### In-run GS-immaturity gate ported from dev6 (19-Sep-2026, late)

Convergence audit of the TypeIIn KNOP_MOM run (11837622) vs Knop et al.
2009 Papers I+II (docs/Knop_et_al_2009_formal_solution.pdf,
docs/aa11793-09.pdf):

- Observed limit cycle: DO_LAM_AUTO releases RD_LAMBDA when MAXCH<50,
  then NUM_LAM=2 alternates 2 LAMBDA + 1 BA.  Every BA step produced
  1e4-1e8% corrections (it87 1.5e5, it102 6.4e5, it119 6.9e7) against
  a GS lagged field still moving by rel. change ~1e3-3e4 per sweep
  (~10-24M cells >50%/sweep).  MAXCH>1e5 then forces LAMBDA again -
  permanent oscillation, never converges.
- dev6 had the matching mechanism (dropped in the dev7 clean restart):
  an in-run gate in SOLVE_FOR_POPS_NONMONO_V1 that forced LAMBDA
  whenever the last GS sweep was still moving >50%, recorded cap
  events, or had significant negative-J resets (GS_UNSAFE_TOL=0.5).
- Ported minimally into solve_for_pops.f (DO_NONMONO_KNOP-gated only):
  GET_FG_GS_MATURITY now also returns the live FG_GS_MAXCHG and
  FG_INEXT_CAPPED (5 args); the MOM getter already returns
  MAXCHG/CAPPED/NEGJ.  When a BA iteration is pending and the GS
  field is immature the iteration is forced to LAMBDA (DIAG,
  MAX_LAM_COR) and FORCE_LAMBDA keeps the mode controller from
  re-requesting BA while the field stays immature.
- Also ported two dev6 audit fixes: NEWMOD resets the SAVE'd scalars
  MAXCH_SUM/COUNT_FULL_LAM_BA, and R_OLD is initialised before the
  REVISE_R_GRID branch (DO_NG_BAND_ACCEL_V4 reads it unconditionally
  when NG fires with REV_RGRID=F - previously uninitialised memory).

Paper-vs-dev7 structural difference (documented, not a bug): the paper
converges the wavelength boundary-value problem to ~1e-8 INSIDE each
formal solve for a fixed source function, then iterates the ALI with
the exact wavelength-tridiagonal Lambda* operator.  dev7 merges the GS
frequency coupling into the SE great-iteration loop (one sweep per
population update) and uses only the diagonal ALO.  The in-run gate
keeps the population solver in the contraction regime (LAMBDA) until
the merged loop matures; whether pure-LAMBDA reaches <0.1% on this
sigma~-70 model is the open question the restarted run must answer.

Build: both files recompiled clean, strings in exe/cmfgen_dev.exe
verified (23:04 build).  Job 11837622 requeued at 23:01:57 just before
the relink, so it picks up the gated binary at its next segment
boundary (~2h); no action needed.

### Frozen-population inner sweeps - Knop Paper-II nested iteration (19-Sep-2026, late)

Second surgical piece closing the Paper-II convergence setup.  The
paper converges the wavelength boundary-value problem internally for
a fixed source function before each ALI/source update; dev7 instead
merged the GS coupling into the SE loop (one sweep per population
update).  Rather than wrapping the fragile ML frequency loop (line-
storage bookkeeping, STEQ/BA accumulators are not re-entrant), the
nested structure is recovered at the great-iteration level:

- In SOLVE_FOR_POPS, when the just-completed sweep's lagged field is
  immature (FG or MOM max rel. change > GS_MATURE_TOL=5%, any capped
  write, or MOM significant-negative-J resets), SOLVEBA_V13 is not
  called - POPS are untouched, so the next great iteration is a pure
  transfer sweep that converges the coupling for FIXED populations.
  The SE update fires automatically on the first mature sweep, which
  is exactly the paper's ordering (converged BVP -> source update).
- Frozen iterations report MAXCH = 100*max(FG,MOM GS change) so
  CONVERGENCE_HISTORY tracks the inner iteration; they stay in
  LAMBDA mode via FORCE_LAMBDA and cannot fire NG/AV/smoothing/
  grid-revision (all gated on .NOT. GS_FREEZE).
- N_FROZEN_MAX=15 caps the freeze: if the lagged map cannot contract
  in 15 sweeps the update proceeds anyway (degrades to the previous
  merged-iteration behaviour, warning printed).
- Checkpoints continue to be written every frozen iteration
  (SCR_RITE_V2 -> POINT1/2 + FGINEXT_GS/MOMJNXT_GS), so the freeze
  state is requeue-safe; the cap counter N_FROZEN resets on NEWMOD.

Evaluation run: TypeIIn_10d_DEV7_KNOP_MOM resubmitted via
run_model_grid.py on the debug partition (FINAL 11853189 / OBS
11853190), cold start.  Watch for 'SOLVE_FOR_POPS: KNOP inner-GS
sweep' + 'frozen iter' lines and whether FG/MOM GS max-change
contracts across frozen sweeps; the first real SE update should
follow the first mature (<5%) sweep.

### Frozen-sweep cap removed (20-Sep-2026)

N_FROZEN_MAX=15 deleted from SOLVE_FOR_POPS; the freeze now persists
until the lagged field is genuinely mature.  Rationale: the first
FIX_T=F run (TypeIIn_10d_DEV7_KNOP_MOM_TFIT, job 11869230) reached
the cap at it16 with the field still moving ~1e6-1e9; the forced SE
update kicked FG/MOM max-change to ~1e16-1e17 and the next formal
solve died in comp_j_blank ("Mean intensity blowing up").  The cap
reintroduced the BA-on-immature-field failure on a timer.  N_FROZEN
is retained purely as the 'frozen iter N' diagnostic counter.
Consequence: a non-contracting GS map now stalls (frozen forever)
rather than degrading - watch consecutive segments for that failure
mode before considering a maturity-metric change (e.g. gating on the
cells>0.5 bulk fraction rather than single-cell max).

### Per-ray/per-frequency Knop refresh + asymmetric spacing (20-Sep-2026, late)

Structural fix for the one-bin-per-great-iteration lag that kept the
GS field immature forever (TypeIIn: FG maxchg ~3e4, MOM ~1e6, frozen
sweeps indefinitely).  Instead of waiting NCF great iterations for
red-to-blue information to cross a compressing zone, each sweep now
captures its transfer inputs and an internal alternating-direction
re-solve converges the lagged coupling before the populations are
updated.  Paper-II ordering is preserved: converged BVP field, then
SE update via the existing maturity gate (unchanged).

FG path (fg_j_cmf_nonmono_v1.f):
- Per-frequency capture during the ML sweep: KN_CHI_EXT/KN_ETA_EXT
  (single precision, extended depth grid) plus scalar/grid data.
- FG_KNOP_REFRESH_V1 (called from cmfgen_sub after 10000 CONTINUE)
  re-solves each ray's (space,frequency) intensity field in
  alternating direction passes on the captured inputs, using the same
  short-characteristics discretisation (DTAU, E-functions, A0-A4,
  monotonic slopes, boundary conditions).  Writes the converged field
  back into I_P/I_M_DEPTH_NEXT and updates the GS maturity metrics.
  Guards: DO_KNOP_REFRESH, HAS_COMPRESSING, complete sweep, INTEGRAL
  method only (TRIDIAG path is not refreshed - known limitation).
- VADAT switch DO_KNOP_REFRESH (default DO_NONMONO_KNOP) via
  control_variable_mod / rd_control_variables / SET_FG_KNOP_REFRESH.
- SET_FG_NU_NEXT(NU_NEXT_FREQ) passes nu(k+1) per call so compressing
  cells use the forward spacing ln(nu_k/nu_{k+1}) instead of the
  backward ln(nu_{k-1}/nu_k); falls back to dLOG_NU when the next
  frequency is unavailable.  MOM solver got the same fix via
  SET_MOM_NU_NEXT and per-depth spacing in mom_j_cmf_nonmono_v1.f
  (midpoint sign held in GAMH_SIGNED_NM).

MOM path (mom_j_cmf_nonmono_v1.f, USE_LAM_ES=F):
- KN_M_* capture store: per-frequency ETA/CHI/ESEC (SM grid), the
  post-GET_MOMS_NON_REL Eddington factors K_ON_J/NMID_ON_J/
  NMID_ON_HMID, boundary data (HBC/NBC/IN_HBC, JPLUS_IB/JMIN_IB,
  HNU_AT_IB), DBB/IC, dLOG_NU, nu(k+1); sweep-fixed args captured once.
- MOM_KNOP_REFRESH_V1 (cmfgen_sub, gated on DO_KNOP_REFRESH) re-drives
  MOM_J_CMF_NONMONO_V1 itself in alternating frequency order, so the
  tridiagonal assembly/THOMAS solve/safeguards are identical by
  construction.  MOM_KN_REFRESHING suppresses GET_MOMS_NON_REL (whose
  EDDFACTOR inputs exist only for the last sweep frequency) and the
  capture writes.  The lagged bluer field is seeded from
  J_NEXT_FRQ(:,L-1); lagged H from a scratch store.  Pass metrics are
  written back into MOM_GS_MAXCHG/MOM_INEXT_CAPPED/MOM_NEGJ_LRG so the
  frozen-iteration gate sees the refreshed field's self-consistency.
- Self-guards: no-op unless a complete KN_M store exists, some depth
  was compressing this sweep, and J_NEXT_VALID.

Pipeline (cmfgen_base/submit_cmfgen.sh):
- On requeue (file_requeue_count>0) with both EDDFACTOR and
  EDDFACTOR_INFO present, the leading COMP_F token in VADAT is flipped
  T->F after prepare_VADAT.py runs, so restarted segments reuse the
  persisted Eddington factors instead of re-Lambda-iterating them
  (previously each requeue reset the scattering J to zero at the
  first frequency).

Build: full `make cmfgen` clean; exe/cmfgen_dev.exe relinked.
Validation: TypeIIn_10d_DEV7_KNOP_MOM_TFIT (interactive, job
11877980 + requeue 11877981) was already running on the pre-refresh
binary; its next segment loads the new build.  Watch for
'FG_KNOP_REFRESH pass' / 'MOM_KNOP_REFRESH pass' lines and whether
FG/MOM GS max-change drops below the 5% maturity gate within the
refresh passes, then whether the first real SE update fires.

FREQ_IDX_ITER off-by-one (20-Sep-2026, later): during frequency k of
the sweep, FREQ_IDX_ITER is k-1 - the INIT branch (FIRST_FREQ=T at
ML=1) resets it to 0 and skips the NEW_FREQ increment.  So the lagged
I_P/I_M_DEPTH_NEXT store actually holds frequency s+1 in slot s
(slot NCF is never written = red boundary zeros; frequency 1 is
never stored).  The first refresh implementation assumed slot s =
freq s, so the capture armed at LK=FREQ_IDX_ITER never reached
KN_NCF, KN_STORE_READY stayed FALSE, and FG_KNOP_REFRESH_V1 returned
silently before its warning checks.  Fixed by: capture at
LK=FREQ_IDX_ITER+1 (KN_* slot = true frequency), SEED_FLD reading
STORE(:,:,l-1), WRITE_STORE writing STORE(:,:,l-1), and the driver
completeness check comparing FREQ_IDX_ITER to KN_NCF-1.  The
FGINEXT_GS checkpoint layout is unchanged (it was always shifted and
is self-consistent).  The MOM path uses FREQ_INDX (=ML) directly and
is unaffected.  Separately fixed the same day: COMP_JCONT_V4.INC
(the OBS-path call site) now calls SET_FG_NU_NEXT(NU(FREQ_INDX+1))
so OBS frequencies do not inherit the stale last-sweep nu(k+1).

### Refresh residual metric + convergence behaviour (20-Sep-2026, late night)

First live refresh runs (jobs 11880567/11881301/11881901, MOM_TFIT
11877980 requeue) exposed three issues the paper does not cover:

1. Faint-cell dominated residual.  |dV|/max(V,1e-25) gives res ~1e21
   when a cell leaves the numerical floor.  Residual is now normalized
   against the per-(ray,frequency) column peak (VMAX_L), so only cells
   that carry a significant fraction of the local intensity count.
   MOM refresh metric got the same treatment (KN_M store column peak).

2. Limit cycle / marginal contraction.  The undamped alternating GS
   converged fast (2e8 -> 0.66 in 6 passes) then either entered a
   period-4 cycle at ~14% residual (cold start) or crawled at ~1%/pass
   (warm restart, lambda ~ 0.99 eigenmode).  Refresh updates are now
   under/over-relaxed per ray: omega adapts between ~0.5 (residual
   growing/oscillating) and ~1.3 (smooth contraction); the fixed point
   is unchanged.  Pass history is printed (FG_KNOP_REFRESH: pass
   residuals).

3. Gate frozen by outliers.  GS_MATURE_TOL on the max-norm meant a
   handful of marginal cells (|Q|->1 at strong compressions) could
   keep the model in transfer-only mode forever.  FG_GS_NBAD/FG_GS_NTOT
   now count significant cells (>1e-4 of column peak) still above
   tolerance; solve_for_pops treats the field as immature only when
   the bad fraction exceeds ~0.1% (or a hard runaway threshold, or
   cap/NaN events).  Counters reset each sweep so a skipped refresh
   cannot open the gate.

Also: FREEZE_I_NEXT (JFNK freeze) now blocks the FG refresh writeback,
matching the in-sweep write guard.

Known open item: MOM_KNOP_REFRESH_V1 printed no output on the
MOM_TFIT job despite the MOM nonmono solver running - a silent guard
return.  A one-shot guard-state diagnostic line
('MOM_KNOP_REFRESH guard: ...') was added at driver entry to identify
which condition fails on the next segment.

### MOM refresh guard + SOLVEBA-on-immature-field fixes (20-Sep-2026, late night)

Guard diagnostic on job 11885512 showed WHY MOM_KNOP_REFRESH_V1 never
ran: KN_M_ANYCMP stayed FALSE all sweep.  HAS_COMPRESSING was only
tested on the midpoint factor (1 + AV_SIGMA_TRUE*NMID_ON_HMID), but
NMID_ON_HMID is floored at 0.01 in get_moms_rel.f, so the midpoint
factor stays ~0.02 > 0 even at sigma_true ~ -98.  The actual
compressing branch in the tridiagonal keys off the DEPTH-point
Q_NM(I)<0 (1 + SIGMA_TRUE*K_ON_J, K_ON_J~0.3), which does go negative
-- the coupling was active but invisible to the flag.  Fixed:
HAS_COMPRESSING now also set on depth-point GAM_SIGNED_NM<0.

Restored original-CMFGEN iteration semantics in solve_for_pops.f:
GS_FREEZE no longer skips SOLVEBA_V13.  Immature-field iterations now
run a forced-LAMBDA update (LAMBDA_ITERATION=T, FIXED_T=T,
COMPUTE_BA=T, TEMP_CHAR='DIAG', T1=MAX_LAM_COR), matching V6's
solve_for_pops_nonmono_v1 which always ran SOLVEBA_NONMONO_V1.  All
standard iteration diagnostics (LAMBDA iteration used, max %
increase/decrease, DEC_VEC/INC_VEC, MINSCALE) print every iteration.
The .NOT. GS_FREEZE gates on AUTO_SMOOTH_POPS / REVISE_R_GRID /
UNDO_IT / AVE_FLIPS remain -- conservative while the field matures.
N_FROZEN renamed in print to 'immature-field iter' (counter kept).

FG refresh adaptive relaxation revised after job 11881901 showed
omega=0.6 under-relaxation HURT (sawtooth envelopes converging ~10%
per direction pair vs undamped reaching 0.13): OMEGA_RAY now starts
at 1.0 and adapts on SAME-DIRECTION residual history (pass N vs N-2,
since alternating sweeps interleave expanding/compressing passes with
intrinsically different residual levels).  Same-direction growth ->
damp toward 0.45; slow same-direction contraction (ratio>0.85) ->
over-relax toward 1.45.  Also added bulk early-exit after pass>=6
when NBAD <= max(50,NTOT/1000), matching the maturity-gate level.

Both running jobs (11881901 timed out mid-it2 at the wall; MOM_TFIT
11877980 timed out at wall) resubmitted warm: FG=11885796,
MOM_TFIT=11886410 (cancelled 11885512/13 mid-sweep to pick up the
HAS_COMPRESSING fix immediately rather than waiting ~1.7h for the
requeue boundary).

### Single-precision capture-store NaN bug + int32 counter overflow (21-Sep-2026)

MOM_KNOP_REFRESH_V1 ran on job 11886410 (CMP=T after the
HAS_COMPRESSING fix) but every re-solve returned NaN rows
('non-finite XM_NM zeroed', SOURCE_NM=NaN at NM point 2) while the
in-sweep solve at the same frequencies was clean.  Root cause:
KN_M_ETASM/KN_M_CHISM/KN_M_ESECSM (and the K_ON_J/NMID stores) were
default REAL = single precision.  Perfectly valid double ETA values
below ~1.2E-38 (cold faint outer cells) underflow to exact 0.0 on
capture; the refresh's LOG(ETA_SM) interp then produces -Inf -> NaN
coefficients -> NaN source rows.  Fix: all six KN_M arrays promoted
to REAL(KIND=LDP); memory cost ~15 MB, irrelevant next to the 8 GB
FG store.  A one-shot store-integrity scan
('MOM_KNOP_REFRESH store scan: bad ETA=.. bad CHI=..') now runs at
first refresh entry to confirm zero non-positive entries.

FG side had the same single-precision store (KN_CHI_EXT/KN_ETA_EXT)
but its LOG(MAX(..,1E-300)) floor prevents NaN; the residual
corruption (values <1e-38 replayed as 1e-300) is physically
negligible at cells emitting ~1e-40, so the single-precision store
stays for its documented ~300 MB memory saving.

FG_GS_NTOT printed negative ('unconverged significant cells=
82210160 of -22911603'): int32 overflow.  With NCF=43475 and ~181
rays the last-pass significant-cell count is ~5e10.  Promoted
FG_GS_NBAD/FG_GS_NTOT/FG_GS_NCONV, the NBD_LS/NTT_LS ray arrays, the
worker NBAD/NTOT/NB_*/NT_*/NBAD_CUR/NTOT_CUR counters, the
GET_FG_GS_MATURITY NBAD/NTOT dummies, and the solve_for_pops
GS_NBAD/GS_NTOT receivers to INTEGER*8; print formats widened to
I12/I14 and MAX(50_8,...) used in the bulk comparisons.

Rebuilt cmfgen_dev.exe at 22:59; MOM_TFIT resubmitted warm as
11887664 (OBS 11887665) to exercise the LDP store + scan.

### MOM refresh GS cascade -> Jacobi snapshot (20-Sep-2026, 23:30)

Job 11887664 (LDP-store binary) ran the refresh cleanly -- no NaN
storm (capped counts ~3k/pass vs 724/call) -- but diverged: descending
passes spiked to 1E5-1E10 residuals with ~9E5 significant negative-J
resets per pass, while ascending passes sat at ~1E3.  Cause: on a
descending pass the routine's compressing read J_NEXT_FRQ(:,L+1) sees
a THIS-PASS write -- a true GS cascade along the coupling direction.
The in-sweep read is pure Jacobi (previous great iteration's field),
so the refresh's descending pass was strictly more aggressive than
the production semantics it is meant to converge.

Fix: JNEXT_RF_READ pass-start snapshot (module array).  While
MOM_KN_REFRESHING is set, all J_NEXT_FRQ(:,FREQ_INDX+1) read sites
(2 interior + 2 boundary + 1 negative-J retry) take JUP_NM from the
snapshot; writes still go through JNEXT_RELAX=0.5 under-relaxation.
This makes the refresh's compressing coupling identical to the
in-sweep lagged coupling.  The expanding side (RSQJNU_PREV chain)
is untouched.

MOM counters promoted to INTEGER*8 (MOM_INEXT_CAPPED, MOM_NEGJ_LRG,
MOM_GS_CAP_FILE, getter dummies, PASS_CAP/NEG, solve_for_pops
receivers): the NaN-storm path could overflow int32 (~3E10 events).
MOMJNXT_GS checkpoint VER bumped 2->3; older files are cleanly
rejected and J_NEXT rebuilt via Lambda warm-up.

Binary relinked 23:27; MOM_TFIT resubmitted warm as 11888720.

### MOM refresh containment + FG bad-cell metric fix (21-Sep-2026, 00:10)

Job 11888720 (Jacobi-snapshot binary) showed the snapshot works:
refresh contracted monotonically for 22 passes (8.1E4 -> ~44), then
an expansive mode blew the tail up to ~1E13.  negJ stayed ~9.7E5
cells/pass throughout -- about 3% of the ND x NCF field sign-flips
on EVERY pass, including in-sweep solves (the 200+ MB
MOM_J_NONMONO_ERRORS dump is mostly in-sweep negative-J rows at
sigma_true ~ -37 to -55 compressions).  The lagged fixed-point map
of the moment solve is genuinely non-contractive on that subspace;
it is a discretization-level instability, not a plumbing bug.

Containment added (mom_j_cmf_nonmono_v1.f):
- XM_NM magnitude cap at 1E40 alongside the non-finite zeroing, so a
  runaway-but-finite value cannot overflow to Inf inside the next
  solve's PSI*JUP/EPS*RSQJNU coupling terms.
- Refresh keeps JNEXT_BEST = field at the lowest-residual pass; on
  exit without convergence J_NEXT_FRQ is restored to it, and
  MOM_GS_MAXCHG reports BEST_CHG (the residual of the field actually
  stored).
- Abort once PASS_CHG > 1E3*BEST_CHG (post pass 3) -- the explosive
  tail has no chance of recovering and only wastes wall time.

FG side (fg_j_cmf_nonmono_v1.f): FG_GS_NBAD counted self-normalized
>5% changes on a ray's last pass.  Rays with no compressing cells
exit after pass 1, where the "change" is seed->exact-solve -- a
transient, not an unconverged state.  CLEAN_EXIT flag now zeroes a
ray's NB_O on RES<TOL or !HAS_CP exits, so FG_GS_NBAD reflects only
cells that are still unconverged.  (The 1.27% bad reported on job
11885796 was dominated by this false positive; the refresh itself
converged 8.4E6 -> 1.5E-2 in 9 passes.)

Binary relinked 00:08; MOM_TFIT resubmitted warm as 11889372.
FG job 11885796 picks this up at its next requeue.

### FG metric fix verified + gate opens (21-Sep-2026, 01:45)

Job 11885796 requeued onto the 00:32 binary and the CLEAN_EXIT/RES
fix confirmed: FG refresh max residual dropped from ~1.8E6 (clean-ray
seed-transient false positive) to 1.08E-1, and unconverged
significant cells from 1.1E8 (1.3%) to 1.98E6 (0.023%).  Pass
residuals converge 1.8E6 -> 0.053 in 10 passes under adaptive omega.
With NBAD/NTOT < 0.1% the bulk gate opens; population updates have
collapsed to physical scale (max Lambda change 4.5E+1 vs 9.5E+94
earlier).

MOM job 11889372: refresh containment stable across iterations --
contracts ~1-2 dex for 9-15 passes, expansive mode explodes, abort
fires, best-pass field restored (~1E4-5E2 residual), no NaN
propagation, ~60-80 s/call vs 350 s before.  The 1E40 XM cap is now
catching ~2.5E5 runaway cells/pass that previously became Inf.
The moment lagged-map instability is intrinsic to this model's
compression strengths (same negative-J rows appear in-sweep);
the MOM path stays in forced-Lambda (V6 gate semantics), which is
the correct conservative state.

### FIX_IMPURITY cleared on KNOP forced-LAMBDA sweeps (22-Sep-2026, 20:40)

Audit of FIX_IMP integration on the non-monotonic path (requested for
the TypeIIn_KNOP_TRIDIAG USE_LAM_ES=F run):

- Wiring verified end-to-end: VADAT FIX_IMP -> RD_FIX_IMP ->
  FIX_IMPURITY (BA mode only, via cmfgen_sub init or the mode
  controller) -> MOD_FIX_IMPURITY (solve_for_pops.f) ->
  FIXPOP_IN_BA_V3 inside GENERATE_FULL_MATRIX_V3 inside
  CMF_BLKBAND_V3 inside SOLVEBA_V13.  Semantics: freeze ions/species
  contributing <1E-15 of the population during the BA update.
  new_ba_up_tst.f and subs/fixpop.f are dead code on this path.
- Gap found: the three KNOP forced-LAMBDA gates (warm-up,
  immature-field, GS_FREEZE) do not clear FIX_IMPURITY, and
  MOD_FIX_IMPURITY is snapshotted before the gates run.  A model
  that reached BA (FIX_IMPURITY=T) and then hit an immature sweep
  would run its forced-LAMBDA relax with trace ions frozen --
  breaking the upstream "FIX_IMPURITY only in non-LAMBDA mode"
  invariant.  Same omission exists in the V6 reference
  (solve_for_pops_nonmono_v1.f).
- Fix (solve_for_pops.f, before CALL SOLVEBA_V13): single guard
  `IF(FORCE_LAMBDA) FIX_IMPURITY=.FALSE.; MOD_FIX_IMPURITY=.FALSE.`
  covering all three gates.  Default path unchanged: the reset only
  fires on forced-LAMBDA iterations where upstream semantics already
  require FIX_IMPURITY=F.
- Config note: with DO_LAM_IT=T + DO_LAM_AUTO=T the option is inert
  in practice -- the MAXCH<50 auto-switch does not restore
  FIX_IMPURITY (upstream soft-off near convergence); it can only
  engage via the CNT_LAM overflow path after the auto-switch.

make cmfgen clean; solve_for_pops.o rebuilt, exe/cmfgen_dev.exe
relinked 20:36.

### Iteration-cycle audit: MOD_FIXED_T + NG guard (23-Sep-2026, 08:45)

Audit of the KNOP convergence cycle in solve_for_pops.f against the
V6 reference (solve_for_pops_nonmono_v1.f) and stock upstream:

Fixed (same snapshot-ordering class as FIX_IMPURITY):
- MOD_FIXED_T/MOD_FIX_T_D_ST/MOD_FIX_T_D_END were computed from the
  incoming FIXED_T/RD_FIX_T before the KNOP gates force FIXED_T=T.
  A forced-LAMBDA iteration entered with FIXED_T=F would leave the
  R.E. row unprotected -- degenerate T equation, singular BA, hard
  STOP in SOLVEBA_V13 (or a garbage dT).  Masked so far by FIX_T=T.
  The FORCE_LAMBDA block now also forces MOD_FIXED_T=T with the
  full-depth range (identical to what the next call computes).
- NG acceleration guard had .NOT. GS_FREEZE but not
  .NOT. FORCE_LAMBDA: NEXT_NG is scheduled off LAST_NG, so a queued
  NG kick could land on a warm-up iteration with a mature restored
  field (GS_FREEZE=F) and extrapolate pops across the I_NEXT reinit.
  AV/UNDO are already protected implicitly via LAST_LAMBDA.

Checked, clean: ITS_SINCE_START/N_FROZEN/MAXCH_SUM/COUNT_FULL_LAM_BA
SAVE'd and NEWMOD-reset; LAST_LAMBDA updates on forced iters;
immature gate + freeze gate coverage has no gap; FIXPOP_IN_BA_V3
DIAG-band mechanics correct; LAPACK IFAIL tested in BLKBAND_V3;
IN_ITS re-read restores FIX_IMPURITY/FIXED_T on the DO_LAM_AUTO flip.

Noted divergence from V6 (not changed): a singular BA is a hard STOP
(SOLVEBA_V13) -- V6's per-depth singularity recovery and
N_BA_FAILURES Lambda-fallback were not ported.  SUCCESS is dead
post-call.  STOP surfaces the failure, which fits the
surface-divergence policy; a fallback would trade that for a silent
Lambda iteration on a bad solve.

Cosmetic: the "LAMBDA iteration used"/"BA matrix computed" headers
print before the KNOP gates, so a frozen BA iteration logs the BA
header without the LAMBDA line (the MAXCH tag still says LAMBDA).

make cmfgen clean; solve_for_pops.o rebuilt, exe/cmfgen_dev.exe
relinked 08:42.

### NG acceleration audit (23-Sep-2026, ~09:00)

Full audit of DO_NG_BAND_ACCEL_V4 + NG_MIT_BAND_OPT_V1 +
DO_TRANS_AND_NG_V1 + NGACCEL, and the call-site in solve_for_pops.f.

Verified correct:
- Iterate ordering: RDPOPS slot 1 = newest record (IREC); NGACCEL's
  D0/D1/D2 differences and the A,B least-squares solve are internally
  consistent with that ordering.  Singular-determinant fallback sets
  one coefficient to zero.
- Per-depth acceptance: LOCINC>10.1x or LOCDEC<0.09x or reversed T
  trend flags a bad depth; >3 bad depths cancels the whole NG and
  restores the newest iterate.  Caller commits POPS/R/V/SIGMA only on
  NG_DONE -- a failed NG leaves the state intact.
- Post-NG control: >50% inc AND dec forces a LAMBDA iteration with
  FIXED_T=T and FIX_IMPURITY=F -- the BA-only impurity invariant holds.
- SIGMA_TRUE is recomputed from the accelerated V after NG success
  (DO_NONMONO_SOB/DO_NONMONO_KNOP branch).
- R/V/SIGMA extrapolation is a no-op on this model: REV_RGRID=F and
  DO_HYDRO=F give WRITE_RVSIG=F, so all 4 history slots carry the
  current R/V/SIGMA and the commit writes back identical values.
  Hydro provenance is preserved.
- Scheduling: NG needs MAXCH < VAL_DO_NG (BEG_NG=10%), >= IBEG_NG=15
  iterations, >=4 iterations since the last NG (ITS_PER_NG=8 spacing),
  re-suspends (NEXT_NG=1500) if MAXCH climbs again, retries +8 after
  failure.  Call-site also excludes LST_ITERATION, GS_FREEZE and now
  FORCE_LAMBDA.

Fixed (upstream bugs in the CHK_NG validity check, dormant under
FIX_T=T because the block only runs when accelerated T differs):
- MEAN was never initialised before MEAN=MEAN+RATIO accumulation --
  the global T-consistency verdict was garbage whenever it ran.
- RATIO=(BIG_POPS-POPS)/T1 divided by T1=RDPOPS(NT,I,1)-RDPOPS(NT,I,2)
  with no guard: stationary T at any depth gives Inf/NaN contaminating
  MEAN/MIN_RATIO/MAX_RATIO.  Depths with T1=0 are now skipped and the
  mean is taken over the counted subset.

Findings (design-level, not patched):
- NG fits a single (A,B) pair over a 60-depth x all-levels band
  assuming one linearly-convergent map.  Under KNOP the per-iteration
  map includes the MOM/FG I_NEXT refresh; the 4-record history can
  span a LAMBDA->BA switch, and on a continuing model
  (NEXT_NG=2000 path) NG can fire on the first iteration of a
  restarted process.  The GS tag mechanism keeps I_NEXT consistent
  across restarts, so this is a narrow edge, not a violation.
- Per-depth bounds (10x/11x) are loose: the observed +87%/-413%
  kicks at iters 57/68 of the ES run passed legally.  With FIX_T=T
  the T-sign and CHK_NG checks are both dormant, so the magnitude
  bound is the only pops guard.  This is why DO_NG=F remains the
  right setting for this model: the kicks are "working as designed",
  not a bug -- the design just tolerates 5x population perturbations.

make cmfgen clean; do_ng_accel_band_v4.o rebuilt, exe/cmfgen_dev.exe
relinked 08:56.

### Full iteration-cycle audit (23-Sep-2026, ~09:30)

Traced the whole great-iteration in cmfgen_sub.f (loop label 20000)
plus solve_for_pops.f, solveba_v13.f, cmf_blkband_v3.f,
generate_full_matrix_v3.f, fixpop_in_ba_v3.f, fiddle_pop_corrections_v2.f,
specify_it_cycle_v3.f, comp_j_blank.f.

Cycle per iteration (verified):
  opacities -> per-freq inner f-loop [FG_J_CMF_NONMONO_V1 (rays->f,Jcmf)
  <-> MOM_J_CMF_NONMONO_V1 (J), or USE_LAM_ES=T: RJ=TC direct from rays]
  until |d f|<ACC_EDD_FAC or 20 iters -> EVALSE_QWVJ_V8 etc. accumulate
  SE residuals + BA Jacobian on the converged RJ -> STEQ_T (RE), STEQ_ED
  (charge cons.) -> FG/MOM_KNOP_REFRESH (mature lagged I_NEXT on this
  sweep's opacities) -> two-photon/chg-exch/advection/adiabatic terms ->
  SOLVE_FOR_POPS [KNOP gates -> FORCE_LAMBDA fixups -> SOLVEBA_V13 ->
  CMF_BLKBAND_V3 -> GENERATE_FULL_MATRIX_V3+FIXPOP -> DGEEQU/DGETRF/
  DGETRS -> POPS *= (1-STEQ*SCALE) MAJOR per-depth limiter] -> POPTOION,
  SUP_TO_FULL, hydro/clump (+SIGMA_TRUE refresh) -> COMP_J_CONV ->
  controller (DO_T_AUTO, MAXCH<EPS, DO_LAM_AUTO, IT_SPECIFIER, IN_ITS).

Maths verified correct:
- Newton step: A.x=F(residual), x converted absolute->fractional,
  multiplicative update POPS*(1-x*SCALE), positivity guaranteed by
  BIG_LIM<1. MAJOR limiter: scale set only by levels >1e-10*Ne,
  T capped at MAX_dT_COR=0.2 and floored at T_MIN=0.3.
- SE block: level SE + ionization eqs + species number conservation
  (EQ_SPECIES rows) + charge conservation (NT-1) + RE (NT). Correct.
- Fixed-point consistency: converged state satisfies SE given J(n),
  J satisfies transfer given n. SE rates use the inner-loop-converged
  RJ, not a stale field.
- KNOP maturity gates measure the self-consistency error of the field
  the SE step actually used (post-sweep refresh residual). Ordering
  sound: gate -> forced-Lambda or BA -> SOLVEBA always runs.
- Convergence (MAXCH<EPS=1%) requires LAST_LAMBDA != MAIN_COUNTER --
  cannot trigger during warm-up/freeze. NUM_ITS>=1 floor, +2 after NG.
- SIGMA_TRUE recomputed after every R/V change (grid revise, NG, hydro).
- SET_POPS_D2_EQ_D1 (T here): STEQ(:,2) zeroed pre-MAXCH, pops slaved
  to depth 1 post-update. Consistent.

Issues found:
1. STALE COMMENT, solve_for_pops.f:84-92 -- claims GS_FREEZE "skips the
   SOLVEBA update entirely ... transfer-only sweep for FIXED pops".
   Actual code forces a LAMBDA iteration and pops DO update (matching
   the later comment and V6 policy). On a divergent field the Lambda
   step still moves pops (TRIDIAG run's 1e24% corrections happened on
   forced-Lambda iters). A true skip-SOLVEBA freeze is NOT implemented.
   Comment is misleading; code matches V6.
2. MAX_CHNG_LIM=1.0D100 in VADAT -- the divergence guard (STOP when
   MAXCH>MAX_CHNG_LIM) is effectively disabled; nothing aborts a
   runaway. Config choice, but explains the TRIDIAG run continuing
   through 1e24% corrections.
3. MAXCH masking: when DO_LEVEL_CHK fires (all but top-10 inc + top-10
   dec < 10%), MAXCH reports the 10th-largest change. Up to 10 levels
   can still move >10% while MAXCH reads <1%; the EPS convergence test
   uses this masked value. Upstream design (ignore few bad levels) --
   flag when interpreting "converged".
4. FIX_NE=T: charge-conservation row zeroed (FIXPOP), Ne frozen at the
   initial value for the whole run (ED read from POPS(NT-1) each iter,
   never updated). Converged ionization is consistent with the frozen
   Ne, not necessarily with sum(z_i n_i). Config choice; matters if
   ionization moves far from the initial ED.
5. FIX_T=T + DO_T_AUTO=F: radiative equilibrium never enforced -- T is
   the input structure for the whole run. SE-converged, not
   energy-conserving. Deliberate per config.
6. BA under KNOP is block-diagonal quasi-Newton: TX/TVX radiation
   response zeroed by design, advection off-diagonals off (INC_ADV=F).
   The "full" iteration has no depth coupling -- convergence relies
   entirely on the outer J<->n map contracting (fails for the MOM
   path under strong compression, works for the ES path). TRIDIAG
   identical to DIAG in this configuration.
7. Inner f-iteration bails at 20 iters with only an OUTGEN warning --
   a non-converged f is used silently for the rest of that frequency.
8. R_GRID_REVISED forces LAMBDA_ITERATION but does not clear
   FIX_IMPURITY (pre-existing upstream edge; REV_RGRID=F here).
9. Jacobian reuse (N_FIX_BA=2, FIX_BA=1%): after MAXCH<1% the BA can
   freeze for >=2 iters; corrections then use the ORIG_POPS fraction
   conversion -- an approximation bounded by the MAXCH_SUM>3*FIX_BA
   recompute trigger.

No source changes made -- all findings are documentation-level or
config-level, not code bugs in the cycle itself.

### KNOP refresh performance pass (24-Sep-2026, commits 41083d9 + d0eb12f)

Audit of FG_KNOP_REFRESH_V1 / KNOP_RAY_SOLVE_V1 found the alternating-
direction fixed-point sweep rebuilt BUILD_FREQ outputs every pass even
though they are pass-invariant (depend only on captured opacities and
static ray geometry, never on the iterated field).  Implemented:

- 41083d9: pass-invariant coefficient pack (DTAU, A0-A4, SR, Q, inner-BC
  ETA/CHI at NI) built once per ray and replayed each pass; guarded by
  KN_RAY_PACK_MB (2 GB/thread), exact fallback to per-pass BUILD_FREQ.
- d0eb12f: DLN_N_SAV hoists the per-frequency LOG(FREQ/FG_NU_NEXT) out of
  three ray-point loops in the in-sweep INTEGRAL path.

Runtime validation (TypeIIn_10d_DEV7_KNOP_ES, FINAL job 11962582, 4 OMP
threads, cold start): refresh t= 749.5/711.3/660.9/646.2/633.2/677.1/
748.7 s over its 1-7 vs 690-773 s pre-change baseline; great-iteration
931-976 s vs 955-1112 s.  Numerics clean (0 diverged rays, 0 clamp
events, 0 cap/bad writes; per-pass residual history matches pre-change
behaviour).  Net refresh gain only ~3-8%, NOT the predicted ~3x:
profiling shows all 4 threads saturated and no compute hotspot -- the
sweeps are memory/latency-bound, so removing BUILD_FREQ flops helps
little.  Kept anyway: verified bit-identical values, small real gain.

Largest remaining lever (not exercised): the refresh is embarrassingly
parallel over KN_NP rays -- raising cpus_per_task (4 -> 16) should cut
it ~3-4x.  Secondary candidates, each needing its own validation: MOM-
path equivalent of the pack (mom_j_cmf_nonmono_v1.f rebuilds+refactors
the tridiagonal per pass); SEED_FLD/WRITE_STORE iterate over the slowest
store dimension (~5 MB strided walk per element group); KN_CHI_COEF/
KN_ETA_COEF (~1 GB double) could be single precision.

### USE_LAM_ES=F maturity fix: quarantine + sweep-drift gate (25-Sep-2026)

The TypeIIn 10d F path (TRIDIAG family) converged in Lambda but the
BA gate could never open: the original test (MOM_GS_MAXCHG<0.05 AND
MOM_NEGJ_LRG==0) is unreachable on this model because the MOM refresh
residual is pinned ~1.3-1.7 by marginally contractive compressing
cells and ~750-2000 significant negJ repairs occur per sweep even in
a settling field.  Three changes, all VADAT-gated (commits 28b0cdd,
7b1e570, 069ed05):

- MOM_RF_KEEPNEG (def F): repaired cells (negJ reset, non-finite,
  cap) keep their prior iterate in J_NEXT_FRQ instead of storing the
  artificial repair.  The repair values are not fixed-point values
  and were reseeding the lagged read with noise: sweep negJ fell
  ~9E5 -> ~2E3.  Quarantined cells are excluded from the residual
  metrics and counted separately (MOM_GS_NQUAR); the write loop
  skips them.
- MOM_RF_GS_DESC (def F): descending refresh passes read the live
  J_NEXT_FRQ(:,L+1) instead of the pass-start JNEXT_RF_READ snapshot
  (true GS along the red-to-blue compressing characteristic).  Tested
  on gsdesc (11966433): it1 best residual 40.8 vs baseline 13.9,
  later iterations track baseline -- neutral-to-worse, left off.
- MOM_KNOP_REFRESH stall exit (unconditional): same-direction
  residual ratio >0.90 twice after pass 8 exits early; ~420 s of
  oscillating passes cut to ~70-100 s/iter.
- MOM_GATE_BULK (def F) + MOM_SW_* snapshot: the bulk counters
  (nbad5/nbad50/ntot/negJ/quar) are now taken at refresh ENTRY, i.e.
  the sweep-level drift |XM_solve - stored J_NEXT| -- the change the
  populations actually see.  The refresh-pass residual measures map
  self-consistency and asymptotes; gating on it can never open BA.
  Freeze while >10% of cells drift >5% per sweep; unsafe on capped
  writes, >0.1% of cells drifting >50%, any cell >10x stored,
  negJ/quar >10%, or refresh residual >1e2.  MOM_SW_* is reset at
  each sweep start (guarded against the refresh's per-frequency
  resets) and snapshotted before the pass loop.

Result (gatebulk 11969828, small partition): freeze released at
it25 when nbad5_sw crossed 10%; first BA computed it28; thereafter
the same Lambda<->BA alternation the T path shows (BA kick ~1e4 ->
1-2 frozen Lambda its -> reopen).  Corrections pre-gate-open were
contracting 17->13%.  Containment throughout: abort guard + full
pre-refresh field restore on unconverged refresh.

Caveat for restarts: resubmitting an immature run replays the
4-iteration warm-up and rebuilds J_NEXT (sweep drift resets to
~2e5); populations carry over via *_IN files.  Only resubmit when
the binary/config change is worth replaying the maturation.

Winning config keys (now in TypeIIn_10d_DEV7_KNOP_TRIDIAG.yaml):
MOM_RF_KEEPNEG=T, MOM_GATE_BULK=T, MOM_RF_GS_DESC=F.

### DO_NONMONO_SOB working setup on Dessart m15 2.82d (27-Sep-2026)

Config `configs/dessart_m15_2p82d_SOB.yaml`, model dir
`cmfgen_validation/m15_2p82d_knop_SOB/cmfgen/2.82days`.
Job 12001805 (interactive, 8 cpu): stable Lambda contraction from a
cold LTE start -- max % increase 1.6e11 (it1) -> 64 (it5) -> 17 (it10)
-> 4.4 (it24), luminosity pinned at L=4.266e8 inner / ~1.04e10 outer
every iteration, DTDR=1.7e-3 nonzero.  Same binary and same
SN_HYDRO_DATA (Luc's 70-point m15 file, sigma=-1.037 outer, ~-157 at
the interface) that diverged 1e52 -> 3e113 -> NaN segfault at it11
under the old config.

The load-bearing fixes, taken from Luc's extracted converged VADAT
(`cmfgen_validation/jnk_extracted/VADAT`) and the Dessart+Audit+
Hillier 2015 setup (arXiv:1503.05463, sec 3.1: inner tau~30 diffusion
boundary, lines treated purely via Sobolev):

- LSTAR = SN_HYDRO_DATA header luminosity (4.265857e8 Lsun).  With
  LSTAR=0 the inner boundary is dead (DTDR=0 every iteration) and the
  statistical-equilibrium iteration has no flux anchor -- this alone
  turned divergence into contraction.
- DIF=T + IB_METH=DIFFUSION.  Inner tau_ross~37 sits in the paper's
  tau~30 diffusion regime.  Side effect: JBAR_W_EDD=T is legal again
  (PERTJFEAU_IBC hard-stops under DIF=F; the earlier JBAR_W_EDD=F /
  NEWPERTJD workaround is no longer needed).
- Solver details matching Luc: ALLOW_OL=F, DC_METH=LTE, THK_LINE=T,
  DO_DIS=T, INC_PEN=T, LAM_VAL=200, FIX_BA=10, N_FIX_BA=2,
  T_INIT_TAU=1.5, SIG_GAU_KMS=3000 (NEG_OPAC_OPT stays forced SRCE_CHK).
- INC_CHG must stay F: CHG_EXCH_DATA is not in the repo and
  prepare_VADAT refuses INC_CHG=T without it.
- Kept VEL_LAW=11 + RD_IN_R_GRID=T (stratified composition preserved).
  Luc used VEL_LAW=12 + RD_IN_R_GRID=F, which resamples onto a
  log-tau_es grid but switches to uniform VADAT abundances -- a
  different physical model; only worth it as a control experiment.
- DO_NONLOC_SOB=F (baseline paper comparison; RH78 coupling is a
  separate experiment and did not cure the TypeIIn divergence anyway).

Forced-by-orchestrator (unchanged): GLOBAL_LINE=SOB, TRANS_*=SOB,
DO_SOB_LINES=T, SOB_FREQ_IN_OBS=T, NON_MONO_SOLVER=T, USE_LAM_ES=T,
SOL_METH=DIAG, NUM_BNDS=1, VEL_LAW=11, RD_IN_R_GRID=T, N_OB_INS=0,
N_IB_INS=0, INC_GRID=F, FG_OPT=INT/INS, FIX_T=T, DO_LAM_IT=T,
DO_LAM_AUTO=T, DO_NG=T, COMP_BA=T, DO_NOMONO_KNOP=F.

Activation still requires BOTH DO_NOMONO_SOB=T and DO_NONMONO_SOB=T in
vadat_overrides; `non_mono_mode: SOB` alone does not open the
non-monotonic block.

Caveat: run still in progress at note time (it24 of up to ~140 inside
the 1:59:00 interactive wall; requeue configured).  The DO_LAM_AUTO BA
transition and final convergence are not yet verified on this run.

### Build provenance + VEL_LAW=12 SIGMA fix (29-Sep-2026)

External review of the dev7 source found a pre-existing numerical bug
and two clean-build failures; all three are fixed and committed:

- `f2f4784` -- `new_main/subs/rv_sn_model_sniin.f`: the interior
  SIGMA estimate divided the weighted average of the two one-sided
  dV/dR slopes by `R(I-1)-R(I)` instead of the full span
  `R(I-1)-R(I+1)`, overestimating `|dV/dR|` by up to 2x on a uniform
  grid and biasing `SIGMA = R/V*dV/dR - 1`.  Reached only under
  `VEL_LAW=12` (cmfgen_sub.f:1044 dispatch); the VEL_LAW=11 path
  (`SET_RV_HYDRO_MODEL_V3` + `adjust_r_grid_v4.f`) uses the correct
  log-difference form and was never affected.  NOT flag-gated --
  the formula is wrong for any smooth V(R).  Every validated run in
  this file used VEL_LAW=11, so no prior result is corrupted; any
  VEL_LAW=12 run predating f2f4784 carries the biased SIGMA.
- `12f1f1a` -- Makefile fixes: all recursive makes now run via
  `$(MAKE) -C dir` (non-zero exit on sub-make failure; the old
  `(cd dir; make)` chains let scripts report success after a failed
  compile).  `d_lte` restored to the `cmfgen` target (link needs
  lte_hydro objects; fc12e0f's removal only worked with stale
  archive members).  `mod_ray_mom_store.f` built as a prerequisite
  of `comp_j_blank.o` in new_main/mod_subs (KNOP dispatch imports
  HNU_AT_OB before plane/ is compiled).
- `22deecd` -- `build.sh` / `build_debug.sh` wrappers
  (set -euo pipefail) + `.gitignore` covers build*.log, mod_dbg/,
  lib_dbg/, exe_dbg/.
- `e3b2b84` -- `d_lte` builds `OBJ_FILES` only.  The debug build
  exposed a real missing-symbol bug in the LTE utility:
  `lte_hydro/lte_sub.f` calls `WR_ABUND_INFO_V2`, but only
  `WR_ABUND_INFO_V3` exists (`subs/wr_x_info.f`) -- an upstream
  V2->V3 rename that lte_sub.f never followed.  Under -O3 the
  optimizer removes the unreachable call and `main_lte.exe` links
  anyway; under -O0 the link fails.  The cmfgen link needs only the
  lte_hydro archive objects, so the utility exes are out of the
  `cmfgen` target (build them manually with `(cd lte_hydro && make)`
  under nodebug if ever needed).

Verified: `./build.sh` builds end to end with zero errors
(`exe/cmfgen_dev.exe` 13:25, contains the SIGMA fix) and
`./build_debug.sh` now completes too (`exe_dbg/cmfgen_dev.exe`
13:22, -O0 -check bounds -traceback).  Both binaries postdate all
source changes.

Reviewed, no change needed in this tree (items were V6/cmfgen_source
artifacts only):

- `set_abund_clump.f` hydro-T clamp: the 0.3..1.0 (3e3..1e4 K)
  clamp exists only in the cmfgen_source variant.  V7's file is
  pristine upstream -- T_HYDRO is read from input_hydro.dat but
  unused; initial T comes from the SN_HYDRO_DATA path
  (RD_SN_DATA/SET_NEW_MODEL_ESTIMATES).  No V7 edit.
- `generate_full_matrix_v3.f` zero-population guard: V7 keeps the
  upstream unconditional `C_MAT(:,J)*POPS(J,L)` column scaling.
  The `ABS(POPS).GT.0` guard that the review flagged as fragile
  exists only in the cmfgen_source variant.  Upstream semantics
  retained per minimal-diff policy.

Runtime validation: a controlled model run on the rebuilt binary is
still required before calling the tree usable for production.
`run_model_grid.py configs/TypeIIn_10d_DEV7_STD.yaml` resolves
install_dir to this tree and dry-runs clean (job graph: one CMFGEN
job).  An interactive submission was attempted but blocked by
AssocMaxJobsLimit and cancelled; resubmit when a slot frees.  A true
VEL_LAW=12 test additionally needs a config carrying `vel_law: 12`
with existing hydro input -- none currently sets both that and
skip_stella.

### Stall-point audit fixes (29-Sep-2026, commits a9ab6b7..5e7d5cd + ed48071)

A static review of the DO_NONMONO_KNOP stall paths flagged six
hazards; all were verified against this tree and fixed as follows.
Runtime validation on a real model is still pending -- these are
compile-verified numerical guards, not a demonstrated converged run.

- `a9ab6b7` -- exit-side hysteresis on the compressing-depth set
  (fg_j_cmf_nonmono_v1.f).  A depth near SIGMA_TRUE=-1 or V=0 flapped
  the set on every NG/grid-revision pass; each flap rebuilt the ray
  geometry and reallocated I_P/I_M_DEPTH_NEXT ZEROED, cold-starting
  the GS coupling and holding GS_FREEZE forever.  Entry test
  unchanged (sigma < -1 or V < 0); exit now requires sigma > -0.98
  or V > +1 km/s.  The stored mask still uses plain thresholds so
  the FGINEXT_GS fingerprint is unaffected.
- `be1213e` -- in-loop bound checks on all five ray-bracket DO WHILE
  searches (J_PNT, both H_PNT loops; J_INDX_NM, H_INDX_NM).  A depth
  that brackets no interval previously marched K past the array end
  -- a silent hang in unchecked builds.  Now STOPs with depth/ray/
  endpoint diagnostics.
- `72641f2` -- quarantine persistence counter (mom_j_cmf_nonmono_v1.f).
  XM_BAD_NM is re-evaluated per call so a structurally failing cell
  retries forever holding a stale iterate; MOM_QUAR_AGE counts
  consecutive quarantined sweeps per cell, MOM_GS_NPQUAR reports
  cells at age 3+, printed via GET_MOM_GS_MATURITY and
  DUMP_MOM_JNEXT_GS.
- `d57d105` -- BA_NREJ_LIM (VADAT, default 20, <=0 disables).  The
  trust radius floors at 1.0, so consecutive rejections could
  alternate LAMBDA/rejected-BA forever and finish un-converged.
  Now a loud STOP.  Also committed `5c198c3` (BA_SE_RES2 rank-based
  second residual channel) which was in the tree uncommitted.
- `5e7d5cd` + `ed48071` (main repo submit_cmfgen.sh) -- GS_STATE
  file written every SOLVE_FOR_POPS call under DO_NONMONO_KNOP;
  generated slurms refuse to requeue or report completion when the
  run ends a segment >= KNOP_STALL_MAX (default 60) consecutive
  frozen iterations; .CMFGEN_KNOP_STALLED marks the run dir.
  Optional VADAT key GS_FREEZE_LIM (default 0 = off) for an
  in-model stop.

Reviewed, no change needed:

- MOM_KNOP_REFRESH_V1 is bounded (KN_M_MAXPASS=60, divergence abort,
  same-direction stall exit, bulk-drift exit) and conditionally
  retains the best-pass field -- the "burns full budget then reverts"
  failure mode is already handled; the compound stall is covered by
  the watchdog above.
- FG_KNOP_REFRESH_V1 is per-ray bounded (KN_RAY_MAXPASS) via
  KNOP_RAY_SOLVE_V1.
- Checkpoint fingerprint (bitwise R/V/SIGMA + COMP_DEPTH_IDX match)
  is deliberately strict; the mask-flap trigger is removed by
  a9ab6b7, and remaining invalidations on genuine grid changes are
  correct behavior.
- Quarantine is not a permanently latched cell state -- XM_BAD_NM is
  rebuilt every call; the persisted-store risk is now visible via
  NPQUAR rather than silently degrading the solution.

### Marginal-release BA trial (30-Sep-2026, commit 7b07fdf)

The single-tier MOM maturity gate kept the FreeNe M15 run frozen
>150 iterations on ~0.12% drifting cells that asymptote but never
converge under the one-bin upwind.  The gate is now two-tier:
MOM_UNSAFE_B is the catastrophic backstop only (capped writes, >1%
of cells drifting >50%, a >10x cell drift, >10% negJ/quarantined);
the strict maturity terms moved to MOM_FREEZE_B.  VADAT key
MOM_TRIAL_MIN (default 30, <=0 off): after that many consecutive
frozen iterations, with DO_BA_RES_ACC armed and the FG side clean,
the iteration is released to a resacc-gated BA step that rolls back
on residual growth.  N_FROZEN keeps counting, so GS_FREEZE_LIM and
KNOP_STALL_MAX still bound a run whose field never matures.
Runtime-validated: NO -- watch for 'MOM field marginal' on the first
released iteration.
(Superseded 30-Sep-2026 PM: accepted trials now refund N_FROZEN to
MOM_TRIAL_MIN -- see "Accepted-trial refund" below.)

### Numerical-robustness audit patch set (30-Sep-2026, commits 5eecc10, 978cd0c, a5d1c20, 06ec8ff, f87e1bf)

External audit of the CMFGEN/CMF_FLUX plumbing around non-monotonic
velocity fields.  All sites verified against the source before
editing; the LAW_TAU1 loop bug is upstream, identical in
../cmfgen_reference/cur_cmf.  Non-flag-gated changes are bug fixes or
degenerate-input guards; each documents why the stock path changes.

- `5eecc10` -- obs/obs_frame_sub_v9.f.  PRIVATE allocatable work
  arrays were unallocated on every non-master thread (writes through
  unallocated/shared descriptors: UB); per-thread ALLOCATE inside the
  region.  NI_MAX now uses max|V| + sum|dV| with CEILING and requires
  MAX_DEL_V_RES_ZONE > 0.  Band hunts bounded at the FREQ_CMF edges;
  sqrt radicands (ray Z, relativistic gamma) guarded; P(LS) > R(NI)
  hunt bounded at NI=1.
- `978cd0c` -- obs/cmf_flux_sub_v5.f.  LAW_TAU1 out-of-bounds read at
  I=ND fixed; TA now the cumulative trapezoidal depth, R ordering
  validated.  LAW_V1 normalized by max|V| with |V| (was V(1)).  SRCE_CHK
  target guarded against zero/subnormal ETA_CONT; ZETA/THETA divide by
  floored CHI; Rosseland accumulator and mean-opacity normalization
  guarded; 1-exp(-x) series branch for x < 1E-4; MU_AT_RMAX radicands
  guarded.
- `a5d1c20` -- plane/fg_j_cmf_nonmono_v1.f.  chi_hat first-failure
  diagnostic made reachable (test now runs on the raw value, catches
  Inf as well as NaN, then clamps).  NEW_R_GRID uses relative
  tolerances for R and VDOP_FRAC instead of exact .NE.
- `06ec8ff` -- new_main/subs/recompute_sigma_true.f +
  new_main/cmfgen_sub.f.  SIGMA_TRUE_CORE returns cleanly for ND < 3;
  clipped cells report locations.  Mean-opacity denominators
  (RLUMST, J_INT, ROSS_MEAN) guarded; 1-exp(-x) and x*coth(x/2)-2
  given stable small-x branches.
- `f87e1bf` -- plane/mom_j_cmf_nonmono_v1.f.  K_ON_J/NMID_ON_J/
  NMID_ON_HMID sanitized after GET_MOMS_NON_REL, after refresh-store
  loads, and after the solve; HUP_NM checks H_NEXT_FRQ/HNEXT_RF_READ
  entries for finite values.  Pre-existing uncommitted 29-30-Sep WIP
  in this file remains uncommitted.

Build: make cmfgen + make -C obs clean (ifort -qopenmp, nodebug);
new diagnostic strings verified in exe/cmfgen_dev.exe, exe/cmf_flux.exe,
exe/obs_frame.exe.  NOTE: `make -C obs` hits a pre-existing missing
`f77` for disp/mod_disp.f; build it with ifort first or run
`(cd disp && make)` workaround as before.

Runtime validation: NOT DONE.  No model has run through
run_model_grid.py against these binaries; the OMP fix especially needs
the 1-thread vs multi-thread cmf_flux reproducibility check.  Audit
item 16 (full +3 a_nu I proportional term) deliberately not
implemented -- proposed as a future runtime flag
(KNOP_USE_FULL_NU_PROP_TERM), requires checkpoint-version gating so GS
state is never mixed between operator definitions.

### Accepted-trial refund + GS_TRIAL_LIM (30-Sep-2026, solve_for_pops.f)

Live evidence from job 12040589 (Dessart m15 dif config): the
marginal-release machinery pulsed correctly -- it30/it32 trials
released, resacc accepted both (nrej=0), lum_ratio moved
0.16 -> 1.03 -- but N_FROZEN counted every frozen iteration
including accepted-release ones, and the strict raw-mature gate
(MOM_SW_MAXCHG <= 10, observed >= ~1.3e3, asymptotes ~1e3 on
marginal cells under one-bin upwind) can never produce the 2
consecutive mature sweeps needed to clear it.  The run stopped at
GS_FREEZE_LIM=60 while making real progress; job 12040591 (dbg
baseline) hit the same wall.

Change: `TRIAL_FIRED` is set when the MOM_TRIAL release arm fires.
After SOLVEBA, an accepted trial (BA_REJ_THIS=.FALSE.) refunds
N_FROZEN to MOM_TRIAL_MIN -- not zero: the post-release field is
re-immature by construction, so the next trial must re-earn marginal
status -- and increments N_TRIAL_ACC.  Rejected trials keep counting
(no refund).  New VADAT key `GS_TRIAL_LIM` (default 15, <=0 restores
the old count-forever behavior) vetoes further trials once
N_TRIAL_ACC reaches it, so a pulsing-forever run is still bounded.
N_TRIAL_ACC is written to `GS_STATE` (`n_trial_acc=`, after the
refund, so the wrapper sees post-refund n_frozen first) and adopted
on warm restart; a genuine 2-sweep mature release resets it.
`MOM_TRIAL terms` prints a 9th `trlim` term plus `ntr=`.
Runtime-validated on itest job 12042808 (see below).

### Paper-fidelity audit fixes (30-Sep-2026, mom_j_cmf_nonmono_v1.f)

Term-by-term check of the discretization against Knop I (A&A 496,
295), Knop II (A&A 501, 813) and Ceulemans+ 2025 (arXiv:2505.08309).
Structure verified correct (chi_hat as |GAM|/dlnv, sign-of-GAM
upwinding, equal-weight compressing coupling, N_CAP=CEILING(T4)
frequency matching, boundary clamp).  Two defects fixed:

- Dead midpoint CHI_HAT loop read GAM_MAG_NM before it was populated
  for the current call (uninitialized on first call, stale after)
  and its output was unconditionally overwritten by the depth-point
  loop.  Removed.
- Expanding-row upwind coefficient PSIPREV used the LAGGED
  GAM_FAC_PREV factor, so a sign flip of 1+sigma*K_on_J across
  adjacent frequencies injected a NEGATIVE upwind source -- Knop I
  Eq. 9 makes that coefficient positive by construction.  Now
  |GAM*GAM_FAC_PREV| at all three assignment sites (the header
  comment already claimed this behavior).

### Diagnostics formatting + frozen= format bug (30-Sep-2026)

Bordered `----` tables added around ITER_SUM, MOM/FG refresh
pass+exit blocks, GS sweep diagnostics, MOM_TRIAL terms, quarantine,
retained-field delta and requeue reports (fix_diagnostics.py).
Runtime caught a real descriptor bug the static audit missed: the
GS/iter line `4(A,A)` was one pair too many, leaving the logical
GS_FREEZE on an A slot -- printed `frozen=` + raw bytes.  Fixed to
`3(A,A)`; all other writes re-verified arg-by-arg.

### misc/ build repairs (30-Sep-2026, unrelated to model exe)

`make` in misc/ failed on the generic link rule missing `-lmod`
(mod_cmfgen symbols), `RETURN` in a main program
(rewrite_old_model_data.f) and an `AES14.4` format typo
(step_init_est_v2.f).  Fixed those three; `plt_steq.f` has further
pre-existing bit-rot left alone -- it does not affect cmfgen_dev.exe.

### Runtime validation: 3-iteration itest PASSED (30-Sep-2026)

`dessart_m15_2p82d_KNOP_MOM_dif_itest.yaml` (submission: debug ->
interactive partition, init_num_its: 3, use_existing_sn_hydro, all
Knop/GS flags as the production config).  Job 12042808 COMPLETED
exit 0 in ~16 min.  Verified live: bordered tables render, 9-term
MOM_TRIAL line with trlim + ntr=, `frozen= T` (format fix),
GS_STATE n_frozen first key + n_trial_acc persisted, ITER_HISTORY
columns, clean RVTJ exit, no fatal/NaN.  (First itest 12042506
passed the same way on the pre-format-fix binary -- its only visible
defect was the frozen= bytes.)

Full model resubmitted through run_model_grid.py:
`dessart_m15_2p82d_KNOP_MOM_dbg_dif.yaml` (+ GS_TRIAL_LIM=15):
FINAL 12043642, OBS 12043643 on small (16 CPUs, 64G, NUM_ITS=500),
cold start on cmfgen_dev.exe rebuilt with the refund, PSIPREV fix,
dead-loop removal and all table formatting.  Production CMFGEN
default raised to 16 CPUs (run_model_grid.py CMFGEN_PROD_CPUS) per
the R1 thread-scaling analysis; debug/interactive stays at 8.
### Non-monotonic mode prerequisites hard-enforced (03-Oct-2026)

The DO_NONMONO_SOB / DO_NONMONO_KNOP line-treatment gates in
rd_control_variables.f were promoted from warnings to hard STOPs, and
rd_cmf_flux_controls.f now applies the same gates to CMF_FLUX_PARAM so
the OBS spectrum cannot be run with a line treatment inconsistent with
the populations:

- DO_NONMONO_SOB: GLOBAL_LINE=SOB, or NONE with every present ion
  TRANS_*=SOB; USE_LAM_ES=T required (the F branch runs the
  clamped-sigma MOM_J_CMF_V11, outside the validated setup);
  LAM_SET=T warns only.  DO_NONLOC_SOB unchanged (warns if SOB off).
- DO_NONMONO_KNOP: GLOBAL_LINE=BLANK, or NONE with every present ion
  TRANS_*=BLANK; LAM_SET=T now STOPs (it rewrites edge lines to SOB);
  USE_FIXED_J / USE_FORMAL_REL / USE_J_REL / USE_DJDT_RTE /
  PLANE_PARALLEL(_NO_V) STOP; INC_GRID/ACCURATE unchanged; both
  USE_LAM_ES branches remain valid (FG formal vs MOM_J_CMF_NONMONO_V1);
  DO_BA_RES_ACC=F warns.
- Both flags T STOP (pre-existing); both F leaves the upstream path.
- GLOBAL_LINE=NONE still requires a TRANS_* key for every species
  slot in VADAT - a missing key aborts in GET_KEY_STRING upstream of
  these checks (pre-existing GET_KEY_STRING behaviour).

Verified on a clone of the jfnk_smoke model dir (no run-model-grid
submission; gate-matrix exercised with the built exe): KNOP+SOB STOP,
KNOP+LAM_SET=T STOP, SOB+KNOP both T STOP, SOB+GLOBAL_LINE=BLANK STOP,
SOB+USE_LAM_ES=F STOP, SOB+SOB+LAM_ES=T proceeded past the reader,
KNOP+BLANK (reference config) proceeded, both flags F proceeded.
control_variable_mod.f flag comments updated; vault notes
[[VADAT keys]], [[Sobolev non-mono solver]], [[Knop formal solver]],
[[Model pipeline]] updated in the same change.

### JFNK order-of-operations audit fixes (03-Oct-2026, commits 9956744, e8abf97, 900d515)

An audit of the JFNK order of operations (admission -> driver ->
residual/JVP -> eval sweep -> fallback) found the architecture sound
but identified defects in state save/restore coverage and the
self-disable path:

- JFNK_DISABLED was a no-op (jfnk_solve_corr_v1.f): after
  JFNK_MAX_STAGNATE consecutive failures the disabled branch set the
  LAMBDA flags but left DO_JFNK_NOW true, so Newton ran anyway under
  mutated module flags.  Now clears DO_JFNK_NOW; dead FORCE_LAMBDA
  local removed.
- Shared helpers JFNK_SAVE_MOD_STEQ_V1 / JFNK_RESTORE_MOD_STEQ_V1
  (jfnk_cmfgen_interface_v1.f) now bracket every
  JFNK_COMPUTE_STEQ_V1 / eval-sweep site and cover SE%STEQ_ADV and
  SE%T_EHB in addition to SE%STEQ/STEQ_ED/STEQ_T/STEQ_T_EHB.  The
  advective and EHB arrays previously leaked perturbed values into
  GENERATE_FULL_MATRIX_V3 (SOLVEBA right-hand side) and the resacc
  metric.
- Frozen mode keeps the module STEQ at the full outer-iteration
  value for the whole activation (the radiative split is now
  bracketed, and the frozen baseline rebuild is bracketed too), so
  a same-activation SOLVEBA_V13 fallback solves the residual the
  outer iteration actually measured.
- Coupled eval admitted without a registered JFNK_SWEEP_CB now
  returns the budget-exhaustion sentinel instead of silently
  evaluating a collisional-only residual.
- JFNK_RAD_SWEEP_V1 guards the ACCESS_F park/restore on
  JFNK_EDD_START_REC > 0 (the capture block is skipped when
  EDD_CONT=F and ACCURATE=F, leaving the -1 sentinel).
- Stale FREEZE_I_NEXT comment in fg_j_cmf_nonmono_v1.f corrected.

Verified: ./build.sh clean (zero errors, only the pre-existing
executable-stack ld warning on cmfgen_sub.o); new helper symbols
and the no-callback sentinel string present in exe/cmfgen_dev.exe;
all changes inside JFNK-gated code paths so the default
(USE_JFNK=F) path is untouched.  Runtime validation through
run_model_grid.py still pending - do not treat as a working
feature until a coupled JFNK activation is exercised on a model.
Vault note [[JFNK solver]] updated in the same change.

Runs resubmitted 03-Oct on the fixed binary (`7221de9`; the stale-
binary jobs were scancelled and the old validation run dirs
deleted): SOB_mask 12073333 (+OBS 12073334,
cmfgen_validation_SOB/m15_2p82d_mask), KNOP_mask 12073335 (+OBS
12073336, cmfgen_validation/m15_2p82d_knop_mask), JFNK_opt
12073580 (+OBS 12073581), JFNK_hyb 12073582 (+OBS 12073583) under
cmfgen_validation_JFNK/.  SOB_v2 (12072559) kept running on the
pre-fix binary - its config exercises no JFNK-gated code.

### JFNK/NG follow-up robustness fixes (05-Oct-2026)

- Newton retry state is restored to the base point before another
  trust-region solve; the GMRES recompute pass is explicit; the
  stagnation test compares preconditioned norms; trust-region norm
  bookkeeping no longer combines coefficients from restarted bases.
- The PTC preconditioner now rescales the actual BA solve by the signed
  effective diagonal ratio. Null rows use the PTC diagonal; the cached
  effective preconditioner diagonal includes LU ridge and identity-row
  repairs (including the accepted banded matrix). The banded solve now
  shares log-x/PTC/mask post-processing. FD retry floors compare the
  perturbation norm h*||v||inf, not scalar h.
- JFNK admission requires USE_LAM_ES=T because the MOM J_NEXT_FRQ /
  H_NEXT_FRQ stores do not yet have an eval freeze. NG history resets on
  held populations; NG_DEFERRED preserves a +1 NG target across an
  auxiliary SCRTEMP counter increment.
- Coupled sweep evals cold-start then restore FOLD, restore EDDINGTON,
  ACCESS_F, and two-photon accumulators, stop the frequency loop after
  JFNK_EVAL_BROKEN, and convert the former eval-time STOP sites into a
  rejected sentinel residual. Outside JFNK_IN_EVAL, the original STOP
  behavior remains.

Verified 05-Oct: `./build.sh` completed with zero errors; all edited
Fortran units recompiled. `exe/cmfgen_dev.exe` contains the broken-eval
fallback and no-callback diagnostics. The build emits the existing
executable-stack linker warning and Makefile duplicate-recipe warning.
The convergence diagnostics tests passed 16/16 (OUTGEN parsing skipped
because no sample OUTGEN was present); targeted `check_convergence.py`
masked/zero/genuine-convergence fixtures passed 3/3. Changes to the
JFNK solver/interface and eval paths are gated; with USE_JFNK disabled
those solver/eval branches do not run. NG history bookkeeping only
changes NG decisions when DO_NG is enabled.

Runtime validation remains NOT DONE: `squeue -u sbrennan` was empty
and the JFNK validation run directories contained no OUTGEN or
ITER_HISTORY. Do not claim coupled JFNK production-ready until a run
through `scripts/run_model_grid.py` shows activation or deterministic
deferral/fallback.
### FG ray-grid V/SIGMA_TRUE fingerprint (05-Oct-2026, shock-robustness review D1)

`FG_J_CMF_NONMONO_V1`'s `NEW_R_GRID` test previously fired only on R
drift, VDOP_FRAC/SOLUTION_OPTIONS change, or compression-set
membership change (hysteresis).  A velocity update at fixed R -- the
DO_V_PARAMS wind revision is the known path; anything that perturbs V
or sigma_true without crossing the set boundaries -- left
V_COEF/SIGMA_COEF, GAM, chi_hat and the ray insertion spacing on the
old field.  `MOM_J_CMF_NONMONO_V1` already fingerprints V (`V_SM_SV`);
the FG solver now does the same, comparing V and SIGMA_TRUE_IN
against the stored main-grid segment of V_EXT/SIGMA_EXT at 1E-6
relative (floor 1.0 in natural units).  No new arrays; the extended
grid is itself the snapshot.

Verified 05-Oct: `./build.sh` clean (zero errors);
fg_j_cmf_nonmono_v1.f recompiled; `V drift`/`SIGMA_TRUE drift`
diagnostic strings present in exe/cmfgen_dev.exe.  The routine only
runs under DO_NONMONO_KNOP, so the default path is unchanged.
Runtime validation pending (run_model_grid.py).

Vault notes updated: [[Knop formal solver]] (rebuild triggers),
[[Failure modes and guards]] (guard table row).
### Shock-robustness review fixes batch 1 (05-Oct-2026)

Audit of non-monotonic velocity handling (shock fronts, velocity
inversions) produced five actionable findings; three landed here:

- D1 `266e06f`: NEW_R_GRID fingerprints V and SIGMA_TRUE (vs the
  V_EXT/SIGMA_EXT snapshot, 1E-6 relative) so a V update at fixed R
  rebuilds the ray geometry -- mirrors V_SM_SV in the MOM solver.
- D3/D4 `04157cc`: GRID_METH=VEL_SPACE rejected inside
  ADJUST_R_GRID_V4 under DO_NONMONO_SOB/DO_NONMONO_KNOP (DO_VEL_
  REGRID_V2 interpolates R(V), assumes monotone V);
  GAMRAY_TRANS=RAD_TRANS rejected under DO_NONMONO_KNOP
  (GAMRAY_SUB_V3/CMF_FORMAL_REL_V4_GAM is a signed-sigma solver fed
  the unclamped gradient).
- THK extension in FG_J_CMF_NONMONO_V1 floors the beta-law exponent
  source at sigma >= -0.999: SIGMA is not re-clamped under KNOP, so
  sigma(1) < -1 made BETA < 0 and V_EXT diverged.  Warns when the
  floor engages (THK_CONT=F is the pipeline default; latent defect).

Verified 05-Oct: `make cmfgen` clean (zero errors), edited files
recompiled, new diagnostic strings in exe/cmfgen_dev.exe.  Gates are
flag-gated; default path unchanged.  Runtime validation pending.
### Shock-robustness review: DO_NONLOC_SOB hardening (05-Oct-2026)

D5 `ff58a30`: sob_rh_couple.f now drops RH78 kernel pairs with
V(I)<=0 or V(J)<=0 (the mu^2 geometry uses V^2 and hides sign flips;
the pair weight carried signed V(J)) and reports the drop count.
GAMH_L uses |V| floored at 1 km/s -- the 1/V table could hold Inf at
a stagnation depth and leak into the net rates.  Infall regions now
get no non-local correction (documented limitation of the kernel
under DO_NONLOC_SOB).

Also note: a concurrent uncommitted solve_for_pops.f edit dropped a
`//` between two WRITE strings (line ~1035); restored so the build
passes.  That WIP (BA_LUM_GATE damping) is left uncommitted for its
owner.

### Shock-safe luminosity diagnostics (06-Oct-2026)

Inner/outer luminosity agreement and drift are no longer convergence
conditions: shocks permit physical energy deposition between the two
boundaries.  `BA_LUM_GATE` / `BA_LUM_DCHG` remain readable no-op keys
for restart compatibility.  `ITER_SUM` reports `lum diag` ratio and
relative drift with `[INFO]`; no luminosity quantity demotes, holds,
caps, or delays SOB or Knop/JFNK updates.  Build and runtime validation
are recorded with the associated submission.

### SOB-to-JFNK staged workflow (06-Oct-2026)

The production sequence keeps the hydro snapshot temperature fixed and
uses a converged SOB population checkpoint to initialize Knop/frozen
JFNK.  The target config sets `warm_restart_from` and
`warm_restart_population_only: true`.  The workflow requires identical
source/target `SN_HYDRO_DATA` and matching MODEL NT/ND, then stages only
POINT1/POINT2/SCRTEMP/MODEL.  It does not copy RVTJ, EDDFACTOR,
FGINEXT_GS, or MOMJNXT_GS across transfer operators.  The target rebuilds
its Knop radiation and lagged-frequency state while preserving the SOB
population basin.  SOB and JFNK configs use the same diffusion inner
boundary, thin outer boundary, dissolution/penetration physics, and
fixed-T controls.

### FIX_IMPURITY disabled in pipeline profiles (06-Oct-2026)

`FIX_IMP=T` was implicated in the SOB stall at job `12101701`: the
generated VADAT confirmed the switch, while CORRECTION_SUM reported
mostly zero corrections during full iterations.  The mechanism is
FIXPOP_IN_BA_V3 zeroing BA equations for low-abundance/low-fraction
stages at each depth; that suppresses exactly the population corrections
the full iteration must apply.  All production SOB/KNOP/JFNK YAML
profiles and the `NON_MONO_SOB_RECOMMENDED` recipe now stamp
`FIX_IMP=F`.  The VADAT key and FIXPOP_IN_BA_V3 remain for explicit
experiments, but are no longer part of the submitted solver recipe.
Fresh SOB submission: CMFGEN `12103504`, OBS `12103505`, convergence
monitor `12103506` on source commit `c13f017`; generated VADAT confirms
`F [FIX_IMP]`.

### Correction-cap defaults for the non-monotonic recipes (06-Oct-2026)

Live SOB job `12103504` is bounded, not divergent: the it21
`BA_PROBE_LAM` trial increased the SE residual and was rolled back;
subsequent Lambda calls reduced `SE_RES` from `6.33e13` to `4.06e13`
by it38.  `MAX_LAM=1e10` is already effectively unclamped and was not
limiting the run.  `MAX_LIN` is now `1e1` in all current SOB, Knop and
JFNK YAML profiles and in the non-monotonic/default recipes; saturated
or probe steps remain separately bounded by `BA_DAMP_CAP=1.3` and
residual acceptance.  Do not tighten `MAX_LAM` as a generic fix: small
Lambda caps made the LTE-to-NLTE cascade impractically slow.  Tune
`BA_DAMP_CAP`, `BA_ENTER_RAW`, or the residual gates only when a live
full-step failure shows that the ordinary `MAX_LIN` cap is binding.

### Fresh frozen-J JFNK submission (06-Oct-2026)

Submitted the frozen-J JFNK validation run
(`DO_NONMONO_KNOP=T`, `USE_JFNK=T`, `JFNK_COUPLED=F`) as a cold start
through `run_model_grid.py` on the standing resolved config
`cmfgen_validation_JFNK/m15_2p82d_jfnk_cpl/config_resolved.yaml`
(model `dessart_m15_2p82d_JFNK_frz_v1`).  The SOB-seeded
`warm_restart_from` keys are absent, so no checkpoint state is staged.
The emptied epoch directory was repopulated with the approved
`SN_HYDRO_DATA` (md5 `ea6c5d5bd2d1c0ae785d70d4c9a56053`, identical to
the SOB/KNOP validation copies) so `use_existing_sn_hydro` applies and
the inner-boundary fingerprint matches the aborted warm-restart design.
Jobs: CMFGEN `12104645` (apu, running), OBS `12104646` (afterok),
convergence monitor `12104647` (30 s cadence).  Source commit
`782c04c`; `cmfgen_dev.exe` (14:52) includes `9213576`.  Generated
VADAT echoes `T [DO_NONMONO_KNOP]`, `T [USE_JFNK]`, `F [JFNK_COUPLED]`.
NOTE: this resolved config predates the `FIX_IMP=F`/`MAX_LIN=1e1`
profile updates and still stamps `FIX_IMP=T`, `MAX_LIN=3.0D0` --
deliberately retained to reproduce the prior recipe; revisit if the
full-iteration corrections stall as in the SOB run.

### VEL_LAW 13/14 + JFNK consistency pass (08-Oct-2026, branch devin/m15-zlaw14-jfnk)

Work on the stalled `m15_2p82d_SOB` 2.82d convergence, all submitted
through `run_model_grid.py`. Commits `c958555`..`6296181`.

- `new_main/subs/rv_sn_model_snhyd.f` (new): `VEL_LAW=13` reads
  `SN_HYDRO_DATA` with the law-12 tau regridding; `VEL_LAW=14` builds
  a blended `dlog(tau_Ross)+1.5*dlog(R)` grid from the kappa column
  (thin layers get radial coverage, thick regions keep tau coverage,
  no 2-dex seam at the zone handoff). Dispatched in `cmfgen_sub.f`,
  flagged SN-model in `rd_control_variables.f`, built via
  `new_main/subs/Makefile`.
- `do_tau_regrid_v2.f`: TAU_OLD(1) uses floored TA, not FLUX_MEAN
  (zero outer opacity put -Inf into the tau grid).
- `cmfgen_sub.f`: DO_LAM_AUTO stable-count guard reverted -- it
  parked non-mono runs in LAMBDA at ~45% demand.
- JFNK consistency: `CMAT_RAD` cached across activations
  (`JFNK_CMAT_MAX_AGE=4`; fixed point is invariant to BA_rad);
  deterministic re-arm `JFNK_REARM_EVERY=10` admitted calls replaces
  the process-lifetime disable that silently reset on requeue;
  SOB zlaw14 config raised `JFNK_MAX_EVALS` 40 -> 200 (one GMRES
  cycle alone can spend ~60 JVPs).
- `solve_for_pops.f`: per-great-iteration `ITER_TYPE` banner naming
  the dispatched solve (LAMBDA/FULL/JFNK-FRZ/JFNK-CPL).
- `Makefile_definitions`: viper `INSTALL_DIR` is self-locating
  (this file's own directory) because ptmp1 and ptmp2 are two
  worktrees of the same repo and each must install into itself.
- `tests/f90` (30 checks: solver core vs stub int_mod, real
  interface routines) + `tests/python/jfnk_outgen_audit.py` (live
  OUTGEN validator + plots under `tests/plots/`). All tests pass.
- `build.sh` now runs the tests/f90 suite after linking; a failed
  check fails the build like a compile error. The suite also writes
  per-iterate |F|/dtau histories that
  `tests/python/plot_solver_tests.py` renders to
  `tests/plots/solver/newton_convergence.png` (plot step non-fatal).

Submitted: FINAL `12131626`, OBS `12131627`, monitor `12131628`
(config `dessart_m15_2p82d_SOB_jfnk_zlaw14_jfnk.yaml`,
`install_dir` ptmp2 worktree). VADAT verified `DO_JFNK=T`,
`JFNK_MAX_EVALS=200`, `VEL_LAW=14`; blended grid executed.
Runtime JFNK admission under SOB: UNVERIFIED -- still deferred
(`LAMBDA-designated iteration`) at it15, normal for the mode
controller while MAXCH is high.

### Status update (2026-10-08, later): null-row projection fix

- First live JFNK admissions (it31/33/35, FINAL 12131626) all
  aborted at the first GMRES matvec: `precond clip -- aborting`,
  info=3, zero useful evals. Root cause: identity-replaced near-null
  BA rows returned z_i=v_i with residual components ~1e12, so the
  |z|>1e8 nonlinear clip fired on every apply. Three consecutive
  aborts disabled JFNK for the rest of that process.
- Fix: `JFNK_ROW_REG(J,K)` mask persisted at factor time; masked
  components projected to zero at apply in all branches (LU, LU-fail
  Jacobi, DGETRS-fail, pure Jacobi, banded, log/PTC post-loop).
  Linear projection keeps M^-1 consistent; no clip event. Detection
  reads raw BA_A (ridge retries cannot un-mark). Counter
  `jfnk_n_reg_proj` in diagnostics.
- Commit 0aa5c4c on main + devin/m15-zlaw14-jfnk, pushed. Fresh
  cmfgen_dev.exe 17:26; 33/33 build-time tests pass.
- Job 12131626 graceful-stopped at it37, now PENDING
  (AssocMaxJobsLimit); requeue picks up the fixed binary and the
  JFNK disable resets on process start.
- Also fixed: .mod-shadowing hazard in tests/f90 (cwd .mod resolves
  ahead of -I; build_and_run.sh now cleans them).

### Status update (2026-10-09): cross-activation state disarm

Numerical audit of module state that survives between
`JFNK_SOLVE_CORR_V1` activations found two leaks; both fixed and
covered by the build-time suite (33/33 checks pass).

- Stale PTC scale/anchor (`jfnk_ptc_scale_ready`,
  `jfnk_ptc_anchor_ready`, `jfnk_ptc_x0`) persisted into the next
  activation's F(0) identity test, which runs before
  `JFNK_SET_PTC_SCALE_V1` re-arms.  Frozen mode would fail the test
  on every reactivation (permanent JFNK disable); coupled mode
  adopted a contaminated baseline and polluted every FD JVP.
  `jfnk_cmfgen_interface_v1.f` now clears both flags at INIT and
  clears `scale_ready` before the rebuild-baseline eval.
- `gmres_m_v1` kept its hookstep basis descriptor (`jj`, `beta`) in
  SAVE locals.  `JFNK_DIAGNOSTIC=T` runs the m-dimension experiment
  (m_test up to 40 > production m), leaving a `jj` that a later
  info==2 recompute indexes past the caller's m-sized h/v/y arrays.
  The pair is now module state (`jfnk_gmres_jj`,
  `jfnk_gmres_beta`), checkpointed per Arnoldi cycle, guarded by
  `jj <= m` on recompute, and saved/restored around the diagnostic
  experiment.
- Detail: `cmfgen_obsidian/solvers/JFNK solver.md`,
  "Cross-activation state disarm (09-Oct-2026)".

### Status update (2026-10-09): coupled JFNK under SOB + cooldown/re-probe

The SOB coupled tier is implemented and runtime-exercised.
`JFNK_SOB_LINE_PASS_V1` (cmfgen_sub.f 6241-6483) replays the
`SUB_SOB_LINE_V3` line loop at the perturbed populations at the end
of `JFNK_RAD_SWEEP_V1`; `JFNK_COUPLED=T` is no longer rejected
under `DO_NONMONO_SOB` (commit `0f6b9f6`).  First runtime test on
`m15_2p82d_SOB_jfnk_cpl` (jobs 12143251/12143252): eval machinery
verified (F_eval vs F_outer ~4-8% rel diff, deterministic evals,
clean state restore), but GMRES stalled on all three admissions --
the residual sits in near-null field-consistency directions of the
coupled Jacobian that the local BA preconditioner cannot reach.
No accepted coupled Newton step yet; admission gates and Lambda
contraction remain the working convergence path.

Lockout is now a cooldown, not permanent: `JFNK_MAX_STAGNATE`
(VADAT, default 3, clamp 1-20) consecutive failed solves park JFNK;
every `JFNK_REARM_EVERY` eligible calls (VADAT, default 10, clamp
1-500) a single-strike re-probe runs; a failed probe re-locks, a
successful one resets the strike count (commits `7e519e6`,
`8d7e068`).  Cooled calls print `JFNK DEFER` at admission.

JFNK OUTGEN prints are now one value per line with expanded labels
(`res` -> `residual`, `it` -> `iteration`, `dtau` ->
`pseudo-timestep`, etc.; commits `9131049`, `ffc18ce`).  Binary
relinked 16:51; 34/34 JFNK tests pass.

### Status update (2026-10-10): SCHUR preconditioner option + git re-init

Research-plan integration: `JFNK_PRECOND=SCHUR` adds a two-level
(n,T) Schur-complement preconditioner.  Routines
`jfnk_schur_form_v1` / `jfnk_schur_apply_v1` live in
`jfnk_solver_v1.f90` so the unit-test driver (`se_driver`) calls the
identical production code via `external` declarations - satisfies the
"exact Fortran in unit tests" requirement.  Per-depth pop blocks use
normal-equation LSQ (Tikhonov `lambda^2 I`, ridge rows, reg-row
normalization) WITHOUT row equilibration: measured bisection on the
shock surrogate shows the `rsf^2` reweighting corrupts `J_nn^-1`
(MZDIFF 0.25 vs 0.0075 raw).  Dense NDxND Schur complement via
DGETRF/DGETRS.

Surrogate verdict (se_driver, systems/shock): SCHUR INFO 3 vs LSQ
INFO 0 NITS 25 - the pop-only elimination drops the dominant T-row
regularization on masked-T specs.  LSQ remains the production
default; SCHUR is an opt-in VADAT value pending a cross-depth-T
surrogate (`tcoup` knob planned for gen_se_system.py).

Git: the worktree parent `cmfgen_nonmono_dev_V7_OLD/.git` is gone
(tree deleted; only a dangling `.git` pointer file remained).  The
tree was re-initialized in place; baseline commit `a15dba1` covers
all session solver work.  Prior history is lost with the parent.
Fresh run submitted against this state: job 12151424
(`m15_2p82d_SOB_jfnk_zlaw14_fresh`, cold start, JFNK_TAUTO=T,
JFNK_PRECOND=LSQ).  Binary relinked 23:07; unit tests pass except
the documented LU-on-shock rank-deficiency blowup (pre-existing).
