# Phase 2.30 - Production capture and baseline replay audit

Date: 22-Nov-2026. Status: PARTIAL. P2 remains disabled.

## Scope and decision gate

The question is whether captured production inputs support a faithful
P2-off replay before testing the midpoint damper. This serves the
KNOP formal / USE_LAM_ES=T / DIAG diagnostic configuration. It does
not validate the moment or TRIDIAG branches. No production numerical
implementation, tolerance, coefficient, or acceptance gate was changed.

The exit criterion is a baseline fidelity verdict, or an exact
missing-state blocker. A manufactured map and a smaller damped step
cannot establish convergence of the production map.

## Source and capture identity

- Requested source baseline: 25f6b52337c2b169846c4db8b682e222a535da05.
- Capture job: 13488131, submitted through scripts/run_model_grid.py.
- Capture source recorded by submission: 92284e61 plus eight dirty
  audit/artifact paths; do not describe the whole tree as clean.
- Frozen executable directory:
  exe_frozen/358814d3e2d34df2_8a62e7ea3fa04408.
- Configuration: configs/dessart_m15_2p82d_VAL_knop_formal_rd.yaml in
  the parent project. LTE_EST=T, warm_restart=false, P2=F, dump rays
  81/82/84. These are diagnostic capture inputs, not a P2-on model.
- Capture directory, relative to the parent project:
  cmfgen_validation/m15_2p82d_VAL_knop_formal_rd/cmfgen/2.82days.

The run reached iteration 250, then returned exit 7 from the wrapper's
GS-stall check: 149 consecutive frozen iterations versus limit 60.
This is evidence of non-convergence, not a segmentation fault.

## Corrections to the earlier handoff

The claim that the capture produced no decks was wrong. Direct inspection
finds KNRDUMP files for all three selected rays through refresh 250,
KNAPPLY files and KNSTAGE files. No new model run is required merely
to obtain those inputs. Completeness and replay fidelity still need
independent checks; presence alone does not prove either.

The claim that a stable committed field is itself a fixed point was also
too strong. Only a verified undamped map defect can establish that.
The previously added mapd prototype has not established arbitrary
full-candidate seeding or raw VNEW capture. It is not acceptance evidence.

## Cold-start field behaviour

ITER_HISTORY records:

| Iteration | Refresh residual | One-sided store drift | Retained-field delta |
|---|---:|---:|---:|
| 1 | 1.6401 | 1.0000 | 1.0000e4 |
| 10 | 2.2162 | 6.2107 | 0.55826 |
| 250 | 1.0429 | 0.10304 | 0.05 |

These are distinct diagnostics. Neither a reduction in retained-field
delta nor population corrections proves refresh convergence. The refresh
residual remains far above KN_RAY_TOL=0.002.

At refresh 250 the per-ray aggregate census gives:

| Ray | Class | Exit | Bad cells | Residual |
|---|---|---:|---:|---:|
| 81 | STAG | 3 | 1371 | 0.0292834254311299 |
| 82 | STAG | 3 | 5829 | 0.0675022582854450 |
| 83 | STAG | 3 | 101307 | 1.04288848798479 |
| 84 | FULL | 1 | 0 | 0.00168456697134306 |
| 85 | INV, reason 5 | -1 | 0 | 0 |

Ray 83, not ray 81, carries the largest reported residual in this
cold-start capture. Treating only rays 81/82 cannot by itself establish
all-ray maturity. No cause for ray 83's residual is claimed here.

## Census telemetry defect found by source review

The per-channel columns in KNRAY_CENSUS_r250.txt contain values such as
4.19955798965059563E-322, 85/85 and 4.94065645841246544E-324, repeated
across rays. They are not usable pass counts or channel residuals.

In plane/fg_j_cmf_nonmono_v1.f:

- KN_PA_RES/OMG/IX/NP are allocated only inside IF(KN_Q_RESID).
- KN_Q_RESID defaults F and is absent from the capture VADAT.
- The independently gated KN_RDUMP census writer passes slices of these
  arrays unconditionally to kn_census_row.

This is an unallocated-array argument defect in diagnostic output.
Do not reconstruct M/P schedule counts from the corrupted columns.
The audit instead requires separate M_SOLVE_END/FLDP_END_DUMP stage
records and checks their P count against the binary deck.
The solver equations have not been patched on the strength of this
telemetry finding. A diagnostic fix needs a regression covering dump
capture with Q both off and on before another production build.

## Reproducible selected-input validation and replay

New utility: unit_tests/python/phase230_capture_audit.py.
Local utility commit: 1e0e02ce; stage-validation test commit: e90044a0
(5/5 PASS). No changes have been pushed. The captured production solver
source is unchanged by these audit helpers.
It uses the independent v2 parser, checks record closure, ray identity,
finite array data and P2-off capture flags for refreshes 1, 10 and 250.
It obtains M/P solve-end counts from stage headers, then uses the existing
phase227_records.normalize_omega bit/order/range checks to create KNOMG2
schedules. It hashes the selected decks and schedules. Original capture
files are opened read-only and are not rewritten.

This is selected-input validation, not a complete 750-deck inventory or
an assertion that every required runtime input was captured. The v2 deck
schema does not directly serialize every PK_* work array: those are
rebuilt by the production solver. Exact coefficient fidelity remains
a replay question.

Allocated audit job 13557300 uses the interactive partition, two requested
CPUs and 30G, with OMP_NUM_THREADS=1 and MKL_NUM_THREADS=1. Its bounded
sequence is selected-input validation, an O3 replay-driver build linked
against existing production objects, then two P2-off recorded-schedule
replays per selected ray at refresh 250 in separate directories.
No P2-on run is in this job.

The first attempt, 13557295, stopped before validation because its Python
environment lacked numpy. The replacement loads python-waterboa/2024.06
explicitly. Runtime artifacts and logs are under:
unit_tests/work/phase230_review_20261122/ (untracked run products).
Compiler output and executable/archive hashes are written by the job.
A matching exit class alone will not be called replay fidelity: full P,
M and IM_NI comparison remains required.

## SOB-population BLANK spectrum

Separate OBS-only job 13557228 was submitted through run_model_grid.py
using configs/dessart_m15_2p82d_VAL_sob_diag_blank.yaml. It stages unchanged
VADAT, RVTJ and POP* from the completed SOB population job 13474708 in
an isolated m15_2p82d_VAL_sob_diag_blank directory. Byte comparisons of
VADAT, RVTJ and every POP* passed after staging.

The CMF parameter preflight selects NON_MONO_SOB_BLANK: GLOBAL_LINE and
TRANS_* are BLANK, DO_SOB_FLUX_OBS=F, NUM_ES=2. The population calculation
remains SOB. The original SOB-rendered spectrum is preserved. The new
spectrum is not yet verified: the job is waiting on QOSGrpCpuLimit.

## Remaining acceptance work

- Complete selected-input validation and inspect any exact failure.
- Compare P2-off replay fields, M-to-P boundary and pass histories with
  production stages; demonstrate repeatability.
- Correct the independent full-candidate raw-map diagnostic before using
  it to accept P2. Do not use mapd's existing source tests as proof.
- Only then run P2 A/B, fixed-cell and significance analyses, clamp-branch
  checks, and thread-count comparisons. None has passed this audit yet.
- Investigate ray 83 separately if the objective is global field maturity.
- Validate the queued BLANK spectrum's controls, finite flux and provenance
  after the job executes. Queue submission is not a completed spectrum.

## Follow-up: gradient singularity versus sweep order

User request: investigate Luc Dessart's SIGMA treatment and adapt the
solver only after identifying the ray failure. The online references
checked on 22-Nov-2026 were:

- Dessart, Audit & Hillier (2015), Sec 3.1 and its angular-optical-depth
  footnote: https://ar5iv.labs.arxiv.org/html/1503.05463 . The substitution
  max(abs(1+sigma*mu^2), -0.2*(1+sigma)) limits the Sobolev optical-depth
  singularity. Here sigma=dln(V)/dln(R)-1. It is not a prescription to
  replace the signed CMF velocity gradient globally.
- Knop, Hauschildt & Baron (2009), large-scale algorithm, introduction
  and Sec 3: https://ar5iv.labs.arxiv.org/html/0909.4264 . A monotonic
  coupling is an initial-value problem in wavelength; changes of coupling
  sign require a wavelength boundary-value solve. The paper discusses
  the lagged C-corner terms of the ordered GS solution.
- Ceulemans, De Ceuster & Decin (2025), Secs 2.3 and 3:
  https://arxiv.org/html/2505.08309v1 . Large Doppler shifts relative to
  frequency spacing can generate oscillatory discrete modes. This is
  relevant background, not proof that their stabilization is needed here.

Source inspection: KN_ORD_SPV_IMPL uses only FLD(r,l-1) when GAM>0,
and no donor when GAM=0. A wholly non-negative-GAM ray therefore has
one-way frequency dependencies. The existing frequency-major SWEEP_FLD
can complete the spatial source stencil at the upstream frequency before
advancing, whereas SWEEP_FLD_ORD previews a not-yet-marched spatial row.
No change to SIGMA, opacity, DS, tolerance or midpoint formula is needed
to compare these algorithms.

The replay-only `sweep ordered|frequency` mode now selects the existing
production algorithms. Both arms keep P2 off and the captured deck
coefficients unchanged. Job 13557757 completed all six arm/repeat pairs
at refresh 250. Artifact directory:
unit_tests/work/ray_gradient_20261122/sweeps/summary.json.

| Ray | GAM minimum | GAM maximum | Ordered residual/class | Frequency residual/class |
|---|---:|---:|---|---|
| 81 | 6.6110e-10 | 1.8551e-9 | 0.01887424 / STAG | 5.805058e-21 / FULL |
| 82 | 8.5514e-10 | 1.6153e-9 | 0.04484838 / STAG | 4.945920e-21 / FULL |
| 84 | 1.3019e-9 | 1.3022e-9 | 0.001684304 / FULL | 6.681912e-23 / FULL |

Frequency-major took two reported P passes in each case, with zero bad
cells. Both repetitions had byte-identical field files in every arm.
The ray-84 max P difference against production remained 3.84775e-8,
within the existing replay discrepancy. This comparison used the native
controller, not recorded-history omf. The ordered replay still differs
from the captured trajectory; this result does not close the omf
fidelity gate. Its M-to-P boundary differences against production are
1.6586e-4 (ray 81) and 8.1129e-5 (ray 82), explicitly unclosed.

Decision: the Sobolev projected-gradient-zero explanation is rejected
for the captured 81/82 failures, since their GAM is strictly positive.
It remains untested for ray 83, whose exact input deck is missing.
The sweep-order comparison supports a narrow non-negative-GAM dispatch
correction, subject to regression, build and live-model checks.
Wholly negative and mixed-sign rays must keep their previous route.
A dispatch regression against the original solver failed 4/8 checks
(job 13557822): the all-positive and zero-gradient ordered requests
selected the ordered path and differed from the frequency-major fields.
The negative and mixed-sign retention checks passed. The new regression
must pass after the correction; no correction has yet been built.

Jobs 13557701 and 13557718 failed before solver execution because the
MKL module did not place the linked libmkl_intel_lp64.so.3 on the runtime
loader path. Job 13557757 explicitly exported MKLROOT/lib and succeeded;
the replay argument tests passed 29/29. These environment failures are
not solver failures.

A one-iteration, P2-off ray-83 continuation configuration has been
prepared in the parent project (dessart_m15_2p82d_KNOP_ray83_probe.yaml).
It is not submitted yet. It preserves the failed cold-start checkpoint
in its original directory, uses an isolated target, disables requeue
and OBS, and requests dump rays 81/83/84. This future capture would be
restart-adjacent, not the identical iteration-250 production trajectory.

Verdict (22-Nov-2026): PARTIAL. Production non-convergence confirmed;
fresh captures found; diagnostic census defect identified; bounded
P2-off replay audit submitted. P2 remains disabled. No numerical cure
or production-operator P2 success is claimed. No GitHub push authorized.

## Follow-on, 22-Nov-2026 (later)

Dispatch candidate applied in plane/fg_j_cmf_nonmono_v1.f: ORD_SWEEP
now requires ANY(GAM(1:NI,LS) .LT. 0), so wholly non-negative-GAM rays
take the existing frequency-major route; mixed-sign and negative-GAM
rays keep the ordered march. No SIGMA clip, no coefficient change,
P2 still off. Build/regression job 13558165 compiled and linked all
production exes cleanly; the aggregate suite reported 9 failures of
which none come from the dispatch edit:

- test_knop_oneway: fixture authored with a mixed-sign profile where
  the case required all-zero GAM; corrected (gzzo -> NCASE 16 with
  GAM zeroed, reference gzzf). Rerun pending in job 13558936.
- test_knop_lagged_field, test_gs_drift_diag,
  test_obs_nonmono_dispatch: stale source pins drifted by earlier
  commits (RES.GT.HUGE -> .NOT.(ABS(RES).LE.HUGE) at 98acf207;
  GS_SW_MAXSYM added to the drift line at 3913cc50; OBS_FG_KNOP
  widened the OBS-grid gate at cmf_flux_sub_v5.f:1754). Pins updated
  to the current semantics and verified against source text.
- test_auto_ion, test_sobew_angular_floor: missing run fixtures
  (VAL_stock dir, FRESH_sob_dieoff RVTJ) - environmental.
- test_par_opac: opac/mod_par_opac.f not a sibling .o and absent
  from lib archives - exact-object rule cannot bind it.
- test_es_conv: undefined set_fg_q_resid_ at link, already recorded
  in CMFGEN_PHASE220_FINAL_ACCEPTANCE_AND_BUILD_CERTIFICATION.md:205
  (GUARD_SRCS link list never updated).
- test_phase227_io: run_all_tests invokes it without its required
  --prepare/--binary/--run argument - harness gap.
- DGELSY parameter-4 error inside the f90 GMRES harness; MKL-version
  artifact under mkl/2026.2, under triage.

Original-binary omf baseline (job 13558443, exe sha256
96ccc2f62336ae7ef1b9557419e31f3fa1cd11d20de0e57640258471acce6595,
linked against the pre-dispatch libsubs): ray 84 replay repeated
byte-identically; rays 81 and 82 abort at pass 3 with a runtime pass
mismatch against the KNOMG2 r250 schedule (np_P 10 and 13, np_M 2).
The omf fidelity gate therefore stays open for the stalled rays:
the recorded-trajectory replay cannot proceed past pass 3 while the
converged control reproduces exactly. Deck gam_min at r250 is
strictly positive for all three armed rays (6.6e-10, 8.6e-10,
1.3e-9), so each is in the population the dispatch change reroutes.

Pending: job 13558936 reruns the corrected batteries and replays
decks 81/82/84 in base mode under the original and candidate
executables for the production-object A/B comparison. OBS BLANK
render resubmitted as job 13558872 with cmf_flux_param_overrides
MIN_CF=2.0D-3 after 13557228 stopped on
GLOBAL_CMF_GRID_UNDERCOVERAGE (red) requiring 2.9313e-3 below the
3.0e-3 grid edge; model VADAT MIN_CF unchanged.

## Candidate A/B result, job 13559005 (22-Nov-2026)

All six focused batteries pass under the candidate build:
test_knop_oneway, test_knop_lagged_field, test_gs_drift_diag,
test_obs_nonmono_dispatch, test_knop_replay_p2, test_knop_ray
(each rc=0). Original replay exe sha256
189d6fccddee6da386b925da925fbd790ba3844469491129093e92177e1f13df;
candidate exe sha256
bb74d5aebad6ff24fa729b612d8775712bfbd76ffa15f2e86d496420b904525f.
Same decks, base mode, p2=0:

| Ray | Original (ordered) | Candidate (dispatched) |
|---:|---|---|
| 81 | res 1.887e-2, IXIT 3, npo 22 | resid_arg 5.805e-21, class 0, frequency_rays 1 |
| 82 | res 4.485e-2, IXIT 3, npo 26 | resid_arg 4.946e-21, class 0, frequency_rays 1 |
| 84 | res 1.684e-3, IXIT 1, npo 3 | resid_arg 6.682e-23, class 0, frequency_rays 1 |

The decks that stall under the ordered march reach machine-precision
residuals under the candidate dispatch with zero bad cells. The
replay-vs-production field difference is unchanged (2.8e-4, 2.6e-4,
3.8e-8): the candidate converges to a fixed point production never
reached, so trajectory agreement is not expected and is not a
fidelity regression. No negative-GAM deck exists in the armed set,
so the ordered-retention path is verified only by the manufactured
mixed-sign regression (gnego/gzb/mixed arms pass); a live ray-83
capture remains the missing production evidence. Dispatch candidate
status: production-object A/B PASS; live-model validation open.
