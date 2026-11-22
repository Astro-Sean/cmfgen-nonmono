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

Verdict (22-Nov-2026): PARTIAL. Production non-convergence confirmed;
fresh captures found; diagnostic census defect identified; bounded
P2-off replay audit submitted. P2 remains disabled. No numerical cure
or production-operator P2 success is claimed. No GitHub push authorized.
