# Future fixes — Knop non-monotonic solver

Status as of 25-Sep-2026 (latest commits 069ed05 + c327942; MOM-path
work validated on TypeIIn_10d_DEV7_KNOP_TRIDIAG_gatebulk job
11969828, FG-path pack work on TypeIIn_10d_DEV7_KNOP_ES 11962582).

The implementation is correct (no bugs found in the audit). The
coefficient-pack + log-hoist changes are in and give ~3-8% on the
refresh. Profiling shows the refresh sweeps are memory/latency-bound
with all threads saturated, so further in-code wins are modest; the
largest remaining lever is hardware-level parallelism.

25-Sep-2026 MOM-path update (USE_LAM_ES=F): the refresh stall exit,
repaired-cell quarantine (MOM_RF_KEEPNEG), and the sweep-drift BA
gate (MOM_GATE_BULK) are in and exercised -- the BA freeze on
TypeIIn 10d now releases (gatebulk opened at it25) and the run
cycles Lambda<->BA like the T path.  Live descending GS reads
(MOM_RF_GS_DESC) were tested and found neutral-to-worse; kept
available, default F.

30-Sep-2026 stall-audit update: the hysteresis, bracket bounds,
quarantine-age, BA_NREJ_LIM, GS_STATE/KNOP_STALL_MAX and
MOM_TRIAL_MIN fixes (commits a9ab6b7..7b07fdf) change control flow
only; none of them alter the R1-R4 cost picture below.

30-Sep-2026 evening update: the accepted-trial refund
(N_FROZEN -> MOM_TRIAL_MIN on accepted releases), GS_TRIAL_LIM /
n_trial_acc, the PSIPREV positivity fix, the dead CHI_HAT loop
removal and the diagnostic-table formatting are likewise
control-flow/diagnostic only -- R1-R4 stand as written.

Items below are ordered by expected impact. Each needs its own
build + controlled model run through scripts/run_model_grid.py before
being called done.

## R1 — cpus_per_task 4 -> 16 (ADOPTED 30-Sep-2026 as the CMFGEN production default)

FG_KNOP_REFRESH_V1 parallelizes over KN_NP rays with SCHEDULE(DYNAMIC);
during the refresh all 4 OMP threads are ~99% busy, so wall time scales
with thread count until the tail (one ray per thread) dominates. The
refresh is ~640-750 s of a ~930-1000 s great iteration, so a 16-CPU job
should save ~400-500 s/iter. KN_MAX_THREADS=16 already caps correctly.

Adopted: `run_model_grid.py` now defaults CMFGEN submission paths
(INIT/edd, FINAL, J-release, OBS) to `CMFGEN_PROD_CPUS=16`;
non-CMFGEN codes keep the generic 4-CPU fallback and debug stays at
8 (interactive QOS caps at 16).  First exercised run: FINAL
12043642.  Watch refresh wall time per pass vs the old ~70-100 s
to confirm the scaling.

## R2 — MOM-path coefficient caching (mom_j_cmf_nonmono_v1.f)

MOM_KNOP_REFRESH_V1 rebuilds and refactors the whole tridiagonal system
and re-interpolates opacities for every (L, pass); only the RHS changes
per pass. Caching the matrix coefficients per frequency (same pattern
as the FG pack) should give a larger relative win on the MOM path since
the factorization is the dominant cost there.  With the stall exit
(25-Sep-2026) the refresh now costs ~70-100 s/iter instead of ~420 s,
so the absolute saving is smaller than originally estimated -- still
worth it if the refresh passes grow again on a converging field.

## R3 — SEED_FLD / WRITE_STORE strided store access
(fg_j_cmf_nonmono_v1.f:4870-4874, 5374-5397)

Both loops index I_*_DEPTH_NEXT over LS2, the slowest store dimension —
a ~5 MB strided walk per element group. Reordering the loops (or
transposing the access) makes the walk contiguous. Small (~few % of
refresh) but trivially safe: pure loop restructure, no math change.

## R4 — KN_CHI_COEF / KN_ETA_COEF precision (fg_j_cmf_nonmono_v1.f:523-524)

The monotone-cubic opacity coefficients are ~1 GB in double precision
(NCF x NI_MAX x 4 x 2). Single precision halves the memory and the
stream traffic in BUILD_FREQ's opacity evals. NOT bit-identical —
changes interpolated chi/eta at the 1e-7 level. Needs a convergence
A/B before adopting; likely harmless but must be shown, not assumed.

## R5 — Frequency-sweep OMP tail (unverified)

One phase-mixed perf window showed ~29% of samples in kmp_wait —
probably the tail of the in-sweep FG_J_CMF frequency loop, not the
refresh (threads are saturated inside the refresh). If confirmed, the
sweep schedule/work split is a candidate; profile a clean sweep-only
window first.

## Deliberately not done

- Tightening the FG refresh exit criteria (passes stall at ~0.1
  residual vs KN_RAY_TOL=2e-3 — the documented -1-eigenvalue limit
  cycle).  Fewer passes would cut wall directly but changes how
  converged the I_NEXT store is; a physics/convergence decision, not
  a safe optimization.  (The MOM refresh did get a same-direction
  stall exit on 25-Sep-2026 — its limit cycle is unconditional, so
  exiting early only saves wall time.)
- Single-precision FLD or the I_*_DEPTH_NEXT store: same reasoning —
  changes the iterated field's rounding.
- Removing the pack replay copies by pointing the SRC_PRIME/SLOPES/
  MARCH consumers at PK_*(:,L) directly: saves ~64 B/cell/pass of
  write traffic (~10-15% of refresh maybe); left as an option if the
  pack is kept long-term.

## N-series: numerical-method candidates (03-Oct-2026 review)

Review of the non-monotonic stack for places where a newer or better
conditioned method could replace the legacy approach.  None are
implemented; each needs its own build + controlled run before being
adopted.  Ordered by expected convergence impact, not cost.

### N1 - Adaptive Anderson/GMRES depth on the SE iteration

NG acceleration and the JFNK GMRES restart (JFNK_M=10 default, hyb 30)
use fixed history depth.  The measured SOB limit cycle (sec 17.4,
[[convergence_issues]]) is exactly the regime where depth-adaptive
Anderson (restart on stagnation of the linearized residual, e.g.
Anderson-acceleration depth control of Walker & Ni 2011) is known to
help.  Touch point: subs/solve_for_pops.f JFNK driver; the machinery
for restart already exists via JFNK_MAX_RESTART.

### N2 - JVP finite-difference step control

The coupled JFNK JVP uses a single noise-scaled FD step (JFNK_EPS,
Pernice & Walker 1998 form).  On cells where the Knop chi_hat
response is near-discontinuous (sigma ~ -1 sign flips in the upwind
selection), a one-sided FD mismeasures dJ/dn.  Options: adaptive
two-sided fallback where the directional derivative fails a
linearity probe (evaluate at eps and eps/2, compare to 2x), or a
complex-step derivative where the code structure allows.  Touch
point: JFNK_RAD_SWEEP_V1 JVP evaluations in cmfgen_sub.f.

### N3 - Frequency remapping on |GAM| instead of fixed Ceulemans count

DO_NONMONO_FREQ_MATCH redistributes frequencies by a fixed N_CAP
count.  An adaptive remap keyed on the local Doppler width
(|GAM|/dln_nu vs line width) would put points where the upwind
source term actually needs them; Ceulemans+ 2025 (arXiv:2505.08309)
discusses the residual structure this targets.  Higher risk: changes
the grid the lagged I_NEXT store indexes.

### N4 - Positivity-preserving population update

BA steps can drive level populations negative on a bad linearization;
the current fix is rejection/demote after the fact (resacc,
BA_SAT_*).  A cheaper guard is a positivity-preserving step-length
(the standard "fraction to boundary" of interior-point methods) on
the BA correction before the residual test.  Touch point:
solveba_v13.f correction application; complements, not replaces,
resacc.

### N5 - Preconditioner reuse diagnostics for JFNK

JFNK_USE_BANDED exists (untested, needs NUM_BNDS>=3).  Before tuning
it, add a cheap diagnostic: report the condition-proxy already
implied by the equilibrated-block LU row/col scales (the scale
min/max are printed - persist them per iteration in ITER_HISTORY) so
preconditioner benefit can be measured rather than guessed.

### N6 - Formal/moment/ALO cross-check mode

do_var_cont.f carries the non-monotonic ALO path.  A diagnostic
switch that, for one iteration, runs both the FG formal solve and
MOM_J_CMF_NONMONO_V1 and reports per-depth |J_FG - J_MOM| would give
an independent consistency probe for the chi_hat discretization -
the Open-issues MOM-path sign-flip question is currently judged only
by convergence behaviour.

Each N-item is speculative; implement singly, verify on the m15
2.82d KNOP config, and document PASS/FAIL here.
