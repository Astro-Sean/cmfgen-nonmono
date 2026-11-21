# JFNK deep audit TODO (baseline c992efd7, branch main)

Task states: NOT STARTED / INVESTIGATING / IMPLEMENTED / TESTED /
VERIFIED. Nothing is VERIFIED on compilation alone.

Baseline: HEAD c992efd7 (post KNOP zero-GAM donor fix). Tree clean;
51 untracked run products left alone.

## T1 JVP finite-difference reliability -- INVESTIGATING
Instrumentation already present in JFNK_JVP_V1
(jfnk_cmfgen_interface_v1.f): EPS0=max(sqrt(eps),JFNK_EPS) floor,
clip-boundary h cap, clip-retry x10 shrink with coupled-mode
fail-closed at the EPS0 floor, per-call NaN guard, one-shot JVP
DIAGNOSTIC (worst rows, 8-point h-sweep, baseline re-eval, -h
direction test, repeat eval), JFNK_HSCAN opt-in sweep, JFNK_HADAPT
plateau selector (adjacent-pair 30% rule, mode-keyed re-arm),
JFNK_LINTEST repeat/homogeneity/additivity probes, JFNK_BA_JVP
analytic band matvec (default F).
Open measurements: effective h per residual mode on a production
model; repeatability/homogeneity/additivity numbers; plateau location
for frozen vs coupled.

## T2 LU preconditioner quality -- NOT STARTED
T21 (per-depth LU backward error) is a tracked XFAIL on mild/hard/
shock surrogate specs (cond ~1e23..1e39).

## T3 Schur failure propagation -- NOT STARTED
jfnk_schur_form_v1 zero-fill on failed block LSQ; JFNK_PRECOND_V1
logs but continues on nonzero LU_INFO from schur_apply.

## T4 BA-JVP coordinate consistency -- NOT STARTED
Analytic band matvec applies BA to raw v; correct vs log-coordinate
derivative only at x=0. Comment claims jfnk_x=0 at admission.

## T5 clip-crossing JVP policy -- INVESTIGATING
Coupled mode fails closed at the EPS0 floor; frozen mode keeps the
last attempt with a note. Check the frozen path.

## T6 non-monotone line-search ring buffer -- NOT STARTED

## T7 residual purity / state restore -- NOT STARTED

## T8 end-to-end Newton/GMRES matrix -- NOT STARTED

## T1 measurements (knop_mom smoke model, it=15, x=0 baseline, x=0)

Probe config: DO_JFNK=T, ANYITER=T, LINTEST/HSCAN/HADAPT=T, LU pc,
RSCALE=EUCLIDEAN, MAX_EVALS=60->400 (VADAT read timing unclear).
OPPROBE battery fires at admission under HSCAN; LINTEST + JVP
DIAGNOSTIC inside first post-factor JVP.

Frozen-J KNOP (admitted it=15):
- REP baseline re-eval: df=0.0 bitwise (deterministic eval).
- heff = 1.4901e-08 = sqrt(eps_mach); no caps/retries/clips on the
  seven resolved directions (rnd1,rnd2,pops,ed,k1,kmid,kworst).
- eFD (FD vs BA_diag*v): pops 7.9e-2, rnd1 3.6e-1, rnd2 2.0e-1,
  ed 3.7e-1, k1 1.66, kworst 8.2e-2, kmid 1.2e-5.
- h-sweep ef FLAT across 4 decades on most dirs -> FD-vs-BA gap is
  h-invariant: the discrepancy is structural (BA reference lacks
  terms the true derivative carries), not FD error. eFC==eFD.
- kworst: ef jumps to ~2.8e11 at mult>=1 and decays ~1/h -- a clip
  crossing at larger h; production h stayed below it.
- DEFECT FOUND+FIXED: probe evals consumed the 60-eval budget;
  battery hit cap mid-flight (blk3..xdep fail=1 = sentinel JVPs;
  ADD d_op=NaN). Added JFNK_EVAL_BYPASS brackets around the OPPROBE
  battery (jfnk_solve_corr_v1) and the in-JVP diagnostic + LINTEST
  eval blocks (jfnk_cmfgen_interface_v1), restored on all exits.
- Residual unreachable-by-population-step gate skipped Newton this
  admission (radiation-pinned rows) -- independent of the budget hit.

Coupled KNOP (admitted it=14):
- REP: df_rel ~2.1-2.6e-16 per repeat -- coupled eval deterministic
  at machine precision (iterative sweep is repeatable).
- Battery ~100 evals at ~4-5 min/eval; harvesting progressively.

SOB modes: knop_mom model is KNOP-based; SOB variants need a flag
flip not yet attempted -- deferred, recorded as a limitation.
