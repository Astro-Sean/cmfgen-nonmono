# JFNK Solver Audit — performance and robustness on non-monotonic m15

Date: 03-Oct-2026.  Tree: `cmfgen_jfnk` @ `545d4c5` (devin/jfnk-bringup).
Runs audited: `m15_2p82d_jfnk_opt` (job 12066716+restarts), `m15_2p82d_jfnk_cpl`,
`m15_2p82d_jfnk_seq/seq2`.  Config: `configs/dessart_m15_2p82d_JFNK_opt.yaml`.

Scope: is the JFNK implementation optimal for difficult non-monotonic systems,
is it hardened enough to solve them correctly, and is the
"Lambda-until-close then JFNK" hybrid the right strategy (user's proposal).
Method: source read of the three JFNK files + host sweep callback, OUTGEN /
JFNK_EVAL.log evidence, literature check (Arramy+ 2024, Janett+ 2024,
Knoll & Keyes 2004, Hauschildt & Storzer 1994, NITSOL/Pernice & Walker,
NKA/Anderson-acceleration literature for radiative transfer).

----------------------------------------------------------------------
## 1. Verdict

The user's proposed strategy is correct and is what the literature
recommends; the code already implements most of it.  JFNK is NOT yet a
net win on this model for three compounding reasons:

1. **Cost asymmetry.**  In coupled mode every residual/JVP eval re-runs the
   full Knop formal sweep (~250-400 s on m15, measured below).  A single
   Newton step needs ~10-40 evals; one activation measured 3.2 h and did
   no useful work (convergence_issues §16.1).  JFNK can only pay for
   itself if one activation replaces >= ~30 Lambda iterations.

2. **Basin of attraction.**  Newton is only worth firing near the root;
   the measured failure at it11 (|F|_w=2.7e13, linearity probe diverges
   at unit perturbation, FD-JVP flat at ~1e32 for eps=1e-4..1e-8, BA-vs-FD
   rel diff ~1e15) is the classic too-early admission signature.
   `JFNK_COR_GATE` (added 162d39e) now defers on demanded correction and
   is working as designed in the opt run (demand decaying ~5x/it:
   1.5e7 -> 4.3e5 -> 7.8e4 at it31).  It is necessary but not sufficient:
   see finding C1.

3. **The inner problem is not the bottleneck.**  The outer loop stalls on
   the J<->n fixed-point map (field drift 4.6e-2, luminosity ratio 0.30
   at it30), not on SE-solve accuracy.  Newton-polishing SE at a moving
   field (frozen mode) or at an immature lagged store (coupled mode)
   cannot converge the real problem faster than the field matures.
   JFNK's value is confined to the endgame, which is exactly the
   Lambda-then-Newton hybrid.

So: the architecture is defensible, the admission logic is the right
idea, but the per-activation cost model, the deterministic-eval problem,
and the single-metric gate are not yet at production quality.

----------------------------------------------------------------------
## 2. Findings

### A. Correctness (no showstoppers found in the solver core; two real issues)

A1. **Coupled-residual JVP has a noise floor, not a truncation floor.**
    `JFNK_JVP_V1` uses forward difference h = sqrt(eps_mach)*max(1,||x||)/||v||
    (jfnk_cmfgen_interface_v1.f:2020-2037).  That is the NITSOL/Knoll-Keyes
    optimal step for a *deterministic, smooth* F.  Here F is produced by an
    *iterative* sweep: FG_J_CMF inner-iteration counts are data-dependent
    (OUTGEN prints "Average number of calls to FG_J_CMF" = 5.6-6.3 varying
    per eval), plus negJ repair, quarantined cells, boundary branch
    switching, Eddington f-loop early bail (convergence_issues §10.2).
    F therefore carries iteration noise at perhaps 1e-4..1e-3 relative.
    JVP noise ~ noise_F / h ~ 1e4-1e5 x signal at h~1e-8.  The FD-vs-eps
    probe in §16.1 (||Jv|| ~ 1e32 at every eps) is consistent with this:
    the difference quotient measures sweep noise/restart effects, not the
    Jacobian.  This is THE structural obstacle for coupled JFNK; a 100x
    better preconditioner will not fix it.  Mitigations (in order of
    preference):
      (a) make the eval sweep deterministic -- fix the FG inner-iteration
          count/tolerance inside evals, freeze repair/quarantine/branch
          decisions to the baseline pass (they are already lagged state),
      (b) use a noise-aware step: h ~ sqrt(noise_F) not sqrt(eps_mach),
          i.e. h*v ~ 1e-2..1e-3 relative (Kelley; Pernice & Walker 2008:
          approximate function values require larger h),
      (c) accept the frozen-J residual for Krylov directions and use
          coupled evals only for trial-step verification.

A2. **Latent state leak: per-frequency quadrature weights are computed at
    perturbed T and not restored.**  In the eval sweep
    (`JFNK_RAD_SWEEP_V1`, cmfgen_sub.f:5047-5945) `QUAD_MULTI_V10` fills
    ATM%WSXzV/dWSXzVdT/WCRXzV at the *eval* T, keyed on `CONT_FREQ`.
    After the eval, `CONT_FREQ` is left at NU(NCF); the next real sweep at
    ML=1 sees NU(1) .NE. CONT_FREQ -> COMPUTE_NEW_CROSS=T -> recomputes,
    so the leak heals at the first frequency.  With JFNK_FIXT=T (current)
    T is never perturbed, so the stale values are identical -- benign now,
    but if JFNK_FIXT is ever set F the first-frequency weights would still
    be recomputed (same key) so this is probably still safe; flagging as
    needs-verification, plus `EDDINGTON` is left at EDD_CONT rather than
    restored to EDD_LINECONT (the outer line section re-sets it; low
    risk).  No action required under current config; note for the free-T
    future.

A3. Things checked and OK:
    - F(0) identity test aborts on NaN and on frozen-mode mismatch;
      coupled mode correctly adopts the eval baseline (lagged-store
      offset is real, not a bug) with a >N/2 spread guard
      (jfnk_solve_corr_v1.f:282-395).
    - Residual restores POPS/ED/T/ATM and all STEQ arrays after each
      eval -- side-effect-free w.r.t. those modules
      (interface:1600-1626).  FREEZE_I_NEXT + ACCESS_F cursor reset +
      JFNK_IN_EVAL prevent lagged-store/record mutation.
    - COMPUTE_BA held false during eval sweeps so SE%BA is not
      double-accumulated (cmfgen_sub.f:5081-5082) -- this was a real bug
      class and is correctly handled.
    - NaN policy: components capped at 1e20 rather than zeroed; JVP
      failure flag propagates to GMRES abort when >10% bad
      (interface:1560-1598, 2168-2192).
    - Log-unknowns (JFNK_LOGX) remove the floor non-differentiability;
      N_NEW_CLIP retry shrinks h when a new constraint boundary is
      crossed (up to 2 retries; each costs a full eval -- see B4).
    - PTC term x*s/dtau with s frozen from baseline |F| and dtau grown
      by tol_prev/tol each Newton step -> dtau -> 1e8 as it converges;
      no fixed-point bias at convergence (solver:1314-1317).
    - EW forcing term is in preconditioned units under left prec
      (solver:535-557); right prec is forcibly disabled because the
      hookstep norm is not valid there -- correct choice.
    - Fallback hierarchy (catastrophic -> conservative safety step;
      stagnation -> MAX_LAM_COR step; JFNK_MAX_STAGNATE self-disable)
      is sound and preserved the contraction in the seq run.

### B. Performance

B1. **Per-eval cost is the dominant term.**  Baseline great iteration
    ~250-400 s; coupled eval = one full continuum formal sweep, similar
    cost.  Measured: cpl run it36 took 8658 s vs ~240 s baseline (the
    F(0)-abort epoch, before the gate).  JFNK_EVAL.log confirms each
    eval runs the full NCF frequency loop with the same inner-iteration
    chatter as a real sweep.

B2. **~13 wasted evals per first activation.**  On the first JFNK fire:
    F(0) identity test (1 eval, needed), NEWTON entry f() re-eval (1
    eval, REDUNDANT -- identical input to the test eval just made),
    COMP_DIAGNOSTIC (4 linearity + 5 eps-sweep + 1 FD-vs-BA = 10
    evals), TEMP_COL_AUDIT (2 evals).  ~55-75 min of m15 walltime before
    GMRES does anything.  The diagnostic block is latched by
    DIAG_DONE/AUDIT_DONE (once per process) but on a multi-restart
    cluster job it re-fires once per resubmit.  These were worth it once
    (they diagnosed the basin problem); they are now pure overhead.
    Recommend gating COMP_DIAGNOSTIC/TEMP_COL_AUDIT behind a new VADAT
    flag (e.g. JFNK_DEEPDIAG, default F) and passing F0 into
    NEWTON_HOOK_V1 to skip the duplicate entry eval.

B3. **Retry storms multiply evals.**  Failure paths each cost full
    sweeps: non-descent retry tightens eta and re-runs GMRES; the clip
    retry runs up to 3 evals per JVP; line search allows JFNK_LS_MAX=20
    evals; the driver allows JFNK_MAX_RESTART=8 baseline/precond
    rebuilds each followed by a fresh Newton (F(0) + GMRES).  Worst case
    per activation is ~100+ sweeps (~8-11 h).  §16.1 measured the
    abortive version of this at 3.2 h.  Recommend a hard eval BUDGET per
    activation (new counter, e.g. JFNK_MAX_EVALS ~ 30) that forces the
    best-iterate fallback path when exhausted; and in coupled mode drop
    JFNK_LS_MAX to ~4 and JFNK_MAX_RESTART to ~2 (the yaml already sets
    MAX_RESTART=2 -- good).

B4. **JVP clip-retry cost.**  Each N_NEW_CLIP>0 retry is another full
    eval (interface:2071-2082).  Under log-unknowns the clip boundary is
    the |x|=40 exp-overflow guard, which should almost never trigger; if
    OUTGEN shows clip retries firing routinely, that is a symptom the
    step is too big (TR should have caught it) -- worth monitoring in
    the counters, no code change.

B5. **GMRES m=10 is small for this operator.**  The preconditioned
    operator is far from identity in coupled mode (depth-block BA
    captures only the local part of a globally-coupled J).  Literature
    uses m=20-50 for NLTE JFNK; restarted GMRES(10) stagnates or wastes
    restarts.  Memory is trivial (N~1.3e5, each Krylov vector ~1 MB).
    Recommend JFNK_M=30-50 for coupled; for frozen mode m=10 is fine.

B6. **EW bounds.**  eta_min=3e-4 forces near-exact inner solves --
    oversolving costs JVPs, and in coupled mode the JVP noise floor
    (A1) makes eta below ~1e-3 meaningless anyway.  Recommend
    eta_min ~ 1e-3..1e-2 coupled.  eta_max=0.1 is reasonable.

B7. **OUTGEN I/O flood mostly fixed**: eval chatter goes to
    JFNK_EVAL.log (cmfgen_sub.f:5092), but the FG_J_CMF "average calls"
    line still lands in OUTGEN per frequency during the *real* sweep;
    fine.

### C. Robustness / admission

C1. **The gate is single-metric and instantaneous.**  JFNK_COR_GATE
    checks only the previous call's raw demanded MAXCH
    (solve_for_pops.f:1012-1019).  It does not require the field to be
    mature (GS_FG_IMMATURE only gates the Knop flag), does not look at
    luminosity balance, and has no hysteresis -- one lucky low-MAXCH
    Lambda step admits JFNK on a field that is still drifting.
    Recommendation: admission = (MAXCH < gate) AND (field drift < tol,
    e.g. 1e-3 or a new VADAT JFNK_FIELD_GATE) AND (sub-gate for K
    consecutive calls, e.g. K=2-3).  All metrics already exist in the
    mode controller; this is a small solve_for_pops.f change.

C2. **Gate value.**  COR_GATE=100 means admit when demanded correction
    < 100%.  The basin evidence (it11) says even ~50% demand was too
    early in the old tests; with the coupled eval noise (A1) the usable
    basin may be smaller.  COR_GATE ~ 20-50 with the hysteresis in C1 is
    a more defensible starting point; the gate should be re-derived from
    the FD-vs-BA diagnostic once JFNK first passes the demand gate --
    the diagnostic already prints the rel-diff that quantifies being
    inside the linear regime (aim: rel diff < ~0.1).

C3. **Self-disable on repeated failure exists** (JFNK_MAX_STAGNATE=3)
    -- good; keep.  Combined with C1/B3 the worst-case cost of a bad
    admission is bounded.

C4. **Preconditioner is the wrong operator for coupled mode -- inherent
    but improvable.**  Depth-block LU of BA is exact for the frozen-J
    residual and only approximate for the coupled one (misses dJ/dn
    nonlocal response).  Options per literature: (i) banded BA
    (JFNK_USE_BANDED=T needs NUM_BNDS>1, i.e. SOL_METH=TRI instead of
    DIAG) -- Hauschildt & Storzer 1994 found optimum ALI bandwidth
    ~5-15, so TRI(3) is a minimum worth testing; note the pipeline
    currently echoes "[JFNK] SOL_METH forced to DIAG (JFNK preconditioner
    is DIAG-only)" which both disables the banded path AND weakens the
    ordinary Lambda-phase BA; (ii) physics-based preconditioner --
    apply one Lambda/operator-split step as M^-1 (the Janett+ 2024
    nested-Krylov / lightweight-model-as-preconditioner idea; also the
    transport-community "source-iteration preconditioner"); that is
    deeper work, see section 4.

C5. **Frozen mode is a wasted asset.**  JFNK_COUPLED=F gives F(x) =
    STEQ_coll(perturbed) + stored radiative + CMAT_RAD*x -- no sweep per
    eval (~seconds, not ~5 min), smooth, and BA is nearly the exact
    Jacobian of that residual.  Newton-solving the frozen-J SE problem
    to convergence is a strict upgrade over one linearized SOLVEBA step
    per great iteration, and it exercises the entire Newton machinery at
    ~zero cost.  It also removes inner-solve noise from the outer
    J<->n map (under-converged SE solves add error to the fixed-point
    iteration).  Recommended sequencing: keep COUPLED=F as the first
    production JFNK stage; promote to coupled only at the endgame.

### D. Minor / hygiene

D1. `JFNK_MAX_RESTART` default 8 at interface layer vs yaml 2 -- make
    the default match the validated value.
D2. Redundant F(0)/entry eval -- B2.
D3. Diagnostics are once-per-process, not once-per-model; on resubmit
    they re-run -- B2.
D4. `jfnk_max_gmres_its` cap exists but is 0 (off) by default -- set a
    real cap for coupled (e.g. 40-60 JVPs per Newton step).

----------------------------------------------------------------------
## 3. Cost model (why JFNK hasn't won yet)

m15 2.82d: N = NT*ND ~ 1.3e5 unknowns; sweep ~250-400 s.

  Lambda iteration:     ~1 sweep + SE solves ~ 400-1000 s per iteration.
  Frozen-J JFNK eval:   ~1-5 s (SE assembly only, no sweep).
  Coupled JFNK eval:    ~250-400 s.

  Frozen-J JFNK activation ~= (1 + nits*(1+m) + LS) * ~3 s ~= ~2-5 min.
    Affordable always.  BA is the right preconditioner -> GMRES should
    converge in < m iterations per Newton step.

  Coupled JFNK activation ~= (1 + nits*(1+m) + LS + retries) * 300 s.
    At nits=10, m=10: ~110 evals ~ 9 h.  Only viable if nits is small
    (basin entry + good preconditioner -> nits ~ 3-6, m-its ~ 5-10
    each -> ~30-50 evals ~ 3-5 h) AND it replaces >>50 Lambda
    iterations -- i.e. a true endgame accelerator only.

Conclusion: the optimal schedule is the user's hybrid, sharpened:

  Phase 0  Lambda/Knop iterations until field mature AND demand small
           (gates: GS_FG_IMMATURE, COR_GATE + drift + hysteresis).
  Phase 1  Frozen-J JFNK (COUPLED=F) -- Newton-polish SE at current
           field; cheap, validates machinery, cleans the outer map.
           Can run for many activations as a SOLVEBA replacement.
  Phase 2  Coupled JFNK only when |F|_w is small AND field drift < eps;
           with deterministic evals (A1a) or noise-aware h (A1b),
           m=30-50, eval budget ~30-50 per activation, best-iterate
           fallback (already implemented).

Literature cross-check: Arramy+ 2024 get <=2x over RH92 in eval count
with a *local* approximate-operator preconditioner on 6-level atoms;
Janett+ 2024 get 10-20x iteration reduction from physics-based
(lightweight-model) preconditioners; Hauschildt & Storzer 1994: optimal
ALI bandwidth 5-15 (ours = diagonal(1), banded(3) available);
Pernice & Walker 2008 / NITSOL: FD step must reflect function noise,
not just machine eps (A1b); transport community uses the fixed-point
sweep itself as the Krylov preconditioner and NKA/Anderson to
accelerate the outer map -- an alternative worth a feasibility note if
coupled JFNK stays expensive (it needs no JVPs at all).

----------------------------------------------------------------------
## 4. Recommendations (ranked by cost/benefit)

R1 (config, immediate):  COR_GATE -> 20-50; add field-drift + consecutive
    sub-gate hysteresis to the admission test (needs small
    solve_for_pops.f edit -- ~20 lines, flag-gated on existing keys).
R2 (config, immediate):  try JFNK_COUPLED=F on the next run -- validates
    the Newton stack cheaply and gives a production-usable Phase-1
    accelerator even if coupled never fires.
R3 (source, small):  JFNK_DEEPDIAG flag (default F) gating
    COMP_DIAGNOSTIC/TEMP_COL_AUDIT; pass F(0) result into NEWTON to
    remove the duplicate eval.
R4 (source, small):  JFNK eval budget counter -> forced best-iterate
    fallback; JFNK_LS_MAX cap for coupled; real jfnk_max_gmres_its cap.
R5 (config):  JFNK_M 10 -> 30-50 for coupled; eta_min -> 1e-3 (needs a
    VADAT key or code change to the 3e-4 module default).
R6 (source, medium):  deterministic eval sweep -- fixed FG inner-count
    or tighter eval tolerance; freeze repair/quarantine/branch
    decisions during evals.  If not feasible, JFNK_EPS_JVP -> ~1e-3 as
    noise-aware FD step (VADAT key already exists but is currently
    ignored in favour of the machine-eps rule -- re-enable the override).
R7 (config/source):  SOL_METH=TRI + JFNK_USE_BANDED=T to exercise the
    banded preconditioner; check why the pipeline forces DIAG.
R8 (research):  physics preconditioner for coupled mode -- one
    Lambda/operator-split step as M^-1 (Janett nested-Krylov analogue),
    or NKA/Anderson acceleration of the outer pops sequence (zero JVPs,
    reuses history; needs safeguarding on non-contractive maps).
R9 (verify):  free-T path state-leak check (A2) before any JFNK_FIXT=F
    experiment.

Nothing above requires restructuring the solver: every change is
flag-gated or config-level except R6/R8.

----------------------------------------------------------------------
## 5. Honest status

- Frozen-J JFNK path: believed functional; never exercised in anger --
  needs a validation run (R2).
- Coupled JFNK: CORRECTLY deferred in the current run; has never
  produced an accepted Newton step on this model.  Not yet validated.
- Lambda fallback after JFNK failure: validated (seq run, res kept
  contracting after JFNK_DONE).
- Gate: functioning (opt run it11-31 all deferred; demand decaying
  ~5x/it).

Status: PASS WITH LIMITATIONS for the admission/fallback architecture;
FAIL for coupled-JFNK production readiness (cost + eval determinism);
BLOCKED for endgame validation (demand still ~7.8e4 > gate at it31).

----------------------------------------------------------------------
## 6. Residual anatomy + partial-coupled-eval audit (12-Oct-2026)

Production reach-gate data (it534-537, zlaw14_lsq jobs) showed the
frozen-J residual is exactly left-null of BA: REACH ~ -1e-28,
|Jz|/|F| ~ 3e-15, cos(F,Jz) ~ -1e-14.  No population step can reduce
it, so every GMRES/JVP eval was wasted (~90 per disabled cycle).

Two changes shipped since Section 5:

- JFNK_REACH_GATE (VADAT, default 0): precheck z=M^-1 F, Jz=BA*z,
  REACH=(2<F,Jz>-<Jz,Jz>)/<F,F>; skip Newton below the gate.  Skip
  sets info=-2 -> Lambda fallback, exempt from the stagnation count.
- JFNK_REARM_EVERY (VADAT, default 10, was hardcoded): after
  JFNK_MAX_STAGNATE failures JFNK pauses for REARM_EVERY x 2^cycle
  (cap x8) admitted calls, then re-arms.  Cycle resets on any
  committed update.  JFNK is never permanently disabled.
- JFNK_ROW_DIAG (VADAT, default 10 rows, 0=off): residual anatomy at
  each admission -- per-ion |F|^2 rollup and top rows labeled
  ion/level/depth with coll vs stored-rad split, Jz_i, BA row norm.
  Row map: J in [EQXzV, EQXzV+NXzV-1] -> ION_ID(ID); J=NT-1 charge,
  J=NT energy.  Classifies pinned rows as radiation-lag
  (|rad|>>|coll|, Jz~0) vs population-driven; decides whether any
  population-space residual can reach them.

Partial-coupled-eval verdict (audit, not implemented):

- A SOB-mode sweep callback replaying only the population-dependent
  Sobolev line-rate update (L_STAR/U_STAR -> CHIL_MAT/ETAL_MAT ->
  SOBJBAR_SIM, continuum frozen via stored EDDFACTOR records /
  GET_J_CHI) is feasible: SOBJBAR_SIM is deterministic local ray
  solves (no iteration noise, clean FD JVPs), ~10-25% of an outer
  iteration per eval.  Mode-dispatched via JFNK_SWEEP_CB; KNOP keeps
  JFNK_RAD_SWEEP_V1.
- BUT it does not fix first-order reachability: the reach check was
  measured against fresh BA, which already contains the linearized
  SOB rate derivatives (VB_SIM/VC_SIM feed BA).  A fresh rate re-eval
  gives Jacobian ~= BA -> same left-null geometry.  What it buys:
  residual consistency (true F has a root; the stored+linearized
  model may not), CMAT_RAD staleness, and finite-step response where
  escape probabilities saturate (d(beta)/dn ~ 0 at the tangent).
- Decision gate: if JFNK_ROW_DIAG shows the pinned rows are
  line-rate rows, implement the SOB partial eval; if continuum-
  driven (photoionization lag), population-space formulation is
  exhausted and the productive direction is accelerating the outer
  J<->n map (NKA/Anderson on the fixed-point residual; no JVPs).
