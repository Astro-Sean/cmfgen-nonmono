# JFNK-CPL development proposal - coupled population-radiation Newton

Date: 11-Nov-2026. Status: PROPOSAL (no implementation in this
task). Trigger: ES_ALI_BA_AUDIT.md decision gate - if the X1-X4
continuum-ES arms do not materially improve the SOB tail
contraction, the slow mode lives in the line-rate / inter-ion
coupling that the continuum ALI stack cannot reach, and the
remaining mathematical option is a Newton method on the coupled
statistical-equilibrium + radiative-transfer system.

This document maps the nine required design considerations onto
what already exists in the tree (much of the machinery is built
and audited) versus what still needs development.

## 0. What already exists

The coupled tier is implemented, not hypothetical:

- `JFNK_COUPLED` (VADAT, default F): residual evaluations
  re-run the full opacity -> COMP_J_BLANK -> radiative-rate
  chain at perturbed populations via `JFNK_SWEEP_CB` /
  `JFNK_RAD_SWEEP_V1` (cmfgen_sub.f ~5755). Each eval is a fresh
  formal solution, so the finite-difference JVP probes the full
  dJ/dn coupling - the Arramy, de la Cruz Rodriguez & Leenaarts
  (2024, A&A 690, A12) formulation: Newton on the coupled
  SE+RTE system.
- `JFNK_COUPLED_SOB` (VADAT, default F): admits the coupled tier
  under DO_NONMONO_SOB by also replaying the SECTION='LINE'
  pass (SUB_SOB_LINE_V3 Sobolev rates + line-group continuum
  solve), so the coupled residual is the same equation the
  outer iteration assembles.
- `JFNK_RESIDUAL_V1` / `JFNK_JVP_V1` (jfnk_cmfgen_interface_v1.f):
  audited 03-Nov-2026 (JFNK_PRODUCTION_OPERATOR_INVESTIGATION):
  depth-major masked SE residual, frozen-J and coupled-J modes,
  state save/restore per call, non-finite capping, PTC anchoring.
- Restarted GMRES + hookstep trust region + preconditioner
  selection (JFNK_M, JFNK_MAXGMRES, JFNK_PRECOND, JFNK_RIDGE,
  JFNK_LSQ_LAM, ...).

The gap is admission policy, not machinery: under SOB the
coupled tier is never demand-admitted - it fires only as the
reach-fail fallback when the frozen tier certifies the residual
unreachable (control_variable_mod.f ~2171-2181). It has zero
convergence commits on record; on the np run the frozen tier
was 322-admitted/0-committed pure overhead before being tuned
off. JFNK-CPL development is therefore: make the coupled
operator cheap enough and reliable enough to be the primary
solver for the stiff sector, not a rarely-fired fallback.

## 1. Coupled radiation-population residual definition

Status: IMPLEMENTED for both frozen and coupled modes.
Development needed:

- Residual scaling audit for the coupled mode. The frozen
  residual is linear in the radiative increment (stored
  CMAT_RAD band product); the coupled residual re-solves J
  iteratively, so its row scaling differs per evaluation.
  Rybicki & Hummer (1991) precondition on the ALI operator;
  the coupled Newton needs the SE equations scaled so that
  radiation-dominated rows (continuum ES, line source
  function) and population-dominated rows (collisional,
  ionisation balance) enter GMRES with comparable weight.
  Proposal: row-equilibrate by the diagonal of the frozen-J
  BA matrix (already stored in BAMAT) - a physics-based
  scaling, not a tuned constant.
- The consistency gate (JFNK_CPL_SOB_TOL=1e-5, noise floor
  ~1e-8 observed) already separates "different equation"
  (O(1) row errors) from solver truncation. Keep it; add a
  per-row max-error report so a partially-missing replay
  term is localised rather than just failed.

## 2. Matrix-free Jacobian-vector products

Status: IMPLEMENTED (forward-difference JVP, h from
JFNK_JVP_H-family keys). Development needed:

- JVP cost is the binding constraint: each matvec = one
  coupled residual eval = one full formal sweep (~97% of an
  outer iteration's cost). Arramy et al. report ~2x fewer
  formal evals than Rybicki-Hummer ALI on small atoms; our
  sweep is 8331 frequencies x NUM_LAM passes, far heavier.
  Options, in order of preference:
  a) Reduce eval cost: the coupled eval need not converge J
     to production tolerance - a fixed small number of
     Lambda/ES passes gives a consistent approximate
     operator (inexact Newton; Eisenstat-Walker forcing
     term can ride the existing trust region).
  b) Reuse the formal sweep across GMRES iterations where
     the perturbation is small (perturbation-size cutoff,
     verified against the consistency check).
  c) Physics-based JVP compression: the dJ/dn coupling is
     dominated by the ES-scattering and strong-line
     transitions identified in SOB_SLOW_MODE_ANALYSIS;
     a reduced-frequency JVP (continuum + flagged line
     frequencies) could cut cost 5-10x if the residual
     stays consistent - must be validated, not assumed.
- Rybicki & Hummer (1992, A&A 262, 209) overlapping
  transitions: line overlap is already treated inside the
  Sobolev replay; the JVP inherits it for free. No extra
  work, but the overlap bookkeeping must be inside the
  perturbed sweep (it is - SET_LINE_OPAC is part of the
  continuum loop).

## 3. Consistent residual scaling

Covered in section 1. Additional: the correction vector is
already fractional populations + log-T + electron density
(control_variable_mod ~2160 area, JFNK_RESIDUAL_V1 docs).
Keep fractional units - they are the positivity-friendly
choice (section 4).

## 4. Positivity-preserving population variables

Status: PARTIAL. x_frac is fractional (n = n0*(1+x)); a
Newton step can still drive 1+x < 0. Existing machinery:
JFNK_CLIP_FRAC_MAX, JFNK_POP_FL floor, non-finite capping,
hookstep trust region. Development needed:

- A hard positivity projection inside the residual eval:
  perturbed populations clipped at a small floor (the
  existing POP_MASK_FLOOR machinery) BEFORE the SE
  evaluation, with the clip counted (JFNK_N_NEW_CLIP
  exists). The residual of a clipped point is then the
  residual of the clipped state - consistent, and the
  Jacobian sees the clip as a kink (acceptable: it is
  measure-zero in practice).
- Alternative considered and rejected for now: log-population
  variables make positivity automatic but destroy the
  conservation-row structure (species sums become
  exponentials of sums) and change the fixed point's
  conditioning. Fractional + clip is the smaller change.

## 5. Particle and charge conservation

Status: STRUCTURAL. Conservation is enforced by the SE
equation set itself (conservation rows are equations in
the same system; the BA matrix includes them). The coupled
residual assembles "the same SE equations the outer
iteration uses" (control_variable_mod ~2176), so a
converged coupled Newton step conserves exactly what the
outer iteration conserves - no separate conservation
repair is needed or wanted. Verify: the coupled residual's
conservation rows must be bit-identical in form to the
production ones (the JFNK_CPL_SOB consistency check
covers this - a missing replay term shows as O(1) row
error).

## 6. Krylov preconditioning

Status: WEAK - this is the main numerical gap.

- Frozen tier: JFNK_PRECOND=LU on the stored BA matrix -
  a direct factorisation preconditioner, appropriate for
  the fixed-J operator.
- Coupled tier: no equivalent. The coupled operator's
  Jacobian is the SE Jacobian plus dRTE/dn terms - never
  formed. Options:
  a) Frozen-J BA LU as a RIGHT preconditioner for the
     coupled operator: captures the stiff collisional +
     stored-radiative part; the Krylov iteration handles
     the dJ/dn correction. This is the MALI analogy -
     Rybicki & Hummer's ALI is exactly a physics-based
     preconditioner for the same splitting. Expected to
     be the highest-value, lowest-risk choice.
  b) ES_ALI_BA-amplified diagonal as a left
     preconditioner - cheaper but weaker; the audit showed
     the continuum channel is not the bottleneck, so do
     not expect this alone to fix the tail.
  c) Multigrid/depth-line preconditioning - research
     grade, defer.
- Implementation note: right-preconditioning keeps the
  residual in physical units (no split residual needed);
  GMRES sees M^{-1} J v. The existing JFNK_LSQ_LAM ridge
  and JFNK_RIDGE machinery carry over.

## 7. Non-monotonic formal-solution consistency

Status: INHERITED. The coupled sweep calls the same
COMP_J_BLANK -> formal solver dispatch as production:
KNOP route uses FG_J_CMF_NONMONO_V1 + GET_FG_ALO_DIAG,
V13 route uses FG_J_CMF_V13, moment route uses the MOM
solver - the ES_ALI_BA_AUDIT section 4 accessor audit
applies unchanged. The perturbed sweep therefore sees
the same multi-resonance, frequency-coupled, sign-aware
velocity-gradient physics as the outer iteration. No
monotonic approximation can enter through the JVP.
Caveat to verify in testing: the KNOP lagged I_NEXT/J_NEXT
stores are frozen during evals (FREEZE_I_NEXT) - the
coupled residual probes dJ/dn at FIXED lagged field, so
the Newton operator is the SE+formal-solve map with the
lagged store as a parameter, not the fully simultaneous
map. That is the correct operator for the outer
iteration's own splitting (the lagged store IS the
iteration's state), but the proposal must state this
explicitly: JFNK-CPL accelerates the population-radiation
fixed point of the lagged-field scheme, it does not
eliminate the lag itself. If the slow mode is the lag
(not the populations), JFNK-CPL will not fix it - the X-arm
data decides this.

## 8. Line-rate and continuum-rate coupling

Status: IMPLEMENTED via JFNK_COUPLED_SOB's SECTION='LINE'
replay. This is the decisive feature vs the continuum
ES_ALI stack: the coupled JVP includes d(Sobolev line
rate)/dn through the same SUB_SOB_LINE_V3 + line-group
continuum solve the outer iteration runs. The
SOB_SLOW_MODE_ANALYSIS species histogram (Al V, O VII,
C V, Cr VI dominating max corrections) points exactly
here. Development needed: none in the replay itself;
validation that the replay is complete (the consistency
check's per-row error report from section 1).

## 9. Newton step globalization and trust regions

Status: IMPLEMENTED (hookstep trust region, JFNK_TR_*
family, correction caps, clip counters, reach gates).
Development needed:

- The coupled eval is ~1 outer-iteration cost; a rejected
  or truncated Newton step wastes a full sweep. The
  existing JFNK_CPL_SOB_MAX_EVALS cap and reach-fail
  denial logic are the right skeleton; add an
  eval-budget-aware forcing term (section 2a) so cheap
  inexact steps are preferred early.
- Watch for the np-run failure mode: 322 admissions,
  0 commits = the trust region never accepted. Before
  production, the coupled tier needs a commit-rate
  diagnostic in ITER_LEDGER (partially present via
  JFNK eval counters) and a defined abort criterion
  (e.g. N consecutive rejected steps -> demote to frozen
  tier for K iterations, already partially present via
  rung logic).

## 10. Development sequence (proposed)

1. X-arm decision gate (ES_ALI_BA_TODO X4): if Case D,
   proceed.
2. Preconditioner spike: frozen-J BA LU right-preconditioner
   for the coupled operator on a single-frequency or
   single-species test problem (unit_tests harness),
   measure GMRES iterations vs unpreconditioned.
3. Consistency instrumentation: per-row error report in
   the JFNK_CPL_SOB check; run on the m15 model at the
   warm-started state, confirm O(1e-8) noise floor holds
   with the Sobolev replay active.
4. Reduced-cost JVP prototype: fixed-pass coupled eval
   (section 2a) behind a new VADAT key, default = current
   behaviour; measure eval cost vs Newton convergence
   on the warm-started DIE-off state.
5. Admission policy: allow demand-admission of the
   coupled-SOB tier when the frozen tier's commit rate
   stalls (the observed plateau signature), keeping the
   reach-fail fallback. Requires ITER_LEDGER evidence
   that frozen-tier stalls correlate with the slow tail.
6. Production candidate: full X-arm-style A/B on the
   DIE-off model vs X1 baseline; same fixed-point
   agreement criterion as the ES_ALI matrix.

Each step lands behind VADAT keys defaulting to current
behaviour; the frozen-J path and stock monotonic CMFGEN
remain bit-identical throughout.

## 11. Literature mapping

- Rybicki & Hummer (1991, A&A 245, 171): ALI =
  physics-based preconditioner for the same
  population-radiation splitting. Applies directly: our
  frozen-J BA LU is the analogous preconditioner; the
  coupled Newton is the "complete operator" limit they
  motivate. Their diagonal-ALI is what ES_ALI_BA already
  implements for the continuum channel.
- Rybicki & Hummer (1992, A&A 262, 209): overlapping
  transitions/continuum - already inside the Sobolev
  replay; no extra development, but it means the coupled
  residual is MORE complete than a naive line+continuum
  split, strengthening the case for the coupled tier
  over further continuum ALI.
- Arramy, de la Cruz Rodriguez & Leenaarts (2024, A&A
  690, A12): JFNK for multilevel NLTE - up to 2x fewer
  formal evals vs RH-ALI on FAL-C/ Ca II / H I. Directly
  the formulation JFNK_COUPLED implements. Their atoms
  are small and plane-parallel; our extension needs the
  items above (preconditioning at 8331-frequency scale,
  non-monotonic consistency, Sobolev replay cost).
- Broyden et al. (2011, A&A): quasi-Newton 4.5x faster
  than MALI - supports Newton-on-SE as the right
  direction; our restarted GMRES is the matrix-free
  analogue.
- MOrad (A&A 2015): Newton-Krylov + PETSc at scale -
  proof the coupled approach survives production
  problem sizes; their exact-Jacobian route is not
  available to us (no assembled coupled Jacobian), so
  JFNK is the correct choice.
- DexRT/Lightweaver (Osborne & Milic 2021; radiance
  cascades 2024): state-of-the-art solar NLTE still uses
  same-preconditioning MALI, "no forms of convergence
  acceleration currently employed" - the field's answer
  is better preconditioning, which is what sections 2/6
  deliver. It also cautions that fully-coupled Newton
  is not yet the community default even where it works.

## 12. Risks

- Eval cost: if the coupled sweep cannot be made cheaper
  (section 2), JFNK-CPL may reduce iteration count but
  not wall time. The reduced-cost JVP is therefore on
  the critical path, not optional.
- The lagged-store caveat (section 7): if the tail is
  radiation-lag-driven rather than population-driven,
  no population-space Newton fixes it. X-arm data must
  distinguish these before this proposal is funded.
- Non-monotonic KNOP route: the coupled sweep inherits
  KNOP's own open drift issues (FG_FIELD_DRIFT_AUDIT);
  JFNK-CPL on KNOP should wait until the field is
  demonstrably converged, or it will Newton-solve a
  moving target.
