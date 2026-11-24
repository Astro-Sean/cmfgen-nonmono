# Phase D: integration design - coupled wavelength-depth GMRES solve

24-Nov-2026. Status: DESIGN for user sign-off (user directive
"implement this plan, same coupling as Baron-Hauschildt" recorded
24-Nov-2026). No integration before sign-off; unit-level pieces
may be built ahead since they touch nothing in production.

## 1. What goes implicit vs what stays outer-iterated

Implicit in the block operator (per direction set):
- depth tridiagonal TA/TB/TC exactly as assembled now, including
  the PSI diagonal absorption (the B&H chi-hat piece, xi=1)
- the frequency coupling coeff*J_upwind as an implicit unknown -
  expanding rows couple k-1, compressing rows couple k+MOM_NM_J;
  coefficients already captured as NU_CUP/NU_KUP

Outer-iterated exactly as now (unchanged machinery):
- opacities, emissivities, source functions (population iteration)
- H reconstruction and the H store (H enters X0/source side; it
  is derived from J algebraically per frequency, so a J-only
  block solve remains consistent - H is recomputed after the solve)
- EPS scatter terms (same argument; sign-split stays in TA/TB/TC)
- boundary closures, LAG_BND/EPS_SCL/XM_POS guards - all fold into
  coefficients and RHS before the solve, unchanged

RHS: X0 minus the frequency-coupling piece (that piece moves from
lagged source to implicit unknown). Everything else in X0 stays.

## 2. The solve itself

Per refresh pass, per direction set:
1. assemble rows exactly as now (single sweep, coefficients and
   RHS captured - the capture machinery from MOM_NUSOLV is reused
   verbatim: it already records TA/TB/TC/X0/coupling per row)
2. GMRES(m~40, tol ~ 1e-3 relative - outer iteration does not
   need tighter), left-preconditioned by per-frequency DGTSV:
   the preconditioner application IS the current sweep solve
3. matrix-free action: per row, TA*x(i-1)+TB*x(i)+TC*x(i+1)
   + NU_CUP*x(NU_KUP) - a few fused loops over the captured
   arrays; never assembled
4. damped non-negative write-back unchanged; negative/nonfinite
   rejections unchanged; H reconstruction pass unchanged

Frequency boundary conditions: expanding rows at k=1 have no
lower neighbour (edge coefficient zeroed, as now); compressing
rows at the red edge likewise - matches current capture.

## 3. What this replaces and preserves

Replaces: the lagged J_next reads inside the linear solve (the
defect). The sweep structure, refresh orchestration, ownership
logic, damping, and census all survive - the GMRES solve slots in
where the per-frequency DGTSV chain runs now.

Preserves: default path untouched. New behaviour gated on a new
VADAT key MOM_BLKGMRES (default F) per AGENTS 9.1 - written
parameter proposal required before the key lands. MOM_NUSOLV
remains as the chain-solve experiment; the block solve supersedes
it if the arm converges.

## 4. Honest risks

1. Outer-iteration feedback (H store scaling with committed J)
   is untouched by this design: the linear defect is removed but
   the nonlinear loop remains. Mitigation already present: damped
   writes, best-pass veto. If the arm still excursions, the H
   channel moves inside the residual evaluation next (bigger
   surgery, scoped only if needed).
2. GMRES convergence at production conditioning is extrapolated
   from a 384-unknown toy - the first real test may need restart
   tuning or a stronger preconditioner (e.g. two frequency
   neighbours folded into the preconditioner blocks).
3. Rings spanning the frequency grid mean the preconditioner is
   not the exact inverse of any submatrix that contains the
   coupling - convergence relies on GMRES seeing coupling in the
   action. Toy says this works; production may differ in degree.
4. Cost: one GMRES solve per pass costs ~ m sweep-passes; if the
   refresh still needs 20-40 passes the cost multiplies. Target:
   contraction should cut pass count sharply - net win expected,
   not guaranteed.

## 5. Implementation stages (phase E preview)

E1. unit_tests/f90/drv_mom_blk.f: dense-reference GMRES test on
    manufactured coupling graphs (extending drv_mom_nu cases);
    python battery pinning convergence counts. No production code.
E2. MOM_BLK_GMRES module: matrix-free action + GMRES shell +
    preconditioner hook, unit-tested standalone against drv cases.
E3. Integration behind MOM_BLKGMRES: capture reuse, write-back
    reuse, census extension (GMRES iterations, residual).
E4. Battery + arm on eps1em3 seed, verdict per checklist rules.

## 6. Verification anchors

- drv_mom_blk: GMRES solution matches dense reference to machine
  precision on every manufactured graph including marginal rings
- iterated probe analogue: the composite refresh map with the
  block solve must show measured growth < 1 in the explosive
  regime where the lagged version measured 1000/iter
- arm exit criteria as for the nusolv arms: converges / stalls
  bounded / detonates, with census evidence either way

## 7. VADAT parameter proposal for MOM_BLKGMRES

- proposed key and default: `MOM_BLKGMRES`, default F.
- user-level decision: whether a research run replaces the lagged
  frequency pass with the coupled wavelength-depth linear solve after
  each MOM refresh pass.
- why an internal default is insufficient: the solve is an unproven
  numerical experiment at production conditioning and costs several
  sweep-like operator applications per pass. It must remain opt-in
  until the eps1em3 arm and wider model matrix establish reliability.
  The default F path is unchanged.
- interactions: `MOM_NUSOLV` remains the per-depth chain experiment.
  If both are true, `MOM_BLKGMRES` takes precedence because both
  consume the same captured rows and the block system contains the
  depth coupling the chain solve leaves lagged. Capture keys, damping,
  write caps and active-set behaviour are unchanged.
- tests: `test_vadat_keys` covers present and absent keys on the
  control store and solver module; `test_mom_blk` covers the matrix-
  free solver against dense references; `test_mom_solve` exercises
  the integrated refresh arm and the precedence path.
- physics or numerical strategy: numerical strategy only. The same
  captured row equations, source terms, guards, boundary closures and
  write-back policy are used; no physical closure or opacity model is
  changed.
