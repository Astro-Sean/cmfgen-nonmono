# Literature audit of the USE_LAM_ES=F moment route

Date: 25-Nov-2026. Read-only audit of the working tree while the
pinned production arms were running. Main files:
`plane/mom_j_cmf_nonmono_v1.f` (7912 lines) and
`plane/mom_blk_gmres_mod.f` (290 lines). Line numbers cite this
working tree, not the 24-Nov audit version.

Question: does the current implementation solve the radiation moment
problem intended by the cited non-monotonic comoving-frame literature,
and is its numerical stabilization mathematically consistent with that
target?

References checked against the original papers:

- Baron and Hauschildt 2004, A&A 427, 987: non-monotonic velocity
  turns the wavelength direction into a boundary-value problem with
  mixed boundary conditions. The wavelength derivative must use a
  local upwind discretization chosen from the local sign of the
  coupling coefficient.
- Knop, Hauschildt and Baron 2009, A&A 496, 295: the positive
  generalized opacity is `chi_hat = chi + xi*a*p^|`. The remaining
  wavelength-coupling terms are retained in the formal solution; the
  paper does not discard or clamp them.
- Knop, Hauschildt and Baron 2009, A&A 501, 813: the physically
  guided Gauss-Seidel ordering follows information propagation in
  space and wavelength. It is a solver ordering for the same coupled
  system, not a different physical system.
- Ceulemans, De Ceuster and Decin 2025, JQSRT 343, 109470: numerical
  oscillations are controlled by following the frequency shift of the
  characteristics and enforcing the extra frequency boundary
  conditions that stabilization introduces.
- Hillier and Miller 1998 and stock `mom_j_cmf_v11.f`: CMFGEN moment
  closure baseline. The stock elimination uses centered midpoint sums
  for the EPS terms.
- Dessart, Audit and Hillier 2015: relevant to the Sobolev non-
  monotonic production route, but not a derivation of this CMF moment
  solver. It does not validate the present `USE_LAM_ES=F` equations.

## 1. Executive verdict

The implementation is numerically careful as engineering, but it is
not yet a literature-correct solve of the full non-monotonic moment
problem. The principal J-frequency advection term uses the right
sign-keyed local-upwind idea, and the GMRES kernel itself is basically
sound. The coupled solve nevertheless stops short of the literature's
wavelength boundary-value problem: the H and EPS wavelength-coupling
terms remain frozen in the right-hand side, H is not reconstructed
after the block J solve, and nonlinear positivity/trust-region edits
are applied after the reported linear convergence without a final
residual check.

The largest correctness concern is not GMRES. It is that several
permanent guards change the equations being solved: depth-upwinding of
the EPS midpoint value, LAG_BND scaling, EPS limiting, XM_POS, and
MOM_H_IMPL. These may be useful continuation or preconditioning
devices, but a converged solution with them active is not yet the
unmodified moment solution. Their activation must be treated as a
homotopy or diagnosed as model error, not as proof of convergence.

Overall verdict: PARTIAL with defects. Stable enough to produce
bounded iterations, not correct enough yet to certify convergence to
the intended equations.

## 2. What matches the literature

### 2.1 Sign-keyed wavelength upwinding: VERIFIED in principle

Baron and Hauschildt require the upwind direction to be local in
space and wavelength because the coefficient of the frequency
derivative changes sign in a non-monotonic flow. The implementation
follows that principle for the principal J advection term:

- expanding rows read the bluer previous frequency;
- compressing rows read the redder stored frequency;
- the diagonal receives the positive magnitude of the local coupling;
- the upwind contribution enters the right-hand side with the matching
  positive magnitude.

Relevant code: `mom_j_cmf_nonmono_v1.f:2121-2148`,
`:2512-2543`, `:2717-2733`, `:2755-2813`. This is algebraically
consistent with first-order upwinding of the principal frequency
advection term. The use of `SIGMA_TRUE` is also consistent with the
project requirement that the physical non-monotonic velocity gradient
not be monotonized.

### 2.2 Positive generalized-opacity split: VERIFIED with a caveat

Knop Paper I puts only the positive part of the frequency-coupling
coefficient into `chi_hat`. The implementation does the same in
`chi_hat = chi + |GAM|/dlnnu` (`:2331-2399`). The caveat is that the
paper keeps all remaining wavelength terms in the transfer problem.
The present implementation keeps some terms explicit and modifies
others, which is the subject of findings F1 and F4.

### 2.3 GMRES kernel: VERIFIED as a linear solver, with reporting limits

`plane/mom_blk_gmres_mod.f` implements a standard restarted GMRES:
Arnoldi basis construction, two classical Gram-Schmidt passes, Givens
rotations, triangular least-squares back-solve, restart, and a final
residual evaluation. The depth-tridiagonal preconditioner indexing is
consistent with the captured row convention. No indexing or sign
defect was found in `BLK_MATVEC` or `BLK_PRECOND_VEC`.

Two limits remain. First, the reported residual is the left-
preconditioned residual, not the raw `b-Ax` residual (`mom_blk_gmres
_mod.f:194-203`). Second, the residual is measured before the wrapper
projects negative or large cells back to the store, so it does not
certify the field actually committed. These are verification defects,
not GMRES arithmetic defects.

## 3. Findings

### F1 - The block solve is not the literature's full wavelength boundary-value solve (DEFECT)

The captured system makes only the principal upwind-J edge implicit:

```text
TA*x(i-1,k) + TB*x(i,k) + TC*x(i+1,k) - CUP*x(i,kup) = X0
```

The capture is explicit that lagged H and EPS scatter terms remain in
`NU_X0` (`mom_j_cmf_nonmono_v1.f:3376-3388`). In the continuous moment
problem those are also frequency-coupled terms. In a compressing
region they must obey the same local-upwind boundary-value structure.
Baron and Hauschildt solve the full boundary-value problem in
wavelength; they do not solve only the largest advection edge and lag
the rest.

This has a second consequence. `MOM_BLK_SOLVE` updates only
`J_NEXT_FRQ` (`:7851-7885`). It never reconstructs `H_NEXT_FRQ` from
the new J field. Yet compressing rows read `H_NEXT_FRQ` as the upwind
H field in the next pass (`:2815-2852`). The post-block moment pair is
therefore internally inconsistent even when GMRES converges. H is not
an independent outer coefficient like the Eddington factors; it is the
first moment reconstructed from the solved system.

Impact: the GMRES residual certifies a partially linearized surrogate,
not the full captured moment equations. This helps explain why a
2-3 iteration linear solve can coexist with a non-contracting outer
loop.

Required direction: derive the eliminated H and EPS frequency coupling
in the same sign-keyed upwind graph and include their J dependence in
the implicit operator, or explicitly label the present solve as a
preconditioned outer splitting and measure the full nonlinear residual
after every commit.

### F2 - The committed block solution is not residually verified (DEFECT)

After GMRES returns, cells are rejected or capped one at a time:

- non-finite values are frozen;
- negative values are rejected;
- large positive values are trust-region capped;
- the damped store is then written (`:7847-7885`).

No residual of the accepted field is evaluated after these edits. The
printed `resid` belongs to the raw GMRES vector before rejection and
projection. In live production these edits are not rare: the fixed-relaxation
arm still rejected roughly 6,200 negative cells per pass at iteration 14,
and the adaptive cold arm rejected hundreds to thousands per pass.

A linear solve followed by a nonlinear projection is a nonlinear
operator. It can be valid, but only if the post-projection residual is
measured and used in acceptance. At present a pass can be called
converged even when the accepted vector does not satisfy the surrogate
linear system at the projected cells.

Required direction: compute the true captured-system residual after
the write-back policy, report raw and preconditioned norms separately,
and make acceptance depend on the post-write residual.

### F3 - The EPS depth-upwind stencil is not derived from the literature (DEFECT)

Stock V11 evaluates the EPS midpoint term with centered sums,
`(J_i+J_{i+1})`, in both the matrix elimination and the H reconstruction.
The non-monotonic solver replaces that midpoint average with a depth-
upwind node chosen from the sign of EPS (`:2612-2641`, `:2882-2987`,
`:4168-4193`).

The justification in the comments is matrix sign control. That is a
numerical reason, not a derivation. EPS is produced by the frequency-
coupling elimination; its sign is a wavelength-coupling direction, not
a spatial propagation direction. The papers justify upwinding in
wavelength. They do not justify upwinding this midpoint value in
radius. Unless a new derivation shows that this term is a physical
radial advection flux, the stencil is first-order interpolation
presented as upwinding and changes the target moment equations.

This is especially important because the H reconstruction uses the
same modified stencil. The implementation is internally consistent
with the modified equations, but those equations are not the stock
CMFGEN moment equations and have not been tied back to the cited
literature.

Required direction: derive the modified stencil from the moment
equations, or restore the centered form and treat the resulting matrix
problem with a solver that does not require an M-matrix.

### F4 - LAG_BND, EPS limiting, and XM_POS are permanent equation modifications (DEFECT)

The three main positivity/stability guards do more than reject failed
iterates:

- `MOM_EPS_LIM` caps signed lagged EPS transfer and rescales the
  implicit coefficient to match (`:2913-2987`);
- `MOM_LAG_BND` scales the row's whole lagged contribution by a
  diagonal-excess bound (`:2990-3030`);
- `MOM_XM_POS` replaces a negative assembled right-hand side by zero
  (`:3032-3052`).

A negative eliminated right-hand side is not automatically unphysical.
The H-divergence and EPS terms are signed coupling terms. Positivity
of the full nonlocal radiation field does not imply positivity of each
algebraically eliminated local right-hand side. These guards therefore
change the discrete equations, not only the iteration path.

LAG_BND is the clearest case: if its scale factor is below 1 at the
accepted fixed point, the fixed point solves the scaled map. The live
capture census shows why this matters: order 1e5 rows per pass have
non-positive diagonal excess (`xs<=0` in all three arms). The bound is
not a rare safety net in this regime.

These mechanisms may still be useful as a continuation method. Then
the continuation must be driven back toward the unmodified equations
and the final residual of that unmodified system reported. They must
not be treated as proof of convergence while permanently active.

### F5 - MOM_H_IMPL is not an implicit H treatment (DEFECT)

`MOM_H_IMPL` adds `|VB|+|VC|` to the diagonal while leaving the same
signed H terms in the right-hand side (`:2648-2661`). That is diagonal
regularization, not implicit coupling. A correct implicit treatment
would express upwind H through its reconstruction equation and move
the corresponding J dependence into the matrix with the correct signs
and frequency edges.

Because the right-hand side is unchanged, even an exactly converged
H_IMPL fixed point solves a different equation. The warm H arm is also
empirically negative: at iteration 697 it still committed a largest
stored-field change of `5.9E10`, with 227,202 cells above 5% and an
outer/inner luminosity ratio of only about 0.039.

Recommendation: retire this experiment. Do not use MOM_H_IMPL as
evidence about the true moment equations.

### F6 - The adaptive same-map relaxation estimator mixes different maps (DEFECT)

The refresh alternates expanding and compressing pass directions
(`:5993-6005`). The V2 estimator nevertheless forms ratios of
consecutive entries in `RES_PASS` (`:7451-7490`), even though
consecutive passes often apply different directional maps. The code
already computes the correct same-direction ratio for stall detection
(`:6630-6632`, `:6794-6824`) but does not use it for omega.

The production signature is exactly the predicted one. In the adaptive
cold arm the estimator repeatedly hit its clamps and moved relaxation
from about 0.81 toward 0.43 after reporting same-map ratios of 0.05
and then 20. Those are not contraction measurements of one affine map.

This is a concrete implementation bug independent of the physics
audit. Fixed relaxation is therefore the correct control, not merely a
convenient test.

### F7 - Final write-back has no uniform acceptance rule (DEFECT)

The 50-times jump guard exists only inside the non-converged best-pass
path (`:6977-7009`). A converged refresh exits before that block and
commits the iterated field without the equivalent guard (`:6714-6754`).
Outer Anderson acceleration is applied after the best-pass decision
(`:7055-7294`) and can modify the field again without re-running the
jump guard or the bulk-regression test.

The live arms demonstrate both holes:

- adaptive cold, iteration 27: a refresh marked converged committed a
  largest stored-field change of `6.8096E10`, with 118,430 cells above
  5% and 7,028 above 50%;
- fixed r01, iteration 14: the raw best-pass field was rejected because
  9,173 cells jumped more than 50x, but outer AA then wrote a correction
  producing a committed change of 7.9121 and 6,463 cells above 50%.

The second case was less damaging, but it shows that the raw veto is
not a final trust region. Acceptance currently depends on path ordering,
not on the final field.

Required direction: after best-pass selection and after any AA update,
evaluate the same final residual, jump, bulk, and positivity policy
once. Either accept the final composite field or reject it as a whole.

### F8 - Frequency boundary conditions are approximate and partly silent (PARTIAL/DEFECT)

Baron and Hauschildt emphasize mixed boundary conditions in
wavelength: at each spatial point the local sign of the coupling
decides which wavelength edge supplies inflow data. In this code,
compressing rows whose upwind slot lies beyond the red edge receive
`JUP=0` (`:2805-2812`). That is a zero-inflow boundary at a finite
grid edge. It may be adequate if the red-edge field is negligible,
but it is not demonstrated from the model physics and is not the
general boundary condition required by the literature.

The off-grid count is now diagnosed, but diagnosis does not make the
boundary correct. Ceulemans et al. also show that stabilized frequency
stencils introduce additional boundary conditions that must be
enforced. The integer stride matching used here does not implement
that boundary machinery.

Required direction: state and implement the physical blue/red frequency
boundary conditions explicitly, then test their influence by extending
the frequency grid or using a controlled boundary model.

### F9 - Frequency matching is not yet the Ceulemans stabilization (PARTIAL)

`MOM_FREQ_MATCH` chooses an integer bin offset with `CEILING` and
widens the log-frequency denominator (`:2302-2395`). This is a useful
local-upwind correction when the characteristic crosses multiple grid
bins, but it is not equivalent to the Ceulemans frequency-connected
stencil. It does not construct a continuous shifted frequency mesh and
does not enforce the additional boundary conditions described there.

Verdict: plausible engineering approximation, currently unverified as
a literature implementation. It needs manufactured cases with a known
characteristic foot and boundary-condition sensitivity tests.

### F10 - Eddington closure and boundary factors are operator-split, with extra closures risk (PARTIAL)

Freezing `K_ON_J`, `NMID_ON_J`, `NMID_ON_HMID`, HBC, and NBC during a
moment refresh is standard operator splitting. That part is acceptable
if the outer iteration converges. The concern is the additional
closure replacement:

- `K_ON_J` is forced into `[1E-4,1]`;
- `NMID_ON_HMID` is capped at magnitude `1E12`;
- boundary HBC/NBC values outside `[0,1]` are replaced by zero
  (`:1929-1955`, `:3075-3137`).

Some bounds prevent division by degenerate measurements, but they can
also change the sign or magnitude of `1+sigma*g` near compression.
The boundary closure assumes an unilluminated outer boundary with
non-negative flux. A non-monotonic interacting ejecta/circumstellar
model can produce inward or sign-changing flux in exactly the region
of interest, so this assumption needs physical verification rather
than a universal floor.

## 4. Live-arm evidence against the current setup

The three running arms are consistent with the static findings:

1. Adaptive cold `13671134`, iteration 27: GMRES residuals are small
   (`2E-4` to `5E-4`), but hundreds to thousands of negative cells are
   rejected per pass; source-inflation ratio `rmax` reaches `1E9`-`1E13`
   in recent passes; H contributes about 15% of the explicit source
   census; outer AA is rejected; a converged refresh still commits a
   `6.8E10` largest stored-field change. Luminosity ratio is about
   0.33.
2. Fixed-relaxation r01 `13671873`, iteration 14: GMRES again converges
   in about three iterations, but about 6,200 negative cells per pass
   are rejected; H/X0 is about 7-8%; the best-pass field is rejected by
   the jump guard; outer AA still writes a smaller correction. The
   luminosity ratio improves to about 0.54, showing that damping helps
   globalization but does not fix the equations.
3. Warm H `13671005`, iteration 697: H/X0 is also about 8%, so the H
   channel is real, but MOM_H_IMPL is not the correct implicit
   treatment. The arm remains near luminosity ratio 0.039 and still
   commits `1E10`-scale field changes.

The main mechanistic reading is unchanged: the captured linear solve
is cheap and well behaved, while the unimplicit H/EPS content and the
nonlinear write-back policy dominate the outer loop.

## 5. Stability assessment

Stable in the limited sense:

- failures are usually counted and contained;
- DGTSV failures do not propagate silently;
- non-finite stores are repaired;
- the GMRES kernel itself is mathematically coherent;
- fixed relaxation bounds the outer correction size.

Not stable in the literature-convergence sense:

- the map being accelerated is repeatedly redefined by clamps,
  damping, projection, ownership suppression, and pass ordering;
- the final accepted field is not the vector whose residual is printed;
- major stabilizers change the fixed point, so convergence of the
  stabilized map need not solve the intended equations;
- the block solve omits part of the wavelength boundary-value coupling
  that motivates the method.

## 6. Ordered recommendations

1. Write down the full discrete moment system, including H elimination
   and EPS terms, as a wavelength-depth operator. Use that derivation
   to decide exactly what belongs in `NU_CUP` and what belongs in the
   depth block. This is the correctness anchor for every later choice.
2. Add a post-write residual evaluation. Report raw and preconditioned
   residuals for the final committed field, after rejection, caps,
   relaxation, and AA.
3. Fix the V2 relaxation estimator to use same-direction or same-
   composite-map residuals. Until then keep `KN_ADAPT_RELAX=F` in all
   controlled arms.
4. Make one final acceptance gate apply after best-pass selection and
   outer AA. Remove the converged-path bypass of the jump guard or
   replace the 50x veto with a clearly defined trust-region update.
5. Treat LAG_BND/EPS_LIM/XM_POS as continuation parameters. Measure
   their activation and drive them off for final acceptance, or replace
   the need for them by implicit treatment of the H/EPS coupling.
6. Re-derive or revert the EPS depth-upwind stencil. Do not keep it
   solely because it improves matrix signs.
7. Retire MOM_H_IMPL. If H coupling matters, implement it through the
   actual H reconstruction equation.
8. Specify physical frequency boundary conditions and test the red-edge
   zero-inflow approximation by boundary sensitivity experiments.

## 7. Dated verdict

25-Nov-2026: PARTIAL with defects. The sign-keyed upwind concept,
positive generalized opacity, tridiagonal depth solve, and GMRES kernel
are sound enough to build on. The current production setup does not yet
solve the full literature wavelength boundary-value moment problem:
H/EPS coupling remains outside the block solve, H is stale after the
block update, post-projection residuals are not measured, several
permanent guards redefine the equations, and final write-back acceptance
is path-dependent. The running arms are useful diagnostics of these
defects, but their output cannot yet be called convergence to the
intended non-monotonic CMF moment equations.
