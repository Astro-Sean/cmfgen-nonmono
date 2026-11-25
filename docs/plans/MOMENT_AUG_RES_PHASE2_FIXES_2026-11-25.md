# Phase-2 fix proposals for the non-monotonic moment solver

Opened 25-Nov-2026. These are the separate implementation proposals
required by the Phase-0 audit task. None of them is part of the
MOM_AUG_RES diagnostic patch; each is a bounded investigation with
its own regression tests and exit criterion. Ordering is by
dependency: FIX 1 and FIX 2 are correctness work on existing
iteration logic; FIX 3-5 restructure the operator and depend on
the augmented-residual evaluator as the acceptance instrument.

Acceptance instrument for all five: `MOM_AUG_RES=T` writes both
RESID_IMPLEMENTED and RESID_TARGET after the last committed
operation. A fix is credited only with a demonstrated reduction
of the target residual on a real case, not with a smaller GMRES
or implemented-equation number.

## FIX 1: adaptive relaxation (RES_PASS ratio misuse)

Defect (verified, audit finding F): MOM_ADAPT_OMEGA_V2 estimates
the spectral radius of the relaxation map from ratios of
residuals measured on ALTERNATING directional passes. A ratio of
residuals from two different maps is not a spectral-radius
estimate of either one; when the compressing and expanding maps
have different contraction factors the estimator reads the beat
note, not the physics.

Proposal:
1. Store the per-pass residual separately for each direction
   class (expanding / compressing / boundary-tagged).
2. Form the ratio only between residuals of the SAME map:
   r_dir(n+2)/r_dir(n), since the sweep alternates. One full
   cycle (two passes) is a valid same-map estimate; single-pass
   ratios are not.
3. Alternatively evaluate the residual of the composite
   two-pass map directly: keep the field at pass start, apply
   both directional maps, and measure ||F(F(x)) - F(x)|| /
   ||F(x) - x|| - the composite-map ratio at the cost of one
   extra field snapshot.
4. Gate the corrected estimator on the same MOM_ADAPT_OMEGA_V2
   key; fixed relaxation remains the controlled default. Never
   blend both estimators in one run.

Regression tests: manufactured alternating maps with known,
different spectral radii per direction in drv_mom_nu style;
verify the estimator returns each map's own rho when queried
same-direction, and the composite rho for the two-pass map.
Unit case where pass residuals oscillate (beat) but both maps
contract - the V2 estimator must not flag divergence.

Exit criterion: estimator reproduces known rho on manufactured
maps to <5%; on a real capture, same-direction ratios are
consistent across cycles. Decision recorded in CHECKLIST.

## FIX 2: uniform final acceptance

Defect (verified, audit finding G): the acceptance sequence -
finite check, positivity check, 50x jump guard, bulk-regression
check, best-pass selection - is applied at pass level, then the
outer Anderson acceleration commits a corrected field AFTER
those guards ran. An AA update can therefore commit a field
that never passed the jump/positivity/regression checks.

Proposal:
1. Move the final acceptance decision to a single routine run
   after best-pass selection AND outer AA: one policy, applied
   once, to the exact field that will be committed.
2. On rejection restore the pre-AA saved field atomically
   (saved state already exists for best-pass selection).
3. Keep the per-pass checks - they feed pass selection - but
   the commit guard is the single post-AA gate.
4. The augmented residual evaluator runs on whatever the gate
   commits, so a rejected AA step is visible in the diagnostic.

Regression tests: drv-level case where AA produces a negative
or jump-violating cell while the pre-AA field passed all checks;
the gate must reject and restore. Case where AA improves the
residual legitimately and commits. Verify the rejection leaves
no partial state (all derived arrays consistent).

Exit criterion: no path from AA to commit bypasses the guard;
manufactured violating-AA case restores cleanly.

## FIX 3: EPS stencil derivation

Defect (verified, audit finding D): the EPS contribution uses a
depth-upwind midpoint donor selected by sign(EPS). Stock
mom_j_cmf_v11.f uses a centred midpoint evaluation. The upwind
choice was introduced to keep the matrix an M-matrix, not
derived from the moment equations; diagonal dominance is not
evidence of physical correctness.

Proposal:
1. Derive the EPS term directly from the CMF moment equations
   at the same midpoint locations the V11 stencil uses; record
   both forms in the derivation doc as equation-numbered steps.
2. Implement the centred stencil behind a diagnostic flag
   (MOM_EPS_CTR, default F) so implemented-vs-target residuals
   can be compared on real captures.
3. Compare three variants on manufactured fields: V11 centred,
   current depth-upwind, and a conservative reformulation if
   one falls out of the derivation (flux-form splitting that
   keeps both accuracy and M-matrix structure).
4. Judge by RESID_TARGET reduction and convergence order on
   smooth and shock profiles - never by matrix conditioning.

Regression tests: manufactured smooth profile (2nd-order
convergence check on grid refinement); manufactured shock
(conservation and positivity); the existing drv_mom_augres
epsmix case extended to flag both stencils.

Exit criterion: the retained stencil is justified by the
derivation and shows equal-or-better target residual; if the
centred stencil wins on accuracy but loses on stability, the
continuation plan goes back to the user, not a hybrid hack.

## FIX 4: implicit H coupling

Defect (verified, audit finding E): MOM_H_IMPL adds diagonal
mass |VB|+|VC| to TB while the true H coupling stays lagged in
X0. This is a stabiliser, not an implicit treatment; the
residual-of-record stays the partial-block residual.

Proposal:
1. Put the H reconstruction equations into the augmented
   operator as explicit rows - z(k) = [x(1),h(1),...] layout
   per docs/plans/MOMENT_H_EPS_IMPLICIT_DERIVATION_2026-11-25.md.
   H remains an independent unknown; no recursive elimination
   (the upwind graph contains chains and rings).
2. First stage: keep the block GMRES on x but evaluate the
   augmented matvec in the residual diagnostics (already done
   via MOM_AUG_RES). Second stage: prototype the augmented
   solve in the unit driver before touching production -
   BLOCK_SOLVE_SCOPING Phase C infrastructure is the vehicle.
3. Retire MOM_H_IMPL only when the augmented operator shows
   equal-or-better contraction without the diagonal mass.

Regression tests: drv_mom_augres him case already isolates the
mass term (implemented exact, target sees the added mass);
extend to measure contraction with and without MOM_H_IMPL on
the iterated-probe analogue of test_mom_blk.

Exit criterion: augmented-operator prototype contracts without
diagonal mass on a case where H_IMPL currently stabilises;
documented before any production change.

## FIX 5: retire equation modifications (continuation)

Defect (verified, audit finding E): MOM_LAG_BND, MOM_EPS_LIM and
MOM_XM_POS each change the equations being solved - they bound
lagged sums, cap boundary EPS flux, and clip the RHS. A field
can converge on the bounded implementation while not solving
the intended equations (RESID_TARGET exposes exactly this).

Proposal:
1. Continuation strategy: scale the modification strengths
   (LAGSCL floor, EPS cap, clip floor) along a homotopy
   parameter, stepping toward the unmodified equations only as
   fast as RESID_TARGET keeps decreasing.
2. Each mechanism retires separately - they do not interact
   through the same code path. Order by observed residual
   damage on real captures (counts now reported per mechanism).
3. The final accepted state must show convergence toward
   RESID_TARGET, not a permanently modified map.

Regression tests: drv_mom_augres lagb/clip/epscap cases verify
each counter and the implemented-vs-target split per mechanism;
a continuation driver that ramps LAGSCL and measures target
residual at each step on a manufactured case.

Exit criterion: a committed field with all three mechanisms
off AND target residual comparable to the implemented residual;
or a documented statement of which mechanism is load-bearing
with the residual cost quantified.

## Status update (25-Nov-2026)

FIX 1 landed. RES_PASS is now slotted per sweep direction
(expanding slot 1, compressing slot 2) and the new
MOM_SAME_DIR_RHO routine ratios only consecutive entries of the
same direction - under ALT mode that is the stride-2, same-map
ratio the proposal asks for. The estimator pools ratios from
both direction slots when both have history and reports LOK
only with three or more ratios, otherwise the old cross-call
residual ratio stands in. Gate stays MOM_ADAPT_OMEGA_V2 with no
blending. Manufactured-map check (drv_mom_gate, test_mom_gate
45/45): rho recovered within 5% on alternating maps with
per-direction rho 0.45 and 0.5, a diverging map reads 1.69,
and the consecutive-ratio beat that the old estimator mistook
for divergence is not flagged. Exit criterion met on the
manufactured suite; the real-capture consistency half is still
owed.

FIX 2 landed. JBEST/HBEST now snapshot the field just before
the outer-AA block, and the new MOM_RF_COMMIT_GATE runs once
after AA and before the maturity census: non-finite or 1e40
entries, cells driven to J<=0, 50x jumps against the retained
store, and bulk regression of the pre-AA n5/n50 counts all
reject, and rejection restores the snapshot in place with a
WARN line per cause. drv_mom_gate verifies accept on a clean
AA step and reject-plus-restore on non-finite, negative,
jump-violating and bulk-regressing commits; the seed-store
case confirms the gate compares against the right reference.
Exit criterion met: every path from AA to commit passes the
single gate.

## Ordering note

FIX 1 and FIX 2 are independent and can proceed in either order.
FIX 3-5 all consume the augmented evaluator; they should follow
the first real-case MOM_AUG_RES runs so the residual damage
ranking is evidence, not guesswork.
