# KNOP direct wavelength-chain solver: matrix derivation (Phase 2.7)

Status: derivation complete; no solver implemented. Per phase
brief, this document reports the matrix structure and the
proposed numerical method before any implementation.
Baseline: 6132c601, 69f5b5e3.

## 1. The assembled per-cell equation (verified in source)

For ray LS, channel c (M/P), ray point i, frequency L, with
frozen coefficients (`fg_j_cmf_nonmono_v1.f`):

- Generalized opacity (diagonal):
  `chi_hat(i,L) = chi(i,L) + GAM_MAG(i,LS)/DLN_C`
  (:4755-4763), where DLN_C = dLOG_NU when GAM>=0 and
  DLN_C = DLN_RAY(i) (the local log-spacing to the donor
  frequency) when GAM<0 (:4752-4753).
- Explicit donor source: S' = ETA/chi_hat + q_c*I_donor with
  q_c = GAM(i,LS)/(DLN_C*chi_hat) for the compressing branch
  (:4734-4737) and the standard Q = (GAM/dLOG_NU)/chi_hat
  for expanding points (:84-85, :4177-4178).
- Donor value: interpolated from the two bracketing sparse
  slots KC1(i), KC2(i) of the same (LS, channel) store with
  weight W_INEXT (:4436 OMP private list; map diag :5926-5937).
- Spatial march: the usual short-characteristics relation
  makes I(i,L) linear in I(i-1,L) plus local source terms.
- KN_COUPLING_OK gate (:4733, :4760): a compressing point
  with no valid bracketing lagged slot drops BOTH sides of
  the coupling (implicit |GAM|/DLN and explicit donor);
  boundary-foot points keep the outflow.

So the equation at unknown x(i,L) = I(i,L) is linear:

   a_i(L) x(i,L) - b_i(L) x(i-1,L) - c_i(L) x(j1(i),L+e)
                                 - d_i(L) x(j2(i),L+e) = s_i(L)

with e = -1 where GAM(i,LS)>0 (donor = previous frequency)
and e = +1 where GAM(i,LS)<0 (donor = next frequency),
j1/j2 the two bracketing depth-slot ray positions, and
a_i(L) = chi_hat + spatial-attenuation diagonal.

Everything on the right except the x terms is frozen:
ETA_RAY, CHI_RAY, GAM, DTAU, boundary intensities. The
system is LINEAR in the field - no approximation needed to
state it. Nonlinear guards in the production code
(intensity floor 1e-30, NONMONO_INEXT_CAP, heal/commit
policy) are post-solve policy, not part of A.

## 2. Matrix structure

Per (LS, channel) pair the system decouples completely:
there is no cross-ray term in the ray solve (rays meet only
through the angular quadrature into JNU_STORE, outside the
inner system), and M/P channels have separate GAM, donors
and boundaries. So the full problem is NP x 2 independent
linear systems.

Inside one (LS, channel) system, the nonzeros of row (i,L)
are columns {(i-1,L), (i,L), (j1(i),L+e), (j2(i),L+e)}.
Order the unknowns by frequency block L, march order inside:

- Block L has edges to block L-1 (from rows where GAM>0)
  AND to block L+1 (from rows where GAM<0). Both directions
  coexist in the same block whenever the ray contains both
  expanding and compressing points - which is always, since
  GAM varies with mu along the ray.
- Inside block L the spatial edges (i-1 -> i) are
  march-ordered: the block diagonal is triangular.

Verdict (Task 1 answer): the matrix is (E) general sparse
with a block-TRIDIAGONAL structure in the L-block ordering,
each off-diagonal block a sparse two-column-per-row map on
the compressing or expanding subset. It is NOT triangular
under any (i,L) lexicographic order, because a compressing
row at frequency L needs a donor in block L+1 while an
expanding row needs L-1.

It IS effectively triangular per row segment if all rows of
a block shared one sign; they do not. A topological-ordering
argument along the transport graph does not yield a global
triangular form either, since the L-edge direction flips
with row index inside each block.

## 3. Size and the Schur-complement reduction

Unknowns per (LS, channel): NI_ray x NCF ~ 150 x 74,704
~ 11M - too large for block-LU directly (150^3 per block).
But the compressing cells are the only ones carrying the
L+1 edge; expanding cells carry L-1 only and can be
eliminated row-by-row in ascending L during assembly. The
only unknowns that must be solved simultaneously are those
feeding forward in L - the compressing-slot values
(N_COMP ~ 12 depths x NCF) plus boundary columns.

Proposed method: Schur complement onto the store slots.

1. March each spatial row; wherever a row's equation needs
   x(j,L-1) with j already eliminated, substitute the
   marching expression - a linear recurrence, exact, O(1)
   per substitution, no fill beyond the two donor columns.
2. The remaining simultaneous system has ~N_COMP unknowns
   per L block with block-tridiagonal coupling; block size
   ~12, so block-LU costs ~NCF x 12^3 ~ 1.3e8 flops per
   (LS, channel) - trivial against the current multi-pass
   iteration.
3. Back-substitute the eliminated expanding rows in reverse
   L order.

This is an EXACT solve of the same discretised equations
(up to roundoff): same chi_hat, same DLN_C, same Q, same
KC1/KC2/W_INEXT donor map, same KN_COUPLING_OK gating, same
boundaries. No previous-generation field enters except as
the documented boundary/edge columns of the formulation.

## 4. Why the lagged iteration stalls (quantitative closure)

Rowwise off-diagonal/diagonal ratio of the wavelength
coupling:

   rho_row = (|GAM|/DLN_C) / chi_hat = 1 - chi/chi_hat

In the stiff region |GAM|/DLN_C >> chi, so rho_row -> 1.
Any pointwise iteration (Jacobi, GS, relaxed writes) is
bounded by max rho_row ~ 0.97 - exactly the code-measured
rho_p90 (Phase 2.6) and the observed ~0.96/refresh outer
decay. The non-contraction is therefore not an ordering or
swap defect: it is the contraction limit of point iteration
on this matrix. The donor-dominated rows carry almost no
diagonal weight, so iteration cannot contract faster than
the coupling allows, and committing a truncated iterate
leaves KN_Q_RESID O(1).

A direct solve removes the iteration entirely for the
frozen-coefficient problem; the remaining outer map
(scattering source feedback, population updates, the
in-sweep/refresh split) is a separate question (Task 6).

## 5. Residual plan (Task 3)

For every direct solve compute r = b - A x with the same
assembled rows; report L2 and L-inf relative residuals,
max absolute and relative residuals, and residuals split by
frequency, depth, GAM sign, and the index bands 59-70,
62, 66-68. KN_Q_RESID stays as the independent map
diagnostic - it must not be substituted for r.

## 6. GMRES assessment (Task 5)

If the Schur construction is not implemented, matrix-free
GMRES on the (i,L) system is viable: matvec = one assembly
sweep per (LS, channel). Unpreconditioned it inherits the
rho~0.97 bottleneck; a block-Jacobi or the L-block-Thomas
preconditioner above is the natural choice. Expected
memory: the Schur matrix needs only ~3 x 12 x NCF nonzero
columns per (LS, channel); dense assembly is neither needed
nor acceptable.

## 7. Validation design (Tasks 4, 7)

- D1 existing lagged run: baseline from caseK2 telemetry
  (already archived).
- D2 direct solve on frozen K2 coefficients: verify
  ||r||/||b|| ~ 1e-12 or better, and that the result is
  independent of the seed field (required for a linear
  solve; any seed dependence is a defect to find).
- D3 existing solver seeded from D2: if the implementation
  is consistent, KN_Q_RESID at the D2 field should collapse
  toward roundoff - the decisive check that both codes
  satisfy the same discrete equations.
- D4 direct solve from perturbed seeds: repeatability.
- Adaptive-omega comparison deferred until D3 validates:
  fixed omega=1 with direct inner solve isolates whether
  the controller was masking inner truncation.

## 8. Constraints honoured

No production default, equation, population, scattering,
RESC_JMAX or J_RESCUE change is made or proposed here. The
direct method solves the same discretised equations; the
documented nonlinear guards remain post-solve policy.

## 9. Driver-level verification: Newton solve of the
   implemented store map (knray p1t1n, 20-Nov-2026)

Before any production-facing implementation, the fixed-point
question was tested directly on the implemented map in the
existing p1t1b finite-difference harness
(unit_tests/f90/drv_knop_ray.f90, new case tag p1t1n;
diagnostic only, no production change).

The state X is the packed compressing store: 3 compressing
depths x 5 active slots x {P,M} channels, NX=30, seeded at
X0V = the 30-iterate fixed point of the refreshed map. For
each map variant the driver forms J2 (FD Jacobian, amplitude
1e-5 of per-element scale), solves
(I - J2) * DXP = F(X) - X by DGESV, and re-evaluates the true
map residual ||F(X*) - X*||. Machine-readable record:
docs/audits/knop_p1t1n_map_solve.json.

### Results by map variant

| variant    | map content                    | J nnz | rho(J)     | Newton RN0 -> RN_final        |
|------------|--------------------------------|-------|------------|-------------------------------|
| fullw05    | sweep + FG_KNOP_REFRESH_V1     |     0 | 0          | 0 -> 0 (X0V already exact)    |
| swp_w05    | sweep only, omega=0.5          |    64 | 0.5 = 1-w  | 5.5e-7 -> 2.4e-14 in 1 step   |
| swp_w10    | sweep only, omega=1.0          |    35 | 0 (nilpotent) | 1.1e-6 -> 0 exactly        |

X*_swp05 = X*_swp10 to 2.7e-14 (L2): relaxation preserves
the fixed point, as required. The true lagged fixed point
sits 9.9e-7 max-element / 2.1e-6 L2 from the committed field
X0V - the same order as the refresh's own truncation
(residual 8.6e-9 relative on |X| ~ 70), so the production
refresh is converging toward the same fixed point the direct
solve lands on. FREEZE_I_PREV ablation (p1t1bf) leaves the
spectrum identical: the packed-store map is not
I_PREV-mediated. Perturbed-iteration probes (p1t1b) contract
at 0.5 for the sweep map and 0 for the refreshed map -
dynamic confirmation of the FD eigenvalues.

### What this establishes

1. The implemented store map is affine: the FD Jacobian is
   amplitude-independent, and one Newton step takes the
   residual to ~1e-14..1e-16. A direct solve of the lagged
   fixed-point problem exists and is reachable in a single
   dense 30x30 solve on this problem.
2. On the p1t1 synthetic problem the unrelaxed lagged
   Jacobian is nilpotent (rho = 0): the wavelength-coupling
   graph is a DAG ending at the unwritten boundary slot.
   Iteration non-contraction cannot occur here; the observed
   eigenvalue 0.5 on the relaxed map is purely the
   under-relaxation factor, not physics.
3. The rho ~ 0.97 production stiffness is NOT reproduced by
   this synthetic problem. Per section 4 it requires either
   cyclic cross-depth donor chains or the |GAM|/DLN >> chi
   regime; the p1t1 geometry has a three-depth acyclic
   chain. A stiffer synthetic case (compression bands with
   cyclic donor mapping, or much smaller DLN) is needed to
   exercise the mechanism the Schur solver targets.
4. Map-level Newton is not the Schur solver: (I - J)^-1
   accelerates the implemented fixed-point iteration. Its
   fixed point coincides with the wavelength-chain solution
   when the map is affine, but it does not expose the
   per-cell matrix r = b - A x residuals of Task 3. The
   Schur-complement derivation in sections 2-3 remains the
   production-relevant formulation; this result validates
   that a direct solve of the same discrete equations is
   both well-posed and computable.

### Status

PASS WITH LIMITATIONS: direct fixed-point solve verified
exact on the implemented map for the synthetic problem;
synthetic problem is acyclic and does not reproduce the
production stiffness regime. No production solver, default,
or dispatch touched; knray p1t1n/p1t1b/p1t1bf all RESULT_OK
(binary rebuilt against libsubs.a at fg_j_cmf_nonmono_v1.f
of 20-Nov-2026, including the concurrent KN_NOCP_PROBE
MODE-1 changes).
