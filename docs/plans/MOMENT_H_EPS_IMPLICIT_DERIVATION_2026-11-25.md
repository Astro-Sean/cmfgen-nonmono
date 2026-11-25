# H and EPS implicit-coupling derivation

Date: 25-Nov-2026. Status: derivation and implementation design only.
No solver behaviour is changed by this document.

Purpose: define the discrete system that `MOM_BLKGMRES` must solve if
the lagged H and EPS terms are to be part of the implicit wavelength
boundary-value problem rather than frozen inside `NU_X0`.

## 1. Notation

Use scaled comoving moments at frequency index k:

- `x(i,k) = R(i)^2 J(i,k)` at depth node `i=1...ND`;
- `h(m,k) = R_mid(m)^2 H(m,k)` at midpoint `m=1...ND-1`;
- `uJ(i,k)` is the local upwind frequency index for the principal J
  frequency-advection coefficient;
- `uH(m,k)` is the local upwind frequency index for the H frequency
  derivative;
- `d(E,m)=m` for `E >= 0` and `d(E,m)=m+1` for `E < 0`, the current
  depth-upwind donor rule;
- coefficients `HU`, `HL`, `HS`, `EPS`, and `EPS_PREV` have the same
  meaning as in `mom_j_cmf_nonmono_v1.f`.

All coefficients are treated as frozen during one linear solve. The
Eddington factors remain outer-iterated quantities.

## 2. Current discrete H recurrence

Ignoring positivity clipping and non-finite repair, the production H
reconstruction at `mom_j_cmf_nonmono_v1.f:4161-4193` has the linear
form

```text
h(m,k) - HS(m,k) h(m,uH(m,k))
  = HU(m,k) x(m+1,k) - HL(m,k) x(m,k)
    + 2 EPS_PREV(m,k) x(d(EPS_PREV,m), uEprev(m,k))
    - 2 EPS(m,k) x(d(EPS,m),k).
```

The signs are carried by the coefficients and donor selectors. This is
the exact linear skeleton of the present implementation before
`MOM_EPS_LIM`, `MOM_LAG_BND`, `MOM_XM_POS`, H clipping, boundary
repair, relaxation, or Anderson acceleration.

The important structural point is that `h(m,k)` depends both on the
same-frequency J field and on the upwind-frequency J field. Therefore
H cannot be treated as a fixed source if the wavelength boundary-value
problem is being solved implicitly.

## 3. Current eliminated moment-row form

For an interior depth row, the present sweep assembles an equation of
the form

```text
TA(i,k) x(i-1,k) + TB(i,k) x(i,k) + TC(i,k) x(i+1,k)
  + CUP(i,k) x(i,uJ(i,k))
  + VB(i,k) h(i-1,uH(i-1,k))
  + VC(i,k) h(i,uH(i,k))
  + (Eup x)(i,k)
  = S(i,k).
```

Here `(Eup x)` denotes the midpoint-scattered lagged EPS contribution
assembled at `mom_j_cmf_nonmono_v1.f:2895-2955`. `TA/TB/TC` already
contain the current-frequency spatial elimination terms and the
current-frequency EPS donor contribution. Boundary rows replace the
first and last depth equations exactly as in the production sweep.

The present block capture puts only the `CUP x(uJ)` term into the
implicit frequency operator. The two H terms and `(Eup x)` remain in
`NU_X0`. That is why the current block solve is only a surrogate of
the full discrete moment system.

## 4. Why substitution is the wrong first implementation

One could substitute the H recurrence into the moment rows. Because
the recurrence contains `h(uH)`, repeated substitution follows the
frequency graph backwards. With alternating compression and expansion
that graph can change direction, contain long chains, or form rings.
Substitution therefore produces dense or solver-dependent frequency
coupling. It also mixes the elimination order with the physical
operator.

The correct sparse formulation is simpler: keep H as an unknown and
solve the moment and H-reconstruction equations together.

## 5. Augmented sparse target system

For each frequency k, use the interleaved unknown vector

```text
z(k) = [ x(1,k), h(1,k), x(2,k), h(2,k), ...,
         x(ND-1,k), h(ND-1,k), x(ND,k) ].
```

There are `2*ND-1` unknowns per frequency. The equations are:

1. existing depth-boundary rows for `x(1,k)` and `x(ND,k)`;
2. existing interior moment rows for `x(i,k)`, including the
   principal upwind-J edge and the H/EPS frequency edges;
3. the H reconstruction equation above for every midpoint `m`.

This system is sparse. It has depth-neighbour coupling, same-frequency
J-H coupling, and local upwind frequency coupling. It does not require
recursive substitution of H. It is the natural moment-equation analogue
of the literature's wavelength boundary-value problem.

## 6. Boundary data

The frequency boundary conditions must be explicit inputs to this
system:

- at a blue inflow edge, prescribe or model incoming upwind `x` and
  the corresponding upwind `h`;
- at a red inflow edge, do the same;
- at an outflow edge, no frequency boundary data are required.

The current implementation silently uses zero for off-grid compressing
reads. That is a boundary model, not a harmless numerical detail. The
augmented implementation should make the boundary model explicit and
separate it from the operator assembly.

## 7. Nonlinear safeguards are outside the linear operator

The following mechanisms make the map nonlinear or change its fixed
point and therefore cannot be hidden inside the linear operator:

- EPS donor limiting;
- LAG_BND scaling;
- XM_POS;
- negative-cell rejection;
- write caps and relaxation;
- H clipping;
- Anderson acceleration.

They may remain as outer nonlinear globalization, but the residual
used for acceptance must be evaluated after them. Ideally the clean
linear operator is first shown to reduce the need for them. If a guard
is still required permanently, its effect on the fixed point must be
reported as model error.

## 8. Implementation sequence

Phase 0, diagnostic only:

- implement a residual evaluator for the augmented system using the
  stored `J_NEXT_FRQ` and `H_NEXT_FRQ`;
- report the moment-row residual, H-row residual, raw norm, scaled
  norm, and worst cells;
- do not change the solver update.

This phase is small and safe. It establishes whether the existing
iteration is converging to the augmented equations or only to the
surrogate captured system.

Phase 1, matrix-free operator:

- implement `z -> A z` for the augmented sparse system;
- test it against finite differences and small manufactured fields;
- test local upwind direction reversal and red-edge inflow boundary
  data explicitly.

Phase 2, block GMRES:

- solve the augmented system matrix-free;
- begin with a simple block preconditioner: the existing depth
  tridiagonal approximation for moment rows and a diagonal or local
  frequency-bidiagonal approximation for H rows;
- report raw and preconditioned residuals separately.

Phase 3, globalization cleanup:

- measure post-projection residuals after damping, caps, and AA;
- apply one final acceptance rule to the composite field;
- retire or phase out equation-changing guards whose need disappears.

## 9. Immediate smallest code change

The first code change should be Phase 0 only: an augmented-residual
diagnostic. It does not alter convergence or the production path. It
creates the measurement needed before any implicit-H solve is worth
building. It also satisfies the project requirement that acceptance be
based on a trustworthy residual rather than solver success flags.

## 10. Open derivation work before Phase 1

1. Write the exact outer and inner spatial boundary rows in augmented
   form. The existing boundary rows are scalar replacements in the
   eliminated system and must be mapped carefully into `z`.
2. Decide whether the depth-upwind EPS donor stencil remains part of
   the target equations. The literature audit has not verified it.
   The augmented operator should initially reproduce the current
   stencil for regression purposes, then compare it with the centered
   V11 stencil as a separate controlled arm.
3. Define the frequency-boundary inflow data for both `x` and `h` at
   blue and red edges.
4. Choose row scaling. The moment and H equations have different
   units and magnitudes; raw GMRES convergence will depend on the
   scaling even when the operator is correct.
5. Decide how matched-frequency strides enter the H rows. The H
   upwind index currently follows the sign of the midpoint coefficient
   and may differ from the principal J stride.
