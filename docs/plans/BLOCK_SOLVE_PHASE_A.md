# Phase A findings: the Baron & Hauschildt coupling setup

24-Nov-2026. Sources: Baron & Hauschildt 2004 (A&A 427, 987,
paper I); Knop, Baron & Hauschildt 2009 (A&A paper II,
arXiv:0909.4264) - the large-scale follow-up that mandates the
fully implicit discretisation. Paper II text extracted to
/tmp/bh2004/paper2.txt for working reference (not committed).

## The coupling setup (the user-requested specification)

The characteristic-form RT equation along path s:

```
dI_l/ds = eta_l - chi_l I_l - 4 a_l I_l - a_l d(lambda I_l)/d lambda   (Eq 1)
```

`a_l` is the velocity-structure coupling coefficient (Mihalas
1980). The wavelength derivative is discretised by LOCAL UPWIND
schemes keyed on the sign of a*lambda:

```
a*lambda >= 0 (expanding): backward difference, reads I_{l-1}
   p- = -lam_{l-1}/(lam_l - lam_{l-1})
   p| = +lam_l  /(lam_l - lam_{l-1})
   p+ = 0
a*lambda < 0 (compressing): forward difference, reads I_{l+1}
   p- = 0
   p| = +lam_l  /(lam_l - lam_{l+1})
   p+ = -lam_{l+1}/(lam_l - lam_{l+1})
```

Crank-Nicholson mixing with parameter xi in [0,1]: the xi-weighted
part is implicit, the (1-xi) part explicit. Generalised opacity
and sources (Knop et al. 2008 form):

```
chi-hat_l = chi_l + xi a_l p|_l                                      (Eq 5)
S-hat_l   = (chi/chi-hat) [ S_l - xi a_l/chi_l (p- I_{l-1} + p+ I_{l+1}) ]
S-tilde_l = -(a_l/chi-hat) [ (1-xi)(p- I_{l-1} + p+ I_{l+1})
             + (4 + (1-xi) p|_l) I_l ]
```

FULLY IMPLICIT = xi = 1: the diagonal upwind weight a*p| is
absorbed into the generalised opacity, and the off-frequency
neighbour terms I_{l+-1} are part of the implicit system - they
are NOT lagged. Paper II states outright why (Sect 2): "for
strong wavelength couplings the interpretation of the entire
wavelength discretisation as a source term is not valid for all
optical depths". That sentence is the literature statement of
the defect our three arms demonstrated numerically.

Also mandated in paper II: a formal solution avoiding negative
generalised opacities (Knop et al. 2008) - chi-hat can go
negative when xi a p| is large negative; this is the same sign
structure as our EPS/generalised-opacity census hits.

## Mapping onto this tree's assembly

Our compressing-row assembly already implements the xi-implicit
DIAGONAL exactly: PSI = |Q_k| on the diagonal is a*p| absorbed
into TB (paper Eq 5 with xi=1 at row level). What is lagged in
our code is only the NEIGHBOUR term coeff*J_upwind (p+ I_{l+1}
compressing / p- I_{l-1} expanding). The block solve is therefore
not a new discretisation - it is the SAME upwind discretisation
with the neighbour term moved from lagged-source to implicit-
unknown. Implementation consequence: the captured NU coupling
(NU_CUP/NU_KUP) already has exactly the right coefficients and
the right graph; the work is assembling them into one simultaneous
system instead of solving frequency-by-frequency with lagged
neighbours.

## Scaling constraint from paper II

Paper I's algorithm: memory scales with (wavelengths) x (layers)^2
entries in the Lambda-operator; practical ceiling ~64 layers x
1000 wavelengths on contemporary hardware. Paper II exists
precisely because of this - large-scale applications need the
iterative/operator-splitting algorithm of its Sect 3. For our
grids (ND ~ 60-100, NCF ~ thousands), a dense assembled operator
is infeasible; the phase-B census must therefore evaluate
matrix-free GMRES with the depth-tridiagonal inverse as
preconditioner (the per-frequency DGTSV we already have) rather
than any assembled direct solve.

## Implementation sketch (for phase C/D)

1. Unknown vector: J(l, i) over the frequency-depth grid, one
   direction set at a time.
2. Matrix-free operator action: for each frequency l, the depth
   tridiagonal (TA/TB/TC as now assembled, PSI in TB); the
   frequency coupling enters the action as coeff*J_upwind read
   from the CURRENT iterate - no store, no lag.
3. Preconditioner: block-Jacobi in frequency = the existing
   per-frequency DGTSV solve. Contraction in the compressing
   case comes from GMRES seeing the coupling in the action.
4. xi = 1 throughout (fully implicit, per the paper-II mandate);
   no explicit wavelength-source piece remains.
5. H closure, EPS scatter, opacities, source functions stay
   outer-iterated exactly as now.

## Verification anchor for phase C

Manufactured compressing case where the lagged map diverges
(drv_mom_nu iter infrastructure): the matrix-free GMRES block
solve must contract there. If it does not, the route dies at
phase C per the scoping plan.
