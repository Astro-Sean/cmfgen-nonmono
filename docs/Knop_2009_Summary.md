# Knop, Hauschildt & Baron (2009) — Summary

**Full citation:** Knop, S., Hauschildt, P. H., & Baron, E. (2009),
"A new formal solution of the radiative transfer in arbitrary velocity fields",
Astronomy & Astrophysics (Research Note), arXiv:0902.1876, DOI:10.1051/0004-6361:200810290

**PDF saved as:** `docs/Knop_et_al_2009_formal_solution.pdf`

## Abstract

A new formal solution of the Lagrangian equation of radiative transfer useful for
solving the equation of radiative transfer in the presence of arbitrary velocity
fields. Normally a term due to the inclusion of the wavelength derivative in the
Lagrangian equation of radiative transfer is associated with a generalised opacity.
In non-monotonic velocity fields, this generalised opacity may become negative. To
ensure that the opacity remains positive, this term of the derivative is included
in the formal solution of the radiative transfer problem.

The new definition of the generalised opacity allows for a new solution of the
equation of radiative transfer in the presence of velocity fields. It is especially
useful for arbitrary velocity fields, where it effectively prevents the occurrences
of negative generalised opacities and still allows the explicit construction of the
Lambda-operator of the system needed for an accelerated Lambda-iteration. Test
calculations showed relative deviations never exceeded 1%.

## Core Equations

### Transfer equation (Eq. 1)

```
dI_l/ds = eta_l - chi_l * I_l - 4*a_l * I_l - a_l * d(lambda*I)/dlambda
```

where `a_l = dlambda/ds` is the wavelength-coupling term.

### Wavelength coupling in spherical symmetry (Eq. 2)

```
a_l = gamma * (gamma^2 * mu*(mu+beta) * dbeta/dr + beta*(1-mu^2)/r)
```

where beta = v/c, gamma = 1/sqrt(1-beta^2), r and mu = cos(theta) are phase-space
coordinates. **It is the sign of the wavelength-coupling term, and not the sign of
the velocity gradient, that determines the upwind sense of the wavelength derivative.**

### Upwind wavelength derivative (Eq. 3)

For sorted wavelength grid lambda_{l-1} < lambda_l < lambda_{l+1}:

```
d(lambda_l * I_l)/dlambda = p_l^- * I_{l-1} + p_l^| * I_l + p_l^+ * I_{l+1}
```

where (for a_l >= 0, upwind = l-1):
  p_l^- = -lambda_{l-1}/(lambda_l - lambda_{l-1})
  p_l^| =  lambda_l/(lambda_l - lambda_{l-1})
  p_l^+ =  0

and (for a_l < 0, upwind = l+1):
  p_l^- =  0
  p_l^| =  lambda_l/(lambda_l - lambda_{l+1})
  p_l^+ = -lambda_{l+1}/(lambda_l - lambda_{l+1})

Key: p_l^| always has the same sign as a_l, so a_l * p_l^| >= 0 always.

### Mixed discretisation (Eq. 4)

With mixing parameter xi in [0,1] (Crank-Nicholson-like):

```
dI_l/ds = eta_l - chi_l*I_l - a_l*(4 + xi*p_l^|)*I_l
          - xi*a_l*(p_l^- * I_{l-1} + p_l^+ * I_{l+1})
          - [1-xi]*a_l*(p_l^- * I_{l-1} + p_l^| * I_l + p_l^+ * I_{l+1})
```

### Standard generalized opacity (Eq. 5) — FAILS for non-monotonic

```
chi_hat_l = chi_l + a_l*(4 + xi*p_l^|)
```

This goes negative when: -4*a_l > chi_l + a_l*xi*p_l^|  (Eq. 6)

### Knop's NEW generalized opacity (Eq. 7) — ALWAYS POSITIVE

```
chi_hat_l = chi_l + xi * a_l * p_l^|
```

This only incorporates the physical opacity and a purely positive contribution from
the velocity field (since a_l * p_l^| >= 0 by construction). The remaining 4*a_l*I_l
term is moved to the explicit source.

### Formal solution (Eq. 8)

```
I_l(tau_2) = I_l(tau_1)*exp(tau_1-tau_2)
           + integral(tau_1 to tau_2) S_hat_l * exp(tau-tau_2) dtau
           - integral(tau_1 to tau_2) S_tilde_l * exp(tau-tau_2) dtau
```

Note the **MINUS sign** before the S_tilde integral.

### Modified source functions (Eqs. 9-10)

```
S_hat_l   = (chi_l/chi_hat_l) * (S_l - xi*(a_l/chi_l)*(p_l^- * I_{l-1} + p_l^+ * I_{l+1}))
S_tilde_l = (a_l/chi_hat_l) * ([1-xi]*(p_l^- * I_{l-1} + p_l^+ * I_{l+1})
                               + [4 + (1-xi)*p_l^|] * I_l)
```

### Discrete integration (Eqs. 11-12)

The S_hat integral uses parabolic interpolation (alpha, beta, gamma coefficients
from Olson & Kunasz 1987). The S_tilde integral uses only LINEAR interpolation
(because it contains the explicit wavelength derivative).

## Key Algorithmic Steps

1. Compute the wavelength-coupling term a_l (Eq. 2) for each ray and depth.
2. Determine the sign of a_l to select the upwind wavelength direction.
3. Compute the upwind coefficients p_l^-, p_l^|, p_l^+ (Eq. 3).
4. Construct the generalized opacity chi_hat = chi + xi*a*p^| (Eq. 7) — always positive.
5. Construct the modified source functions S_hat and S_tilde (Eqs. 9-10).
6. Integrate along the characteristic using Eq. 8 with piecewise interpolation (Eqs. 11-12).
7. The formal solution is now a boundary-value problem (not initial-value) because
   the upwind direction can reverse.
8. The Lambda-operator construction from Baron & Hauschildt (2004) works without
   modification.

## Test Results

- Monotonic velocity fields: new solution matches old solution to < 1%.
- Scattering continua: differences ~1e-3 %.
- Non-monotonic (alternating) velocity field: required xi=1 (fully implicit) because
  xi=0 was ill-conditioned. The solution was robust and stable.
- The removal of the 4a*I term from opacity and addition as explicit source had
  "apparently little influence on the solution."

## Companion Paper

Knop, Hauschildt & Baron (2009b), "Comoving-frame radiative transfer in arbitrary
velocity fields - II. Large scale applications", A&A 501, 813 (arXiv:0909.4264).
Saved as `docs/aa11793-09.pdf`. This paper adds:
- Domain decomposition for memory reduction
- Gauss-Seidel type iterative solvers for speed
- Quasi-analytic formal solution
- Parallel computing support

## References Cited

- Baron & Hauschildt (2004), A&A, 427, 987 — original boundary-value algorithm
- Chen et al. (2007), MNRAS, 380, 104 — complete CMF transfer description
- Hauschildt (1992), JQSRT, 47, 433 — generalized opacity formulation
- Hauschildt & Baron (2004), A&A, 417, 317 — wavelength discretisation schemes
- Knop et al. (2007), A&A, 463, 315 — earlier formal solution
- Mihalas (1980), ApJ, 237, 574 — Lagrangian CMF transfer equation
- Olson & Kunasz (1987), JQSRT, 38, 325 — interpolation coefficients
- Scharmer (1981), ApJ, 249, 720 — ALI method
