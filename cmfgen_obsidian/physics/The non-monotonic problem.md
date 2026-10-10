# The non-monotonic problem

Why this tree exists. Physics of `sigma < -1` and what breaks in stock
CMFGEN.

## The velocity-coupling coefficient

CMFGEN solves the $O(v/c)$ comoving-frame transfer equation in nu-form:

$$ \frac{dI}{ds} = \eta - \chi I + a_\nu \frac{dI}{d\ln\nu} $$

with the velocity-coupling coefficient, per ray point,

$$ a_\nu = \Gamma = \frac{V}{R\,c}\left(1 + \sigma\mu^2\right),
\qquad \sigma = \frac{d\ln V}{d\ln R} - 1 . $$

The $dI/d\ln\nu$ term is Doppler shift along the ray: photons move in
comoving frequency as they propagate through a velocity gradient.

- **Monotonic expansion** (homology, $V \propto R$): $\sigma = 0$,
  $\Gamma > 0$. The characteristic in (space, frequency) marches
  redward. The equation is an **initial-value problem** in frequency:
  upwind is always the bluer neighbour, already computed on the
  downward frequency sweep.
- **Compression** ($V$ decreasing outward, or any region with
  $\sigma < -1$): $\Gamma$ **changes sign** at compressing cells
  ($1 + \sigma\mu^2 < 0$ for some $\mu$). Photons are transported
  *blueward* there. Upwind is now the *redder* neighbour - which on a
  blue->red sweep has **not been computed yet**. The IVP has become a
  **boundary-value problem** in wavelength.

> [!quote] Baron & Hauschildt 2004, A&A 427, 987 ([arXiv:astro-ph/0408212](https://arxiv.org/abs/astro-ph/0408212))
> "Since the equation is first order in wavelength at each spatial
> point there is only one boundary condition, whose wavelength sense
> depends on the local sign of the coefficient."

Analogy: the blue->red frequency sweep is like driving a one-way
street. In expansion all traffic flows one way, so each cell can read
off the car that already passed. In a compression zone the traffic
**reverses** locally - the cell now needs the car that has not arrived
yet. You either stop and wait (the full banded BVP solve of
Baron & Hauschildt 2004, which does not scale) or you drive past with
yesterday's photo of that spot and keep re-taking the picture until it
stops changing (Knop's lagged `I_NEXT` store + refresh iteration).

![[fig_sigma_velocity.png]]

## What stock CMFGEN does wrong

`FG_J_CMF_V13` builds (stock form)

```
chi_hat  = chi + GAM/dLOG_NU
Q        = GAM/(dLOG_NU*chi_hat)
S'       = eta/chi_hat + Q*I_prev
```

With $\Gamma < 0$ in a compression, `chi_hat` can go **negative**
(absorption into a "negative opacity"), `Q` points the wrong way, and
the solver returns negative or NaN intensities. Stock code also clamps
`SIGMA >= -1` in places (e.g. `MOM_J_CMF_V11`'s `W`/`EPS` terms), so
feeding it the true sigma corrupts those equations.

Sobolev line transfer has its own version: $|1+\sigma\mu^2|$ in the
escape probabilities goes through zero at $\sigma \to -1$, and the line
resonates at **three** Sobolev points on one ray instead of one.

![[fig_sobolev_points.png]]

![[knop_solver.png]]

## The physical setting: Type IIn

Fast ejecta run into slow, dense circumstellar material; a reverse
shock stalls the shell into a cold dense shell (Dessart, Audit &
Hillier 2015). Real measured numbers on the test models:

- `TypeIIn_10d` snapshot (~9.84 d): ND=169, sigma_min = -98, ~39% of
  zones with sigma < -1.
- `dessart_m15_2p82d` (Luc's 70-point hydro): sigma ~ -1.04 outer,
  ~ -157 at the interface; **all 70 depths have sigma < -1** and
  per-cell |dV| reaches ~1650 km/s - up to ~330 frequency bins of
  Doppler drift per cell at 5 km/s bin spacing (convergence_issues
  Sect. 14.1).

That drift/bins ratio matters: Ceulemans et al. (2025) show the
standard discretization admits oscillatory eigenmodes
($-1 < \lambda < 0$) when the per-cell drift exceeds the bin spacing -
deep extrapolation, not interpolation, on the upwind read.

## The three solver responses in this tree

| Approach | Flag | Idea | Status |
|---|---|---|---|
| Sobolev patch | `DO_NONMONO_SOB` | floor the vanishing factor; carry `SIGMA_TRUE` unclamped into line rates | works on m15 SOB runs; diverged on TypeIIn (see [[Sobolev non-mono solver]]) |
| Knop formal solver | `DO_NONMONO_KNOP` | sign-aware upwinding + lagged `I_NEXT` store (Knop 2009) | stable & contracting; see [[Knop formal solver]] |
| JFNK | `USE_JFNK` | replace the population solve with Newton-Krylov near the root | machinery complete; endgame admission not yet demonstrated, see [[JFNK admission]] |

## Key quantities the code carries

- `SIGMA` (`mod_cmfgen.f`): keeps its upstream infrastructure clamp -
  the moment equations and continuum formal solver are untouched.
- `SIGMA_TRUE` (`mod_cmfgen.f`): **unclamped** $d\ln V/d\ln R - 1$ from
  centred differences on R,V. Filled by `RECOMPUTE_SIGMA_TRUE`
  (`new_main/subs/recompute_sigma_true.f`) at setup, after
  `DO_WIND_VEL_V1`, after NG acceleration, and after grid revision.
- `GAM(I,LS)`: per-ray-point coupling coefficient; its **sign** selects
  the upwind direction everywhere in the Knop solver.
- `COMP_DEPTH_IDX` / `HAS_COMPRESSING`: which depths are on the
  compressing branch (sigma_true < -1 or V < 0, with exit-side
  hysteresis sigma > -0.98 / V > +1 km/s so the set cannot flap
  iteration to iteration).

## References

- [[Literature]] - Knop Papers I & II, Dessart+2015, Ceulemans+2025,
  Baron & Hauschildt 2004, Feldmeier & Nikutta 2006.
