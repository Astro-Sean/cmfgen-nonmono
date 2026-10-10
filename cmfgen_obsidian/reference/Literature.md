# Literature

Papers the design cites, and what each contributes.

## Non-monotonic comoving-frame transfer

- **Knop, Hauschildt & Baron 2009, A&A 496, 295 (Paper I)** -
  positive-generalized-opacity formal solution:
  `chi_hat = chi + |GAM|/dln nu`. The discretization this tree
  implements (audited in [[KNOP_MATH_AUDIT]]).
  [arXiv:0902.1876](https://arxiv.org/abs/0902.1876) ·
  [DOI](https://doi.org/10.1051/0004-6361:200810290)
- **Knop, Hauschildt & Baron 2009, A&A 501, 813 (Paper II)** -
  large-scale implementation: physically ordered Gauss-Seidel march,
  converged-BVP-before-SE ordering, ALI with wavelength-tridiagonal
  Lambda*. The ordering rule the whole control layer enforces.
  Key claim: "the solution is then direct and its speed is optimal"
  (Sec 3.2, on the ordered march).
  [ADS](https://ui.adsabs.harvard.edu/abs/2009A%26A...501..813K) ·
  [DOI](https://doi.org/10.1051/0004-6361/200911793)
- **Baron & Hauschildt 2004, A&A 427, 987** - the original non-mono
  CMF BVP algorithm: assemble the full wavelength-space system and
  solve exactly (SuperLU). The reference dev7's iterative march
  approximates.
  [arXiv:astro-ph/0408212](https://arxiv.org/abs/astro-ph/0408212) ·
  [DOI](https://doi.org/10.1051/0004-6361:20040660)
- **Hauschildt & Baron 2004, A&A 417, 317** - wavelength
  discretizations; xi=1 fully implicit is required for strong
  coupling ("the interpretation of the entire wavelength
  discretization as a source term is not valid for all optical
  depths").
  [astro-ph/0401164](https://arxiv.org/abs/astro-ph/0401164) ·
  [DOI](https://doi.org/10.1051/0004-6361:20034473)
- **Ceulemans, De Ceuster & Decin 2025, JQSRT 343, 109470** -
  oscillatory eigenmode analysis (-1 < lambda < 0) when per-cell
  Doppler drift exceeds bin spacing; the frequency-matching fix
  (`DO_NONMONO_FREQ_MATCH`). Also: storing/reusing computed
  intensities at re-encountered frequencies is the correct boundary
  treatment - our lagged I_NEXT store is "the correct variant".
  [arXiv:2505.08309](https://arxiv.org/abs/2505.08309)
- **Dessart, Audit & Hillier 2015, MNRAS 449, 4304** - the Type IIn
  interaction model and Sobolev tau-floor treatment
  (`DO_NONMONO_SOB`).
  [arXiv:1503.05463](https://arxiv.org/abs/1503.05463) ·
  [DOI](https://doi.org/10.1093/mnras/stv609)
- **Rybicki & Hummer 1978, ApJ 219, 654** + **Feldmeier & Nikutta
  2006, A&A 446, 661** - non-local Sobolev coupling across distant
  resonance points (`DO_NONLOC_SOB` / `sob_rh_couple.f`).
- **Olson & Kunasz 1987** - tridiagonal/pentadiagonal approximate
  operators give 2-7x fewer ALI iterations (motivates the untested
  `JFNK_USE_BANDED` path).
- **Trujillo Bueno & Fabiani Bendicho 1995, ApJ 455, 646** - GS/SOR
  inside the formal march at source-function level (a different object
  than this tree's "GS" store lag; possible future inner-map
  acceleration).

## JFNK / Newton-Krylov

- **Knoll & Keyes 2004, JCP 193, 357** - JFNK + pseudo-transient
  continuation (`JFNK_PTC`). Newton "most effective near the
  solution"; Krylov needs preconditioning - the two facts the
  admission gates exist to exploit.
  [DOI](https://doi.org/10.1016/j.jcp.2003.08.010)
- **Kelley & Keyes 1998, SIAM J. Optim. 9, 339** - PTC convergence
  theory (succeeds where line search/trust region stagnate).
  [DOI](https://doi.org/10.1137/S1052623496303961)
- **Arramy, de la Cruz Rodriguez & Leenaarts 2024, A&A 690, A12** -
  JFNK for multilevel SE: GMRES rtol sweet spot 3e-4..3e-2 (our eta
  floor), block-Jacobi preconditioner inadequacy for many-frequency /
  strong-scattering problems, "few operator-split iterations before
  Newton" prescription (= our admission strategy), line search +
  continuation as missing robustness. Abstract: JFNK "converges
  faster than our reference case by up to a factor 2" in
  formal-solution evaluations - the reason JFNK here is an endgame
  accelerator, not a replacement.
  [arXiv:2406.17234](https://arxiv.org/abs/2406.17234) ·
  [DOI](https://doi.org/10.1051/0004-6361/202449963)
- **Janett et al. 2024, A&A 682, A48** - physics-based (lightweight
  model) preconditioners cut Krylov counts 10-20x; candidate direction
  for the coupled-mode preconditioner (source-iteration
  preconditioner).
  [DOI](https://doi.org/10.1051/0004-6361/202348048)
- **Pernice & Walker 1998 (NITSOL)** - FD step sized to sqrt of the
  function noise level (`JFNK_EPS` rationale).
- **Hauschildt & Storzer 1994** - optimal ALI bandwidth ~5-15
  (preconditioner bandwidth note).
- **Chen & Shen 2006** - preconditioner update from Krylov products
  (deliberately not implemented).
- **Lambert et al. 2015** - false-convergence stalls of stationary
  iterations (the plateau detector's motivation).

## Verification status of the mapping

- Structure verified against the papers: chi_hat/upwinding (Paper I),
  ordering rule + ordered march (Paper II), PTC, EW band, noise-aware
  FD step - see [[KNOP_MATH_AUDIT]] and
  [[KNOP_CONVERGENCE_STRATEGY]] claim tags.
- Empirical (this tree's own calibration): all thresholds (5%/50%,
  ~0.1% bulk, churn fraction, gate values 1e3/1e2, NREQ=2) - the
  implementation's expression of "sustained contraction", not derived
  numbers.
