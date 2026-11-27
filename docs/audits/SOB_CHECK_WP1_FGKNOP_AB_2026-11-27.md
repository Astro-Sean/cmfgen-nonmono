# WP1 audit: DO_SOB_FG_KNOP vs clamped-V13 on m15 @ 2.82d

Dated 27-Nov-2026.  Status: PROVISIONAL - fgknop arm at 90%
convergence (iteration budget exhausted at it448); leg-2 warm
restart submitted (FINAL 13774145 -> OBS 13774146).  All numbers
below are the budget-exhausted state unless re-extracted after
leg 2; the verdict row will be updated when the restart lands.

## Question (SOB_check.md WP1)

Does solving the formal field with the sign-aware Knop solver
(`DO_SOB_FG_KNOP=T`, GS_REL_STABLE=2) change the CSM-line excess
seen under the clamped-V13 formal solver on the m15 interacting
model, and by what mechanism?

## Arms

| arm | run dir | state |
|---|---|---|
| ref (clamped V13) | `cmfgen_validation/m15_2p82d_VAL_sob/cmfgen/2.82days/` | converged 100% at it414 (24-Nov job) |
| fgknop | `cmfgen_validation/m15_2p82d_VAL_sob_fgknop/cmfgen/2.82days/` | it448, 90% conv, NUM_ITS spent; pops demand ~3% declining, worst margin 0.78 dex |

Luc reference: `cmfgen_validation/luc_bundle/m15_z1p4em2_ekin1p2_mdot0p01_rcsm6e14_nb10a_vt50.fl`.

## Resolved-control diff (OUTGEN banners)

The fgknop arm adds `DO_SOB_FG_KNOP=T`; auto-sets
`GS_REL_STABLE=2` and `DO_NONMONO_FREQ_MATCH=T`.  The
`DO_FG13_ALO=T`/`DO_ES_ALI=T` absent-key auto-sets do NOT fire
(DO_ES_ALI suppressed "when COHERENT_ES and not
DO_SOB_FG_KNOP"): 0 ES_ALI lines in the fgknop OUTGEN vs 402 in
ref.  All other resolved controls identical (VEL_LAW=14,
DO_NONMONO_SOB=T, USE_LAM_ES=T).

## (a) Per-depth J_nu in the CSM (JH_AT_CURRENT_TIME)

Extraction: `cmfgen_validation/wp1_jdepth_bfac.py`; grids
verified identical (max |dV| = 0).  Records are `RSQ_J(70) +
RSQ_H(69) + H_INBC + H_OUTBC + NU`; J = RSQ_J/R^2.  8331
frequencies, nu 0.003..200.0 (1e15 Hz).

At the O V -> O VI gs edge (nu = 27.54e15 Hz, 113.9 eV):

- mid/inner CSM (depths ~18-45, V = 73-468 km/s): fgk/ref =
  0.79-0.85 - the Knop solve carries ~15-20% LESS EUV J.  Same
  sign at 129.9 eV (0.63-0.81) and 150.1 eV (0.63-0.80).
- outermost cells (depths 1-12, V = 49-55 km/s, the near-CSM
  edge): fgk/ref = 10^3-10^5 - the clamped solve leaves the
  outer CSM nearly dark at EUV (ref J ~1e-18 at 113.9 eV, ~1e-25
  at 150 eV); Knop transports ~1e-13/~1e-19 there.  Boundary
  records confirm it is transport, not normalization:
  RSQ_J[0] 3.5e-8 -> 3.9e-3 (nu 27.5), H_OUTBC 0.89 -> 0.98 and
  0.94 -> 0.16 (nu 36.3).
- crossover at depth ~13-16 (V ~57-68 km/s).

## (b) Departure coefficients (POP* `XzV_F` blocks)

O V 2s3s ^3Se / ^1Se, 2s3p ^1P^o / ^3P^o (levs 11-16): fgk/ref
0.42-0.95 across the CSM - deepest deficit ~0.43-0.60 at the
outer edge, rising through ~0.85-0.95 mid-CSM and crossing 1.0
at the deepest CSM cells (V >~ 380 km/s).  Follows the EUV J
deficit: less photo-pumping -> lower excited-state populations.

C IV 3s ^2Se (lev 4) and 3p ^2P^o (levs 5-6): fgk/ref =
1.05-1.35 - OPPOSITE SIGN to O V (higher populations under
Knop), largest at the outer edge and ~1.06-1.13 through the
mid-CSM.

## (c) Nine line probes vs Luc .fl (OBS renders, 90% state)

`plot_fgknop_ab_wp1.py`, +-8 A window maxima; figure
`cmfgen_validation/plots/validation/fgknop_ab_wp1.png`.

| probe | ref/luc | fgk/luc | fgk/ref |
|---|---|---|---|
| O V 3145 | 1.248 | 1.219 | 0.977 |
| O V 5115 | 0.918 | 0.909 | 0.991 |
| He II 5411 | 1.043 | 1.042 | 0.998 |
| O V 5573 | 1.288 | 1.254 | 0.974 |
| O V 5584 | 1.371 | 1.334 | 0.973 |
| O V 5605 | 1.425 | 1.387 | 0.973 |
| O V 5616 | 1.435 | 1.397 | 0.974 |
| C IV 5801 | 1.274 | 1.318 | 1.035 |
| C IV 5812 | 1.416 | 1.466 | 1.035 |

## J-ratio spectrum by band (mid-CSM depths 20-45)

Median fgk/ref over frequency bands:

| band (eV) | depth 20 | depth 30 | depth 45 |
|---|---|---|---|
| 3-10 (FUV) | 0.93 | 0.93 | 0.93 |
| 10-20 | 0.89 | 0.89 | 0.89 |
| 20-40 | 0.88 | 0.88 | 0.89 |
| 40-60 | 0.92 | 0.92 | 0.93 |
| 60-90 | 0.74 | 0.77 | 0.90 |
| 90-130 | 0.78 | 0.78 | 0.83 |
| 130-250 | ~295* | 1.85 | 0.81 |
| >250 | >>1 | >>1 | >>1 |

*outer-depth ratios at 130-500 eV are inflated by ref J ~ 0 there.
Above ~250 eV the clamped solve transported essentially no
field into the CSM (J underflow); Knop carries real X-ray flux -
ratio is formally huge but absolute J is small in both.

## Interpretation (provisional)

The excess MOVES, but only by a few percent at line flux: the
Knop formal field carries ~15-20% less EUV J through the
bulk of the CSM (sign-aware transfer treats the compression
shell, sigma ~ -82, as genuinely optically thick to redshifted
photons), which lowers O V n=3 departure coefficients ~10-15%
and weakens the O V multiplet excess by ~2.5%.  C IV goes the
other way (~+7-13% DC, +3.5% line flux) - the 5801/5812 excess
is not reduced by the Knop solve.  Line-flux response is damped
relative to the DC change because the emitting weight sits in
cells where the ratio is nearer unity and the lines are partly
density-squared weighted.

The ~10^3-10^5x J boost at depths 1-12 is a real transport
difference (boundary H_OUTBC differs), but the absolute J there
is small in both arms.

The opposite signs of the O V and C IV DC shifts have a coherent
explanation.  The C IV boost is uniform across 3s and both 3p
levels (~1.05-1.35), which is the signature of an ionization-
fraction change, not an excitation change: the C IV -> C V
edge sits at 64.5 eV, inside the 60-90 eV window where the Knop
J is ~10-25% lower, so less C IV is destroyed and every C IV
level scales up together.  The O V n=3 suppression is an
excitation deficit: those levels are photo-pumped through the
114-150 eV edge region (deepest part of the J deficit), so the
excited states lose population while the underlying O V fraction
may itself rise.

## Open items before verdict

- leg-2 restart (13774145/13774146) queued on QOSGrpCpuLimit;
  re-extract (a)-(c) at 100% convergence and finalize the
  verdict row.
