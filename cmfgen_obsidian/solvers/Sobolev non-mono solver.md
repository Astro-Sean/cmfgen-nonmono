# Sobolev non-mono solver (`DO_NONMONO_SOB`)

The Dessart, Audit & Hillier (2015, arXiv:1503.05463, Sec 3.1)
treatment of line transfer under a non-monotonic velocity field.
Flag: `T [DO_NONMONO_SOB]` in VADAT (optional, default F). Requires
`GLOBAL_LINE='SOB'`; `rd_control_variables.f` **hard-stops** otherwise
(enforced 03-Oct-2026, previously a warning). The `GLOBAL_LINE=NONE`
form is accepted only when every ion present has `TRANS_*=SOB` (note:
NONE requires a `TRANS_*` key for every species slot - a missing key
aborts in GET_KEY_STRING first). `USE_LAM_ES=F` is also rejected: it
routes J through the clamped-sigma `MOM_J_CMF_V11`, outside the
validated SOB setup. `LAM_SET=T` only warns (redundant; it resets
GLOBAL_LINE_SWITCH to NONE). The main VADAT reader keeps SOB and `DO_NONMONO_KNOP` mutually
exclusive. In the OBS reader, both flags may be T: `DO_NONMONO_KNOP`
selects Knop formal J while SOB lines and `SIGMA_TRUE` EW estimators
preserve the population line treatment. `prepare_CMF_PARAM.py` carries
`COH_ES=F` into OBS when requested, so non-coherent electron-scattering
wings are iterated.
`GRID_METH='VEL_SPACE'` is rejected inside `ADJUST_R_GRID_V4`
(05-Oct-2026): `DO_VEL_REGRID_V2` interpolates R(V), which assumes a
monotonic velocity field. TAU_SPACE, FULL_R and OB_ONLY remain
allowed.

Two further REV_RGRID hazards were closed (13-Oct-2026). The
interpolated-sigma check in `ADJUST_R_GRID_V4` clamped sigma < -1
only for V < 0.2 km/s and STOPped otherwise - on a shocked field
sigma < -1 sits at thousands of km/s, so any revision would abort.
The check is now skipped under the non-mono flags: the interpolated
SIGMA is transitional anyway since `RECOMPUTE_SIGMA_TRUE` reclamps
it under SOB right after the routine returns, and KNOP handles
sigma < -1 natively. Second, `BA_MATRIX_STORED` is now cleared on
`R_GRID_REVISED`: a stored BAMAT's depth-column indices no longer
map to the same radii, so a `FIX_BA`/`N_FIX_BA` recall would apply
corrections to the wrong cells. That fix is not flag-gated - a
stale matrix is wrong on every path that revises the grid.

## What it changes

- **`SIGMA_TRUE`** (`mod_cmfgen.f`): unclamped $d\ln V/d\ln R - 1$ from
  centred differences on R,V. `SIGMA` keeps its infrastructure clamp
  so the moment equations and continuum formal solver are untouched.
  `RECOMPUTE_SIGMA_TRUE` fills it at setup and after every R/V change;
  under `DO_NONMONO_SOB` it also re-clamps `SIGMA >= -0.999` because
  the model readers leave `SIGMA` unclamped and `sigma < -1` reaching
  stock routines (e.g. `MOM_J_CMF_V11` `W`/`EPS`) breaks them. Not
  applied under `DO_NONMONO_KNOP` (commit `5c85085`).
- **Tau floor** in `sobjbar_sim.f`, `sobew.f`, `sobew_grad_v2.f`,
  `sobew_grad_v3.f`: for sigma < 0,

  $$ |1+\sigma\mu^2| \to \max\left(|1+\sigma\mu^2|,\; -0.2(1+\sigma),\; 10^{-4}\right) $$

  Bit-identical to upstream for $\sigma \ge -1$.
- `SUB_SOB_LINE_V3` and the dielectronic `SOBEW` call in
  `cmfgen_sub.f` pass `SIGMA_TRUE` when the flag is set.
- **Transient-iterate guards** in `SUB_SOB_LINE_V3` (commit `b107b4b`):
  `VB_SIM`/`VC_SIM`/`SRAT` zero when `ETAL_MAT`, `BB_COR`, `CHIL` go
  non-positive; self-line derivatives use expanded forms avoiding
  division by `POPS(NUP)` and `U_STAR_RATIO`. Without them a bad
  Lambda iterate puts Inf/NaN or wrong-sign entries into the BA matrix
  and the correction runaways. Identical results for valid iterates.

## Non-local Sobolev coupling (`DO_NONLOC_SOB`, default F)

`subs/sob_rh_couple.f` - Rybicki & Hummer (1978) coupling between
distant Sobolev resonance points, in the Feldmeier & Nikutta (2006)
geometry-free form. Photons trapped at sigma ~ -1 get a real
escape/reprocessing channel at remote resonance points.

$$ \delta J_l(r) = \tfrac12 \int dr'\, r' v'\,
  \Gamma^H_l(r')\, S_l(r')\, K(r,r') $$

with the two-point kernel $K$ built from unclamped $|1+\mu^2\sigma|$
(the Dessart floor applies only to the local beta integral). Applied
as a lagged rate correction `DRC_SIM = CHIL*dJ/BB_COR` subtracted in
the STEQ residual only - **not** folded into `ZNET_SIM`/`JBAR`, which
was the Jacobian-contaminating bug of revision 1 (folding a lagged
remote term into `ZNET_SIM` gave the BA matrix a derivative it does
not possess).

> [!warning] Expansion-side coupling only (05-Oct-2026)
> RH78 assumes a monotone expanding flow. `SOB_RH_GEOM` now drops any
> pair with `V(I)<=0` or `V(J)<=0` (the mu^2 geometry uses V^2 and
> the pair weight carries signed V(J)) and reports the drop count;
> `GAMH_L` uses `|V|` with a 1 km/s floor so a stagnation depth can
> no longer leave Inf in the table. Infall regions get no non-local
> correction - a documented limitation, not the full non-monotonic
> kernel.

## Runtime record

- **m15 2.82d** (`dessart_m15_2p82d_SOB.yaml`, job 12001805): stable
  Lambda contraction from cold LTE start - max % increase
  1.6e11 -> 64 -> 17 -> 4.4 (it24); luminosity pinned
  L=4.266e8 inner / ~1.04e10 outer. The load-bearing setup, from
  Luc's converged VADAT: `LSTAR` = the SN_HYDRO_DATA header luminosity
  (4.265857e8 Lsun - with `LSTAR=0` the inner boundary is dead and the
  SE iteration has no flux anchor), `DIF=T` + `IB_METH=DIFFUSION`
  (inner tau ~37, the paper's tau~30 regime), `ALLOW_OL=F`,
  `DC_METH=LTE`, `THK_LINE=T`, `DO_DIS=T`, `INC_PEN=T`, `LAM_VAL=200`,
  `FIX_BA=10`, `N_FIX_BA=2`, `INC_CHG=F` (no CHG_EXCH_DATA in repo).
- **TypeIIn 10d**: diverged both with the plain floor and with
  `DO_NONLOC_SOB` - abort ~1e107 corrections at it11. The Sobolev
  line-rate map itself is anti-contractive there: remote absorption
  demands populations that cannot balance it, each Lambda iterate
  amplifies. Conclusion: the escape-probability line formulation is
  the diverger for this sigma structure; Knop+BLANK is the path for
  TypeIIn.
- **Late m15 SOB behaviour** (jobs 12065584/12067329): a period-7-ish
  `FULL -> LAMBDA -> LAMBDA` limit cycle. Each full step demanded the
  saturated 1e7% decrease (H II lev1 ~ depth 61-70) and rebounded the
  residual ~2.4x; `NUM_LAM: 50` let the `LAM_N_STABLE` merit gate
  replace the `cnt_lam` escape leash; `BA_ENTER_RAW` blocks the
  Lambda->full transition while raw demand stays saturated. Residual
  floor declining ~10%/cycle.

## Related

- [[The non-monotonic problem]]
- [[BA matrix and SOLVEBA]] - `BA_ENTER_RAW`, resacc, trust radius
- [[Run log highlights]]
- [[Literature]] - Dessart+2015, RH78, Feldmeier & Nikutta 2006
