# Sobolev non-mono solver (`DO_NONMONO_SOB`)

The Dessart, Audit & Hillier (2015, arXiv:1503.05463, Sec 3.1)
treatment of line transfer under a non-monotonic velocity field.
Flag: `T [DO_NONMONO_SOB]` in VADAT (optional, default F). Requires
`GLOBAL_LINE='SOB'`; `rd_control_variables.f` **hard-stops** otherwise
(enforced 03-Oct-2026, previously a warning). The `GLOBAL_LINE=NONE`
form is accepted only when every ion present has `TRANS_*=SOB` (note:
NONE requires a `TRANS_*` key for every species slot - a missing key
aborts in GET_KEY_STRING first). `USE_LAM_ES=F` was hard-rejected
until 12-Oct-2026; it is now an opt-in variant that routes the
continuum J through the stock `MOM_J_CMF_V11` on clamped `SIGMA` and
restores the `VAR_MOM_J_CMF_V12` dJ/dn, dJ/dT linearization in BA
(see "J-update mode" below). The Sobolev line path is unchanged.
`LAM_SET=T` only warns (redundant; it resets
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

## J-update mode: USE_LAM_ES=T vs F (12-Oct-2026)

Under `USE_LAM_ES=T` (pipeline default) the continuum J comes from the
FG ray solve and `do_var_cont.f` zeroes `TX/TVX/dJ_DIF_d_T/dRSQH` for
`DO_NONMONO_SOB` - no ALO store exists on that path, so the SE update
is a pure Lambda iteration for every continuum term. Measured on
`m15_2p82d_SOB_jfnk_zlaw14_fresh`: the endgame at depths 65-70 contracts
at ~0.9986/iteration, with ~850 of 1839 levels at d65-66 still demanding
>1% corrections - a broad optically-thick Lambda mode, not a trace-row
artifact (`CORRECTION_SUM` histogram, it ~780).

With `USE_LAM_ES=F` the J update now runs `MOM_J_CMF_NONMONO_V1`
(12-Oct-2026 wiring), the Knop Paper II sign-split moment solver
already used under `DO_NONMONO_KNOP`: `chi_hat = chi + |GAM|/dLOG_NU`
stays positive in compression zones and the upwind direction follows
the sign of GAM (`J_PREV` expanding, `J_NEXT` compressing). The GS
lagged stores `J_NEXT_FRQ`/`H_NEXT_FRQ` are read only under
`J_NEXT_VALID`, so a cold start simply runs expanding-coupled for the
first sweep; `MOMJNXT_GS` checkpoints survive warm restarts.
`do_var_cont.f` pairs it with fixed-J + `GET_MOM_ALO_DIAG` +
`GET_MOM_CELL_QUAR` - identical to the KNOP+LAM_ES=F recipe. The
Eddington factors (`N_ON_J` etc.) are computed live inside the solver
by `GET_MOMS_NON_REL` from `MOD_RAY_MOM_STORE`, which `FG_J_CMF_V13`
fills every frequency (all four moments, same accumulations as the
Knop formal solver). `SIGMA_TRUE` still drives the Sobolev lines and
feeds the solver's GAM term. JFNK admission is gated on `USE_LAM_ES`
(`solve_for_pops.f`) so the JFNK tiers stay inactive in this mode.
SN-model consistency check in `check_param_consistency.f` prints a
`USE_JREL or USE_DJDT_RTE` warning for `USE_LAM_ES=F`; benign.

**Refresh is mandatory but not sufficient in this mode.** The first
moment probe (job 12160057) ran without `MOM_KNOP_REFRESH_V1` - the
refresh driver was gated on `DO_NONMONO_KNOP`, so the compressing-zone
`J_NEXT`/`H_NEXT` map only matured one great iteration at a time.
Result: emergent luminosity ran away 3.9E8 -> 2.2E11 over 12
iterations with quarantined cells growing to ~7500/sweep.
`DO_KNOP_REFRESH` now defaults on for `DO_NONMONO_SOB .OR.
DO_NONMONO_KNOP` and `cmfgen_sub.f` calls both `FG_KNOP_REFRESH_V1`
and `MOM_KNOP_REFRESH_V1` under either flag (each self-guards when
its solver did not run).

With refresh enabled (job 12160690, commit bfa3291) the field stayed
bounded for ~5 iterations (~5-7E9 against a converged level ~1E10,
demands down to ~200%), then diverged anyway once the lagged map
matured: luminosity 1.7E10 (it446) -> 1.6E17 (it457), job cancelled.
GS-sweep diagnostics at the onset: `drift 4.2` (lagged map moving 4x
per sweep), `cmp_frac=0.945` (94.5% of cells compressing - the whole
field rides on the J_NEXT upwind), `max|EPS|/HU=5E6` with ~5M cap
events, `max|sigT|=113`. The refresh delta oscillated 0-1E5 without
converging. All population steps stayed LAMBDA-mode (demands under
LAM_ENTER=400), so the BA/ALO linearization never acted - the
divergence is a pure J-field runaway, not a Jacobian failure.

**KNOP control (job 12161004): the divergence is not SOB wiring.**
The same probe under `DO_NONMONO_KNOP+LAM_ES=F` (BLANK lines, Knop
formal ray moments, SOL_METH=TRI, same probe_snap warm restart)
diverged identically: bounded ~5e9 for 4 iterations, then 7.7E10
(it446) -> 9.6E12 (it448), steeper than the SOB run. Two differences
are still worth noting: under KNOP `max|EPS|/HU` stays ~11 (vs 5E6
under SOB - V13's ray-moment factors inflate the advection coupling
at shock cells), and the KNOP refresh intermittently reached
conv=T while the field still ran away, i.e. the map's fixed point
itself is a runaway state on this field from this restart.

**Cold-start probes (jobs 12161660 SOB / 12161668 KNOP): both
diverged.** SOB: bounded ~5E9 for ~3 iterations then 2.9E11 (it8)
-> 1.4E13 (it11). KNOP: it1=2.79E10 (vs September's 2.64E10),
then 3.1E10 -> 4.2E10 -> 4.95E11 (it4) with demands pinned at
~9.5E5 - where the September cold run DECLINED 2.64E10 -> 1.72E10.
So on the current build the moment path diverges on m15 from cold
AND warm starts under BOTH pathways. `mom_j_cmf_nonmono_v1.f` is
byte-identical to the stale pre-1ec74d5 snapshot, so there is no
intra-tree regression to bisect; the September convergence ran on
the `cmfgen_nonmono_dev_V7` worktree (gone), so config/tree
differences cannot be ruled out, but the practical conclusion is
unchanged: **the lagged-map moment solve is not a viable J-update
on m15 in the current build**. The converging pathways are the
formal-solve modes: SOB+LAM_ES=T (+DO_FG13_ALO) and KNOP+LAM_ES=T.

**Code audit of the USE_LAM_ES=F blow-up (13-Oct-2026).** The
earlier attribution of the SOB-vs-KNOP `max|EPS|/HU` gap (8.8E6 vs
~15) to clamped-SIGMA Eddington factors is CONTRADICTED by the
frequency grid: the ratio scales as 1/dLOG_NU, and the SOB
continuum grid (9665 points, read back from the EDDFACTOR records
of job 12161660) carries 1774 bound-free edge doublets at
dLOG_NU = 2.0E-11 (18% of the grid), while the KNOP line grid's
minimum spacing is 1.67E-5. 15 x (1.67E-5/2E-11) = 1.3E7 matches
the observed 9.1E6. The clamped-SIGMA inconsistency in V13's
factors is real but second order; it shows up as the sign
classification of the near-resonant outer region (cmp_frac 0.94
under SOB vs 0.71 under KNOP), not as the 1E6 coupling gap.

Defects found in `plane/mom_j_cmf_nonmono_v1.f` (all VERIFIED by
reading the executed path; their quantitative share of the runaway
is SUPPORTED, not measured):

1. *Edge doublets.* Every frequency-advection coefficient (GAM,
   GAMH, CHI_HAT, PSI, EPS) carries 1/dLOG_NU. At a 2E-11 doublet
   they are ~1E8x the physical terms. Stock V11 survives because
   its implicit and explicit terms are built from the same
   coefficients and the row degenerates to J_k = J_{k-1}. The
   nonmono solver does not degenerate cleanly: at compressing rows
   the J-equation reads the lagged J_{k+1} while the EPS term reads
   the bluer J_{k-1}, and the limiter in (2) breaks the
   implicit/explicit balance. A relative O(1) error on a term 1E7x
   the physics is the SOB energy source. KNOP never saw doublets.
2. *MOM_EPS_LIM is non-conservative.* The donor-cell limiter (lines
   ~2319-2347) caps only the explicit lagged pair (+2E*J_up on row
   m, -2E*J_up on row m+1); the implicit pair on TA/TB/TC (lines
   ~2101-2110) is uncapped. The two pairs together are the
   frequency derivative of sigma*n*J; capping one side turns a
   derivative into a net source/sink on rows where |EPS| exceeds
   the local RHS - 4.1M cap events per SOB sweep, 2.2M under KNOP.
3. *N_ON_J disables the midpoint sign split.* `GET_MOMS_NON_REL`
   sets NMID_ON_HMID = 0 under N_ON_J, so GAMH_FAC = 1 + sigma*0 = 1
   and GAMH_SIGNED >= 0 at every midpoint. The compressing-midpoint
   branches (H_NEXT_FRQ read, HS = W/(1+W), EPS_PREV = EPS, JUP read
   for the EPS source) documented in convergence_issues Sect. 7 are
   dead code in every N_ON_J run. The sigma*n*J term - the whole
   velocity-gradient part of the H equation in this mode, and the
   dominant term at shock cells (|EPS|/HU >> 1) - is always
   differenced toward the bluer frequency, which is the downwind
   direction wherever its coefficient sigma*n < 0 (negEPS = 1.85M of
   1.95M midpoints under SOB). Downwind differencing of a dominant
   advection term is anti-diffusive. Both pathways run N_ON_J.
4. *Zero-GAM resonance over the outer envelope.* At depths 1-8
   sigma_true ~ -1.04 and K_ON_J ~ 0.92, so 1 + sigma*K ~ +-0.04:
   GAM_FAC sits at the 0.01 floor (max|GAM| = 1.00E2 exactly in
   every sweep summary) and the upwind direction of these rows
   flips with small changes in f between frequencies and
   iterations. The negative-J dump (MOM_J_NONMONO_ERRORS) is
   dominated by these rows, with TB ~ |TA|+|TC| (no diagonal
   excess) so their J is set entirely by neighbours and lagged
   terms.
5. *Sweep J vs stored J.* The J returned to the population solver is
   the raw solve, while J_NEXT_FRQ receives a relaxed, suppressed
   (94% of rows on the SOB grid under MOM_DOM_CMP) or frozen value.
   Until the refresh converges (it did not: delta 8.7E4) the
   populations equilibrate to a field that is not the lagged map's
   input.

The September cold run (convergence_issues Sect. 8) ran 3
iterations from zero lagged stores; the compressing coupling only
matures at it2-3, so it never tested the regime where today's runs
diverge. Treat it as consistent with, not contrary to, these
findings.

Fix candidates, in order of expected leverage: (a) short-circuit
doublet steps - for dLOG_NU below a floor (~1E-8) copy J, H and the
stores from the previous frequency instead of solving; (b) make the
limiter conservative by scaling EPS_NM and EPS_PREV_NM of a midpoint
by the same factor so the implicit and explicit halves stay a
derivative; (c) give N_ON_J a midpoint direction indicator (sign of
1 + sigma*N/H from the stores where |H| is resolved, else sign of
the adjacent depth-point Q) so the Sect. 7 compressing branches
engage; (d) revisit the resonance floor so near-zero-GAM rows are
treated as non-advective (coefficient -> 0) rather than floored to
a direction. None implemented yet.

**V11 probe result (12-Oct-2026): stock moment solver not viable on
m15.** The earlier wiring ran `MOM_J_CMF_V11` + `VAR_MOM_J_CMF_V12`
under `USE_LAM_ES=F`; a warm restart from `probe_snap` died in
`comp_j_blank.f` after ~20 s - V11 fails to converge f at the bluest
frequencies (log nu ~ 11.4-11.6, J ~ 1E-60) and the `|RJ(1)| > 1E30`
guard trips. Both `N_TYPE=N_ON_J` (job 12158938) and `G_ONLY` (job
12159013) failed identically: the stock discretization assumes
expanding-flow upwind and is unstable on this field. That failure is
what motivated routing `USE_LAM_ES=F` to the sign-split solver
instead; pairing it with `VAR_MOM_J_CMF_V12` would also have been
inconsistent (V12 linearizes the V11 system), so the fixed-J + MOM
diagonal pairing is used, as under KNOP.

## Opt-in DO_FG13_ALO (12-Oct-2026, commit 336e899)

`FG_J_CMF_V13` now stores the same diagonal Olson-Auer-Buchler
approximate Lambda operator `FG_J_CMF_NONMONO_V1` uses on the Knop
path: per-ray captures of `chi_hat = CHI + GAM/dLOG_NU` and the
direction-dependent modified sources `S'_M/S'_P`, accumulated in the
serial J_LOOP into `dJ/dETA`, `dJ/dCHI` and `LAMSTAR` diagonals with
the same `JQW` weights and `J_PNT` mapping as `JNU_STORE`. Exposed to
`do_var_cont.f` via `GET_FG13_ALO_DIAG`; under `DO_NONMONO_SOB` +
`DO_FG13_ALO=T` the continuum `TX(K,K,1:2)` diagonals are filled, so
the BA step regains continuum self-coupling while J still comes from
the formal solve (`USE_LAM_ES=T` unchanged). Default F keeps the
validated all-zero continuum Jacobian. The construction assumes
`dS'/dCHI = -S'/chi_hat` (Q and I_prev frozen), skips rows where
`chi_hat <= 0` (possible in compression zones where the unsigned GAM
term dominates), and runs only when `SOLUTION_METHOD='INTEGRAL'` -
otherwise `ALO_VALID=F` and that frequency falls back to the plain
Lambda update.

Probe `dessart_m15_2p82d_SOB_jfnk_alo_probe` (job 12159823,
interactive partition) warm restarts from `probe_snap` with
`DO_JFNK=F`; `SOB ALO (FG)` activation confirmed at it444/446 with
the expected signs (dJ/dCHI < 0, dJ/dETA > 0). Note `probe_snap`
predates the EDDFACTOR warm-restart fix, so the first iteration
recomputes factors (one bounded cold start); the measurement window
is the plateau after recovery.

Companion pipeline defect found in the same audit: every SINGLE_FINAL
requeue rewrote `VADAT` leaving `COMP_F=T`, so CMFGEN recomputed
EDDFACTOR and zeroed `RJ` at the first frequency of each restarted
segment - the first post-requeue iteration then demanded 100-300%
corrections at depths 9-18 and ~30% of all iterations were spent
re-converging. Fixed in `cmfgen_base/submit_cmfgen.sh` (requeue block
now flips `COMP_F` T->F when EDDFACTOR/EDDFACTOR_INFO exist, matching
the FINAL template); see [[Run log highlights]].

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
