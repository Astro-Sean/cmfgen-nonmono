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

Fixes implemented (13-Oct-2026, this commit):

- (a) Spacing floor instead of a doublet short-circuit: MOM_DLN_J /
  MOM_DLN_H are floored at MOM_DLN_MIN=1E-6 where they feed
  CHI_HAT / GAMH / GAM.  The doublet rows then act as a small
  physical step (bounded coefficients) rather than being skipped,
  so all store/ALO/refresh bookkeeping is unchanged.  Real
  spacings (min ~1.7E-5 on the KNOP grid, ~1E-5 on SOB) are
  untouched.
- (b) Conservative EPS limiter: the donor-cell cap now records a
  per-midpoint factor EPS_SCL = capped/raw, applies (EPS_SCL-1)*2*EPS
  to the implicit TA/TB/TC stencil, then rescales EPS_NM and
  EPS_PREV_NM so the boundary row, negJ rebuild and H reconstruction
  all see the capped coefficient.  Implicit and explicit halves stay
  a derivative.  Inactive when MOM_EPS_LIM=F.
- (c) N_ON_J midpoint direction: GET_MOMS_NON_REL fills
  NMID_ON_HMID = N/H under N_ON_J where |H| resolves (|H| > 1E-3|N|),
  else a signed 1E10 proxy giving 1+sigma*g the sign of sigma*N.
  Opt-in via GMNR_NMJ_G in MOD_RAY_MOM_STORE (set by
  MOM_J_CMF_NONMONO_V1; stock V8-V11 callers keep the zero fill).
  The sigma*n*J frequency upwind is now keyed on sign(EPS) itself
  (scatter, boundary row, negJ rebuild, H reconstruction) and
  EPS_PREV takes the current weights at EPS<0 midpoints -- a flux
  split with one direction per term.
- (d) not implemented: the zero-GAM resonance (1+sigma*K ~ +-0.04
  across depths 1-8) still floors GAM_FAC at 0.01; direction there
  remains data-driven.

Residual: the EPS magnitude still uses the H-part spacing
MOM_DLN_H when its own direction is compressing while
GAMH_SIGNED >= 0 (resolved-H mixed-direction midpoints) -- a
second-order coefficient inconsistency, documented for the next
iteration if needed.

**Lagged-map gain fixes (13-Oct-2026, commits 66756a3, 6373327)
#status/unverified.** Two changes target the observed expansive
channels of the lagged map; both are VADAT-gated inside
`MOM_J_CMF_NONMONO_V1` only and leave the monotonic path untouched:

- `MOM_GAM_PIN` (default T): the lagged comoving-gradient factors
  `GAM_FAC_PREV`, `GAMH_FAC_PREV`, the Eddington boundary
  `PSIPREV`, and the `EPS_PREV` angle factor are evaluated at the
  CURRENT frequency's Eddington factors instead of the bluer
  neighbour's.  Near the `1+sigma*K=0` resonance the two
  evaluations drift in sign and magnitude, which let
  `PSIPREV/PSI > 1` -- the explicit upwind weight then exceeds the
  implicit diagonal and the lagged map amplifies.  With the pin
  `PSIPREV=PSI`, `WPREV=W`, `EPS_PREV=EPS` at every point, the
  compressing branch's existing convention extended to all
  channels.  F restores the independent lagged evaluation.
- `MOM_LAG_BND` (default T): per-row bound on the lagged frequency
  source.  A row's lagged coupling weight (upwind-J coefficient +
  `|VB|+|VC|` H terms + `2|EPS_PREV|` midpoint transfers) may not
  exceed its diagonal excess `TB-|TA|-|TC|`; where it does, the
  whole lagged contribution (accumulated in `XLAG_NM`) is scaled
  by `min(1, excess/weight)`.  Applied to interior rows, both
  outer-boundary branches, and the negative-J rebuild (via
  `LAG_SCL_NM`).  Rows whose unbounded map was already contractive
  keep factor 1.  The H reconstruction stays unscaled - it solves
  the eliminated H equation, not the bounded J-row map.
  Diagnostics: `lagbnd scaled=<n> mins=<s>` line in the sweep
  report.

Run-testing in progress (13-Oct-2026): four probes on binary frozen
at `bdf2f10` (sha 98cc0b41fccec88c) -- SOB cold job 12167354, KNOP
cold 12167426 (resubmit of diverged 12165640, which was pinned to
pre-fix `32b0d2a` and cancelled), SOB warm 12167431, KNOP warm
12167439.  Early SOB cold signal (it1-6): outer luminosity bounded
and falling 8.3E6 -> 7.6E6 (previously ramped x1.15/it past
physical to 4E21), zero negJ/quar repairs, refresh bulk-converged
from it3.  `lagbnd scaled` ~1.4M entries/sweep with `mins=0` --
rows with no diagonal excess get their whole lagged source nulled;
whether the clip fraction contracts as the field settles is the key
consistency question.  NOT yet classified converged: the unfixed
SOB cold divergence was also slow over its first tens of
iterations.

(b) implicit H coupling and (d) the
KNOP `ESEC>CHI` pure-scattering channel remain candidate follow-ups.

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

## OBS render: DO_SOB_FLUX_OBS (12-Oct-2026)

`TRANS_TYPE='SOB'` lines are excluded from the cmf_flux formal
transfer (`VEC_TRANS_TYPE .NE. 'BLAN'` skip in `cmf_flux_sub_v5.f`), so
without a deposition path `obs_fin_10` contains continuum only - a
`DO_SOB_LINES=F` + `SOB_FREQ_IN_OBS=F` OBS run produced 815
continuum frequencies and no lines (job 12162809). `DO_SOB_LINES=T`
alone only writes `EWDATA`; it never touches `OBS_FLUX`.

`CMF_FLUX_PARAM` key `DO_SOB_FLUX_OBS` (default F) enables the
resolved Sobolev render. Inside the `DO_SOBOLEV_LINES` blend-group
loop, `ADD_SOB_LINE_PROF` writes each `TRANS_TYPE='SOB'` member's
profile-integrated opacity/emissivity (`CHIL_MAT`/`ETAL_MAT`, scaled
by `CLUMP_FAC`) into the stored `CHI_CMF_ST`/`ETA_CMF_ST` as a
normalized Doppler core of width `FL*VDOP_MIN/C_KMS` at the line's own
frequency. The core weights are normalized by local CMF cell width and
carry the `1.0D-15` conversion the per-Hz `CHIL_MAT`/`ETAL_MAT` units
require (same convention as `LINE_PROF_SIM` in `SET_PROF_V6`), so the
stored opacity integrated over Hz equals `CHIL_MAT*CLUMP_FAC`
regardless of how coarsely the CMF grid resolves the line. Masering
cells (`CHIL_MAT<=0`) contribute neither opacity nor emissivity; the
dropped-cell count is reported per blend group. After the line loop
`OBS_FRAME_SUB_V9` is re-run on the augmented
arrays, so the observer's-frame ray integral draws each profile from
the isovelocity geometry (Dessart et al. 2015, Sec. 3.1) and
`OBSFRAME` is rewritten; `OBSFLUX` keeps continuum+BLANK only as a
comparison product. `ETA_DATA`/`CHI_DATA` on disk are written before
the EW block, so a standalone `obs_frame.exe` render still lacks SOB
lines.

Two bugs were found by comparing against a BLANK control render on the
same checkpoint (job 12166771, max 1.6e7 Jy - sane): (a) the first
version floored the *total* `CHI_ST` at 1e-30 when a masering line's
negative opacity overwhelmed the continuum, driving `S=ETA/CHI` to
~1e30 and producing 1e18 Jy spikes (job 12166626); fixed by dropping
maser cells entirely (2db3338). (b) The profile weight omitted the
1.0D-15 Hz conversion, injecting 1e15x the correct opacity so every
line saturated at its source function - a ~200x elevated
pseudo-continuum remained after fix (a) (job 12167387); fixed in
60c707d.

This replaces the first-cut EW deposition (`DEP_SOB_OBS_FLUX`,
superseded 12-Oct-2026), which spread each line's `SOBEW_GRAD_V2` EW
over an invented flat/one-sided box in wavelength; it produced boxy
emission and comb-like UV artifacts (job 12165745) because the
profile shape was assumed, not computed by the formal integral.

Motivation: on the non-monotonic m15 field the CMF/BLANK render emits
~8-19x the Sobolev line luminosity (He II 4686: 838 A vs 84 A EWDATA;
Luc's .fl ~100 A), so BLANK is not a consistent render for SOB
populations. BLANK remains the render for KNOP runs.

Requirements and constraints: `DO_SOB_FLUX_OBS=T` requires
`DO_SOB_LINES=T` (the EW block owns CHIL_MAT/ETAL_MAT; the reader
STOPs otherwise). `SOB_FREQ_IN_OBS=T` is needed so the obs grid is
refined across each resonance zone - without it the rendered profile
lands on the continuum-spaced grid. `prepare_CMF_PARAM.py`
stamps the combination from `cmfgen.obs_line_treatment='sob'` (the
SOB default); `'sob_ew'` gives T/T flags with injection off
(EWDATA diagnostic only).

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
### Probe results, 12-Oct-2026 build (commit 06480eb + gate fix 32b0d2a)

Four-cell moment grid resubmitted with fixes (a)-(c) (jobs 12163155,
12163156, 12163159, 12163169, interactive partition). CORRECTED
READING: an earlier draft of this section misread the `lum diag`
column (outer/inner ratio drift, ~1E-3) as evidence of a stable
field. The actual outer luminosity tells the true story -- all four
cells diverge:

- **SOB warm** (job 12163155): L_outer ramps ~x1.15/it for ~40 its
  from ~1E7 (3 dex below physical ~1E10), crosses the physical level
  WITHOUT pausing, then accelerates ~x2/it to 4E21 before T went NaN
  and the run SIGSEGV'd (exit 174) at it550. The ramp restarted
  identically after the requeue boundary (lagged J_NEXT/H_NEXT
  stores are process-local) -- same growth law regardless of
  restart. `max|EPS|/HU` did collapse 9.1E6 -> ~96 (the dln floor
  works as designed) and the refresh reports `conv=T`, but the
  converged map is expansive: `delta` (retained field vs the field
  the populations saw) stays ~1E2 and L doubles anyway. The refresh
  conv flag can come from the bulk exit (<0.1% cells moving >50%
  per pass) -- it does not certify the map's fixed point.
- **SOB cold** (it99+, running): identical signature --
  L_outer 8E6 -> 3E13 over it1-47, reset to ~2E7 at the requeue,
  ramping again. Same map, same gain, populations irrelevant.
- **KNOP warm** (job 12163167): L_outer 4E8 -> 1.2E32 over
  it441-454 (~x10-100/it, accelerating) with populations HELD the
  whole time -- frozen matter, so the runaway is purely the
  lagged-J map plus the coefficient-level feedback (Eddington
  factors rebuilt each iteration on the growing J through the
  formal solve). `refresh conv=F delta=1E6-1E7`, eps_cap ~40M/sweep
  (~10x per cell). Exited "COMPLETED" via the false converged exit
  below (old binary), n_frozen=14 < watchdog 60.
- **KNOP cold** (job 12163169 then 12165640): first run false-exited
  at it6; resubmitted on the gate-fixed binary -- now correctly
  reports demand (pops 1E7 not 0), continues iterating, still
  diverging (L 4E12 -> 7.8E17 by it6, refresh delta 1.2E15).

Failure locus in both pathways is the frequency-advection chain at
shock cells: refresh argmax sits at R~1.4-2.1E4, V~850-5000 km/s,
sigma_t ~ -5 to -64 (cmp_frac ~0.8 = 80% of cells read a lagged
frequency neighbour). The refresh's per-pass residuals floor at
~1E4 with same-direction ratio ~1.0 under KNOP -- the map is
non-contractive, and under SOB its "converged" iterate still grows.
Also present under KNOP: ESEC/CHI-1 up to 15.6 (negative-absorption
cells floored to pure scattering by MAX(0,1-COH)), and the formal
I_NEXT/fg store churns (drift 5-11, cap count overflows) because
the formal solve re-derives factors on the inflated J each
iteration.

False converged exit (now fixed, commit 32b0d2a): with
GS_POPS_HELD the populations are frozen while the field is
immature; BA_RES_ACC masked all 70 columns (|STEQ|>BA_TRUST),
zeroed STEQ, and DO_LEVEL_CHK read the all-zero post-mask state as
"all corrections <10%", overwriting MAXCH=0 -- discarding the
STEQ_MASKED_MAX boost added 04-Oct. `MAXCH < EPS` in CMFGEN_SUB
then set NUM_ITS_TO_DO=1 -> RVTJ + clean exit -> wrapper
`READY=0 + valid RVTJ` marked it COMPLETED. Fixes: `L_ANY_MASK`
disarms DO_LEVEL_CHK whenever any column was masked (armed-step and
Lambda mask paths), and the MAXCH<EPS gate now requires
`.NOT. GS_POPS_HELD`. Remaining gap: the wrapper's completion check
never inspects GS_STATE (frozen=T, mom_refresh_conv=F, pops_held=T
in a "completed" run is accepted); the n_frozen watchdog only fires
at 60 consecutive frozen iterations.

Conclusion: LAM_ES=F diverges on m15 under both pathways from both
starts. The divergence is not population feedback (KNOP warm held
populations throughout). It is the lagged frequency map's own
spectral radius >1 at compressing shock cells, plus -- under KNOP --
the J -> formal factors -> J coefficient-level coupling via the
ES emissivity and Eddington stores.
