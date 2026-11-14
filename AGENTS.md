# AGENTS.md - cmfgen_jfnk (tree rules)

Subordinate to the project-root `../AGENTS.md`, which controls the whole
project: mission, the five pillars, clamp policy, VADAT rules, the
testing ladder, and the working posture (root section 2). This file adds
rules specific to this tree. Where both apply, the stricter rule wins;
tree rules may tighten project policy but never weaken it. Standing user
directives recorded here (e.g. the push-to-origin rule in section 8) take
precedence over project defaults within this tree only.

You are an expert in radiative transfer and numerical methods, with
specific expertise in non-monotonic velocity fields, Sobolev
approximations, and Jacobian-Free Newton-Krylov (JFNK) solvers applied to
stellar atmosphere modeling (CMFGEN).

Outcome for this tree (see `../AGENTS.md` section 2): the roadmap below,
carried through to verified model runs submitted via
`scripts/run_model_grid.py`. A flag that compiles, a clean diff, or a
passing unit test is progress toward that outcome, not the outcome
itself. The outside judges are the verification checklist in section 7.

## 1. Purpose

This tree is the **clean restart** of the CMFGEN non-monotonic velocity
development pipeline. It was extracted directly from
the pristine upstream archive `cur_cmf_18jun25.tar.gz`
(md5 `1d5a7a168e77329dd97c3019bcad3c09`) and carries a **minimal diff**.

## 1.1 Development roadmap (user directive 18-Sep-2026)

1. `DO_NONMONO_SOB` on a non-monotonic model; produce a spectrum.
2. Non-mono SOB on an interacting model (type IIn-like spectrum).
3. Knop arbitrary-velocity solver as a plug-in; touch base CMFGEN as
   little as possible.

Run area: `/ptmp/sbrennan/SN2024pba/cmfgen_validation` (config
`configs/dessart_m15_2p82d_****.yaml`).

## 1.2 Pipeline defaults (implemented in scripts/, not this tree)

- `ADJUST_CORRECTIONS`, `IT_SPECIFIER`, `table*.dat` excluded from
  run directories by default (opt-in only).
- Grid adjustments off: `REV_RGRID=F`, `N_OB_INS=0`; the initial
  `SN_HYDRO_DATA` arrives pre-optimized.
- `IT_ON_T=F` - no initial T iteration; hydro T used as-is.
- Pinned T is the default T policy (user directive 13-Nov-2026):
  `cmfgen.fix_temperature` defaults to true in `prepare_VADAT.py`, so
  `FIX_T=T` + `FIX_T_AUTO/DO_T_AUTO/DO_GT_AUTO=F` and T stays on the
  `SN_HYDRO_DATA` profile at all depths. `fix_temperature: false`
  selects the previous thin-region release; `solve_temperature: true`
  is fully free T.
- `THK_CONT=F` (10-Oct-2026) - thin outer continuum boundary; the
  m15 wind models are optically thin at the outer edge. Thick stays
  available via `vadat_overrides: THK_CONT: T`.
- Interactive (debug) partition jobs use the QOS max of 16 CPUs.
- Test/debug runs (debug exe, FP-trap or A/B diagnostics) go to the
  `interactive` partition via sbatch (2 h wall, 16-CPU QOS cap,
  30G mem cap), never on the login node (user directive 03-Nov-2026).
- Production CMFGEN jobs default to 4 CPUs (`CMFGEN_PROD_CPUS` in
  scripts/run_model_grid.py, 07-Nov-2026; was 16 then 8, both queued on
  QOSGrpCpuLimit).  `cmfgen.cpus_per_task` overrides per model.  The
  generated slurm keeps the small request for schedulability but bumps
  OMP_NUM_THREADS to the scheduler-allocated count (`squeue %C`) when
  viper grants more than requested; `SLURM_CPUS_ON_NODE` stays banned
  (whole node on general).
- KNOP test/diagnostic runs request a minimum of 8 CPUs
  (`cmfgen.cpus_per_task >= 8` in configs/*knop*/*KNOP*.yaml,
  user directive 14-Nov-2026). Applies to future KNOP runs; the
  production 4-CPU default above is unchanged for non-KNOP work.
- `cmfgen.do_gt_auto` YAML key controls `DO_GT_AUTO` in IN_ITS
  (grey-T update); set `false` when it hurts convergence.

## 1.3 The cmfgen_obsidian/ vault (living docs - keep updated)

`cmfgen_obsidian/` is an Obsidian vault: markdown notes documenting CMFGEN,
the non-monotonic problem, the Knop/GS/JFNK solver stack, control
gates, VADAT keys, diagnostics, and run history. Index:
`cmfgen_obsidian/Home.md`.

Mandatory: any change to solver behaviour, gating logic, VADAT keys,
diagnostic output, or run procedure must update the affected notes in
the same edit - a stale note is the same class of defect as a stale
comment.  Status entries logged in `AGENTS_HISTORY.md` or
`cmfgen_obsidian/history/convergence_issues.md` are mirrored into the
relevant vault note as architecture-level facts.

Vault conventions: plain ASCII; `[[wikilinks]]` between notes;
`#status/verified|partial|unverified` tags.

## 2. Prime directive: stability first

Every change to this tree must satisfy all of the following before it
is considered done:

1. **Default path unchanged.** With `DO_NONMONO_SOB` absent or `F` in
   VADAT, the code path must be identical to pristine upstream. New
   behaviour is gated on explicit flags, never on by default.
2. **Minimal diff.** Keep the difference vs the pristine reference
   (`../cmfgen_reference/cur_cmf`) as small as possible. No wholesale
   ports, no refactoring of working upstream code, no reformatting.
3. **Surgical edits.** Make small, effective changes that do not
   drastically alter the overall structure of the CMFGEN code unless
   warranted. Prefer: a flag-gated branch in an existing routine over
   a forked routine; a new standalone subroutine over surgery inside
   a large upstream routine; a few lines of interpolation/guard code
   over reworking an algorithm. If a change cannot be made small,
   that is evidence the approach is wrong - stop and reconsider
   before restructuring. Every edit should be reviewable as a short
   diff a reader can verify against the reference in minutes.
4. **One change at a time.** Do not batch independent modifications.
   Rebuild and verify after each change.
5. **No silent physics changes.** No velocity monotonisation, no new
   floors or clamps outside the documented Dessart tau floor, no
   population/opacity substitutes. Any numerical guard must be local,
   documented, and physically justified (see the clamp policy in
   `../AGENTS.md`, section 6).
6. **Build must stay clean.** `./build.sh` (the full production target)
   must complete with zero errors after every edit. Treat new warnings
   on edited files as failures until understood.
7. **Honest status.** Compiling is not validated. A feature is only
   "working" once a model run through `run_model_grid.py` demonstrates
   it. Report PASS / PASS WITH LIMITATIONS / FAIL / BLOCKED.

## 2.1 Comment conventions (mandatory)

Comments are part of the reviewable diff. Every comment added or
touched in this tree must satisfy all of the following.

### Content

1. Explain the "why", not the "what". Do not restate the code. If a
   line needs a comment to be readable, prefer rewriting the code.
2. Be detailed but concise. A good dev comment answers: what physical
   or numerical assumption justifies this, what breaks without it,
   and what upstream behaviour it replaces.
3. Cite the source of the method where one exists:
   `Knop et al. (2009, Paper I Eq 7)`, `Dessart, Audit & Hillier
   (2015, arXiv:1503.05463)`. Use the paper's own equation/section
   numbers.
4. State the gating flag explicitly: which VADAT key controls the
   branch (`[DO_NONMONO_KNOP]`), its default, and that the default
   path is unchanged.
5. Name the exact routines/keys that exist in THIS tree. Do not cite
   V6-only machinery (`SOLVE_FOR_POPS_NONMONO_V1`, `NON_MONO_SOLVER`,
   `SOLVE_FOR_POPS_JFNK_V1`, `SIGMA_TRUE_EXT`, `EXTEND_SIGMA_TRUE`,
   `VAR_MOM_J_CMF_NONMONO_V1`) as if present. If a comment must
   reference them, mark it: `(V6, not ported to this tree)`.
6. Label temporary diagnostics as temporary, with what they are meant
   to catch (e.g. the `DBG-T` blocks in `cmfgen_sub.f`).
7. When code is removed rather than replaced, leave a short comment
   saying what was removed and why (see the `S_tilde` removal block in
   `fg_j_cmf_nonmono_v1.f`).
8. A bug fix that is NOT flag-gated must say so and justify why the
   stock path changes (see the `'ZERO FLUX'` fix comments in
   `comp_j_blank.f`).

### Form

- Plain ASCII only. No em-dashes or en-dashes; use `-` or `--`. No
  smart quotes. No Greek letters or math symbols; write `sigma`,
  `mu`, `chi_hat`, `nu`, `GAM`. No arrows other than `->`.
- Simple direct English. Banned words and phrases: delve, tapestry,
  testament, leverage, utilize, seamless, robust, crucial, paramount,
  landscape, navigate, intricate, "it is important to note".
- Match the file's existing comment style: `!` in the free-form
  `.f`/`.INC` files; `C` in column 1 only where the file already uses
  fixed-form `C` comments (e.g. `subs/sob_rh_couple.f`, `subs/gff.f`).
- Keep comment lines within 132 columns (files compile with
  `-extend-source`); match the local indentation of the code the
  comment describes.
- Follow the upstream layout conventions: `! Altered : DD-Mon-YYYY :`
  history lines at the top of a routine, `! edited by <name>` on new
  blocks, `! Dimensioned ...` notes on module arrays, `!` on its own
  line to separate comment blocks.

### Discipline

- Update the comment in the same edit as the code it describes; a
  stale comment is a defect, not a shortcut.
- Comment-only edits do not require a rebuild, but they must still be
  reviewed for factual accuracy (routine names, flag names, equation
  numbers) before commit.
- Do not add comments to untouched upstream lines. The diff vs
  `../cmfgen_reference/cur_cmf` must stay minimal in comments too.

## 3. Current feature: DO_NONMONO_SOB

Implements the Dessart, Audit & Hillier (2015, arXiv:1503.05463, Sec 3.1)
Sobolev treatment for non-monotonic velocity fields:

- `SIGMA_TRUE` (`MOD_CMFGEN`): unclamped `dlnV/dlnR - 1` from centred
  differences on R,V. `SIGMA` keeps its infrastructure clamp, so the
  moment equations and continuum formal solver are untouched.
- `RECOMPUTE_SIGMA_TRUE` (`new_main/subs/recompute_sigma_true.f`):
  fills `SIGMA_TRUE`; called at model setup, after `DO_WIND_VEL_V1`,
  after NG acceleration, and after `ADJUST_R_GRID_V4` grid revision.
  Under `DO_NONMONO_SOB` it also re-clamps `SIGMA` to `>= -0.999`:
  the model readers leave `SIGMA` unclamped, which let `sigma < -1`
  reach stock routines that assume `sigma >= -1` (`MOM_J_CMF_V11`'s
  `W`/`EPS` terms).  Not applied under `DO_NONMONO_KNOP` (26-Sep-2026,
  commit 5c85085).
- Tau floor in `sobjbar_sim.f`, `sobew.f`, `sobew_grad_v2.f`,
  `sobew_grad_v3.f`: for sigma < 0, `|1+sigma*mu^2|` is replaced by
  `max(|1+sigma*mu^2|, -0.2*(1+sigma), 1e-4)`. For sigma >= -1 the
  expression is bit-identical to upstream.
- `SUB_SOB_LINE_V3` and the dielectronic `SOBEW` call in `cmfgen_sub.f`
  pass `SIGMA_TRUE` when the flag is set.
- `SUB_SOB_LINE_V3` carries transient-iterate guards on the BA/Jacobian
  terms (26-Sep-2026, commit b107b4b): `VB_SIM`/`VC_SIM` and `SRAT` are
  set to zero when `ETAL_MAT`, `BB_COR` or `CHIL` go non-positive, and
  the self-line derivatives avoid dividing by population-proportional
  denominators (`POPS(NUP)`, `U_STAR_RATIO`) that can vanish on a bad
  Lambda iterate.  Identical results for valid iterates.
- VADAT: `T [DO_NONMONO_SOB]` (optional, default F). Requires
  `GLOBAL_LINE='SOB'` (or `NONE` with every present ion `TRANS_*=SOB`).
  `USE_LAM_ES` is un-forced over explicit values: an absent key
  auto-sets T since 10-Nov-2026 (the `FG_J_CMF_V13` formal J on
  clamped SIGMA is the solver the SOB populations are validated
  against); an explicit F selects the opt-in moment path (next
  entry).
  `LAM_SET` is normalized to F by the SOB auto-set (hard STOP
  backstop) since 22-Oct-2026 - was warn-only before.  Solver modes
  that bypass the COMP_J_BLANK non-mono branch (`USE_FIXED_J`,
  `USE_FORMAL_REL`, `USE_J_REL`, `USE_DJDT_RTE`, `PP_MOD`, `PP_NOV`)
  are rejected under `DO_NONMONO_SOB` since 22-Oct-2026.
- `USE_LAM_ES=T` under `DO_NONMONO_SOB` (DEFAULT since 10-Nov-2026,
  superseding the 21-Oct-2026 F default; the m15 recipe, job
  12001805 / VAL_sob_jfnk): continuum J comes from the
  `FG_J_CMF_V13` formal solve on clamped `SIGMA` -- the clamped
  Sobolev approximation of Dessart, Audit & Hillier (2015,
  arXiv:1503.05463, Sec. 3.1). `SIGMA_TRUE` still drives the Sobolev
  line optical depths and their angular-factor floor. This is the
  intended SOB product: less accurate than `DO_NONMONO_KNOP` but a
  converged initial guess for it (`warm_restart_from`).
- `USE_LAM_ES=F` under `DO_NONMONO_SOB` (opt-in experiment since
  10-Nov-2026; was DEFAULT 21-Oct -> 10-Nov-2026, rewired
  12-Oct-2026):
  continuum J comes from `MOM_J_CMF_NONMONO_V1` (the Knop sign-split
  moment solver native to non-monotonic fields), paired with fixed-J +
  `GET_MOM_ALO_DIAG` + `GET_MOM_CELL_QUAR` in `do_var_cont.f` - the
  same recipe as KNOP+LAM_ES=F. `FG_J_CMF_V13` still runs each
  frequency, so `N_ON_J` factors come live from its ray moments via
  `GET_MOMS_NON_REL`; `SIGMA_TRUE` drives the Sobolev lines and the
  solver GAM term. JFNK stays inactive in
  this mode (admission is gated on `USE_LAM_ES`). `DO_KNOP_REFRESH`
  defaults on here as well as under KNOP (`MOM_KNOP_REFRESH_V1`
  converges the lagged J_NEXT/H_NEXT map within each sweep).
  Known caveat: the V13 ray moments feeding `K_ON_J` are built on
  clamped `SIGMA` while the solver's GAM term uses `SIGMA_TRUE`
  (real but second-order per the 13-Oct audit); the dominant
  advection-coefficient inflation comes from bound-free edge
  doublets on the continuum frequency grid (~2E-11 dLOG_NU).
  Runtime status on the current build is FAIL - bounded
  non-contraction (job 12316397, `m15_2p82d_SOB_l14_lte_np_mom`,
  cancelled at it742, 22-Oct-2026): the reworked containment holds
  (no luminosity runaway) but the lagged map never contracts -
  sweep drift ~1.0 every sweep, ~200K cells quarantined/sweep
  (~183K persistently frozen), refresh stalls at residual ~1.1,
  model convergence pinned at 63-75%.  The `USE_LAM_ES=T` partner
  (`m15_2p82d_SOB_l14_lte_np`) converged 100% at it718, so the
  formal-solve J update remains the working recipe on m15; the
  moment path is not viable there yet.  See
  `cmfgen_obsidian/history/convergence_issues.md` 22-Oct-2026 (d).
  The earlier V11+V12 variant is not restorable by design - V11
  assumes expanding flow and V12 linearizes the V11 system.
  The GS freeze/maturity/watchdog stack stays KNOP-scoped
  (`GS_FREEZE` is forced F outside `DO_NONMONO_KNOP`), so the moment
  path under SOB runs bare - `MAX_CHNG` and `DO_BA_RES_ACC` are the
  only rollbacks.
- SOB auto-sets (absent-key only, 21-Oct-2026; explicit VADAT values
  always win): `DO_NONMONO_FREQ_MATCH=T` under `USE_LAM_ES=F`,
  `DO_FG13_ALO=T` under `USE_LAM_ES=T`, `DO_BA_RES_ACC=T`,
  `DO_BA_RES_REUSE=T` (when `DO_BA_RES_ACC` ends up T), `NG_ABSTAIN=T`,
  `DO_ES_ALI=T` (11-Nov-2026, when `COHERENT_ES` and not
  `DO_SOB_FG_KNOP` - both `USE_LAM_ES` routes carry a stored
  diagonal).
  `LAM_SET=F` is an override auto-set (22-Oct-2026) - applied even to
  explicit T, like KNOP, because the key is structurally
  incompatible, not a knob.
  Overrides of explicit user values (e.g. `GLOBAL_LINE`) print a
  banner warning to OUTGEN; auto-sets print ` SOB auto-set:` lines.
- Solver preference (user directive 11-Nov-2026): new production and
  experiment configs use `DO_NONMONO_SOB` + `USE_LAM_ES=T` (the V13
  formal-J route). `DO_NONMONO_KNOP` is retained for reference and
  diagnostic comparisons only - do not select it for new configs.
  `DO_ES_ALI_BA` stays opt-in: the X3 warm arm matched X2
  iteration-for-iteration (`docs/audits/SOB_ACCELERATION_RESULTS.md`),
  so it is not part of the default accelerator set.
- `BA_DEM_RELFRAC` (VADAT key, removed 17-Oct-2026, commit b398571a):
  was a relative demand floor in `SOLVEBA_V13` - rows below
  `RELFRAC*max(POPS)` at a depth excluded from MAXCH/MAXCH_APP.  It hid
  demand at populated depths; superseded by `BA_DEM_PFLOOR` and the
  exit-only `BA_EX_RELFRAC` floor.  Do not re-add.
- `DO_FG13_ALO` (VADAT key; absent-key auto-set T under
  `DO_NONMONO_SOB` + `USE_LAM_ES=T` since 21-Oct-2026 - under the
  `LAM_ES=F` moment path the matching diagonal comes from
  `GET_MOM_ALO_DIAG` instead): `FG_J_CMF_V13`
  stores the same diagonal OAB approximate-Lambda operator the Knop
  solver uses (`GET_FG13_ALO_DIAG`); `do_var_cont.f` fills the
  continuum dJ/dCHI, dJ/dETA diagonals under `DO_NONMONO_SOB`. J
  stays the `USE_LAM_ES=T` formal solution.
- `DO_ES_ALI` (VADAT key, default F; absent-key auto-set T under
  `DO_NONMONO_SOB` since 11-Nov-2026 when `COHERENT_ES` and not
  `DO_SOB_FG_KNOP`; `ES_ALI_TOL` 1e-4,
  `ES_ALI_MAX` 20-clamped): under `DO_NONMONO_SOB` + `USE_LAM_ES=T`
  + `COHERENT_ES` on the V13 path (not `DO_SOB_FG_KNOP`), the
  `RJ=TC` Lambda update in `comp_j_blank.f` becomes an
  Olson-Auer-Buchler Jacobi step on the stored V13 diagonal
  (`GET_FG13_ALO_DIAG`): `RJNEW = RJ + (TC-RJ)/max(1 -
  DJDETA*chi_scat, 0.01)`, and the J loop keeps passing until the
  last relative update is below `ES_ALI_TOL` or `ES_ALI_MAX`.
  Targets the ES-scattering-dominated m15 convergence mode (~2
  Lambda steps/great iteration under the Eddington-only exit).
  Added 27-Oct-2026; verified in `adp2_esali` (J loop converges to
  1e-4 in ~2 passes, 0 unconverged freqs).
  Warm-arm evidence 11-Nov-2026 (`docs/audits/SOB_ACCELERATION_RESULTS.md`):
  X2 (`DO_ES_ALI=T` alone, warm restart from the converged FRESH_sob
  checkpoint) converged in ~30 new iterations where the unaccelerated
  X1 twin was still pinned at maxapp ~0.5 after ~80; X3 adding
  `DO_ES_ALI_BA` matched X2, so the BA factor stays opt-in.
  The `USE_LAM_ES=T` + `DO_NONMONO_KNOP` branch was inert until
  02-Nov-2026: the update always read the V13 diagonal via
  `GET_FG13_ALO_DIAG`, but the KNOP formal solver is
  `FG_J_CMF_NONMONO_V1`, which stores its diagonal in
  `FG_J_CMF_NONMONO_MOD_V1`, so every frequency fell back to plain
  Lambda (`KNOP_fresh`: 56377/56377 `ES_ALI` fallbacks each
  iteration).  The accessor is now selected on the same predicate
  as the solver dispatch (`GET_FG_ALO_DIAG` under KNOP /
  `DO_SOB_FG_KNOP`); `KN_NO_GAIN` is forced T under KNOP since
  08-Nov-2026, so the DEN floor (0.01) albedo <= 1 assumption holds
  there automatically.  Runtime status
  on the KNOP branch: pending the `KNOP_esB` A/B run.
  Extended 29-Oct-2026 (`ES_ALI_MOM_ON`): same OAB Jacobi update on
  the `USE_LAM_ES=F` branch under `DO_NONMONO_SOB` (the KNOP route
  is closed -- `USE_LAM_ES` is forced T there since 08-Nov-2026),
  where `MOM_J_CMF_NONMONO_V1` solves J with the
  scattering source lagged (TA = ETA + chi_scat*RJ_prev).  Uses
  `GET_MOM_ALO_DIAG`, the mom solver's stored bare local response.
  Gated on `DO_ES_ALI` + `COHERENT_ES` + `.NOT.USE_LAM_ES` +
  `.NOT.DOB_FG_KNOP`; shares ES_DJ/ES_ALI_TOL/ES_ALI_MAX and the
  fallback/print accounting.  Runtime status unverified pending
  `KNOP_warmS`.
- `DO_ES_ALI_BA` (VADAT key, default F; `ES_ALI_BA_MIN` 0.1):
  same Jacobi factor as DO_ES_ALI applied to the stored ALO when
  `do_var_cont.f` fills TX -- the BA linearization then sees the
  ES-converged J response `dJ ~ ALO/(1-DJDETA*chi_scat)` instead of
  the bare-Lambda one.  Without it the BA under-corrects a
  scattering-coupled mode every iteration (fixed-ratio ~0.97/it
  decay on the m15 2.82d model).  Applied per depth-row to every
  filled band; composes with DO_FG13_ALO_TRI.  Added 28-Oct-2026;
  requires `DO_ES_ALI=T` since 11-Nov-2026
  (CHECK_PARAM_CONSISTENCY hard-stops the pairing; the X3 warm arm
  showed no gain over X2, so it stays opt-in).  The J-loop DEN floor
  (0.01) and `ES_ALI_BA_MIN` (0.1) differ by design; activations are
  counted in `ES_ALI_NFLOOR`/`ES_ALI_BA_NFLOOR`.
- `DO_FG13_ALO_TRI` (VADAT key, default F): `FG_J_CMF_V13` stores
  the off-diagonal Lambda* entries `dJ(K)/dS'(K+/-1)` from the
  direct A1/A2 (+ propagated A0*A2) ray coefficients
  (`GET_FG13_ALO_TRI`); `do_var_cont.f` fills `TX(K,K+/-1,1:2)` so
  the banded BA matrix couples adjacent depths (Olson & Kunasz
  1987).  Effective only with `SOL_METH=TRI`/`NUM_BNDS=3`
  (pipeline: `depth_independent_solver: false`); with NUM_BNDS<3
  the fill is skipped.  Ray-node neighbors map to depth neighbors -
  exact on the non-inserted grid, approximate on inserted/tangent
  nodes; source-gradient (dS) couplings ignored, matching the
  diagonal.  Added 28-Oct-2026; runtime status unverified pending
  `adp2_alo3`.
- `DO_SOB_DEP_WINGS` (OBS CMF_FLUX_PARAM key, default F;
  `SOB_DEP_WING_EXT` 2000 km/s): in `ADD_SOB_LINE_PROF`, continue
  the deposited Gaussian core as 1/x^2 damping wings out to
  SOB_DEP_WING_EXT instead of truncating at 3 Doppler widths.  The
  full CMF solve carries Voigt/damping wings, so the pure-Gaussian
  deposit renders optical line wings systematically too narrow
  (m15 2.82d SOB-render comparison vs the Dessart .fl: HeI 5876
  FWHM 18A vs blank-render ~30A; with wings on, 32A).  Opt-in via
  `cmfgen.cmf_flux_param_overrides`.  Wing shape uses the true
  per-line Voigt (`VOIGT` with A_VOIGT=(ARAD+1.55E4*|C4|^0.667*Ne)/
  (4*pi*dnu_D), same as SET_PROF_V6) since 28-Oct-2026.
- `DO_SOB_DEP_FLOOR` (OBS CMF_FLUX_PARAM key, default F;
  `SOB_DEP_CHI_FLOOR` 1.0): masing/negative-opacity cells in
  ADD_SOB_LINE_PROF keep their emissivity with opacity floored
  (optically-thin emitters) instead of being dropped.  The dropped
  cells sit in the dense inner shell and carry the P-Cygni wing
  flux; with the floor on, injected SOB emissivity is 3.35E17 vs
  ~2/3 of that dropped before (m15 2.82d esali OBS).  Added
  28-Oct-2026.
- `LUM_DIV_FAC3`/`LUM_DIV_N3` (VADAT keys, KNOP only; default
  N3=0 = off): third, drift tier of the `LUM_DIV_*` runaway guard
  in SOLVE_FOR_POPS.  Trips `LUM_BLOCKED` on N3 consecutive calls
  each with `L_out` growth above factor FAC3 (default 1.005).
  The FAC=1E2/N=3 and FAC2=5/N2 tiers both require per-call
  growth >= 5x, so a sustained ~2%/call drift is invisible to
  them -- the m15_2p82d_KNOP_l14_freet restart climbed 1.6E9 ->
  1.2E10 L_sun over ~300 iterations under `hold_pops` (a cold
  `I_NEXT` rebuild after a checkpoint-less restart; `FGINEXT_GS`
  absent because the parent run predated the dump code) and only
  stopped on non-finite `L_out`.  Opt-in: a legitimate field
  re-equilibration after a cold rebuild produces the same slow
  climb signature, so enabling is a per-run policy.  Set
  `LUM_DIV_CAP` to a model-physical level (m15 ~3E9) alongside.
  Added 29-Oct-2026.
- `GS_FOREIGN_SEED` (VADAT key, KNOP only, default F): on restart,
  accept a fingerprint-valid but tag-mismatched `FGINEXT_GS` (e.g.
  from a converged SOB run of the same model) as the lagged-field
  seed instead of a cold `I_NEXT=0` rebuild.  Seed stays immature
  (`I_NEXT_GS_OK=F`) so `N_LAM_RWARM` re-converges the map first.
- `KN_LUM_TRACE` (VADAT key, default F): first-detection scan of
  `STEQ_T` after each energy-equation accumulation site -- reports
  site, freq, depth, operands.  Diagnostic only.
- Empty-level net-bracket guards (29-Oct-2026, not flag-gated):
  `ZNET_SIM`, `ZNET`, `BETAC` all divided by an emissivity that is
  exactly 0 where the upper level is empty (thin outer depths on
  checkpoint restarts).  `RJ*CHIL/0=-Inf` then `ETAL*ZNET=0*Inf=NaN`
  poisoned `STEQ_T` -> `DEP_RAD_EQ` -> `RLUMST(ND)`, tripping the
  SOLVE_FOR_POPS non-finite stop and killing every warmS/l14_freet
  requeue at iteration 864 (trace: SIM=16, K=2, ML=1204).  Now the
  `ETAL=0` resonance contributes no bracket term (ZNET_SIM/ZNET
  consumers all multiply by ~0 anyway); `BETAC=0` likewise.
  Verified: warmS it864/865 complete (~315-470 s, STEQ_T finite,
  max 8.5E12 at K=70).  warmS run status: STABLE but not yet
  converged -- outer/inner L ratio ~10 (4.2E9 vs 4.3E8), KNOP
  refresh residual ~4E-4 frac unconverged.
- `adp2_esali` (m15 2.82d SOB+LAM_ES+DO_ES_ALI): clean convergence
  at it851, Model Convergence 100% (MAXCH<EPS=1e-3), ES_ALI J loop
  converging all freqs to dJ/J~1E-4 in ~2 passes.
- Lagged-store fixes (02-Nov-2026, from esB job 12586072 review):
  (a) `FG_SEED_COMMIT` (VADAT, default T) - per-ray seed-fill
  waiver: a ray that has never committed via WRITE_STORE skips the
  FG_COMMIT_JUMP veto once, so a cold/foreign store can take its
  converged first fill instead of being rejected ray-by-ray (esB:
  67/85 rays jump-rejected at fill gaps 1E2-6E2). Validated
  FGINEXT_GS restores count all rays seeded. (b) `SNP_P`/`SNP_M`
  per-frequency snapshots: the in-sweep store write now relaxes
  against the frozen pre-frequency column on EVERY solver call, so
  under DO_ES_ALI (~12 comp_j_blank re-entries/freq) the store
  lands on the last (ES-converged) pass's field rather than the
  first; churn metrics keep their NEW_FREQ one-sample census.
  Both changes touch `fg_j_cmf_nonmono_v1.f` only; default path
  unchanged (KNOP-scoped).

## 4. Recent corrections and runtime status notes

### 4.1 KNOP precondition-phase corrections (21-Oct-2026, m15_2p82d_KNOP_np cold-start stall, job 12298049)

The cold KNOP start pinned the `GS_PRECOND` gate forever: a handful
of ill-conditioned rows demanded +1E15..+1E72% every LAMBDA call
(documented m15 stall class - degenerate equations whose SE target
moves faster than the applied clamp can chase it).  Three surgical
fixes:

- `fiddle_pop_corrections_v2.f` MAJOR scale: under DO_NONMONO_SOB or
  DO_NONMONO_KNOP each row's contribution to the per-depth SCALE is
  capped at the bound the per-element apply clamp already enforces
  (DPTH_BIG_LIM/DPTH_LIT_LIM).  One unsatisfiable row used to veto
  the whole depth's corrections (SCALE ~1E-69); the bad row still
  takes its identical clamped step, the rest of the depth moves.
  Stock path unchanged.
- `solveba_v13.f` LAMBDA-phase quarantine: a row demanding beyond
  `BA_LQ_DEM`% (default 1E12) for `BA_LQ_PAT` (default 8)
  consecutive LAMBDA calls parks into `BA_QUAR_MASK`; its correction
  is then suppressed on every call, so it leaves MAXCH and stops
  churning.  Armed only with `DO_BA_RES_ACC`+`DO_BA_QUAR`; bounded
  by `BA_QUAR_MAX`; each parking prints its EQ_DESC and the
  suppressed max demand prints per call.  Floor-level rows
  (BA_DEM_PFLOOR) and held calls are not counted.  This reverses
  the deliberate "parked rows still update on LAMBDA calls" choice
  for LAMBDA-parked rows only: under the precond gate every call is
  LAMBDA, and these rows demonstrably do not contract.
- `solve_for_pops.f`: `N_PRECOND_TOT` now counts only calls that
  fail to improve on the stretch's best demand by 10% (`GS_PC_BEST`
  tracker), matching the documented "no-progress" semantics - a
  still-contracting cold start was previously stopped at
  `GS_PRECOND_LIM` regardless of progress.

### 4.2 SOB nonloc it64 NaN localisation and solver hardening (05-Nov-2026, m15_2p82d_VAL_sob_nonloc job 12689657)

The `sob_nonloc` model segfaulted at it65 entry: `GFF_VEC` indexed
the 12x11 Gaunt table with `INT(NaN)` after `SCRTEMP` carried an
all-NaN `POPS(:,11)` column.  `STEQ_VALS` showed it64's final
solve returned `STEQ(:,11)` entirely NaN while every residual dump
before it was finite - the defect is inside the BA matrix assembly
or the `CMF_BLKBAND_V3` block solve, not in the `DRC` injection
into the residual (that path was audited and cleared).  Changes:

- `solveba_v13.f`: after `CMF_BLKBAND_V3` returns, `STEQ` is
  scanned for non-finite entries; the first bad (level, depth,
  value, count) is printed and the run stops.  Previously NaN
  corrections sailed through every `.GT.`/`.LT.` scale guard into
  `POPS` and `SCRTEMP`.  Unconditional - a non-finite solve output
  is always a defect.
- `sob_rh_couple.f`: the `DRATE` finiteness stop now reports DJ,
  CHIL, BB_COR and partner count; the dominance census counts
  cells where remote `dJ` exceeds the local mean intensity
  (`dJ>J_local`, the `RH_MAXDJ` normalisation) instead of the
  degenerate `|DRC|>|ETAL*ZNET|` test that was trivially true at
  detailed balance.  New VADAT key `SOB_RH_TEST_S1` (default F)
  accumulates the partner sum with unit source function; remote
  `dJ(S=1)` above ~1 means the two-point kernel overcounts.
- `mom_j_cmf_nonmono_v1.f`: `MOM_RF_DIV_BULK` tolerated excursions
  are now bounded - max streak 3 consecutive and excursion
  magnitude capped at 1E6 x MIN_RES (the 1E3 abort bound stands);
  committed-field finite/>1E40 and 50x-jump guards already
  existed.  `MOM_XMWORST` reports the largest negative RHS value
  clipped by `MOM_XM_POS` each sweep so the safeguard's cost is
  measured, not assumed.

### 4.3 KNOP stability settings forced over explicit VADAT (08-Nov-2026)

Under `DO_NONMONO_KNOP=T`, `RD_CONTROL_VARIABLES` Rule 1 now forces
three keys to T even when VADAT sets them F (explicit F prints a
` KNOP override:` banner; absent keys print ` KNOP auto-set:`):

- `KN_NO_GAIN=T`: with F the stock `NEG_OPAC` floors can leave
  `CHI < CHI_SCAT`, so the source coefficient sum exceeds 1 and the
  cell amplifies photons; the lagged `I_P/I_M_DEPTH_NEXT` store
  compounds that gain each great iteration (~8.6x/iter measured;
  `KNOP_fresh` 02-Nov-2026 diverged to `L_out` ~1e29). The floor
  makes the lagged map non-expanding in sup norm.
- `USE_LAM_ES=T`: the `MOM_J_CMF_NONMONO_V1` moment path never
  converged under KNOP -- its lagged `J_NEXT` map is expansive
  there (m15_2p82d job 12747893, 06-Nov-2026: adaptive
  `JNEXT_RELAX` contraction factor 1.24, every refresh
  non-converged, NaN in SOLVEBA at it1478; TypeIIn same verdict).
  The MOM code stays in the tree and stays reachable under
  `DO_NONMONO_SOB`; only the KNOP route is closed.
- `DO_NONMONO_FREQ_MATCH=T`: with F every compressing upwind read
  uses the adjacent bin; on m15 the per-cell Doppler drift reaches
  ~330 bins, so the read is at the wrong frequency in the
  oscillatory-eigenvalue regime (Ceulemans et al. 2025, JQSRT,
  arXiv:2505.08309; Knop et al. 2009 II, A&A 501, 813).
  Correctness requirement backed by theory -- no failed run is yet
  attributed to F.

The removed `DO_ES_ALI`+KNOP-without-`KN_NO_GAIN` read-time warning
is dead code under the force. Under `DO_NONMONO_SOB` the
`USE_LAM_ES` default flipped back to T on 10-Nov-2026 (absent-key
auto-set; explicit F still honoured as the opt-in moment
experiment), and `KN_NO_GAIN`/`FREQ_MATCH` keep their read values.

## 5. Deliberately NOT ported from V6

Do not re-add these without an explicit, tested reason:

- `solveba_nonmono_v1` / `generate_full_matrix_nonmono_v1` BA damping
- `solve_for_pops_nonmono_v1` population-solver fork
- `SIGMA_TRUE_EXT` / `EXTEND_SIGMA_TRUE` (extended-grid variant used
  only by the Knop formal solver)
- Knop `chi_hat` formal solver (`DO_NOMONO_KNOP` path)
- JFNK machinery and all V6 debugging edits

The Dessart SOB path uses the standard CMF formal solver; none of the
above are required by the paper. If convergence at compression depths
proves unstable in a real model, revisit `solveba_nonmono` first -
and only after documenting the instability.

## 6. Build

```bash
module load intel/2024.0 mkl
cd /ptmp/sbrennan/SN2024pba/cmfgen_jfnk
./build.sh           # make KEEP_EXE=1 clean && make DEBUG=nodebug all;
                     # non-zero exit. Old production exes remain installed
                     # during compile and are replaced by temp-link + mv.
                     # -> exe/cmfgen_dev.exe, exe/cmf_flux.exe,
                     #    exe/obs_frame.exe and the utility tree
./build_debug.sh     # make dbg_clean && make DEBUG=debug cmfgen
                     # -> exe_dbg/cmfgen_dev.exe (-O0 -check bounds)
```

Build config: `viper_intel_nodebug`, ifort `-O3` + AVX-512 flags,
`-align array64byte`, MKL. PGPLOT is stubbed (`pgplot_stubs/`).

The OBS executables are part of `./build.sh` (`make all`). To rebuild
only `cmf_flux.exe` after OBS source changes, run `make -C dust`,
`make -C disp mod_disp.o` (parent Makefile first - the stale
`disp/subs` rule would invoke unavailable `f77`), then `make -C obs`.

Build-order constraints (each enforced in the Makefiles, all reached on
a clean tree since 29-Sep-2026):

- `d_nm_subs` depends on `d_mod_subs` (new_main/subs needs the
  CONTROL_VARIABLE_MOD module file; fc12e0f).
- `d_lte` is in the `cmfgen` target but builds `OBJ_FILES` only
  (12f1f1a, e3b2b84): `do_cmf_hydro_v2.o` in libsubs.a needs
  `lte_hydro/` objects or the final link fails on a clean tree.  The
  utility exes are deliberately excluded: `main_lte.exe` cannot link
  under -O0 (`lte_sub.f` calls `WR_ABUND_INFO_V2`, absent here; the
  -O3 link dead-code-eliminates the call).  To build them anyway:
  `(cd lte_hydro && make)`.
- `plane/mod_ray_mom_store.f` is compiled as a prerequisite of
  `comp_j_blank.o` inside `new_main/mod_subs` (12f1f1a): the
  DO_NONMONO_KNOP dispatch imports `HNU_AT_OB` from
  `MOD_RAY_MOM_STORE`, which otherwise is not built until `d_plane`.
- All recursive makes run via `$(MAKE) -C dir` (12f1f1a) so a
  sub-make failure returns a non-zero exit status.  The previous
  `(cd dir; make)` form let build scripts print success after a
  failed compile; use `./build.sh`, not a bare `make | tee` with an
  unchecked exit code.

## 7. Verification checklist for any source change

1. `./build.sh` (`make KEEP_EXE=1 clean && make DEBUG=nodebug all`)
   completes, zero errors, and leaves `exe/cmfgen_dev.exe`,
   `exe/cmf_flux.exe`, and `exe/obs_frame.exe` present.
2. `diff -rq ../cmfgen_reference/cur_cmf .` shows only intended files.
3. Changed files actually recompiled (check build log for the ifort
   lines).
4. Symbols/strings present in `exe/cmfgen_dev.exe` where relevant.
5. Default-path equivalence: with the flag off, no edited code is
   reached (verify gating condition, not just intent).
6. Before claiming a feature works: a model run through
   `scripts/run_model_grid.py` (never raw sbatch) with `OUTGEN` showing
   the option echoed.

## 8. Git tracking (mandatory)

This tree is a git repository and **must stay tracked**. The working
tree state at any submitted run must map to a commit.

This section is a standing user directive for this tree: it takes
precedence over the project default in `../AGENTS.md` (do not push
unless requested) within `cmfgen_jfnk/` only.

Rules:

1. **Commit every verified source change.** One logical change per
   commit, matching the one-change-at-a-time rule. Commit after the
   build is clean and the change is verified - not before.
2. **Clean tree before submission.** A model submitted through
   `run_model_grid.py` must be built from a committed state. If the
   working tree is dirty at submit time, either commit the change or
   state explicitly in the run notes that the binary does not match
   any commit.
3. **Record the commit.** Note the short hash (`git rev-parse
   --short HEAD`) in the dated status entry (in `AGENTS_HISTORY.md`)
   for any run that depends on a code change.
4. **Never commit build artifacts.** `mod/`, `mod_dbg/`, `*.mod`,
   `*.o`, `*.a`, `exe/`, `lib/` and `build_*.out`/`build_*.err` are
   covered by `.gitignore` - do not force-add them.
5. **Push every commit to `origin/main`, never to
   `upstream`.** `origin` is the private dev repo
   `git@github.com:Astro-Sean/cmfgen-nonmono.git` (private -
   CMFGEN is licensed code). `upstream` is the pristine CMFGEN repo
   (`pitt-phyast/cmfsource`), used only for fetching and diffing.
   The dev repo has a single working branch, `main` (renamed from
   `development` 03-Nov-2026). Local `main` tracks `origin/main`
   with `push.default=upstream`, so a bare `git push` lands on
   `main`. Push after every commit - `main` must always mirror
   local work; do not leave work unpushed. Commit messages should
   be short and concise. Never force-push or rewrite history.
6. **Do not commit run products.** Checkpoints (`FGINEXT_GS`,
   `SCRTEMP`, `POINT*`), `OUTGEN`, and model-directory files do not
   belong in this repo - they live in the run directories under
   `/ptmp`.

Baseline: the full dev7 diff was committed as `019e1a7`
(21-Sep-2026); later changes are individual commits on top.

## 9. Validation and testing protocol

Validation and unit tests must be conducted using the Python scripts
located in the `cmfgen_jfnk/unit_tests` directory. These scripts execute
small-scale examples with plots saved to `cmfgen_jfnk/unit_tests/plots`,
designed to mimic the specific data structures, edge cases, and physical
regimes CMFGEN is expected to encounter in production runs. Any code
modification must be verified against these tests.

## 10. Development checklist (living document - mandatory)

`docs/CHECKLIST.md` is the running ledger for this development line
and must be maintained continuously:

- Update it in the same edit as the change it records - a solver
  behaviour change, a new or retired VADAT key, a diagnostic, an exit
  class, a validation result, or a changed run status all land there.
  A stale checklist is the same class of defect as a stale comment.
- Boxes are marked `[x]` only when the item has evidence - an audit
  section, a commit hash, a job id, or a test count. `[~]` marks in
  progress or done-with-limitations; `[ ]` stays open until the
  evidence exists. Blocked or failed items are labelled as such, never
  silently removed.
- The file is the handoff document for resuming this work (including
  by pasting it into an external review): keep items short, cite the
  audit/test/job that backs each one, and keep open items honest
  about what is still unmeasured.
