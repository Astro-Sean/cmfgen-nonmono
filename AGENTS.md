# AGENTS.md - cmfgen_jfnk

You are an expert in radiative transfer and numerical methods, with specific expertise in non-monotonic velocity fields, Sobolev approximations, and Jacobian-Free Newton-Krylov (JFNK) solvers applied to stellar atmosphere modeling (CMFGEN).


## 1. Purpose

This tree is the **clean restart** of the CMFGEN non-monotonic velocity
development pipline. It was extracted directly from
the pristine upstream archive `cur_cmf_18jun25.tar.gz`
(md5 `1d5a7a168e77329dd97c3019bcad3c09`) and carries a **minimal diff**.

## 1.1 Development roadmap (user directive 18-Sep-2026)

1. `DO_NONMONO_SOB` on a non-monotonic model; produce a spectrum.
2. Non-mono SOB on an interacting model (type IIn-like spectrum).
3. Knop arbitrary-velocity solver as a plug-in; touch base CMFGEN as
   little as possible.



Run area: `/ptmp/sbrennan/SN2024pba/cmfgen_validation` (config
`configs/dessart_m15_2p82d_****.yaml`).

Pipeline defaults (implemented in scripts/, not this tree):

- `ADJUST_CORRECTIONS`, `IT_SPECIFIER`, `table*.dat` excluded from
  run directories by default (opt-in only).
- Grid adjustments off: `REV_RGRID=F`, `N_OB_INS=0`; the initial
  `SN_HYDRO_DATA` arrives pre-optimized.
- `IT_ON_T=F` - no initial T iteration; hydro T used as-is.
- `THK_CONT=F` (10-Oct-2026) - thin outer continuum boundary; the
  m15 wind models are optically thin at the outer edge. Thick stays
  available via `vadat_overrides: THK_CONT: T`.
- Interactive (debug) partition jobs use the QOS max of 16 CPUs.
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
   documented, and physically justified (see the repository
   `../AGENTS.md`).
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
  `USE_LAM_ES=T` remains the validated path; `=F` warns but is
  permitted as an opt-in diagnostic (see next entry). `LAM_SET=T`
  warns only.
- `USE_LAM_ES=F` under `DO_NONMONO_SOB` (opt-in, rewired 12-Oct-2026):
  continuum J comes from `MOM_J_CMF_NONMONO_V1` (the Knop sign-split
  moment solver native to non-monotonic fields), paired with fixed-J +
  `GET_MOM_ALO_DIAG` + `GET_MOM_CELL_QUAR` in `do_var_cont.f` - the
  same recipe as KNOP+LAM_ES=F. `FG_J_CMF_V13` still runs each
  frequency, so `N_ON_J` factors come live from its ray moments via
  `GET_MOMS_NON_REL`; `SIGMA_TRUE` drives the Sobolev lines and the
  solver GAM term. Warns, no longer STOPs. JFNK stays inactive in
  this mode (admission is gated on `USE_LAM_ES`). `DO_KNOP_REFRESH`
  defaults on here as well as under KNOP (`MOM_KNOP_REFRESH_V1`
  converges the lagged J_NEXT/H_NEXT map within each sweep). Status
  on m15: refresh enabled is necessary but not sufficient - job
  12160690 diverged 3.9E8->1.6E17 once the map matured (lagged-map
  drift ~4x/sweep, 94.5% compressing cells, |EPS|/HU~5E6). Whether
  the expansive map is model-inherent or SOB-specific is under test
  in KNOP_MOM_probe (job 12161004). The earlier V11+V12 variant
  diverged on m15 (jobs 12158938, 12159013); do not restore it -
  V11 assumes expanding flow and V12 linearizes the V11 system.
- `BA_DEM_RELFRAC` (VADAT key, default 0=off, 12-Oct-2026): relative
  demand floor in `SOLVEBA_V13` - rows below `RELFRAC*max(POPS)` at a
  depth are excluded from the MAXCH/MAXCH_APP metric (their corrections
  still commit). Catches physically empty ions at populated depths that
  `BA_DEM_PFLOOR` cannot reach.
- `DO_FG13_ALO` (VADAT key, default F, 12-Oct-2026): `FG_J_CMF_V13`
  stores the same diagonal OAB approximate-Lambda operator the Knop
  solver uses (`GET_FG13_ALO_DIAG`); `do_var_cont.f` fills the
  continuum dJ/dCHI, dJ/dETA diagonals under `DO_NONMONO_SOB`. J
  stays the `USE_LAM_ES=T` formal solution. Opt-in; under test in
  `dessart_m15_2p82d_SOB_jfnk_alo_probe` (job 12159725).

## 4. Deliberately NOT ported from V6

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

## 5. Build

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

## 6. Verification checklist for any source change

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

## 6.1 Git tracking (mandatory)

This tree is a git repository and **must stay tracked**. The working
tree state at any submitted run must map to a commit.

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
5. **Push to `origin`, never to `upstream`.** `origin` is the private
   dev repo `git@github.com:Astro-Sean/cmfgen-nonmono.git` (private -
   CMFGEN is licensed code). `upstream` is the pristine CMFGEN repo
   (`pitt-phyast/cmfsource`), used only for fetching and diffing;
   `remote.pushDefault=origin` routes a bare `git push` to the dev
   repo. Never force-push or rewrite history.
6. **Do not commit run products.** Checkpoints (`FGINEXT_GS`,
   `SCRTEMP`, `POINT*`), `OUTGEN`, and model-directory files do not
   belong in this repo - they live in the run directories under
   `/ptmp`.

Baseline: the full dev7 diff was committed as `019e1a7`
(21-Sep-2026); later changes are individual commits on top.


VALIDATION AND TESTING PROTOCOL:
Validation and unit tests must be conducted using the Python scripts located in the `cmfgen_jfnk/unit_test` directory. These scripts execute small-scale examples  with plots saved to `cmfgen_jfnk/unit_test/plots`  designed to mimic the specific data structures, edge cases, and physical regimes CMFGEN is expected to encounter in production runs. Any code modification must be verified against these tests.
