# Build and git workflow

## Build

```bash
module load intel/2024.0 mkl
cd /viper/ptmp2/sbrennan/SN2024pba/cmfgen_jfnk
./build.sh           # make clean && make cmfgen; non-zero exit on failure
                     # -> exe/cmfgen_dev.exe
./build_debug.sh     # make dbg_clean && make DEBUG=debug cmfgen
                     # -> exe_dbg/cmfgen_dev.exe (-O0 -check bounds -traceback)
```

- Config `viper_intel_nodebug`: ifort -O3 + AVX-512,
  `-align array64byte`, MKL. PGPLOT stubbed (`pgplot_stubs/`).
- Spectrum executables: `make cmf_flux` (exe/cmf_flux.exe,
  obs_frame.exe). `make -C obs` hits a missing-`f77` for
  `disp/mod_disp.f` - build with ifort first or `(cd disp && make)`.
- On NFS/lustre, if a change does not appear in the binary, `rm` the
  stale `.o` before make (clock-skew suspicion is a known failure).

### Build-order constraints (enforced in the Makefiles)

- `d_nm_subs` depends on `d_mod_subs` (CONTROL_VARIABLE_MOD module
  file).
- `d_lte` is in the `cmfgen` target but builds `OBJ_FILES` only:
  `do_cmf_hydro_v2.o` references lte_hydro objects needed at link.
  The utility exes are excluded - `main_lte.exe` cannot link under
  -O0 (`lte_sub.f` calls `WR_ABUND_INFO_V2`, only V3 exists; -O3
  dead-code-eliminates it).
- `plane/mod_ray_mom_store.f` compiles as a prerequisite of
  `comp_j_blank.o` inside `new_main/mod_subs`: the KNOP dispatch
  imports `HNU_AT_OB` before `d_plane` is reached.
- All recursive makes go via `$(MAKE) -C dir` so sub-make failure
  propagates. Use `./build.sh`, not bare `make | tee`.

## Git rules (mandatory per AGENTS.md)

- Commit every verified source change, one logical change per commit;
  clean tree before submitting a run; record `git rev-parse --short
  HEAD` in the dated status note.
- `origin` = private dev repo `git@github.com:Astro-Sean/cmfgen-nonmono.git`;
  `upstream` = pristine `pitt-phyast/cmfsource` (fetch/diff only;
  `main` tracks `upstream/main`). Never push to upstream, never
  force-push. CMFGEN is licensed - keep the dev repo private.
- Never commit build artifacts (`mod/`, `*.o`, `*.a`, `exe/`, `lib/`,
  `build_*.out/.err`) or run products (`OUTGEN`, POINT*, FGINEXT_GS,
  SCRTEMP).
- Diff-vs-pristine hygiene: `diff -rq ../cmfgen_reference/cur_cmf .`
  must show only intended files; comments are part of the reviewable
  diff (see AGENTS.md Sect. 2.1 comment conventions - plain ASCII, no
  Greek symbols, banned-word list, `! Altered :` history lines).

## Verification checklist for any source change

1. `make cmfgen` clean, zero errors.
2. `diff -rq` vs pristine reference shows only intended files.
3. Changed files actually recompiled (check the ifort lines).
4. New symbols/strings present in `exe/cmfgen_dev.exe`.
5. Default-path equivalence: flag off => edited code unreachable.
6. A model run through `scripts/run_model_grid.py` (never raw sbatch)
   with OUTGEN showing the option echoed.

Status labels used in the tree docs: PASS / PASS WITH LIMITATIONS /
FAIL / BLOCKED; claims tagged VERIFIED / SUPPORTED / UNVERIFIED /
CONTRADICTED / BLOCKED.

## Related

- [[Model pipeline]]
- [[Checkpoints and state files]]
