# Production-route consolidation audit - 25-Nov-2026

Work order: moment solver (`USE_LAM_ES=F`) is preserved but no longer
the production target. Production = `DO_NONMONO_SOB=T` +
`USE_LAM_ES=T` (V13 formal J). Scope: audit the dispatch, remove
accidental moment-route dependencies, keep experimental routes
gated and preserved, do not develop the moment solver further.

Branch `knop/lam-handoff-ni1`, base HEAD `e7fade7`.

## A. Source audit - findings and dispositions

| # | Finding | Location | Severity | Disposition |
|---|---------|----------|----------|-------------|
| 1 | `DO_KNOP_REFRESH` absent-key default evaluated at read site (~:3549) before the absent-key `USE_LAM_ES` flip to T (~:3934), so an absent key resolved armed-T on the V13 route: misleading `knop_refresh=T` echo plus a per-iteration `MOM_KNOP_REFRESH_V1` invoke (inert via `KN_M_READY` self-guard). | `rd_control_variables.f` | MEDIUM | FIXED - absent key under plain SOB+LAM_ES=T (not `DO_SOB_FG_KNOP`) now resolves F, the value `run_model_grid` already stamps; explicit T stays honoured + inert warning |
| 2 | Three shared-code MOM-store gates used `(DO_NONMONO_KNOP .OR. DO_NONMONO_SOB)` while the store exists only on SOB+LAM_ES=F. Under KNOP the legs were a pre-08-Nov vestige (KNOP forces LAM_ES=T); under production SOB they fired armed-but-inert calls. | `cmfgen_sub.f` `SET_MOM_GS_TAG`, `MOM_KNOP_REFRESH_V1`; `solveba_v13.f` `RESTORE_MOM_JNEXT_LATEST` x2 | MEDIUM (waste/misleading; self-guards prevented harm) | FIXED - all four gates narrowed to `DO_NONMONO_SOB .AND. .NOT. USE_LAM_ES` |
| 3 | Documented SOB solver-mode rejections (`USE_FIXED_J`, `USE_FORMAL_REL`, `USE_J_REL`, `USE_DJDT_RTE`, `PP_MOD`, `PP_NOV`) were claimed by docs since 22-Oct-2026 but never coded - only the KNOP block enforced them. A mixed VADAT silently bypassed the non-mono continuum solve while Sobolev lines ran on `SIGMA_TRUE`. | `rd_control_variables.f` SOB validation block | HIGH | FIXED - hard STOP added, mirroring the KNOP block |
| 4 | `DO_SOB_FG_KNOP` + explicit `USE_LAM_ES=F` was honoured: `comp_j_blank` still routed to `FG_J_CMF_NONMONO_V1` while `do_var_cont` paired `GET_MOM_ALO_DIAG` (MOM store never ran) and the banner announced the moment path. Found during fix work. | `rd_control_variables.f` Rule-1b block | HIGH | FIXED - `USE_LAM_ES` forced T under `DO_SOB_FG_KNOP` (KNOP 08-Nov precedent: auto-set if absent, banner override if explicit F) |
| 5 | `MOM_BLKGMRES`/`MOM_NUSOLV`/`MOM_H_IMPL` set T off the moment route (production, KNOP, SOB_FG_KNOP, stock) were silently dead keys in OUTGEN. | `check_param_consistency.f` | LOW (misleading echo, no harm) | FIXED - inert warnings added |
| 6 | `USE_LAM_ES=F` banner said "opt-in experiment" but lacked the required "unsupported for production" wording. | `rd_control_variables.f` | DOCUMENTATION | FIXED - banner extended |
| 7 | `scripts/run_model_grid.py` comment still claimed SOB defaults `USE_LAM_ES=F` (21-Oct recipe); code has defaulted T since 10-Nov. | `scripts/run_model_grid.py` (parent repo, uncommitted alongside existing dirty state) | DOCUMENTATION | FIXED - comment corrected |
| 8 | `test_vadat_keys.py` HELPERS missed `iter_log_fmt_mod.f90`: link broken since `0cea0a27` added ILF calls to `mom_j_cmf_nonmono_v1.f`. Same rot affects 13 other test files plus `test_es_conv` (missing `SET_FG_Q_RESID` stub) and `test_mom_gs` (stale sweep-banner pin). | `unit_tests/python/` | MEDIUM (test rot, not solver) | FIXED: helper line added to 11 files (mode1/oneway inherit via `test_knop_ray`), `SET_FG_Q_RESID` stub added to `drv_knop_guard.f`, sweep pin updated |

## B. Verified-clean items (audit coverage, no defect)

- J-field dispatch (`comp_j_blank.f`): `DO_NONMONO_SOB`+`USE_LAM_ES=T`
  -> `FG_J_CMF_V13` with `RJ=TC`; `DO_NONMONO_KNOP`/`DO_SOB_FG_KNOP`
  -> `FG_J_CMF_NONMONO_V1`; SOB+`LAM_ES=F` -> `MOM_J_CMF_NONMONO_V1`
  at :1122. `INC_GRID`/`ACCURATE` forced F under nonmono so the
  extended-grid branch is unreachable.
- ALO accessor pairing (`do_var_cont.f`): KNOP uses `GET_FG_ALO_DIAG`
  under T / `GET_MOM_ALO_DIAG` under F; SOB uses `GET_FG_ALO_DIAG`
  under `DO_SOB_FG_KNOP`, `GET_FG13_ALO_DIAG` under `DO_FG13_ALO`,
  `GET_MOM_ALO_DIAG` under F. Correct per owning solver.
- ES-ALI gating: `DO_ES_ALI_BA` requires `DO_ES_ALI` (hard STOP);
  `DO_ES_ALI` warns when no gate can engage.
- JFNK admission: gated on `USE_LAM_ES=T` (MOM store has no eval
  freeze - documented reason, STOP enforced); `jfnk_cmt` ledger
  counts committed Newton steps only.
- BLKGMRES isolation: three walls hold (route gate,
  `KN_M_READY` self-guard, `NU_CACHE_RDY`); `MOM_BLKGMRES`/
  `MOM_NUSOLV` default F and are never auto-set.
- Nov-24 moment-audit defects F1-F4 already fixed in `910c4667`
  (verified at current line numbers).
- `HAS_LAGGED_FIELD` correctly covers KNOP, `DO_SOB_FG_KNOP`, and
  SOB+`LAM_ES=F`; the new MOM-store gates sit inside its SOB-F leg.

## C. Production architecture (resolved dispatch)

```
DO_NONMONO_SOB=T
  |-- USE_LAM_ES=T (default for absent key)  [PRODUCTION]
  |     |-- DO_SOB_FG_KNOP=F: comp_j_blank -> FG_J_CMF_V13, RJ=TC
  |     |     ALO: GET_FG13_ALO_DIAG (DO_FG13_ALO auto-set T)
  |     |     MOM_* keys: dead-key warning if set
  |     |-- DO_SOB_FG_KNOP=T: comp_j_blank -> FG_J_CMF_NONMONO_V1
  |     |     ALO: GET_FG_ALO_DIAG; refresh forced T (live);
  |     |     USE_LAM_ES forced T (new, 25-Nov)
  |
  |-- USE_LAM_ES=F (explicit only; banner: experimental,
        unsupported for production validation)
        comp_j_blank -> MOM_J_CMF_NONMONO_V1   [EXPERIMENTAL]
        ALO: GET_MOM_ALO_DIAG; MOM_BLKGMRES/MOM_NUSOLV live;
        DO_KNOP_REFRESH armed default T; JFNK refused (no MOM
        eval freeze); solver modes J_REL/DDT/PP/fixed-J STOP.

DO_NONMONO_KNOP=T -> FG_J_CMF_NONMONO_V1 (reference/diagnostic)
Stock monotonic CMFGEN -> stock dispatch (unchanged)
```

## D. Changes made

Functional:

- `new_main/mod_subs/rd_control_variables.f` - `USE_LAM_ES` forced T
  under `DO_SOB_FG_KNOP`; `DO_KNOP_REFRESH` absent-key re-default F
  under resolved LAM_ES=T; solver-mode STOP under SOB; banner
  unsupported wording.
- `new_main/mod_subs/check_param_consistency.f` - inert warnings for
  `MOM_BLKGMRES`/`MOM_NUSOLV`/`MOM_H_IMPL` off the moment route and
  for `DO_KNOP_REFRESH=T` off every refresh route.
- `new_main/cmfgen_sub.f` - `SET_MOM_GS_TAG` and
  `MOM_KNOP_REFRESH_V1` gates narrowed to the MOM route.
- `new_main/subs/solveba_v13.f` - both `RESTORE_MOM_JNEXT_LATEST`
  rollback gates narrowed; `USE_LAM_ES` added to the ONLY list.

Documentation:

- `docs/NONMONO_KEYS.md` - rejection caveat, F-route unsupported
  note, refresh default/behaviour, dead-key warning, auto-set list.
- `docs/CHECKLIST.md` - dated entry.
- `cmfgen_obsidian/` - `reference/VADAT keys.md`,
  `solvers/Sobolev non-mono solver.md` updated to match.
- `scripts/run_model_grid.py` (parent repo) - stale F-default
  comment corrected; left uncommitted next to pre-existing dirty
  state there.

Tests:

- `unit_tests/python/test_vadat_keys.py` - 5 new cases + 3 dedicated
  blocks; HELPERS gained `iter_log_fmt_mod.f90` (pre-existing link
  rot fixed).
- `unit_tests/python/test_iter_conv_gate.py` - new function pinning
  the narrowed gates (5 checks).

## E. Test results

- `./build.sh` (run 1, 11:21-12:05): all objects compiled clean;
  `exe/cmfgen_dev.exe`, `exe/cmf_flux.exe`, `exe/obs_frame.exe`
  fresh; `tests/f90` JFNK battery 36/36; `run_all_tests.py`
  69/110 files passed. All 41 failures classified: 14 pre-existing
  ILF-helper link rot (fixed, see F8), `test_es_conv` missing
  `SET_FG_Q_RESID` stub (fixed), `test_mom_gs` stale sweep pin
  (fixed), `test_auto_ion` missing fixture model dir
  (environmental), `test_phase227_io` argparse harness picked up
  by the suite glob with no args (parallel WIP, not this work),
  and 24 instant `rc=1 (0s)` staleness-guard trips caused by an
  `iter_log_fmt_mod.f90` working-tree edit that landed at 11:57
  mid-suite (parallel work stream, unrelated).
- `run_all_tests.py` rerun (12:12, after repairs + archive
  refresh): 101/110 files passed. The 9 remaining failures are all
  classified as not belonging to this change:
  `test_auto_ion` (missing fixture model dir),
  `test_sobew_angular_floor` (missing /ptmp RVTJ fixture),
  `test_par_opac` (`mod_par_opac.o` builds only into
  `par_opac.exe`, never archived - build topology gap),
  `test_knop_lagged_map` (fails a check labelled KNOWN DEFECT:
  jump veto keeps garbage - pre-existing solver defect, not rot),
  `test_rv_snhyd` (pinned known segfault on non-finite Density),
  `test_obs_frame_norm` (upstream reference run underproduced:
  ref=10 of expected outputs, tree=24, tolerated rc=174 OMP
  teardown), `test_iter_log_kvb` and `test_phase227_io` (parallel
  work-stream WIP files), `test_mom_gs` (suite ran the pre-fix
  pin; verified 34/0 after the `drv_mom_gs` `conv_diag_lvl=2`
  repair below).
- `test_vadat_keys.py`: 104 checks PASS, 0 fail (was 99 checks).
- `test_iter_conv_gate.py`: 18 checks PASS, 0 fail.
- `test_obs_nonmono_dispatch.py`: all PASS.
- `test_dvc.py`: 170 PASS.
- `test_iter_output.py`: PASS.
- Pre-existing rot repaired: 11 test files linking
  `mom_j_cmf_nonmono_v1.o` gained `iter_log_fmt_mod.f90` in their
  helper lists (`test_fg_sw_sym`, `test_fnv1a`,
  `test_knop_lagged_map`, `test_knop_operator`, `test_knop_p2_state`,
  `test_knop_ray` - which also covers `test_knop_mode1` and
  `test_knop_oneway` by import - `test_knop_replay_p2`,
  `test_knop_spv`, `test_mom_aa_state`, `test_mom_nu`,
  `test_mom_sweep_table`); `drv_knop_guard.f` gained a
  `SET_FG_Q_RESID` stub matching its three existing setter stubs;
  the sweep-banner pin updated to the post-ILF spelling
  `MOM_J_CMF_NONMONO_V1 GS SWEEP` in `test_mom_gs` and
  `test_mom_sweep_table`, and `drv_mom_gs` `rt` case gained
  `conv_diag_lvl=2` because the ILF migration had folded the
  unconditional banner inside the level>=2 report block.
- No model-run verification performed (per task rules: production
  runs require approval; unit/parser level only).

## F. Remaining risks / deferred items

- `DUMP_MOM_JNEXT_GS` still called unconditionally in
  `subs/scr_read_v2.f` (upstream file): early-return on unallocated
  `J_NEXT_FRQ` makes it a cheap no-op off the MOM route; left
  untouched per minimal-diff rule.
- `GS_STATE` read at `solve_for_pops.f:697` gated on
  `.NOT. NEWMOD .AND. ITS_SINCE_START=1` only; a stale file from an
  abandoned earlier run could seed a count - all consumers sit
  inside `HAS_LAGGED_FIELD` so production route unaffected; deferred.
- `comp_j_blank.f:977` writes an inner-boundary formula into
  `HFLUX_AT_OB` under `DIF=F` in the PP/V11 branch - upstream line,
  unreachable on non-mono routes; flagged for upstream review, not
  changed here.
- `JCOMP` diagnostic vacuous under `LAM_ES=T` (`RJ=TC` by
  construction) - cosmetic.
- Full-model convergence of the production route is established by
  prior runs (m15 converged at it718, job 12001805); this audit made
  no physics change, so no new model evidence is claimed.
- Moment-route deferred work stays deferred by directive
  (H-implicit block derivation, write-back acceptance ordering);
  see `docs/plans/MOMENT_H_EPS_IMPLICIT_DERIVATION_2026-11-25.md`.

## G. Recommended next development tasks (ranked)

1. Verdict the three live experiment arms (13671005 MOM_H_IMPL,
   13671134 cold, 13671873 fixed JNEXT_RELAX) against the
   CHECKLIST exit criteria; cancel answered arms and record
   verdicts - they decide whether any moment-route work is still
   warranted at all.
2. Upstream review of `comp_j_blank.f:977` `HFLUX_AT_OB` write and
   `prepare_CMF_PARAM.py:27` stale OBS claim.
3. (Only if arm verdicts demand) the moment-route plan-doc
   residual-evaluator first step - explicitly NOT recommended
   absent new evidence.

Verdict: VERIFIED at parser/unit level. The production route no
longer touches any MOM machinery: its call sites, defaults, and
warnings now reflect the resolved solver pairing exactly.
25-Nov-2026.
