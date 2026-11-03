# P1-003 audit: shared module state and SAVE variables with new flags off

Date: 2026-11-03. Tree: `cmfgen_jfnk` @ `b12a9c27` (main).
Scope: every module-level or SAVE state introduced or touched by the
fork, classified by whether it can leak into the stock path
(`DO_NONMONO_SOB=F`, `DO_NONMONO_KNOP=F`, `DO_JFNK=F`).

Method: enumerated new SAVE declarations across the 109-file baseline
diff (`P1_BASELINE_DIFF_INVENTORY_2026-11-03.txt`); traced each new
variable's write sites and read sites for flag gating; checked the
unconditional call sites (`RECOMPUTE_SIGMA_TRUE`, `SET_FG_NU_NEXT`,
`FG_J_CMF_V13` ALO store, `SOLVEBA_V13` snapshot re-arm,
`SOLVE_FOR_POPS` counters) in executable code.

Classes: (a) written only under new flags - cannot leak;
(b) written unconditionally but read only under flags or for
printing - overhead only, no result effect;
(c) live in the stock path - needs justification.

## Class (a) - flag-gated writes only

| State | Writer gate |
|---|---|
| JFNK eval buffers in `cmfgen_sub.f` (`S_*`, `SV_*`, `UP_TWO_SAV`, `DOWN_TWO_SAV`, `FOLD_SAVE`, `SCL_SAV`, `NOSCL_SAV`, `QAGE_SAV`, `XMBAD_SAV`, `CELLQ_SAV`, `ACCESS_F_SAVE`, `COMPUTE_BA_SAVE`, `EDDINGTON_SAVE`) | DO_JFNK eval bracket only |
| `MOD_HOLD_LEV_DIS` (`mod_lev_dis_blk.f`) | written only in `jfnk_cmfgen_interface_v1.f:734-769` (eval bracket); read in `sup_to_full_v4.f` |
| `ES_ALI_NUNC/NFALL/DJMAX` (`comp_j_blank.f`) | DO_ES_ALI branch only (SOB+LAM_ES+COHERENT_ES or KNOP) |
| GS machinery in `solve_for_pops.f` (`N_CHURN_STREAK`, `HIST_*`, `PLAT_CNT`, `N_FROZEN*`, `N_MATURE`, `N_TRIAL_*`, `N_RELEASES*`, `GS_*`, `GS_DIV_*`, `N_REL_*`, `N_PRECOND*`) | `GS_FREEZE` forced F outside DO_NONMONO_KNOP |
| `JFNK_*` counters (`JFNK_FRZ_HITS`, `JFNK_CPL_HITS`, `JFNK_FRZ_HOLD`, `JFNK_LAST_DEM`) | JFNK admission path |
| `LUM_DIV_*` latch state (`LUM_DIV_LATCH`, `LUM_DIV_BASE/B2/B3`, `LUM_DIV_NCNT*`) | evaluated only under `DO_NONMONO_KNOP .AND. LUM_DIV_N>0` (solve_for_pops:781,805); `LUM_BLOCKED` therefore always F in stock mode |
| `POPS_GS_SAVE`, `GS_POPS_HELD` | hold armed only under `(GS_HOLD_POPS .AND. GS_FREEZE) .OR. LUM_POPS_HELD` (solve_for_pops:2092) - both terms KNOP-scoped |
| Nonmono module state (`SET_FG_NU_NEXT`, `mod_ray_mom_store`, `fg_j_cmf_nonmono`/`mom_j_cmf_nonmono` module arrays, `jfnk_*` modules, `knop_status_mod`) | written only inside the flag-gated dispatch (`comp_j_blank.f:835-862`) |
| `BA_QUAR_MASK`, `BAACC_PENDING`, `BA_TRUST*`, `BA_NREJ*`, `BA_HOLD_RESTORE_REVIVAL` | `DO_BA_RES_ACC`/`DO_BA_QUAR` machinery only |
| `KN_REPAIR_T` T-repair | `IF(KN_REPAIR_T)` gate, VADAT default F (`cmfgen_sub.f:1501`) |
| `GMNR_NMJ_G` N/H midpoint block (`get_moms_non_rel.f:210`) | opt-in flag |

## Class (b) - unconditional writes, flag-scoped or print-only reads

| State | Assessment |
|---|---|
| `ALO_*` store + `ALO_VALID`/`ALO_ND` (`fg_j_cmf_v13.f:435-471,1870-1953`) | allocated, zeroed and filled every call; consumed only via `GET_FG13_ALO_*` under `DO_FG13_ALO`/`DO_FG13_ALO_TRI`. Extra flops only. |
| `POPS_BA_SAVE`, `BAACC_COL_PREV`, `BAACC_COL_MIN`, `BA_RBK_MASK` (`solveba_v13.f:925-995`) | refreshed on every non-LAMBDA call (`IF(.NOT. LAMBDA_IT .OR. BA_SNAPSHOT_MISSING)`); read only inside `BAACC_PENDING` branches that never arm without DO_BA_RES_ACC. |
| `ITS_SINCE_START` | incremented every call (sfp:612); readers are KNOP/JFNK-gated sites plus the first-call GS_STATE restart adopt (sfp:631-700) whose adopted counters feed only the GS machinery above. |
| f-iteration-cap census counters (`comp_j_blank.f`), `RLUMST_NAN_CNT` | write-only telemetry, printed at end. |
| `BA_DEM_EXM`, `BA_DEM_EXR`, `BA_TR_ABS` | exported telemetry feeding `BA_DEM_EXIT_OK` (`cmfgen_sub.f:5009`); `BA_TRACE_EVAL` ANDs every output with `L_NM=DO_NONMONO_SOB .OR. DO_NONMONO_KNOP`, so with both flags off `BA_DEM_EXIT_OK` is identically T - the conjunct in the exit test is a no-op. |
| `MONO_WARNED`, `VNEG_WARNED`, `BAACC_R2GATE_PRINTED`, `GS_RAW_NUDGE_DONE`, `ITER_LEG_DONE`, `BA_MATRIX_STORED`, `PREV_COMPUTE_ENT`, `BA_REUSE_STALL` | warn-once / bookkeeping flags for gated machinery. |
| `K_ON_J_SAVE`, `_PREV` Eddington state (`get_moms_non_rel.f`) | pre-existing upstream state, not new. |

## Class (c) - non-gated edits that can change stock results

Each is a deliberate upstream-defect fix, documented in place; none is
a silent physics change.

| Edit | Stock-path effect |
|---|---|
| `get_moms_non_rel.f` degenerate-J guards (4 sites): `JNU_STORE<1e-30` -> K/J=1/3, H/J=0, N/J=0, mid N/J=0 | Changes Eddington factors only where upstream divided by a sub-tiny J and produced garbage ratios. Justified: undefined-ratio defect fix; documented "dev7" comments. |
| `get_moms_non_rel.f` realloc-on-ND-change | Fixes out-of-bounds writes when callers pass differing ND; upstream callers pass a fixed ND so benign. |
| `subs/ins_line*`/`cmfgen_sub` empty-level ZNET_SIM/ZNET/BETAC guards (29-Oct) | `RJ*CHIL/0=-Inf` then `ETAL*ZNET=NaN` poisoned `STEQ_T` on restarts with empty upper levels; now contributes nothing where `ETAL=0` exactly. Changes only runs that would have produced NaN. |
| `new_main/subs/cmf_blkband_pin_null_v1.f` | pins exactly-disconnected zero equations (row, column and RHS all zero); upstream `DGETRF` failed the whole depth block. Test: `test_cmf_blkband_pin.py`. |
| strict interp guards (`lin_interp`, `mon_interp`, `mon_interp_fast*`, `mon_int_funs_v2`) | loud abort on non-monotonic/out-of-range inputs vs silent garbage upstream. Bit-identity on valid paths: `test_interp_guards.py`. |
| W-series fail-loud IOSTAT/bounds guards (`get_ibound`, `read_bcd_mat`, `scr_read_v2`, `rd_cont_j`, `add_dust_eta`, `get_jh_at_prev_time_step`, `characteristics_v2*`) | abort/report where upstream read stale/garbage data (`ROBUSTNESS_AUDIT_2026-10.md` W5,W10,W11,W12). |
| `comp_j_blank.f` `'ZERO FLUX'` inner-boundary fix | documented non-gated bug fix per AGENTS.md comment conventions (upstream defect). |
| `solveba_v13.f` `BAACC_RES_MIN` update moved inside `BAACC_PENDING` gate (W1) | the ungated update itself was the defect; `BAACC_RES_MIN` is only read under DO_BA_RES_ACC. |

## Findings

- No unexplained class-(c) shared-state leakage found. Every new
  SAVE/module variable is flag-gated (a) or overhead-only (b); every
  non-gated edit in (c) is a documented defect fix whose stock-path
  effect is confined to inputs where upstream produced NaN, garbage,
  or a failed solve.
- `fg_j_cmf_v13.f:1069` `SIGMA_RAY <= -1 -> -0.999` is UPSTREAM code
  (cur_hmi, Aug-2015), not a fork addition - do not attribute it to
  this development.
- `RECOMPUTE_SIGMA_TRUE` under flags-off evaluates sigma on a scratch
  copy and warns once; writes nothing shared.
- Overhead note (class b): the ALO store and POPS_BA_SAVE snapshot
  machinery add per-call zeroing/copying of ND- and NT*ND-sized
  arrays in stock mode. Result-neutral; worth a gating tidy only if
  profiling shows cost.

## P1-004 addendum - formal-transfer dispatch, flags off

Stock chain verified per modified formal-path file:

| File | Flags-off behaviour |
|---|---|
| `comp_j_blank.f` | `FG_J_CMF_V13` called in the ELSE branch (:869) exactly as upstream; `MOM_J_CMF_V11` called in the nonmono `IF(.NOT.USE_LAM_ES)` chain only under `DO_NONMONO_*`; `dLOG_NU=LOG(FL/NU_NEXT_FREQ)` gated on nonmono flags; `MAXCH` exclusions `.AND. .NOT. DO_NONMONO_KNOP` are identity in stock. Non-gated edits: 'ZERO FLUX'->'ZERO_FLUX' typo fix (upstream dead branch - VADAT value carries the underscore; documented dev7), FCONV_* counters (telemetry), JFNK_IN_EVAL eval branches (JFNK-scoped). |
| `fg_j_cmf_v13.f` | Only additions are the ALO diagonal/TRI store + accessors (class b) and upstream `SIGMA_RAY` clamp (Aug-2015, not ours). Solver algorithm untouched. |
| `mom_j_cmf_v11.f` | 0 diff lines - unmodified. |
| `get_moms_non_rel.f` | class-(c) degenerate-J guards + realloc fix (see table above). |
| `do_var_cont.f` | two hunks: declarations + a single `ELSE IF(DO_NONMONO_KNOP .OR. DO_NONMONO_SOB)` branch holding all ALO/TX/quarantine fills. Stock ELSE path byte-identical. |
| `characteristics_v2*.f` | W12 `MIN(P,0.9999*R_GRID)` clamp (fires only where upstream built NaN Z on a tangent ray); `ABS(dV)` point-insertion criterion - inserts on steep deceleration that upstream missed (grid-accuracy fix; changes stock results only on decelerating profiles, where it adds resolution). |
| `cmf_formal_rel_v4*.f` | `MAX(0,.)` SQRT guards + `N_STORE>=1` floor (degenerate-input guards; USE_FORMAL_REL is opt-in anyway). |
| `add_dust_eta.f`, `get_ibound.f`, `get_jh_at_prev_time_step.f` | W5/W12 bounded-walk/IOSTAT/NaN guards - fail-loud vs upstream silent bad reads. |
| `fg_j_cmf_v12.f` | dead code - no callers in either tree; its I_CORE/negative-intensity repairs are unreachable. |
| `mod_cmfgen.f` | `SIGMA_TRUE(:)` + `NU_NEXT_FREQ` declarations; all consumers flag-gated. |

Conclusion: the flags-off formal-transfer path executes the upstream
`FG_J_CMF_V13` + `MOM_J_CMF_V11` pair; every fork change on this path
is flag-gated, telemetry, or a documented defect fix that only alters
inputs where upstream produced NaN/garbage (or decelerating profiles
for the ABS(dV) insertion).

## P1-005 addendum - Sobolev and line-transfer dispatch, flags off

| File | Flags-off behaviour |
|---|---|
| `sub_sob_line_v3.f` | `SOBJBAR_SIM`/`SOBEW` route to `SIGMA_TRUE` only under `DO_NONMONO_SOB .AND. .NOT. DO_SOB_CLAMPED_LINES`; ELSE calls with `SIGMA` as upstream. `DRC_SIM`/`SOB_RH_COUPLE` gated `DO_NONMONO_SOB .AND. DO_NONLOC_SOB`; `DRC_SIM=0` keeps the `-DRC_SIM` term an exact no-op. Non-gated: `ETAL_MAT>0` guards on `VB_SIM`/`VC_SIM`/`SRAT` (upstream divided unconditionally -> Inf/NaN on a bad Lambda iterate; guard returns 0), `dRATE_dUP` rewrite (algebraically identical: `RATE_FAC=EINA*POPS(NUP)*U_STAR_RATIO`, so `RATE_FAC*ZNET/POPS(NUP)` == `EINA*U_STAR_RATIO*ZNET`; avoids POPS(NUP)/U_STAR_RATIO denominators on transient iterates). |
| `sobjbar_sim.f`, `sobew.f`, `sobew_grad.f`, `sobew_grad_v2.f`, `sobew_grad_v3.f` | Dessart et al. tau floor `MAX(|T1|, -0.2(1+sigma), 1e-4)` applied only when `SIGMA<0`; `sigma>=0` path is `ABS(T1)` identical to upstream. `test_sob_tau_floor.py` asserts bit-identity vs reference at sigma>=0 and sigma=-0.5 (14/14 pass 2026-11-03); floor binds only at sigma<=-1 geometry. |
| `ins_line_v6.f`, `ins_line_obs_v5.f` | `.AND.` loop conditions rewritten as `DO WHILE(idx<=N)` + explicit `IF 'BLA' EXIT` - Fortran does not guarantee short-circuit evaluation; upstream could read `TRANS_TYPE(N_LINES+1)` when no blank transition exists. Same valid-path semantics; documented non-gated fix. `ins_line_obs_v5` also has comment-only edits. |
| `cmfgen_sub.f` (SIM_LINES) | `SIGMA_TRUE` computed unconditionally by `RECOMPUTE_SIGMA_TRUE` on a scratch copy (warn-once, no writes to `SIGMA`); consumers flag-gated. |

Test evidence: `test_sob_tau_floor.py` 14/14, `test_sob_rh_couple.py`
9/9 (homology -> zero pairs, V<=0 and NEG_OPACITY rows never couple).
Remaining limitation: line insertion and dRATE derivative fixes are
verified by inspection + compile, not by a dedicated unit test; a
full model A/B comparison remains open under P1-007/008.

## P1-006 addendum - population/statistical-equilibrium solver, flags off

Non-flag-gated stock edits in the BA/Lambda chain:

| Edit | Classification |
|---|---|
| `cmf_blkband_v3.f`: `POP_DIV=MAX(POPS, 1e-40*max(|ORIG_POPS|,1e-300))` floor on the `STEQ*(ORIG_POPS/POPS)` rescale (2 sites) | Defect fix: upstream divided by a transiently non-positive POPS -> Inf/NaN. Counted via `N_DIV_FLOOR` and sets `FLAG`; only alters results where upstream produced non-finite STEQ. Limitation: N_DIV_FLOOR is counted but never printed to LUER - activation is invisible unless debug dumps run. |
| `cmf_blkband_v3.f`: `CMF_BLKBAND_PIN_NULL_V1(C_MAT,STEQ(1,K),N,N_PIN)` before DGETRF | Defect fix: pins only index whose row, column, AND RHS are all exactly zero (0=0 subproblem); ordinary matrices take the upstream LU path. Prints `pinned N uncoupled zero rows at depth K` when active. Covered by pin-null unit test. |
| `cmf_blkband_v3.f`: BCD CSCRATCH/BCSCRATCH/DSCRATCH recall | Opt-in machinery: only entered under `WR_BA_INV .AND. .NOT. BA_COMPUTED`; missing cache files fall back to recomputing LU with a LUER warning. Stock solves unaffected. |
| `fiddle_pop_corrections_v2.f`: `DPTH_BIG/LIT_LIM` | Flag-gated: stock branch sets `+/-HUGE` - upstream behaviour (no bound). `DO_SOB_FG_KNOP` POPS floor block also gated. |
| `generate_full_matrix_v3.f`: `.AND. .NOT. JFNK_QUIET_IO` on debug prints | Internal runtime flag only set during JFNK evals; no VADAT key; stock output unchanged. |
| `fixpop_in_ba_v3.f` | Message text only ("depths this call"). |
| `solve_for_pops.f` / `solveba_v13.f` | `POPS_BA_SAVE`/`POPS_GS_SAVE` snapshots (class b, result-neutral), `BAACC_RES2_MIN` gate fix (W1 documented), `LUM_BLOCKED=F` sites dead in stock, BA_DEM telemetry prints, `ITS_SINCE_START` counter. JFNK branches gated `DO_JFNK`/`L_NM`. |

Conclusion: no unconditional physics change enters the stock BA
iteration; every non-gated edit is a documented defect fix confined
to inputs where upstream produced Inf/NaN or read out of bounds, or
is print/telemetry-only.
