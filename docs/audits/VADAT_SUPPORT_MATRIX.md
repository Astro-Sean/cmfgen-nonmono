# VADAT configuration support matrix (INT-001)

- Date: 2026-11-03
- Axes: transfer {stock, SOB, SOB+FG_KNOP, KNOP} x J update
  {USE_LAM_ES=T formal, USE_LAM_ES=F moment} x solver
  {conventional, JFNK} x matrix {DIAG, TRI}.
- Status codes: PROD = exercised on production m15 runs;
  TEST = covered by unit/systems tests; SRC = gate verified at
  source, no model evidence; REJ = loudly rejected; UNSUP =
  allowed but documented unsupported/unstable.

## Hard gates (all verified at source)

| Rule | Enforcement | Ref |
|---|---|---|
| DO_NONMONO_KNOP and DO_NONMONO_SOB mutually exclusive | STOP | rd_control_variables.f:858-865 |
| DO_JFNK requires DO_NONMONO_KNOP or DO_NONMONO_SOB | STOP | :1424-1433 |
| NUM_BNDS<3 forces METH_SOL=DIAG | loud banner + rewrite | :2951-2957; cmf_blkband_v3.f:437 |
| SOL_METH not in {DIAG,TRI} (PEN/GSIT/MIN) | STOP at solve | solveba_v13.f:1066-1071 |
| DO_NONMONO_SOB overrides GLOBAL_LINE to SOB | loud banner auto-set | :3473-3485 |
| JFNK_COUPLED rejected under DO_NONMONO_SOB | rejected at read | :3464-3473 comment |
| DO_CO_MOV_DDT / PURE_HUBBLE_FLOW / DO_HYDRO vs nonmono | STOP | :873-902 |
| stock + USE_LAM_ES=F on SN model needs USE_J_REL/USE_DJDT_RTE | STOP (nonmono exempt) | check_param_consistency.f:43-47 |
| JFNK banded preconditioner needs NUM_BNDS>1 | silent DIAG fallback | jfnk_cmfgen_int_mod.f:66-74 |

## Matrix

### Conventional solver, DIAG (production default)

| transfer \ J update | LAM_ES=T (formal) | LAM_ES=F (moment) |
|---|---|---|
| stock | PROD - upstream default path | REJ on SN unless USE_J_REL/USE_DJDT_RTE (check_param :43) |
| SOB | PROD - m15 runs (jobs 12667827/12667984) | UNSUP - default J path under SOB; P3-012 non-contraction documented |
| SOB+FG_KNOP | SRC - supported route; auto-sets N_LAM_RWARM=4, GS_REL_STABLE=2 (rd_control_variables.f:3426-3446) | inert flag - moment path runs anyway |
| KNOP | TEST - 10/10 drv_knop_ray + vault model history | UNSUP - 93/93 mom tests pass but P3-012 non-contraction on m15 |

### Conventional solver, TRI (NUM_BNDS=3)

| transfer \ J update | LAM_ES=T | LAM_ES=F |
|---|---|---|
| stock | SRC - upstream banded path; assembly audited (P5-003) | as DIAG cell |
| SOB | SRC - NUM_BNDS=3 fills ALO_TRI off-bands; m15 untested | UNSUP - same P3-012 caveat |
| SOB+FG_KNOP | SRC - same as DIAG cell | inert flag |
| KNOP | TEST - knop_ray driver solves the tridiagonal sweep | UNSUP - same caveat |

### JFNK solver (requires a nonmono flag)

| transfer \ J update | LAM_ES=T | LAM_ES=F |
|---|---|---|
| stock | REJ - STOP "DO_JFNK requires DO_NONMONO_KNOP or DO_NONMONO_SOB" | same |
| SOB | PROD - frozen tier admitted on m15; every admission defers on the FD non-descent guard (BUG-005). Coupled tier REJ at read | UNSUP - inherits the LAM_ES=F caveat; tier untested |
| SOB+FG_KNOP | SRC - same frozen-tier machinery as SOB cell | inert flag |
| KNOP | TEST - coupled tier ran historically (true/Hessenberg ratios 1.09-1.95 in vault logs); BUG-001 tres path reachable here | UNSUP - same LAM_ES=F caveat |

### JFNK solver, TRI (NUM_BNDS=3)

Same admission gates as the DIAG column; additionally
JFNK_USE_BANDED=T enables the banded-LU preconditioner
(DGBTRF/DGBTRS). Under NUM_BNDS=1 the banded preconditioner is
silently disabled (comment-documented). No production run has
exercised the banded preconditioner; unit coverage = the
preconditioner-factor tests in test_jfnk_regressions.py.

## Notes

- SOB+FG_KNOP only changes behaviour under USE_LAM_ES=T (it swaps
  the formal solver V13 -> FG_J_CMF_NONMONO_V1). Under
  USE_LAM_ES=F the moment path already runs; the flag is inert.
- "UNSUP" cells are not hard-rejected: they run, but the
  documented non-contraction (P3-012) means convergence is not
  guaranteed on strong-compression models.
- The conventional DIAG cell is the only PROD-verified column
  besides SOB+LAM_ES=T JFNK admission; every other cell derives
  from unit drivers, source gates, or vault history.
