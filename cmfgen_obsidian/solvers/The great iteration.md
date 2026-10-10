# The great iteration

What one CMFGEN iteration does in this tree, in call order. Verified
against `cmfgen_sub.f` (loop label 20000) in the 23-Sep-2026
iteration-cycle audit (AGENTS.md).

Analogy for the nesting: each iteration is a negotiation between two
parties who update at different speeds. The radiation field ("the
factory floor") must finish its own internal production run - every
frequency, the f/Eddington-factor loop - before the populations
("management") see the finished output and adjust policy once.
Management moving policy mid-production is exactly what the freeze
controller forbids.

![[iteration_cycle.png]]

## Call sequence

```text
per great iteration (label 20000 in cmfgen_sub.f):
  opacities (COMP_OPAC etc.)
  -> per-frequency inner loop over ML = 1..NCF:
       inner f-loop: FG_J_CMF_NONMONO_V1 (rays -> f, Jcmf)
                     <-> MOM_J_CMF_NONMONO_V1 (J)        [USE_LAM_ES=F]
                     or USE_LAM_ES=T: RJ = TC directly from rays
       until |d f| < ACC_EDD_FAC or 20 iters (bail = warning only)
       EVALSE_QWVJ_V8 etc. accumulate SE residuals + BA Jacobian
       on the converged RJ
       STEQ_T (radiative equilibrium row), STEQ_ED (charge cons.)
  -> FG_KNOP_REFRESH_V1 / MOM_KNOP_REFRESH_V1   [KNOP only]
       converge the lagged I_NEXT store on this sweep's opacities
  -> two-photon / charge-exchange / advection / adiabatic terms
  -> SOLVE_FOR_POPS:
       KNOP maturity gates -> FORCE_LAMBDA fixups
       -> JFNK gate -> maybe JFNK_SOLVE_CORR_V1
       -> SOLVEBA_V13 -> CMF_BLKBAND_V3 -> GENERATE_FULL_MATRIX_V3
           (+FIXPOP_IN_BA_V3) -> DGEEQU/DGETRF/DGETRS
       -> POPS *= (1 - STEQ*SCALE) with the MAJOR per-depth limiter
       mode-controller bookkeeping for the NEXT call
  -> POPTOION, SUP_TO_FULL, hydro/clump (+ SIGMA_TRUE refresh)
  -> COMP_J_CONV -> iteration controller
       (DO_T_AUTO, MAXCH<EPS convergence test, DO_LAM_AUTO, IN_ITS)
  -> SCR_RITE_V2 checkpoints (POINT1/2 + FGINEXT_GS + MOMJNXT_GS)
```

## Invariants verified in the audit

- **Newton sign**: `A.x = F(residual)`, `x` converted absolute ->
  fractional, multiplicative update `POPS*(1 - x*SCALE)`; positivity
  guaranteed by `BIG_LIM < 1`. `MAX_dT_COR=0.2` caps dT; `T_MIN=0.3`
  floors it.
- **MAJOR limiter** (`SCALE_OPT=MAJOR`): the step scale is set only by
  levels above `1e-10*Ne` - trace levels cannot throttle the model.
- **Fixed-point consistency**: at convergence, SE holds given `J(n)`
  and the transfer solve is self-consistent at `n` - the SE rates use
  the inner-loop-converged `RJ`, not a stale field.
- **Convergence test**: `MAXCH < EPS` (1%) requires
  `LAST_LAMBDA != MAIN_COUNTER` - cannot trigger inside warmup/freeze.
- `SET_POPS_D2_EQ_D1` (T in these configs): `STEQ(:,2)` zeroed
  pre-MAXCH, pops slaved to depth 1 post-update - the outer-BC helper.
- `SIGMA_TRUE` is recomputed after every R/V change (grid revision, NG
  acceleration, hydro).

## Config-level caveats recorded in the audit

> [!warning] Not code bugs, but they change what "converged" means
> - `MAX_CHNG_LIM=1e100` disables the runaway-abort guard - nothing
>   stops a divergence by itself.
> - `DO_LEVEL_CHK` masks MAXCH: it reports the 10th-largest change, so
>   up to ~10 levels may still move >10% while MAXCH reads <1%.
> - `FIX_T=T` + `DO_T_AUTO=F` => radiative equilibrium is never
>   enforced; T is the input structure (deliberate for snapshot runs).
> - `FIX_NE=T` (if ever set) freezes Ne at its initial value.
> - Inner f-iteration bails at 20 its with only an OUTGEN warning.
> - Jacobian reuse (`N_FIX_BA`, `FIX_BA`): below `FIX_BA` the BA can be
>   frozen; corrections then use the `ORIG_POPS` conversion, rebounded
>   by the `MAXCH_SUM > 3*FIX_BA` recompute trigger.

## Where the new machinery plugs in

| Stage | New code | Flag |
|---|---|---|
| Formal solve | `FG_J_CMF_NONMONO_V1` chi_hat/upwind | `DO_NONMONO_KNOP` |
| Post-sweep | `FG_KNOP_REFRESH_V1` converges `I_NEXT` | `DO_KNOP_REFRESH` |
| Line rates | `SUB_SOB_LINE_V3` + `SIGMA_TRUE` + tau floor | `DO_NONMONO_SOB` |
| Population update | GS freeze / mode controller gates | (always on under KNOP) |
| Population update | `JFNK_SOLVE_CORR_V1` replaces SOLVEBA | `USE_JFNK` + gates |

## Related

- [[BA matrix and SOLVEBA]]
- [[Knop formal solver]]
- [[Field maturity and GS freeze]]
- [[JFNK admission]]
