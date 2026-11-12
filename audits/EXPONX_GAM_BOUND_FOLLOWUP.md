# Follow-up: the claimed GAM >= -5000 lower bound

Work commit: `cd0f4374f1be07645a6a79f205a48f4285ea864b` (HEAD at time of writing: `247c9ad6`, knop-only commits on top).
Read-only analysis; no implementation change. Corrects and refines the
bound discussion in `EXPONX_NEGATIVE_OPACITY_AUDIT.md` /
`EXPONX_RECOMMENDATIONS.md`.

## 1. The claim under review

Earlier audit text used the argument

    GAM = GAMH/T1,  T1 floored at 1.0e-4  =>  GAM >= -5000

As stated in the review prompt, the denominator floor alone does **not**
establish that bound. `T1 >= 1e-4` only gives

    GAM >= GAMH/1e-4  (for GAMH < 0),

which is unbounded below unless `GAMH` itself has a demonstrated lower
bound. The bound needs *two* ingredients, not one.

## 2. The traced quantity chain

- `sobjbar_sim.f` (and identically `sobew.f:55`):
  `GAMH(I) = CHIL(I)*3.0E-10*R(I)/V(I)/FL`.
- `CHIL` at that point is the **SIM-group-summed** line opacity
  (`CHIL = sum of CHIL_MAT` over the blend group; per-line `CHIL_MAT`
  never reaches `EXPONX` - see `EXPONX_CALLSITE_TABLE.md`).
- `CHECK_LINE_OPAC` (`new_main/cmfgen_sub.f:4113-4141`) evaluates the
  *same algebraic expression* `CHIL*3e-10*R/V/FL` - i.e. `GAMH`
  itself - and when it is `< -0.5` flags `NEG_OPACITY`, logs to
  `NEG_OPAC`, resets `CHIL=1.0` and every `CHIL_MAT` in the group to
  `+1/NUM_SIM_LINES`. The JFNK sweep replays the same crush
  (`cmfgen_sub.f` ~7294-7306); the obs path mirrors it
  (`obs/cmf_flux_sub_v5.f:2880`); `sub_sob_line_v3.f:210,297` zeroes
  variation terms for flagged cells.

So `GAMH` *is* the group-summed angle-free Sobolev optical depth:
the `R/V` factor, the `3e-10` constant, and the frequency factor `FL`
are identical in both expressions, and the crush threshold is applied
to exactly the quantity that later forms `GAMH`.

## 3. Does `-0.5 < group_tau < 0` bound GAM?

With `CHECK_LINE_OPAC` active, any group reaching `SOBJBAR_SIM` with
negative summed opacity has `-0.5 < GAMH < 0`. Combined with
`T1 >= 1e-4` (the floor applies for `sigma<0`; for `sigma>=0`,
`T1=|1+sigma*mu^2| >= 1`):

    -5000 < GAM < 0      (strictly, GAM in (GAMH/1e-4, 0])

The `-5000` figure is therefore **justified - but only through the
combination** of (a) the `-0.5` crush threshold bounding `GAMH` and
(b) the `1e-4` denominator floor. The denominator floor alone does not
produce it; earlier wording that implied otherwise is corrected here.

Conditions under which the bound fails:

- `CHECK_LINE_OPAC` disabled (control variable; it was active in the
  tested model - `NEG_OPAC` contains 159 logged cells);
- any future path that computes `GAM` from a `CHIL` that did not pass
  through the group-sum check.

## 4. The bound does not remove the overflow channel

`EXPONX(X)` overflows below `X ~ -709.78` (double precision).
`-5000 < -709.78`, so inside the justified bound there remains the
interval

    GAMH in (-0.5, -0.071),  T1 ~ 1e-4   =>   GAM in (-5000, -710)
    =>  EXPONX(GAM) = +Inf.

i.e. the crush leaves a factor ~1000 of headroom in `GAMH` below the
overflow point when the angular denominator sits at its floor. The
bound narrows where the hazard lives but does not eliminate it.

## 5. Mathematical vs reachable vs populated

- Mathematically possible overflow: YES - `EXPONX` accepts any `X`,
  overflows for `X < -709.78`.
- Reachable code path: YES - an uncrushed group with
  `-0.5 < GAMH < -0.071` traversing a near-singular direction
  (`T1 ~ 1e-4`, requires `sigma < -1` locally and `mu` near the
  DAH15-floor regime) reaches `EXPONX` in `SOBJBAR_SIM` with
  `GAM < -710`. `sobew.f`/`sobew_grad*.f` reach the same
  `EXPONX(GAM)` for the EW diagnostic.
- Actually populated in the tested model: NO evidence. The run
  completed; `NEG_OPAC` shows only crushed cells; no NaN/Inf in the
  outputs examined. This is a latent vulnerability, not an observed
  failure.

## 6. Corrections applied to prior wording

- Earlier text "the floor bounds GAM >= -5000" is imprecise: replace
  with "with `CHECK_LINE_OPAC` active, surviving negative-opacity
  groups have `-0.5 < GAMH < 0`; together with `T1 >= 1e-4` this gives
  `-5000 < GAM < 0`. The denominator floor alone bounds nothing."
- The recommendation stands: if diagnostics ever show the
  `(-0.5, 0)` group-tau window populated at `sigma < -1` depths, the
  fix belongs in the callers (a saturation-aware treatment of
  `EXPONX`/`d_EXPONX_dX` arguments), not inside `EXPONX` itself.
