# SOBEW_GRAD family short-ray stale-GAM audit

Base commit: `cd0f4374f1be07645a6a79f205a48f4285ea864b` (HEAD now
`247c9ad6`, knop-only) plus the uncommitted working-tree fix. See
SOBEW_SHORT_RAY_GAM_AUDIT.md section 0 for the repository-state
reconciliation.

## 1. Scope

Three sibling routines of `SOBEW`, each combining the EW diagnostic
with a radiative force multiplier (`GLINE(I) += TB(I)*EXPONX(GAM(I))
*HAQW(I,LS)` -> `FORCE_MULT`):

- `subs/sobew_grad.f`   (v1, `DIF` logical inner boundary)
- `subs/sobew_grad_v2.f` (`INNER_BND_METH`, `THOMAS_RH` solver)
- `subs/sobew_grad_v3.f` (adds `EMISS_LINE_FLUX`/`EMISS_EW`)

All three carried the same defect as `SOBEW`: `ZALONGP` + `GAM`
computed only inside the `NI>2` branch, while the common post-branch
loop evaluated `EXPONX(GAM(I))` for `I=1,NI` - so `NI=1` and `NI=2`
rays reused the preceding `NI=3` ray's `GAM`.

## 2. Call-site and reachability verification

Searched the whole tree (not relying on comments):

| routine | production callers | reachable? |
|---|---|---|
| `SOBEW_GRAD` (v1) | **none** | dead code in this tree. Fixed anyway for consistency with the family fix; fork comment already said "not currently called, but should not carry the unguarded form". |
| `SOBEW_GRAD_V2` | `obs/cmf_flux_sub_v5.f:2980` (DO_NONMONO_SOB, `SIGMA_TRUE`), `:2986` (default, `SIGMA`) under `IF(.NOT. DO_CMF_EW)` -> `ABS_EW`; `disp/maingen.f:3956` and `:6486` (force-multiplier + EW in the disp analysis tool) | YES - obs/disp analysis path, called with the model `NC,NP,ND` (signature `CMF_FLUX_SUB_V5(ND,NC,NP,...)`, `AQW(ND,NP)`, `P(NP)`) |
| `SOBEW_GRAD_V3` | **none** | dead code in this tree; fixed anyway |

With `NP=ND+NC` (the tested model: `NP=85=70+15`) the last two rays
have `NI=2` and `NI=1`, so the V2 stale values were live in the
obs/disp `ABS_EW` and `FORCE_MULT` outputs.

## 3. Denominator semantics (Part 6 check)

Fork (all four files, pre- and post-fix): `T1 = ABS(1+sigma*mu^2)`
with `T1 = MAX(ABS(T1), -0.2*(1+sigma), 1e-4)` for `sigma<0` — the
adopted `tau_S = GAMH/|1+sigma*mu^2|` prescription (Lamers,
Cerruti-Sola & Perinotto 1987; DAH15 floor retained verbatim).

Reference tree:

- `sobew.f`, `sobew_grad_v2.f`, `sobew_grad_v3.f`: `GAMH/ABS(1+s*mu^2)`
  (unfloored).
- `sobew_grad.f`: `GAMH/(1+s*mu^2)` — **signed, no ABS**. For
  `sigma<-1` the reference v1 produces sign-flipped/huge GAM and NaN or
  +Inf bursts — a separate upstream defect unrelated to the stale-ray
  bug. The fork's v1 already used ABS+floor; preserved.

## 4. Correction applied

Identical minimal relocation in each file:

1. `CALL ZALONGP(R,Z,P(LS),NI)` moved before the `THICK`/`IBOUND` block
   and the NI-specific branches (immediately after `AV(1:NI)=0.0`).
2. The complete `T1`/`GAM` block moved with it, computing `GAM(1:NI)`
   for every ray length.
3. Duplicate `ZALONGP` removed from `NI=2`; `ZALONGP`+`GAM` block
   removed from `NI>2`.

Preserved per routine: the distinct transfer solvers (`NORDTAU`+2-pt
formulas for `NI=1,2`; `TCOMPD`/`THOMAS` or `THOMAS_RH` for `NI>2`),
inner-boundary treatments (`DIF`/`INNER_BND_METH`, `TOR`/`IBOUND`),
`GAMH` sign, the DAH15 floor, `EXPONX`, `FORCE_MULT`/`EW`/`EMISS_EW`
formulas. `SOBJBAR_SIM` unchanged.

Post-fix anchors: `sobew_grad.f` ZALONGP:73, GAM:94, branches:113/116,
EXPONX:181; `sobew_grad_v2.f` 88/108/127/130/233; `sobew_grad_v3.f`
91/111/130/133/236.

## 5. Test results

`python/test_sobew_short_ray.py` — 8 checks per routine (see
SOBEW_SHORT_RAY_GAM_TESTS.md):

- `sobewg2`: 8/8 PASS (all via new-vs-ref ratio discrimination).
- `sobewg3`: 8/8 PASS.
- `sobewg1`: 8/8 PASS — `sigma=+2`, `sigma=-1`, negative-GAMH checks
  use the ref ratio; `sigma=-1.05`/`sigma=-2` checks use the absolute
  current-ray prediction because the reference v1's signed denominator
  makes its output nonfinite/meaningless there (meas-vs-pred agreement
  ~2e-9, e.g. T1: pred 1.7308527831516e12 vs meas
  1.7308527831476e12).

Force diagnostics: `GLINE` accumulates `TB*EXPONX(GAM)*HAQW`, so the
same stale-GAM values fed `FORCE_MULT`; the fix corrects both EW and
force outputs on the two outermost rays.

## 6. Output impact

For the reachable `SOBEW_GRAD_V2` (`obs/cmf_flux_sub_v5.f`,
`disp/maingen.f`): `ABS_EW`, `EW`, and `FORCE_MULT` change wherever
the two outermost rays carry weight — order-percent level for
uniform-weight quadratures (see the synthetic EW table in
SOBEW_SHORT_RAY_GAM_AUDIT.md section 6). No impact on populations or
the statistical-equilibrium solver.

## 7. Limitations

- v1 and v3 remain uncalled; the fix is consistency, not a live-path
  repair.
- The `sigma<-1` sobewg1 checks cannot compare against the reference
  (signed denominator) and rely on the absolute prediction; a shared
  bug in the prediction itself is excluded by the cross-checks passing
  for sobew/sobewg2/sobewg3 in ratio mode.
- No permanent production diagnostics were added; the GAM dumps used
  throwaway `/tmp` instrumentation only.
