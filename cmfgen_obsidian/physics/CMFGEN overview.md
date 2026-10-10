# CMFGEN overview

#status/verified

CMFGEN (Hillier & Miller) is a 1-D spherical, non-LTE
radiative-transfer / model-atmosphere code. For this project it is run
in **SN snapshot mode**: it reads a hydrodynamic structure (radius,
velocity, density, temperature, composition) from `SN_HYDRO_DATA` at a
fixed post-explosion time, then iterates the coupled problem

```
statistical equilibrium (level populations, T, Ne)
    <->  comoving-frame radiative transfer
```

until the populations stop changing, then produces a spectrum with the
`cmf_flux` / `obs_frame` codes (`obs/`).

Upstream: `pitt-phyast/cmfsource`, baseline tarball
`cur_cmf_18jun25`. This tree (`cmfgen_jfnk`) = `cmfgen_v7_ord_run`
(carrying the resacc fixes `567d774`/`a085c47`) + the V4/V6 JFNK
machinery re-ported on top (port commit `10272f6`).

## What the code solves

| Piece | Equation | Unknowns |
|---|---|---|
| Statistical equilibrium | rate equations per level per depth | `POPS(NT,ND)` |
| Charge conservation | sum z_i n_i = Ne | electron density `ED` |
| Radiative equilibrium | net heating/cooling = 0 (when `FIX_T=F`) | temperature `T` |
| Comoving-frame RT | $dI/ds = \eta - \chi I + a_\nu\, dI/d\ln\nu$ | intensities on rays x frequencies |
| Moment equations (optional J source) | variable-Eddington J,H,K system | `J` per freq x depth |

- $a_\nu = \Gamma = (V/R)(1 + \sigma\mu^2)/c$ is the velocity-coupling
  coefficient; see [[The non-monotonic problem]].
- Unknowns vector per depth block: level populations, then `ED`, then
  `T`; `NT` rows per depth, `ND` depths. `N ~ NT*ND ~ 1.3e5` on the
  Dessart m15 2.82d model.

## Iteration styles

CMFGEN has two classic population updates, switched by the mode
controller in `SOLVE_FOR_POPS`:

- **Lambda iteration** (`LAMBDA_ITERATION=T`): fixed-point step, cheap,
  bounded by `MAX_LAM_COR`. Safe far from the solution, slow near it.
- **Full / BA step** (`LAMBDA_ITERATION=F`): Newton-like solve of
  `A x = F` where `A` is the "BA" matrix - a block-diagonal-per-depth
  approximate Jacobian built in `GENERATE_FULL_MATRIX_V3`, factored
  with LAPACK (`DGEEQU`/`DGETRF`/`DGETRS`) inside `CMF_BLKBAND_V3`,
  driven by `SOLVEBA_V13`. See [[BA matrix and SOLVEBA]].

The layered additions in this tree sit on top of that: the
[[Field maturity and GS freeze|GS-freeze controller]] forces Lambda
while the lagged radiation field is still moving, and the
[[JFNK solver]] may replace the BA linear solve near the root.

## Where things live

```
new_main/cmfgen.f              main program
new_main/cmfgen_sub.f          CMFGEN_SUB: model setup + the great-iteration loop
                               (label 20000); also hosts JFNK_RAD_SWEEP_V1
new_main/mod_subs/             modules + solver-level routines:
  solve_for_pops.f             mode controller + all the gates (GS, BA, JFNK)
  solveba_v13.f (new_main/subs) BA solve + residual-acceptance rollback
  generate_full_matrix_v3.f    builds the BA matrix
  cmf_blkband_v3.f             banded/block solver (DIAG/TRIDIAG)
  do_var_cont.f                dJ/dchi,dJ/deta ALO assembly (TX/TVX)
  sub_sob_line_v3.f            Sobolev line rates + BA terms
  control_variable_mod.f       VADAT key definitions
  rd_control_variables.f       VADAT reader + flag guards
  mod_cmfgen.f                 global state (R,V,SIGMA,SIGMA_TRUE,...)
  jfnk_*.f/.f90                JFNK machinery (see [[JFNK solver]])
plane/                         formal solvers along rays:
  fg_j_cmf_v13.f               stock monotonic CMF solver
  fg_j_cmf_nonmono_v1.f        Knop solver: chi_hat, I_NEXT store, refresh
  mom_j_cmf_nonmono_v1.f       moment-equation variant (USE_LAM_ES=F)
subs/                          Sobolev kernels, sobjbar_sim.f, sobew*.f,
                               sob_rh_couple.f (RH78 non-local coupling),
                               scr_read_v2.f (checkpoints)
obs/                           cmf_flux.exe / obs_frame.exe spectrum stage
lte_hydro/                     LTE hydro utilities (objects only in `cmfgen` target)
knop_solver_files/             flat export of the touched files (reference set)
cmfgen_obsidian/  docs/        project docs vault; diagrams (copied into attachments/)
```

## The two radiation paths (`USE_LAM_ES`)

- `USE_LAM_ES=T` - J comes from the **formal solution** on rays
  (`FG_J_CMF_*`); this is the only path OBS supports and the
  paper-conform Knop route.
- `USE_LAM_ES=F` - J from the **moment solver** (`MOM_J_CMF_*`); a
  moment analogue, globally tridiagonal in depth per frequency. Its
  lagged map proved non-contractive on strongly compressing models
  (see [[Open issues]]).

> [!note] TRIDIAG is a no-op in this configuration
> With `USE_LAM_ES=T` the BA matrix has exactly zero off-diagonal depth
> blocks (VJ=0 branch in `do_var_cont.f`), so `SOL_METH=TRIDIAG` is
> bit-identical to `DIAG`. Verified empirically (AGENTS.md, 18-Sep-2026).

## Next

- [[The non-monotonic problem]] - why the stock solver fails
- [[The great iteration]] - the per-iteration call sequence
- [[Model pipeline]] - how a run is actually submitted
