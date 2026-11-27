# escatter.md checklist diagnosis (27-Nov-2026)

Diagnosis pass on `checklists/escatter.md` (43 lines). The file names it
cites (`src/cmf_flux.f90`, `src/rates.f90`, `src/ionization.f90`,
`src/velocity_grid.f90`) do not exist in this fork. The real modules:
`cmfgen/subs/eval_ed.f`, `cmfgen/new_main/mod_subs/set_new_model_estimates.f`,
`cmfgen/new_main/mod_subs/grey_t_iterate.f`, `cmfgen/plane/fg_j_cmf_v13.f`,
`cmfgen/plane/fg_j_cmf_nonmono_v1.f`,
`cmfgen/plane/mom_j_cmf_nonmono_v1.f`, `cmfgen/obs/cmf_flux_sub_v5.f`.
HEAD: 08efa7f7 (origin/main).

## Verdict per item

### 1. Electron density summation - CONFIRMED CORRECT, no fix needed

`n_e = sum_i sum_j n_ij Z_ij` is computed per element by
`EVAL_ED` (cmfgen/subs/eval_ed.f): the J=1..NSPEC-1 loop accumulates
`LOC_CH = sum_j T1*ED*(Z_j + ...)` over every ionization stage in the
element's partition-function table, closed by `LOC_CH += Z(NSPEC)`
(the fully stripped stage). The driver loops
`DO ISPEC=1,NUM_SPECIES` over all `SPECIES_PRES` elements
(`set_new_model_estimates.f:687-695`, `grey_t_iterate.f:366`) and
Newton-iterates ED to 1e-4. Iron and every other modelled element's
high stages are inside their own EVAL_ED call by construction. The
prescription's failure mode (dropped trace-element electrons) does
not exist here.

### 2. Thermalization epsilon floor - REJECTED / N-A

There is no `EPS = kappa_a/(kappa_a+kappa_e)` variable to floor:
the continuum source function is evaluated directly as
`SOURCE_RAY = ETA_RAY/CHI_RAY` (fg_j_cmf_v13.f:1247), and the
scattering-vs-absorption split enters the moment/FEAUTRIER solve
through the opacities, not a fixed-point on epsilon. `EPS_NM` in
mom_j_cmf_nonmono_v1.f is the comoving-frame epsilon
(GAMH*sigma_true term, :2557) - a different object entirely.
A Mihalas epsilon floor would be a hidden-physics clamp on the
wrong quantity (AGENTS.md section 6 forbids it). Matrix
singularity under pure scattering is handled by the solver
structure, not a source-function floor.

### 3. Non-monotonic CMF frequency mapping - CONFIRMED, already handled

The `dv/dr<0` fold is exactly what the Knop non-monotonic
formal/moment solvers exist for: sign-split characteristic pairs
selected per local SIGMA_TRUE sign (KN_GEO_REPS in
fg_j_cmf_nonmono_v1.f; the upwind node is chosen by the sign of
EPS_NM at mom_j_cmf_nonmono_v1.f:2531). Signed optical-depth steps
with signed-upwind weights are the paper method (Knop et al.
2009/2010), not the prescription's `|ds|` absolute-path hack -
which would break the sign bookkeeping the solver depends on.
The stock V13 path clamps the gradient instead
(`fg_j_cmf_v13.f:1070` SIGMA_RAY pin), and ESEC itself is
`6.65E-15*ED` (cmf_flux_sub_v5.f:929,1014) with a cumulative-ES-depth
bug already fixed at :917-929.

### 4. Smooth Eddington factor - REJECTED / N-A

CMFGEN does not use an analytic f(tau) closure. Variable
Eddington factors are computed from the formal-solution ray
moments and stored on EDDFACTOR records (referenced by
mom_j_cmf_nonmono_v1.f:631,1991,5645; the mom_jrel_v*
family). `f = 1/3 + (2/3)(1-e^{-tau})^2` would replace an exact
ray-derived closure with a fitted one - worse, and it ignores
that electron scattering is already inside the opacity/ray
optics the factors are built from.

## Bottom line

All four prescriptions are either already correct in the fork or
target a different code structure. No source changes warranted.
The real, live electron-scattering defect is a different one: the
ES_ALI OAB update amplifying compressing-cell noise on the
USE_LAM_ES=F route (docs/CHECKLIST.md, `DO_ES_ALI=F` isolation arm,
job 13749064) - tracked there, not in this checklist.
