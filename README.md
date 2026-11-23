# CMFGEN with non-monotonic velocity fields

> **Warning: this code is under active development.**
> Solver
> behaviour, input keys and outputs change between commits. Nothing
> here should be treated as production-ready or science-validated.

Development fork of [CMFGEN](https://sites.pitt.edu/~hillier/web/CMFGEN.html)
(baseline `cur_cmf_18jun25`, `pitt-phyast/cmfsource`) for supernova
ejecta where the velocity law is not monotonic. The case we care
about is a Type IIn: fast ejecta hit slow, dense circumstellar
material and a reverse shock stalls the shell. There
`sigma = d ln V / d ln r` drops below -1, which the stock code
assumes cannot happen, so the transfer inside the compression zone
returns negative or NaN intensities.

## What this fork adds

The non-monotonic transfer and JFNK solver paths are opt-in.
The observed-frame continuum rescue described below is enabled by default.

- **Knop formal solver** (Knop, Hauschildt & Baron 2009): sign-aware
  comoving-frame frequency coupling for velocity reversals. This route
  remains diagnostic; full-model convergence and physical correctness
  are not established. Lagged `I_NEXT`/`J_NEXT` maps are checkpointed
  across restarts and refreshed before each population update.
- **Non-monotonic Sobolev treatment** (Dessart, Audit & Hillier
  2015): the unclamped velocity gradient drives the line opacity in
  the compression zone and the optical-depth denominator is floored
  so it cannot cross zero.
- **Newton-Krylov population solve** (Knoll & Keyes 2004; Arramy et
  al. 2024): an optional GMRES population update with residual-based
  acceptance and bounded budgets. Admission requires the formal-J path
  (`USE_LAM_ES=T`); input checks reject `DO_JFNK` on the
  moment-store route. Enabling the flag does not guarantee admission.
  Frozen-J updates use stored radiative terms; coupled transfer requires
  separate controls. Full-model effectiveness remains unverified.
- **Trace-species handling**: demand floors plus a Stage-B polish
  phase stop physically empty levels from stalling the convergence
  exit.
- **Population audit diagnostics** (`POP_AUDIT` VADAT key): a
  per-iteration census of every species and ionisation stage, with
  conservation residuals and departure coefficients, for tracking
  down suspect ionisation balances.
- **Iteration diagnostics**: `ITER_HISTORY` plus a machine-readable
  `ITER_LEDGER` row per great iteration on the non-monotonic paths,
  recording solver call counts, refresh and stage timings,
  freeze/release events and memory use.
- **KNOP replay telemetry** (`KN_RDUMP_LS`, `KN_RDUMP_LS2`,
  `KN_RDUMP_LS3` VADAT keys, default off): writes a binary
  per-refresh deck for each armed ray (march coefficients, donor
  topology, both store columns, seeded and final fields, applied
  outcome) plus a per-pass applied-omega record, so a ray solve can
  be re-run offline under the same solver. The replay driver lives
  in `unit_tests/f90/drv_knop_replay.f90`; the audit trail is in
  `docs/audits/knop_phase2*`. Diagnostic only; the solve is
  unchanged.
- **ES accelerated Lambda iteration** (`DO_ES_ALI`, `DO_ES_ALI_BA`):
  an Olson-Auer-Buchler Jacobi step on the stored approximate-Lambda
  diagonal converges the coherent electron-scattering J loop in a
  few passes instead of the ~2 Lambda steps per iteration the stock
  update allows, with an optional matching factor for the continuum
  BA Jacobian. `DO_ES_ALI_BA` requires `DO_ES_ALI`; denominator
  floor activations are counted in `ES_ALI_NFLOOR` and
  `ES_ALI_BA_NFLOOR`.
- **ES redistribution convergence diagnostic** (`ES_RESID_DIAG` OBS
  key, default off): writes per-sweep residuals of the non-coherent
  electron-scattering field (`RJ_ES`, the `ES_J_CONV` generations)
  to `ES_RESID_DIAG`, plus a final render-generation audit.
  Diagnostic only; the sweep and render are unchanged.
- **Continuum J rescue** (`J_RESCUE`, `RESC_JMAX` OBS keys, on by
  default): on dense interacting-SN structures the stock continuum
  J/f iteration can ratchet to an unphysical amplitude in the
  pure-scattering infrared, where the formal-solve source is nearly
  homogeneous in J and the `FEDD` convergence test is
  amplitude-blind. A frequency that exits with a non-finite or
  over-cap field restores the previous frequency's converged
  field and the render always returns a bounded spectrum; rescued
  frequencies are logged. Explicit `J_RESCUE=F` restores stock
  behaviour. See `docs/audits/BLANK_RENDER_J_RESCUE.md`.
- **Hydro-structure velocity grids** (`VEL_LAW` 12, 13 and 14):
  `RV_SN_MODEL_SNHYD` builds the depth grid from `SN_HYDRO_DATA`,
  uniform in `log tau_es` (VEL_LAW=13) or in a blended
  tau-Rosseland and log-radius coordinate with structure-aware
  insertion at velocity reversals and density spikes (VEL_LAW=14).
  The `RV_SN_MODEL_SNIIN` reader was repaired so VEL_LAW=12 again
  runs on the standard seven-column `input_hydro.dat`.
- General correctness fixes in the population solve, the
  observed-frame tools and the build.

## Layout

| Path | Contents |
|---|---|
| `new_main/` | Main program and the VADAT input readers |
| `new_main/subs/` | New-solver support routines |
| `plane/` | The Knop formal solver |
| `subs/`, `newsubs/`, `com/`, `disp/` | Upstream source, surgical edits only |
| `obs/` | Observed-frame tools (`cmf_flux`, `obs_frame`) |
| `unit_tests/` | Test drivers; `run_all_tests.py` runs the suite |
| `docs/` | `audits/`, `plans/`, `figures/`, and the VADAT key reference `NONMONO_KEYS.md` |

## Building

```bash
module load intel mkl
./build.sh          # production build
./build_debug.sh    # bounds-checked debug build
```

VADAT input keys select the solver paths and diagnostics. Defaults can
also depend on the selected route; see `docs/NONMONO_KEYS.md` for details
and `AGENTS.md` for the development rules and verification checklist.

## Validation status

A clean build or a generated spectrum does not establish that the
new paths converge or produce a valid interacting-supernova
spectrum. Current runtime status is tracked in `docs/CHECKLIST.md`
and the audit notes under `docs/audits/`.

See [`docs/plans/NUMERICAL_METHODS_DEVELOPMENT_PLAN.md`](docs/plans/NUMERICAL_METHODS_DEVELOPMENT_PLAN.md)
for the staged investigation plan and acceptance tests.

## Attribution

All credit for CMFGEN itself goes to its original authors, led by
D. John Hillier (Hillier & Miller 1998, *ApJ* 496, 407,
[DOI](https://doi.org/10.1086/305350) |
[ADS](https://ui.adsabs.harvard.edu/abs/1998ApJ...496..407H)), and
to the authors of the methods this fork builds on (see References).
The maintainer takes no credit for any of the code in this repository.
All credit belongs to the original authors. Please cite the upstream
CMFGEN papers and the method papers listed below, not this fork.

## References

- Knop, Hauschildt & Baron 2009, *A&A* 496, 295,
  [DOI](https://doi.org/10.1051/0004-6361:200810290) |
  [arXiv:0902.1876](https://arxiv.org/abs/0902.1876)
- Knop, Hauschildt & Baron 2009, *A&A* 501, 813,
  [DOI](https://doi.org/10.1051/0004-6361/200911793)
- Dessart, Audit & Hillier 2015, *MNRAS* 449, 4304,
  [DOI](https://doi.org/10.1093/mnras/stv609) |
  [arXiv:1503.05463](https://arxiv.org/abs/1503.05463)
- Feldmeier & Nikutta 2006, *A&A* 446, 661,
  [DOI](https://doi.org/10.1051/0004-6361:20053380)
- Knoll & Keyes 2004, *J. Comput. Phys.* 193, 357,
  [DOI](https://doi.org/10.1016/j.jcp.2003.08.010)
- Arramy, de la Cruz Rodriguez & Leenaarts 2024, *A&A* 690, A12,
  [DOI](https://doi.org/10.1051/0004-6361/202449963) |
  [arXiv:2406.17234](https://arxiv.org/abs/2406.17234)
