# CMFGEN with non-monotonic velocity fields

> **Warning: this code is under active development.** Solver
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

All new solver paths are opt-in. With their flags unset the code
follows the stock upstream path.

- **Knop formal solver** (Knop, Hauschildt & Baron 2009): sign-aware
  comoving-frame frequency coupling, so `sigma < -1` zones no longer
  produce negative intensities. Lagged `I_NEXT`/`J_NEXT` maps are
  checkpointed across restarts and refreshed before each population
  update.
- **Non-monotonic Sobolev treatment** (Dessart, Audit & Hillier
  2015): the unclamped velocity gradient drives the line opacity in
  the compression zone and the optical-depth denominator is floored
  so it cannot cross zero.
- **Newton-Krylov population solve** (Knoll & Keyes 2004; Arramy et
  al. 2024): an opt-in GMRES update for stalled models where the
  stock Newton step overshoots, with residual-based acceptance and
  bounded budgets.
- **Trace-species handling**: demand floors plus a Stage-B polish
  phase stop physically empty levels from stalling the convergence
  exit.
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
| `docs/` | README figures and audit notes |

## Building

```bash
module load intel mkl
./build.sh          # production build
./build_debug.sh    # bounds-checked debug build
```

All new behaviour is off by default and switched on through VADAT
input keys. See `AGENTS.md` for the development rules and
verification checklist.

## Validation status

A clean build or a generated spectrum does not establish that the
new paths converge or produce a valid interacting-supernova
spectrum. Current runtime status is tracked in `AGENTS_HISTORY.md`
and the audit notes under `docs/audits/`.

## Attribution

All credit for CMFGEN itself goes to its original authors, led by
D. John Hillier (Hillier & Miller 1998, *ApJ* 496, 407,
[DOI](https://doi.org/10.1086/305350) |
[ADS](https://ui.adsabs.harvard.edu/abs/1998ApJ...496..407H)), and
to the authors of the methods this fork builds on (see References).
This repository is a development fork; the maintainer takes no
credit for the original code. Please cite the upstream CMFGEN
papers and the method papers listed below, not this fork.

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
