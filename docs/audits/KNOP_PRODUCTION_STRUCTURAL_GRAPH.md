# KNOP Production Structural-Graph Audit (P1-T1f, PARTS 1-4)

Status: **BLOCKED on a production checkpoint.** PART 1 (checkpoint
survey) is complete and its conclusion is negative: the existing
`FGINEXT_GS` file cannot reconstruct a production fixed-opacity sweep
offline. PARTS 2-4 (graph construction, Tarjan SCC, targeted FD) are
specified below but **NOT executed** - no production nilpotency claim
is made or implied.

Scope: `plane/fg_j_cmf_nonmono_v1.f`. No production numerical
equation was modified. This report records the minimum data a
checkpoint must carry and the graph/SCC procedure to run once one
exists.

Evidence classes: **[MATH]** code structure, **[MANUF]** manufactured
`drv_knop_ray` model, **[PROD]** real non-monotonic model (not yet run).

The node/edge methodology and the full read/write trace of the
radiation stores are in `KNOP_PRODUCTION_DEPENDENCY_GRAPH.md` (P1-T1c);
this report does not repeat them.

---

## PART 1 - Checkpoint survey (COMPLETE)

### Existing checkpoint: `FGINEXT_GS`

Written by `DUMP_FG_INEXT_GS` (~line 5759), read back by
`RESTORE_FG_INEXT_GS`. Current format `VER=3`. It stores:

- magic/version/tag, `N_COMP`, `NP`, `NCF`, `ND`
- grid fingerprint `R`, `V`, `SIGMA`, `COMP_DEPTH_IDX`
- convergence/cap metrics
- frequency grid `KN_FREQ`
- the two lagged stores `I_P_DEPTH_NEXT`, `I_M_DEPTH_NEXT`
- end marker

It does **not** store the arrays the fixed-opacity operator is built
from. On restart those come from live model state, not the file.

### Minimum data to reconstruct a production fixed-opacity sweep

To rebuild the operator offline (and hence the store map `T`), a
diagnostic checkpoint must additionally carry, per the solver inputs
at `FG_KNOP_REFRESH_V1` / `KNOP_RAY_SOLVE_V1`:

| Array / value | Used for |
|---|---|
| `CHI`, `CHI_EXT` | absorption coefficient per (depth,freq) |
| `ETA`, `ETA_EXT` | emissivity per (depth,freq) |
| `GAM` | advection/velocity-gradient term (drives expanding vs compressing branch) |
| `DTAU` | radial optical-depth increment |
| `Q` | frequency-coupling/source coefficient (`SOURCE_PRIME = SOURCE_RAY + Q*I_*_PREV`) |
| `R_RAY`, `Z` | ray geometry / impact-parameter coordinate |
| `J_PNT`, `RAY_PNT` | ray-to-depth and point mapping |
| `COMP_DEPTH_IDX` | sparse compressing-depth index (already in FGINEXT_GS) |
| `NM_RAY`, `DLMN`, `NFM`, `NNF` | frequency matching / donor offsets |
| `I_P_DEPTH_NEXT`, `I_M_DEPTH_NEXT` | lagged store seed (already in FGINEXT_GS) |
| inner-boundary config (`KN_INNER_METH`, `KN_IC`, `IM_NI`) | `I_P(NI)=I_M(NI)` closure |
| `NP`, `NCF`, `N_COMP`, `ND` | dimensions (partially in FGINEXT_GS) |
| `FG_FREQ_MATCH`, `INEXT_RELAX`, `KN_RAY_*` tolerances | operator config |
| great-iteration identifier | provenance / repeatability |

`FGINEXT_GS` supplies only the store seed, the sparse-depth index,
the frequency grid, and the grid fingerprint. It is **insufficient**.

### Decision

Add a **diagnostic-only checkpoint writer** (e.g. `DUMP_KNOP_SWEEP_STATE`)
that serialises the full table above at a chosen great iteration,
gated on an explicit diagnostic flag and off by default. Normal
restart semantics (`FGINEXT_GS`) are untouched. Alternative: an
in-situ perturbation run that builds the graph from live arrays
without an intermediate file - cheaper to implement but cannot be
re-run or shared. The file route is preferred because PART 4's
targeted FD needs repeatable state restore.

---

## PART 2 - Structural dependency graph (SPECIFIED, NOT RUN)

Node set: `(CHANNEL, KC, LS, SLOT)` over active store elements;
CHANNEL in {P inward, M outward}; `KC` a compressing depth index;
`LS` the ray; `SLOT` the frequency-bin index within the matching
window. Both channels included.

Edges traced conservatively (every possible dependency, no numerical
threshold, no `GAM`-sign-only ordering) through:

1. compressing `I_NEXT` donor interpolation (`JPR1/JPR2`, `WINX`),
2. expanding `I_PREV` (`SOURCE_PRIME = SOURCE_RAY + Q*I_*_PREV`),
3. `SOURCE_PRIME` -> `S`, `dS` (`SLOPES`),
4. monotonicity clamps (`VNEW` floor/cap in `MARCH`),
5. inward + outward short-characteristic recurrence (`MARCH`),
6. inner boundary `I_P(NI)=I_M(NI)` (`KN_INNER_METH`),
7. frequency matching `NM_RAY`/`NFM`/`NNF`/`DLMN`,
8. final store writes (`WRITE_STORE`/`WRITE_STORE_CAP`),
9. `I_*_STORE -> I_*_PREV` copy on `NEW_FREQ` (the expanding channel
   confirmed `I_PREV`-mediated in `KNOP_PREV_ABLATION_AUDIT.md`).

Ray independence must be verified from the implementation before
per-ray analysis: the sweep state (`FLD`, coefficient pack `PK_*`,
`I_*_PREV`) is per-ray and the lagged stores are per-ray columns, so
fixed-opacity transport is expected ray-independent **[MATH, to be
confirmed against `FG_KNOP_REFRESH_V1`'s ray loop]**.

The conservative rule set is fixed by P1-T1c; this audit adds the
requirement that the graph be built from the **checkpoint's actual
geometry and branch conditions** (`HAS_COMPRESSING`, active
`KN_INNER_METH`, `FG_FREQ_MATCH`, ordered vs unordered sweep path),
not from the manufactured profile.

## PART 3 - Graph analysis (SPECIFIED, NOT RUN)

Tarjan SCC on the conservative graph. Report: node count, edge count,
self-loops, SCC count, non-trivial SCC count, largest SCC, maximum
dependency-chain length, topological order if acyclic, cross-channel
dependencies, same-slot dependencies, backward-frequency
dependencies. No numerical thresholding - the graph is structural.

## PART 4 - Targeted FD validation (SPECIFIED, NOT RUN)

For each candidate cycle edge: perturb the input store element, run
the complete ordered production sweep at fixed opacity, measure the
output store element, two central-difference amplitudes, restore all
state between evaluations, record repeatability and derivative
uncertainty. Validate at least one edge of each mechanism:
compressing interpolation, expanding `I_PREV`, parabolic source
stencil, inner-boundary M-to-P coupling, multi-bin frequency
matching. If the conservative graph is acyclic, nilpotency is
established **for that checkpoint only**. If a cycle exists, extract
the SCC submatrix and compute eigenvalues with LAPACK `dgeev`; report
spectral radius and dominant coupling mechanism.

---

## What is and is not established

- **Established [MANUF]:** the manufactured store map is nilpotent
  (P1-T1b); same-slot edges are `I_PREV`-mediated (P1-T1e ablation).
- **Established [MATH]:** `FGINEXT_GS` lacks the operator inputs; a
  diagnostic checkpoint writer is required.
- **NOT established [PROD]:** whether production geometry creates a
  directed cycle. No claim is made.

## Blocking item

A real non-monotonic production checkpoint carrying the PART 1 table.
Until it exists, PARTS 2-4 cannot run and production nilpotency
remains open.
