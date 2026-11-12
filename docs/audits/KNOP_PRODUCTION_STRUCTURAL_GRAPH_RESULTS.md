# KNOP Production Structural Graph Results (P1-T1h)

Status: **PARTIAL - PART 1 complete, Layer-A structural graph built and
analysed on a manufactured checkpoint. PARTS 3-8 (Layer-B1/B2 graphs,
in-situ finite-difference probe, eigen analysis, full test matrix) are
NOT yet implemented.** No production equation, convergence criterion,
population equation, commit behaviour, or restart behaviour was
modified. All new functionality is diagnostic-only and disabled by
default.

Evidence classes: **[SRC]** source trace; **[RUN]** measured on the
`knray` unit-test driver (production `*.o`, `ifort`, `-extend-source`);
**[GRAPH]** computed by `unit_tests/python/knop_graph.py` from
`unit_tests/python/knop_topo.py`.

Scope: `plane/fg_j_cmf_nonmono_v1.f`, module
`FG_J_CMF_NONMONO_MOD_V1`, refresh path `FG_KNOP_REFRESH_V1` ->
`KNOP_RAY_SOLVE_V1` (OMP) -> `SWEEP_FLD`/`MARCH`/`WRITE_STORE`/
`WRITE_STORE_CAP`/`HEAL_STORE`.

Prior docs: `docs/audits/KNOP_CHECKPOINT_DESIGN_REVIEW.md` (the
approved two-tier plan this phase implements),
`docs/audits/KNOP_PRODUCTION_DEPENDENCY_GRAPH.md` (the Layer-A edge
rules), `docs/audits/KNOP_PRODUCTION_STRUCTURAL_GRAPH.md` (the
specified-but-blocked graph algorithm).

Scientific references: Knop et al. (2009), A&A 496, 295-298,
DOI 10.1051/0004-6361:200810290; Knop et al. (2009), A&A 501, 813-820,
DOI 10.1051/0004-6361/200911793; Trefethen & Bau (1997), Numerical
Linear Algebra.

---

## PART 1 - Compact topology checkpoint (COMPLETE)

A versioned, atomic, diagnostic-only topology checkpoint was added to
`FG_KNOP_REFRESH_V1`. It is gated on the new module logical
`KNOP_DUMP_TOPO` (default `.FALSE.`), declared last in the module
variable list so existing module variable offsets are preserved for
precompiled archive objects. The writer is `DUMP_KNOP_TOPO_V1`, invoked
immediately before the OpenMP ray loop (the serial point chosen in the
P1-T1g read-set audit), and writes through a `NEWUNIT` handle to
`KNOP_TOPO` via a `KNOP_TOPO.tmp` temporary plus `RENAMEFILEQQ`
atomic rename, matching the existing diagnostic-dump IOSTAT
conventions.

The checkpoint is explicitly a **topology** checkpoint, not a portable
offline transport checkpoint. It deliberately excludes the ~300 MB
`KN_CHI_COEF`/`KN_ETA_COEF`/`KN_CHI_EXT`/`KN_ETA_EXT` opacity packs:
the transport sparsity pattern is opacity-independent (verified -
`KN_RF_COUPLING` and `MATCH_FREQ` read only geometry/frequency, no
opacity), so the packs set Jacobian *weights*, not the nonzero
*pattern*. They remain required only for the later in-situ
finite-difference probe, which runs on the live frozen `KN_*` arrays
rather than a serialised copy.

### Checkpoint schema

Unformatted sequential Fortran records, little-endian:

- Record 1 (header): magic `20261120`, version `1`, kind tag `1`.
- Record 2 (kind params): precision-kind metadata.
- Record 3 (dims): `KN_ND`, `KN_NP`, `KN_NC`, `KN_NCF`, `KN_ND_EXT`,
  `N_COMP`, `NRAY_MAX` (written as `UBOUND(RAY_PNT,1)` so a stale
  module scalar cannot corrupt the bound).
- Record 4 (integer config `ICNF`): `KN_RAY_MAXPASS`,
  `KN_MAX_THREADS`, `KN_RAY_MAX_INS`, `FG_ORD_AA_M`, `FG_ORD_P2N`,
  `FG_ORD_P2REL`, `KN_RAY_HOTPASSES`, `KN_HOT_TOTAL_PASSES`,
  `FREQ_IDX_ITER`, `ND_ADD`, `N_TH`.
- Record 5 (real config `RCF`): `INEXT_RELAX`, `KN_RAY_TOL`,
  `KN_RAY_BULK_TOL`, `KN_RAY_OMEGA`, `KN_RAY_STALL`, `KN_RF_TOT_MB`,
  `KN_RAY_PACK_MB`, `FG_COMMIT_JUMP`, `FG_COMMIT_NC`, `FG_DIV_HEAL_CAP`,
  `NONMONO_INEXT_CAP`.
- Record 6 (logical config `LCFG`): `DO_KNOP_REFRESH`,
  `FREEZE_I_NEXT`, `FREEZE_I_PREV`, `HAS_COMPRESSING`, `KN_STORE_READY`,
  `KN_ADAPT_RELAX`, `FG_FREQ_MATCH`, acceleration switches, commit
  switches, seed state, plateau state, `PRINT_KNOP_DEBUG`.
- Record 7 (character config): method/boundary strings.
- Record 8 (array-presence masks): integer and real presence masks.
- Integer topology arrays, then real topology arrays (below).
- Store fingerprints: per channel/slot, sum, sum-abs, sum-squares,
  max-abs, nonzero count; `-1` sentinels for absent stores.
- End marker record.

Arrays serialised (present masks permit optional absence):
`COMP_DEPTH_IDX`, `J_PNT`, `RAY_PNT`, `NI_RAY`, optional `MAX_LS`,
optional `NM_RAY`, `GAM`, `GAM_MAG`, `dGAMdR`, `R_RAY`, `Z`, optional
`DLN_RAY`, `KN_FREQ`, `KN_DLN_P`, `KN_DLN_N`, `KN_DBB`, `KN_IC`,
`KN_R`, `KN_P`. Derived per-ray locals (`NFM_L`, `DLNM_L`, `DRF_L`,
`NNF`) are deliberately **not** serialised - they are rebuilt by
`MATCH_FREQ` inside the solve.

### Reader

`unit_tests/python/knop_topo.py` decodes the unformatted-sequential
framing manually (mixed-type records defeat `scipy.io.FortranFile`'s
homogeneous `read_record`), validates magic/version/dims/record
structure/array presence/fingerprints, and reports shapes, types,
extrema, configuration, and fingerprints. Optional arrays absent at
the checkpoint (`MAX_LS`, `NM_RAY`, `DLN_RAY` in the manufactured
case) are handled via the presence masks.

### Driver case and verification

Driver case `topo` (case 63) in `unit_tests/f90/drv_knop_ray.f90`
reuses the refresh setup with `KNOP_DUMP_TOPO=.TRUE.` and reports
`TOPO_FILE=T` when `unit_tests/work/KNOP_TOPO` exists. Verified
**[RUN]**: `RF_FAIL=F`, `TOPO_FILE=T`, `RESULT_OK`; the manufactured
checkpoint is 2116 bytes. Full `test_knop_ray.py` regression:
**43 passed, 0 failed** after the change.

### Manufactured checkpoint contents (the analysed state)

`unit_tests/work/KNOP_TOPO` from the `topo` case:

- Dims: `KN_ND=8`, `KN_NP=3`, `KN_NCF=6`, `N_COMP=3`, `NRAY_MAX=8`.
- `icnf`: `KN_RAY_MAXPASS=14`, `KN_MAX_THREADS=16`, `N_TH=16`,
  `ND_ADD=0`, `KN_RAY_HOTPASSES=0`, `KN_HOT_TOTAL_PASSES=0`,
  `FREQ_IDX_ITER=5`, `FG_ORD_AA_M=2`, `FG_ORD_P2N=2`,
  `FG_ORD_P2REL=3`.
- `rcf`: `INEXT_RELAX=0.5`, `KN_RAY_OMEGA=1.0`, `KN_RAY_STALL=1.1`,
  `FG_COMMIT_JUMP=50.0`, `FG_COMMIT_NC=0.05`,
  `FG_DIV_HEAL_CAP=1e6`, `NONMONO_INEXT_CAP=1e20`.
- `lcfg` true: `DO_KNOP_REFRESH`, `HAS_COMPRESSING`, `KN_STORE_READY`,
  `FG_COMMIT_CONV`, `FG_SEED_COMMIT`. `FG_FREQ_MATCH` is **false** in
  this checkpoint, so `NM_RAY=1` is implied even though `NM_RAY` is
  not serialised; `FREEZE_I_NEXT`/`FREEZE_I_PREV` are false.
- Present arrays: `COMP_DEPTH_IDX`, `GAM`, `GAM_MAG`, `J_PNT`,
  `KN_DBB`, `KN_DLN_N`, `KN_DLN_P`, `KN_FREQ`, `KN_IC`, `KN_P`,
  `KN_R`, `NI_RAY`, `RAY_PNT`, `R_RAY`, `Z`, `dGAMdR`. Absent:
  `MAX_LS`, `NM_RAY`, `DLN_RAY` (per-frequency derived work arrays).

**This is a manufactured checkpoint, not a production checkpoint.** No
production `KNOP_TOPO` has been captured; every result below is for
the manufactured `topo` driver state only.

---

## PART 2 - Layer-A transport graph (built, conservative)

`unit_tests/python/knop_graph.py` builds the directed transport graph
from the checkpoint. Vertex = a writable radiation-store element
`(channel, KC, LS, frequency_slot)` where `channel` is 0 for `I_M`
(inward) and 1 for `I_P` (outward), `KC` indexes the compressing depth
`COMP_DEPTH_IDX`, `LS` indexes the ray pencil `1..KN_NP`, and
`frequency_slot` is the store slot `s` corresponding to frequency
`s+1`. Phantom/unwritable elements (non-compressing depths, out-of-range
slots, invalid `J_PNT`) are excluded. On the manufactured checkpoint
this gives **90 vertices** = 2 channels x 3 rays x 3 compressing depths
x 5 writable slots (`NCF-1 = 5`).

Edges are built from the source-level dependency rules in
`KNOP_PRODUCTION_DEPENDENCY_GRAPH.md` section 2.1, traced against
`KNOP_RAY_SOLVE_V1`/`SWEEP_FLD`/`MARCH`/`WRITE_STORE`:

- **Compressing-cell donor read.** For a compressing ray-point `I`,
  `K=RAY_PNT(I,LS)-ND_ADD`, `KC1=COMP_DEPTH_IDX(K)`,
  `KC2=COMP_DEPTH_IDX(K+1)`; the donor store slot is
  `FREQ_IDX_ITER+NM_RAY(I)` interpolated in `R_RAY`/`KN_R` between
  `KC1`/`KC2`. With `FG_FREQ_MATCH=F`, `NM_RAY=1`. Both bracketing
  depths contribute when valid; a single valid depth contributes alone;
  neither contributes nothing; out-of-bounds slots drop the coupling;
  `GAM=0` reads no donor.
- **Spatial stencil.** `I_M` marches with increasing ray index
  (`I_M(I+1)=I_M(I)A0+SOURCE_PRIME(I)A1+SOURCE_PRIME(I+1)A2`);
  `I_P` marches with decreasing ray index
  (`I_P(I)=I_P(I+1)A0+SOURCE_PRIME(I+1)A1+SOURCE_PRIME(I)A2`). A
  compressing store point therefore depends on `SOURCE_PRIME` over its
  upwind stencil, propagated through `J_PNT`/`WINX` interpolation.
- **Expanding-cell `I_PREV` chain.** `SOURCE_PRIME(I)=SOURCE_RAY(I)+
  Q(I)*I_*_PREV(I,LS)` for `GAM>0`; `I_*_PREV` is the same frequency's
  previous-sweep intensity, i.e. the same ray-point at frequency `s`,
  which itself writes store slot `s-1`. This produces backward and
  same-slot frequency edges; it is **not** inferred from `GAM` sign
  alone - the exact donor indexing and active branch are traced.
- **Cross-channel boundary.** `I_P(NI)=I_M(NI)` for core rays couples
  the two channels only where the implementation genuinely does so;
  under `KN_INNER_METH='ZERO_FLUX'` the inner boundary uses `I_CORE`,
  so no M<->P edge is generated for this checkpoint.
- **Per-ray independence.** Rays are separate OMP work items; no
  transport edge crosses `LS`.

The builder retains the **conservative structural graph** (all
possible dependencies) and is designed to be intersected later with
the **numerically validated graph** (FD-confirmed edges). It does not
use a single arbitrary numerical threshold to define structure.

### Manufactured Layer-A result [GRAPH]

```
vertices = 90
edges    = 789
self-loops = 75
weak components = 6, each of size 15
SCCs = 9 total, 6 nontrivial, largest = 15
is_dag = False
channel pairs present: (0,0):366, (1,1):423   (no cross-channel edges)
slot diffs (dst_slot - src_slot):
  -1:204   (forward compressing donor reads, input slot = output+1)
   0:195   (same-slot, expanding I_PREV reach)
  +1:156, +2:117, +3:78, +4:39  (backward expanding-chain reach)
```

The six weak components are exactly the six (channel x ray) pairs -
transport keeps rays and, under `ZERO_FLUX`, channels disjoint. Each
15-vertex block is internally one nontrivial SCC in the **conservative**
graph: the expanding `I_PREV` chain reaches back across frequency
slots, closing cycles that the earlier manufactured `p1t1`/`p1t1b`
Jacobians showed carry ~zero weight. This is the expected conservative
superset; whether any of these cycles carry nonzero Jacobian weight is
a Layer-A numerical question that the in-situ FD probe (PART 5) is
specified to answer. **The conservative graph is cyclic; it does not
by itself establish or refute nilpotency.**

---

## PARTS 3-8 - Status

- **PART 3 (Layer-B1 iteration/acceleration dependencies):** NOT
  implemented. The per-ray dispatch path (`ORD_SWEEP`, `HAVE_PACK`,
  `NI`, `KN_RAY_PACK_MB`, `FG_ORD_SWEEP`, `KN_RAY_OMEGA`,
  `INEXT_RELAX`, `KN_ADAPT_RELAX`, `FG_ORD_AA`, `FG_ORD_AITKEN`,
  `FG_ORD_P2DAMP`) is configured in the checkpoint (`icnf`/`rcf`/
  `lcfg`) but the per-ray dispatch outcome and the relaxed/accelerated
  inner-iteration dependency recording are not yet built.
- **PART 4 (Layer-B2 conditional commit graph):** NOT implemented.
  `KN_JMP_SIBREF_V1` sibling medians, `FG_JUMP_MAX` veto,
  `FG_RAY_SEEDED` first-commit waiver, `FG_COMMIT_CONV` verdicts,
  `WRITE_STORE_CAP` caps, `HEAL_STORE` healing, and seed dependence
  are identified as conditional branches in the design review but no
  commit-graph edges have been recorded.
- **PART 5 (in-situ central-FD probe):** NOT implemented. Requires a
  Fortran driver that snapshots both stores + `FG_RAY_SEEDED` + mutable
  state, runs the baseline map, restores, applies +/- perturbations on
  the live frozen `KN_*`, recomputes sibling references, and measures
  the central difference at two amplitudes under `KN_MAX_THREADS=1`.
  `DRV_P1T1B` already contains reusable snapshot/restore + two-amplitude
  FD + `dgeev` machinery for the manufactured store map.
- **PART 6 (Tarjan SCC + FD validation + eigen):** Tarjan SCC is
  implemented in `knop_graph.py` (`graph_stats`). FD edge validation
  and the SCC-Jacobian eigen/spectral-radius analysis are blocked on
  PART 5.
- **PART 7 (separate map-A/map-B1/map-B2 reporting):** NOT implemented;
  depends on PARTS 3-5.
- **PART 8 (manufactured test matrix + 41-vs-43 count):** PARTIAL. The
  `topo` case (63) and `repl`/`replo`/`replc`/`replh` cases exist;
  `test_knop_ray.py` reports 43 passed / 0 failed. The full listed
  matrix (alternating bands, multi-bin freq match, both channels,
  cross-channel BC, sibling-median perturbation, jump-veto activation,
  capped commit, healed store, first-commit vs seeded, ordered vs
  unordered, per-ray pack fallback, two FD amplitudes, `N_TH=1`
  repeatability) is not yet built. The 41-vs-43 discrepancy is resolved:
  the design review reported 41/41 before the `repl*` cases were added;
  the current suite is 43/43 including T38/T38b.

---

## What is established vs unproven

**Established [RUN]/[GRAPH] on the manufactured `topo` checkpoint only:**

- A compact, versioned, atomic, gated topology checkpoint exists and is
  readable; it carries dims, integer/real/logical/character config,
  `N_TH`, array-presence masks, topology arrays, and store fingerprints,
  and excludes the opacity packs by design.
- The Layer-A conservative transport graph builds correctly: 90
  vertices, per-(channel x ray) disjoint components, correct writable
  slot set, and dependency edges following the traced source rules
  (compressing donor slot `FREQ_IDX_ITER+NM_RAY`, `J_PNT`/`WINX`
  interpolation, expanding `I_PREV` backward reach, `ZERO_FLUX`
  boundary giving no cross-channel edge here).
- The conservative graph is **cyclic** (6 nontrivial SCCs, largest 15,
  not a DAG) because the expanding `I_PREV` chain can close frequency
  cycles. This is a statement about *possible* dependencies, not about
  Jacobian weights.

**NOT established:**

- **No production checkpoint was used.** Every result is for the
  manufactured `topo` driver state. Nothing here is a production result.
- **Nilpotency is NOT established** for the selected state or any
  production state. The conservative graph is cyclic; only an
  FD-validated (weight-thresholded) graph plus eigen analysis could
  speak to nilpotency, and that probe is not yet implemented.
- The FD-validated edge set, the numerical Jacobian, spectral radius,
  and eigenvalue sensitivity are all unmeasured.
- Layer-B1 (relaxed/accelerated inner iteration) and Layer-B2
  (conditional commit wrapper: sibling medians, jump veto, caps,
  healing, seed state) graphs are unbuilt.
- Whether the conditional commit wrapper is locally smooth or
  nonsmooth at a production checkpoint is unknown.
- Population-coupled CMFGEN convergence, refresh fixed-point
  convergence, and the actual numerical iteration are distinct claims
  and none are addressed here.

**Remaining unproven for production:** transport-graph acyclicity or
nilpotency on real geometry, the weight of the conservative cycle
edges, the commit-wrapper's active-set structure, and any
production-representative spectral result. A single checkpoint -
manufactured or production - could not establish universal production
nilpotency regardless.

---

## Files

- `plane/fg_j_cmf_nonmono_v1.f`: `KNOP_DUMP_TOPO` gate +
  `DUMP_KNOP_TOPO_V1` + call site before the OMP ray loop.
- `unit_tests/f90/drv_knop_ray.f90`: `topo` case (63).
- `unit_tests/python/knop_topo.py`: checkpoint reader/validator.
- `unit_tests/python/knop_graph.py`: Layer-A graph builder + Tarjan
  SCC (`graph_stats`).
- Commits: `90d41165` (writer), `b63345c3` (driver case),
  `355d2f4b` (reader), `0dd58b02` (graph builder).
