# KNOP Checkpoint Design Review (P1-T1g)

Status: **DESIGN REVIEW COMPLETE.** No serializer, reader, or replay
infrastructure was implemented in this phase. This report verifies the
proposed checkpoint contents against source, resolves the `DLMN`
question, characterises OMP reproducibility, runs a repeatability
experiment, and recommends a two-tier architecture. No production
equation, convergence criterion, population equation, or restart
behaviour was modified.

Evidence classes: **[MATH]** code structure only; **[RUN]** measured on
the `knray` unit-test driver (production `*.o`, `ifort`).

Scope: `plane/fg_j_cmf_nonmono_v1.f`, module `FG_J_CMF_NONMONO_MOD_V1`,
refresh path `FG_KNOP_REFRESH_V1` -> `KNOP_RAY_SOLVE_V1` (OMP) ->
`SEED_FLD`/`SWEEP_FLD`/`SWEEP_FLD_ORD`/`BUILD_FREQ`/`MATCH_FREQ`/
`SRC_PRIME`/`SLOPES`/`MARCH`/`WRITE_STORE`/`WRITE_STORE_CAP`/
`HEAL_STORE`.

Prior audit: `docs/audits/KNOP_PRODUCTION_CHECKPOINT_AUDIT.md`
(commit `14ebd1de`). This review supersedes its read-set where the
new evidence differs (notably: `DLMN` does not exist; `SNP_P`/`SNP_M`
and the in-sweep `I_*_PREV`/`I_*_STORE`/`A0`-`A4` arrays are **not** in
the refresh read set).

---

## PART 1 - Two distinct operators

The refresh path implements **two different mathematical objects** that
must not be conflated. They do **not** have identical Jacobians.

### A. Fixed-opacity transport operator

Per-ray map from the lagged input stores to the pre-commit transported
fields:

```
FLD,  FLD_P  =  T_ray(  I_P_DEPTH_NEXT, I_M_DEPTH_NEXT ;  frozen KN_* )
```

- Frozen opacities (`KN_CHI_COEF`/`KN_ETA_COEF`/`KN_CHI_EXT`/
  `KN_ETA_EXT`), fixed geometry (`GAM`/`GAM_MAG`/`dGAMdR`/`R_RAY`/`Z`),
  fixed frequency matching (`KN_FREQ`/`KN_DLN_*`/`FG_FREQ_MATCH`),
  fixed boundary (`KN_R`/`KN_P`/method flags).
- `FLD`/`FLD_P` are per-ray locals inside `KNOP_RAY_SOLVE_V1`
  (allocated there, freed on return). They are **not externally
  observable** without production instrumentation.
- This is the object whose sparsity/Jacobian the structural graph is
  meant to analyse: it is (piecewise) smooth in the store inputs.

### B. Complete implemented refresh map

```
I_*_DEPTH_NEXT(next) = R( I_*_DEPTH_NEXT(cur) )
```

`FG_KNOP_REFRESH_V1` adds, on top of `T_ray`:

- adaptive inner iteration and per-ray convergence exits (`IXIT`),
- the shared hot-pass budget `KN_HOT_SPENT`/`KN_HOT_TOTAL_PSES`,
- ordered vs unordered sweep selection (`ORD_SWEEP`),
- relaxation (`INEXT_RELAX`, `KN_ADAPT_RELAX`),
- acceleration/damping (`FG_ORD_AA`, `FG_ORD_AITKEN`, `FG_ORD_P2*`),
- commit jump vetoes (`FG_COMMIT_JUMP`, `FG_JUMP_MAX`, sibling refs
  `JMP_SIB_M/P`), convergence vetoes (`FG_COMMIT_CONV`),
- capped commits (`WRITE_STORE_CAP`, `NONMONO_INEXT_CAP`),
- divergence healing (`HEAL_STORE`, `FG_DIV_HEAL_CAP`),
- store write-back and telemetry (`FG_RAY_VC`, `KN_RF_RESID`, ...).

`R` is the map actually iterated to a fixed point and the map whose
**committed stores** are the observable output.

### Which object does the structural graph analyse?

The production structural dependency graph / SCC / nilpotency
investigation targets the **fixed-point operator `R` (map B)**, because
that is what is iterated. However:

- The *transport edges* (which store cells a given cell reads during
  the formal solution) come from map A and are opacity-independent in
  pattern (see PART 6).
- Map B adds **data-dependent edges** through the commit vetoes and
  healing: whether a given store element influences another's committed
  value depends on runtime magnitudes (jump thresholds, sibling medians,
  caps), not only on topology. These must be represented as
  conditional/active-set edges, not unconditional arcs.

`FLD`/`FLD_P` were **not** instrumented: they are per-ray locals with
no external hook, and capturing them would require modifying
`KNOP_RAY_SOLVE_V1`. For the repeatability question (PART 5) the
committed store `I_*_DEPTH_NEXT` is the correct observable of map B;
the pre-commit `FLD` distinction is documented here for the graph work
but deferred as instrumentation.

---

## PART 2 - Module-state read/write audit

All identifiers below are module-associated in
`FG_J_CMF_NONMONO_MOD_V1` unless noted. Line numbers are
`plane/fg_j_cmf_nonmono_v1.f` at the audited revision.

### 2.1 Entry conditions (gates, not data)

Read at `FG_KNOP_REFRESH_V1` entry (6811-6842); the routine returns
silently unless all hold:

| Variable | Decl | Read at | Class |
|---|---|---|---|
| `DO_KNOP_REFRESH` | `LOGICAL` (871) | 6811 | config gate |
| `FREEZE_I_NEXT` | logical | 6812 | config gate |
| `HAS_COMPRESSING` | `LOGICAL,SAVE` (312) | 6813 | persistent state (set at alloc, 2566) |
| `I_P_DEPTH_NEXT` allocated? | `REAL(LDP),ALLOCATABLE` (251) | 6814 | state presence |
| `KN_STORE_READY` | `LOGICAL` (959) | 6815 | persistent state (armed by sweep) |
| `FREQ_IDX_ITER` | integer | 6826 (`== KN_NCF-1`) | sweep-position state |
| `SOLUTION_METHOD` | `CHARACTER(10)` (769) | 6835 (`=='INTEGRAL'`) | config |

`KN_RF_WARNED` (logical) is a one-shot warning latch - telemetry, not
operator state.

### 2.2 Physical / frozen-opacity inputs

| Variable | Decl (line) | Shape | Role |
|---|---|---|---|
| `KN_CHI_EXT` | `REAL` (889) | `(KN_ND_EXT,KN_NCF)` | frozen opacity, extended grid |
| `KN_ETA_EXT` | `REAL` (890) | `(KN_ND_EXT,KN_NCF)` | frozen emissivity, extended grid |
| `KN_CHI_COEF` | `REAL(LDP)` (901) | `(KN_ND_EXT,4,KN_NCF)` | cubic log-coeffs, rebuilt 6854 |
| `KN_ETA_COEF` | `REAL(LDP)` (902) | `(KN_ND_EXT,4,KN_NCF)` | cubic log-coeffs, rebuilt 6859 |
| `KN_FREQ` | `REAL(LDP)` (891) | `(KN_NCF)` | frequency grid |
| `KN_DLN_P` | `REAL(LDP)` (892) | `(KN_NCF)` | backward log-span |
| `KN_DLN_N` | `REAL(LDP)` (893) | `(KN_NCF)` | forward log-span |
| `KN_DBB` | `REAL(LDP)` (894) | `(KN_NCF)` | bandwidth/spacing aux |
| `KN_IC` | `REAL(LDP)` (895) | `(KN_NCF)` | intensity-coupling aux |

`KN_CHI_COEF`/`KN_ETA_COEF` are **rebuilt every refresh** from
`KN_CHI_EXT`/`KN_ETA_EXT` via `MON_INT_FUNS_V2` (6850-6860), so the
independent inputs are the `*_EXT` arrays; the `*_COEF` arrays are
derived work arrays (rebuild, do not serialise as independent state).
Note `KN_CHI_EXT`/`KN_ETA_EXT` are **single-precision `REAL`**, the
`*_COEF`/grid arrays are `LDP`.

### 2.3 Geometry / ray mapping

| Variable | Decl (line) | Shape | Notes |
|---|---|---|---|
| `R_RAY` | real array | `(?,?)` | radial coord per ray-point |
| `Z` | real array | `(?,?)` | depth coord per ray-point |
| `GAM` | `REAL(LDP)` | `(?,KN_NP)` | velocity gradient - **read AND written** (sanitise) |
| `GAM_MAG` | `REAL(LDP)` | `(?,KN_NP)` | |GAM| - **read AND written** |
| `dGAMdR` | `REAL(LDP)` | `(?,KN_NP)` | dGAM/dr - **read AND written** |
| `RAY_PNT` | `INTEGER,ALLOCATABLE` (55) | `(:,:)` | point->ray map |
| `J_PNT` | `INTEGER,ALLOCATABLE` (38) | `(:,:)` | ray->point map |
| `NI_RAY` | `INTEGER,ALLOCATABLE` (57) | `(:)` | per-ray active depth count |
| `MAX_LS` | `INTEGER,ALLOCATABLE` (58) | `(:)` | per-ray max slot |
| `NM_RAY` | `INTEGER,ALLOCATABLE` (863) | `(:)` | matching structure |
| `KN_R` | real | `(KN_ND)` | radial grid (alloc 2661) |
| `KN_P` | real | `(KN_NP)` | impact-parameter grid (alloc 2661) |

**Mutable error-handling state:** `GAM`, `GAM_MAG`, `dGAMdR` are written
in place by the non-finite sanitisation branches in `BUILD_FREQ`
(9573-9584 and 9691-9702): if `ABS(x) .LE. HUGE(x)` fails (NaN/Inf), the
module element is zeroed and the per-ray `DIVERGED` flag is set. Each
ray writes only its own `LS` column (thread-safe), but the mutation is
**sticky across refresh calls**. On finite data it is a no-op (the
repeatability fingerprint in PART 5 confirms inputs are unchanged). A
checkpoint must treat these as potentially-mutated geometry, not
immutable physical input.

### 2.4 Radiation stores (the state vector)

| Variable | Decl (line) | Shape | Role |
|---|---|---|---|
| `I_P_DEPTH_NEXT` | `REAL(LDP),ALLOCATABLE` (251) | `(N_COMP,KN_NP,KN_NCF)` | lagged P store - input AND output |
| `I_M_DEPTH_NEXT` | `REAL(LDP),ALLOCATABLE` (252) | `(N_COMP,KN_NP,KN_NCF)` | lagged M store - input AND output |
| `COMP_DEPTH_IDX` | `INTEGER,ALLOCATABLE` (274) | `(KN_ND)` | depth->compressing-slot map (0=skip) |

Store slot `s` holds frequency `s+1`; slot `KN_NCF` is never written
(freq-1 has no retained slot). Writable elements are
`(KC=1..N_COMP, LS=1..KN_NP, L=1..KN_NCF-1)`.

**Not in the refresh read set** (verified: zero hits in 6771-10400):
`I_P_PREV`, `I_P_STORE`, `I_M_PREV`, `I_M_STORE`, `SNP_P`, `SNP_M`,
`SOURCE_PRIME`, `A0`-`A4`, `S`, `dS`, `EE`, `E0`-`E3`,
`ALO_*_RAY`, `FG_NOSLOT_RAY`. These are in-sweep (`FG_J_CMF_NONMONO_V1`)
arrays only.

### 2.5 Persistent algorithm state

| Variable | Decl | Read | Write | Class |
|---|---|---|---|---|
| `FG_RAY_SEEDED` | `INTEGER,ALLOCATABLE,SAVE` (586) | 7744 (`L_UNSEED`) | 7817 (post-commit) | persistent state |
| `KN_HOT_SPENT` | `INTEGER*8,SAVE` (619) | 8013,8669,10359 | 8012,8668 (CRITICAL), reset 6882 | shared mutable state |
| `KN_HOT_TOTAL_PSES` | integer | 8013,8669,10358 | `SET_FG_HOTPASS` (6733) | config |
| `KN_RAY_HOTPASSES` | integer (617) | 8015,8669,10357 | `SET_FG_HOTPASS` | config |
| `JMP_SIB_M`,`JMP_SIB_P` | allocatable | `FG_JUMP_MAX` | rebuilt 6940-6941 | derived work array |

`FG_RAY_SEEDED(LS)` is the "ray has committed a full store" flag: read
by the first-commit jump-veto waiver (7743-7744) and set to 1 after a
successful commit (7817). It is genuine persistent algorithm state and
**must** be in any complete snapshot; omitting it changes the veto
behaviour on the first commit after restore.

`JMP_SIB_M/P` are allocated and rebuilt each refresh by
`KN_JMP_SIBREF_V1(I_*_DEPTH_NEXT, ...)` (6940-6941) from the **input**
stores, then read by `FG_JUMP_MAX`. They are deterministic functions of
the input stores - derived work arrays, not independent state (but see
PART 2.6 for the cross-ray consequence).

### 2.6 Cross-ray dependency (sibling reference)

`KN_JMP_SIBREF_V1` (1175-1244): for each frequency slot `LS2=2..KN_NCF`,
`REF(LS2)` = **median over ALL rays `RS=1..KN_NP`** of that ray's peak
`|STORE(KC2,RS,LS2-1)|` over compressing depths `KC2` (finite,
`>1e-28`, `<=HUGE` only). `REF(1)=-1` sentinel.

Consequence: within one refresh call `JMP_SIB` is fixed (built before
the OMP loop), but ray `LS`'s **commit verdict depends on every other
ray's input store** through the sibling-median jump threshold. Across
refresh iterations this is a genuine **all-to-all cross-ray edge** in map
B. The per-ray transport (map A) has no such edge. The structural graph
must include these sibling edges; they are unconditional in structure
(the median is always computed) but the *magnitude* of the veto effect
is data-dependent.

### 2.7 Ordered-sweep dispatch (per-ray, not global)

- `FG_ORD_SWEEP` (`LOGICAL`, 1049, default `.TRUE.` 1133) is the VADAT
  `[KN_RAY_ORDERED]` setter target (`SET_FG_ORD_SWEEP`, 6645).
- `KN_RAY_PACK_MB` (`REAL(LDP)`, 1027, default 2048.0, 1145) is the
  per-thread pack memory cap.
- Inside `KNOP_RAY_SOLVE_V1`: `HAVE_PACK` (local logical, 7446) is set
  at 7579-7585 - the pack is allocated only if
  `64*NI*NCF <= 1e6*KN_RAY_PACK_MB`.
- `ORD_SWEEP = FG_ORD_SWEEP .AND. HAVE_PACK .AND. NI .GE. 3` (7632).

`HAVE_PACK`/`ORD_SWEEP` are **per-ray locals** (depend on that ray's
`NI`), so ordered vs unordered is decided **per ray**, and a ray can
silently fall back to unordered if its pack does not fit or `NI<3`.
`KN_RAY_ORDERED` (VADAT) -> `FG_ORD_SWEEP` is only the first of three
gates.

### 2.8 Acceleration / damping controls (config)

`FG_ORD_AA`, `FG_ORD_AA_M`, `FG_ORD_AITKEN`, `FG_ORD_P2DAMP`,
`FG_ORD_P2N`, `FG_ORD_P2REL`, `KN_RAY_OMEGA`, `INEXT_RELAX`,
`KN_ADAPT_RELAX`, `KN_INNER_METH`, `KN_METHOD_SAV`,
`KN_INSERT_SAV`, `FG_FREQ_MATCH`, `FG_COMMIT_JUMP`,
`FG_COMMIT_CONV`, `FG_COMMIT_NC`, `FG_SEED_COMMIT`, `FG_DIV_HEAL_CAP`,
`NONMONO_INEXT_CAP`, `FG_BULK_CONTR`, `FG_PLATEAU_DET`,
`KN_RAY_TOL`, `KN_RAY_BULK_TOL`, `KN_RAY_MAXPASS`, `KN_RAY_HOTPASSES`,
`KN_RAY_STALL`, `KN_HOT_TOTAL_PSES`, `KN_MAX_THREADS`, `KN_RF_TOT_MB`,
`KN_RAY_PACK_MB`.

These are algorithm **configuration** (set once, read-only during the
refresh). They select branches but are not state. Any replay must hold
them fixed (the work order already requires this).

### 2.9 Output / telemetry (written, not read as input)

`FG_RAY_VC` (per-ray commit verdict), `KN_RES_PASS`, `KN_HOT_CNT`,
`KN_RF_RESID`, `FG_RF_DCHG`, `FG_RF_PEND`, `FG_RF_FAILED`,
`FG_DIV_HEAL`, `FG_COMMIT_REJ`, `FG_NC_REJ`, `FG_NC_CAP`,
`FG_RAY_DCH`, `FG_RAY_WI`, `FG_RAY_WL`, `FG_RAY_WS`, `KN_RAY_ORDERED`
(reporting). Written under `!$OMP CRITICAL` or per-ray slots; not
operator inputs.

### 2.10 Classification summary

- **Physical input:** `KN_CHI_EXT`, `KN_ETA_EXT`, `KN_FREQ`,
  `KN_DLN_P`, `KN_DLN_N`, `KN_DBB`, `KN_IC`, `R_RAY`, `Z`, `GAM`,
  `GAM_MAG`, `dGAMdR`, `RAY_PNT`, `J_PNT`, `NI_RAY`, `MAX_LS`,
  `NM_RAY`, `KN_R`, `KN_P`, `COMP_DEPTH_IDX`.
- **Algorithm configuration:** all `FG_*`/`KN_*` scalars in 2.8 plus
  entry gates in 2.1.
- **Persistent algorithm state:** `FG_RAY_SEEDED`, `KN_STORE_READY`,
  `HAS_COMPRESSING`, `FREQ_IDX_ITER`, `KN_HOT_SPENT` (reset per call).
- **Rebuilt local work array:** `KN_CHI_COEF`, `KN_ETA_COEF`,
  `JMP_SIB_M`, `JMP_SIB_P`, and all per-ray locals (`IKC*`, `JPR*`,
  `WINX`, `IGRD`, `IM_NI`, `NFM_L`, `DLNM_L`, `DRF_L`, `FLD`, `FLD_P`,
  `PK_*`, `HAVE_PACK`, `ORD_SWEEP`).
- **Output/telemetry:** section 2.9.
- **Mutable error-handling state:** `GAM`, `GAM_MAG`, `dGAMdR`
  (sanitisation writes, 2.3).

---

## PART 3 - `DLMN` resolution

**`DLMN` does not exist.** `grep -c '\bDLMN\b'` over
`plane/fg_j_cmf_nonmono_v1.f` returns 0; there is no declaration and no
read/write site. The prior audit's "unresolved `DLMN` shape" was a
naming artefact.

The actual frequency-matching quantities, all **per-ray locals** inside
`KNOP_RAY_SOLVE_V1` (allocated there, rebuilt each call):

| Name | Kind | Built by | Depends on |
|---|---|---|---|
| `DRF_L` | `REAL(LDP)(NI)` | ray setup (~7563) | `ABS(GAM)*` grid spacing (`Z`) |
| `NFM_L` | `INTEGER(NI)` | `MATCH_FREQ` (9893+) | `DRF_L`, `KN_DLN_N`, `KN_FREQ`, `FG_FREQ_MATCH`, `NCF` |
| `DLNM_L` | `REAL(LDP)(NI)` | `MATCH_FREQ` | `KN_DLN_P`/`KN_DLN_N`, `KN_FREQ`, `NFM_L` |
| `NNF` | `INTEGER` (local scalar) | `MATCH_FREQ` (9891) | `DRF_L/KN_DLN_N`, capped `NCF-L` |

`MATCH_FREQ(L)` (9881+): defaults `NFM_L=1`, `DLNM_L=KN_DLN_P(L)`;
for compressing cells (`GAM<0`) with `FG_FREQ_MATCH` on, computes
`NNF=CEILING(DRF_L/KN_DLN_N)` capped at `NCF-L`, refines `DLNM_L` from
the `KN_FREQ` grid. With the flag off, compressing cells keep offset 1
and `KN_DLN_N(L)`.

**Independent inputs vs derived:** `KN_FREQ`, `KN_DLN_P`, `KN_DLN_N`,
`KN_DBB`, `KN_IC`, `FG_FREQ_MATCH`, and the geometry feeding `DRF_L`
are independent inputs. `NFM_L`, `DLNM_L`, `DRF_L`, `NNF`, `IKC*`,
`JPR*`, `WINX`, `IGRD`, `IM_NI`, `PK_*` are **deterministically
rebuilt** inside the ray solve and must **not** be serialised as
independent state. `DLMN` should be dropped from the inventory.

---

## PART 4 - OMP reproducibility

### Mechanism (source-confirmed)

- Ray loop: `!$OMP PARALLEL DO NUM_THREADS(N_TH) SCHEDULE(DYNAMIC)`
  over `LS=1..KN_NP`; `N_TH=min(OMP_GET_MAX_THREADS(), KN_MAX_THREADS)`
  further bounded by `KN_RF_TOT_MB`/`MAXVAL(NI_RAY)`.
- Hot-pass budget: `TRY_HOT()` (10354) grants a hot pass only while
  `KN_HOT_SPENT < KN_HOT_TOTAL_PSES` (read **outside** any lock, 10359).
  Each hot pass then does `KN_HOT_SPENT=KN_HOT_SPENT+1` inside
  `!$OMP CRITICAL (KN_HOT_SPD)` (8011-8014 unordered, 8667-8670
  ordered) and exits `IXIT=4` if `KN_HOT_SPENT > KN_HOT_TOTAL_PSES`.

**Two identical refresh calls CAN assign hot-pass extensions to
different rays.** The `TRY_HOT` gate-check and the `CRITICAL`
increment are not atomic together: under `SCHEDULE(DYNAMIC)`, several
threads can read `KN_HOT_SPENT < cap` simultaneously, all enter HOTMODE,
then each increments - so `KN_HOT_SPENT` can **overshoot** the cap by up
to `N_TH-1`, and which rays win the finite budget is scheduling-order
dependent.

### Demonstrated [RUN]

`replh` case (DRV_REPL, `KN_RAY_TOL=1e-30`, `KN_RAY_MAXPASS=4`,
`KN_RAY_STALL=0`, `KN_RAY_HOTPASSES=60`, `KN_HOT_TOTAL_PSES=2`,
`KN_ADAPT_RELAX=F`, 3 rays), `OMP_NUM_THREADS=4`, 6 repeats:

- Pass 1 hot-pass spend alternated `2` and `3` (overshoot of cap=2
  observed: `hot-pass budget spent=3 cap=2`); pass 2 always `2`.
- **Committed stores bit-identical every run:** `REPL_MAXABS=0`,
  `REPL_MAXREL=0`, `REPL_RMS=0`, `REPL_VCDIF=0`, `NMIS_P=NMIS_M=0`.

So on this manufactured model the shared-budget race perturbs the
*iteration path* (hot-pass count) but **not the committed fixed point** -
the extra hot passes refine the same converged store. This does **not**
prove production parallel reproducibility: on a harder model where a
hot pass decides converge-vs-diverge, the race could change outputs.

### Deterministic diagnostic configuration (recommended)

For any replay/nilpotency probe: force `N_TH=1`
(`OMP_NUM_THREADS=1` or `KN_MAX_THREADS=1`), record actual
`OMP_GET_MAX_THREADS()`, `KN_HOT_SPENT`, whether
`KN_HOT_TOTAL_PSES` was exhausted, and per-ray `NPASS`/`IXIT`/verdict
(`FG_RAY_VC`). Single-thread replay establishes **serial**
repeatability only; it is necessary but not sufficient for parallel
production claims.

---

## PART 5 - Repeatability experiment [RUN]

New diagnostic subroutine `DRV_REPL` in
`unit_tests/f90/drv_knop_ray.f90` (cases `repl`=59 unordered,
`replo`=60 ordered, `replc`=61 converging baseline, `replh`=62
hot-stress). It snapshots `I_P_DEPTH_NEXT`, `I_M_DEPTH_NEXT`,
`FG_RAY_SEEDED`, arms the production entry gates
(`DO_KNOP_REFRESH=T`, `FREEZE_I_NEXT=F`, `HAS_COMPRESSING=T`,
`KN_STORE_READY=T`, `FREQ_IDX_ITER=KN_NCF-1`,
`SOLUTION_METHOD='INTEGRAL'`), computes an input fingerprint over
`KN_FREQ`/`KN_DLN_*`/`KN_DBB`/`KN_IC`/`KN_CHI_EXT`/`KN_ETA_EXT`/
`GAM_MAG`/`dGAMdR`/`GAM`, runs `FG_KNOP_REFRESH_V1`, restores the full
input snapshot, recomputes the fingerprint, re-runs, and compares every
writable store element plus `FG_RAY_VC`.

Results (`OMP_NUM_THREADS=1`, `knray` binary, production objects):

| Case | Config | MAXABS | MAXREL | RMS | NMIS | VCDIF | FRQFP |
|---|---|---|---|---|---|---|---|
| `replc` | converging baseline | 0 | 0 | 0 | 0/0 | 0 | identical |
| `repl` | unordered, hot budget | 0 | 0 | 0 | 0/0 | 0 | identical |
| `replo` | ordered sweep | 0 | 0 | 0 | 0/0 | 0 | identical |
| `replh` | hot-stress, 1 thr | 0 | 0 | 0 | 0/0 | 0 | identical |
| `replh` | hot-stress, 4 thr x6 | 0 | 0 | 0 | 0/0 | 0 | identical |

`RES1==RES2` in every case (`KN_RF_RESID` identical across the two
passes). The input fingerprint is unchanged after pass 1, confirming
the refresh leaves its `KN_*`/geometry inputs untouched on finite data
(the `GAM_MAG`/`dGAMdR` sanitisation is a no-op here).

**Conclusion:** at `N_TH=1` the refresh is a **bit-identical
deterministic function** of its captured inputs on all tested
configurations, including the ordered path and a contended hot-pass
budget. The only schedule-sensitive element is the shared hot-pass
*path* (PART 4), which did not alter committed output on these models.
No missing state was exposed - the snapshot (`I_*_DEPTH_NEXT` +
`FG_RAY_SEEDED` + held-fixed config) is sufficient for serial replay on
this model set. Regression: `test_knop_ray.py` 41/41 PASS.

Caveat: `FLD`/`FLD_P` (pre-commit transport) were not captured; the
comparison is on committed stores (map B). A transport-level (map A)
repeatability check would need the deferred `FLD` instrumentation.

---

## PART 6 - Compact checkpoint + in-situ FD probe

### Does the structural graph need the opacity pack?

**No - for the sparsity pattern.** Verified:

- `KN_RF_COUPLING(I,L)` (9865-9872) reads only `GAM` sign, `NFM_L`,
  `IKC1/2` - no opacity.
- `MATCH_FREQ`/`DRF_L` depend on `GAM`, grid spacing (`Z`), `KN_DLN_*`,
  `KN_FREQ`, `FG_FREQ_MATCH` - no opacity.
- Upwind direction, frequency-coupling reach, depth ordering, and the
  store read/write index set are set by geometry, velocity sign,
  frequency matching, and boundary flags.

Opacity values (`KN_CHI_*`/`KN_ETA_*`) set the **weights** of the
Jacobian entries, not which entries are nonzero (assuming generic
nonzero opacities). The commit-veto edges (jump/convergence/cap/heal)
are **data-dependent** and are not captured by topology alone.

### Compact topology checkpoint (sufficient for the structural graph)

| Item | Why needed |
|---|---|
| `COMP_DEPTH_IDX(ND)` | compressing-depth -> store-slot map |
| `J_PNT`,`RAY_PNT` | ray<->point index maps |
| `NI_RAY`,`MAX_LS`,`NM_RAY` | per-ray depth/matching structure |
| `GAM`,`GAM_MAG`,`dGAMdR` | velocity sign (direction) + matching span |
| `R_RAY`,`Z` | spatial coords -> `DRF_L` spacing, ordering |
| `KN_FREQ`,`KN_DLN_P`,`KN_DLN_N`,`KN_DBB`,`KN_IC` | frequency grid/matching |
| `KN_R`,`KN_P` | boundary points |
| dims `KN_ND,KN_NP,KN_NC,KN_NCF,KN_ND_EXT,N_COMP` | array bounds |
| config scalars (PART 2.8) + entry gates | branch selection |
| store **index set** (N_COMP,NP,NCF-1 + maps) | which elements exist |

This is order ~MB (geometry + grids + integer maps), **not** the ~300MB
opacity pack. It determines the structural nonzero pattern of map A and
the unconditional transport edges of map B.

### In-situ fixed-opacity FD probe (for numerical Jacobian/nilpotency)

Candidate-cycle finite-difference validation does **not** require
serialising the coefficient pack: run it **in situ** where the frozen
live `KN_*` arrays are already populated. Perturb one store element,
re-enter `FG_KNOP_REFRESH_V1` (or `KNOP_RAY_SOLVE_V1` for a single ray),
measure the output delta -> numerical Jacobian column. This uses the
production coefficients directly and needs only the compact topology
checkpoint to know *which* elements to probe.

### Recommendation: two-tier architecture

- **A. Compact production topology checkpoint** - for the structural
  dependency graph and SCC.
- **B. In-situ fixed-opacity finite-difference probe** - for numerical
  Jacobian entries, spectral-radius/nilpotency checks on candidate
  cycles.

This is **not** a fully portable offline operator checkpoint: it cannot
reproduce the transport operator outside the live run (the opacity pack
stays in memory). If a portable checkpoint is later required (e.g.
offline replay on another node), the full ~40-array serializer remains
the fallback - but it is not needed for the stated graph/nilpotency
goal.

### Data-dependent edges caveat

Map B's veto/cap/heal edges depend on runtime store magnitudes. The
compact checkpoint gives the *transport* skeleton; veto edges must be
reconstructed from a live (or FD-probed) run's active set, or treated
as conditional arcs in the graph.

---

## PART 7 - Recommendation and remaining uncertainties

### Recommended sequence

1. Keep `DRV_REPL` as a permanent diagnostic (it is opt-in, neutral when
   unused, and already regression-clean). Add a pytest entry asserting
   `REPL_MAXABS==0` for `replc`/`repl`/`replo` at `OMP_NUM_THREADS=1`.
2. Build the **compact topology checkpoint** (PART 6-A) for the
   structural graph - much smaller than the originally scoped full
   serializer.
3. Implement the **in-situ FD probe** (PART 6-B) for Jacobian/
   nilpotency validation on candidate cycles.
4. Defer the full ~40-array serializer unless a portable offline
   operator is demonstrably required.

### Remaining uncertainties

- **Parallel production reproducibility is unproven.** The hot-pass
  budget race is real (overshoot observed); on models where hot passes
  change convergence outcomes, multi-thread committed stores could
  differ run-to-run. Recommend `KN_MAX_THREADS=1` (or a fixed thread
  count + recorded `KN_HOT_SPENT`) for any nilpotency/FD measurement.
- **`FLD`/`FLD_P` transport-level repeatability** not yet measured
  (needs instrumentation of `KNOP_RAY_SOLVE_V1` locals).
- **Map B veto edges** are data-dependent; the structural graph captures
  transport edges only.
- **`KN_CHI_EXT`/`KN_ETA_EXT` are single-precision `REAL`** while the
  coefficient/grid arrays are `LDP` - a fingerprint must respect the
  mixed precision.
- Ordered-vs-unordered is **per-ray** (`HAVE_PACK`/`NI>=3`); a "ordered"
  configuration can silently mix ordered and unordered rays.

### Explicit non-claims

- No production nilpotency is claimed.
- No Tarjan/SCC/graph analysis was performed (out of scope).
- No serializer/reader/replay was implemented (deferred per work order).
- P1-T1h not started.
- The repeatability result is for the manufactured `knray` models, not
  a production supernova model.

### References

- Knop et al. (2009), A&A 496, 295-298. DOI:10.1051/0004-6361:200810290
- Knop et al. (2009), A&A 501, 813-820. DOI:10.1051/0004-6361/200911793
- Trefethen & Bau (1997), Numerical Linear Algebra.
