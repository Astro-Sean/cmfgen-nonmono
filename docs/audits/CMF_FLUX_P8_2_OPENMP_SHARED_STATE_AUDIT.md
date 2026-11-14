# P8.2 follow-up: OpenMP ray-loop shared-state audit

Scope: source-level audit of the OpenMP ray loops flagged in the
P8.2 convergence report (cross-run field drift ~1e-4 relative
from generation 3 between runs, including same-node cases).
Companion to `CMF_FLUX_P8_2_ES_CONVERGENCE_AUDIT.md`.

## 1. Parallel regions audited

`plane/fg_j_cmf_nonmono_v1.f` (gfortran/ifort -fopenmp/-qopenmp
build; LAPACK link uses the *threaded* MKL
`-lmkl_gnu_thread` / `-lmkl_intel_thread`):

| region | lines | schedule | role |
|---|---|---|---|
| `!$OMP PARALLEL DO SCHEDULE(DYNAMIC)` | 4213-5431 | dynamic | FG_J_CMF_NONMONO_V1 per-frequency ray solve, `DO LS=1,NP` |
| `!$OMP PARALLEL DO PRIVATE(I)` | 5436-5455 | default | AV/CV mean-intensity assembly, own-column writes |
| `!$OMP PARALLEL DO ... SCHEDULE(STATIC)` | 7225-7283 | static | RESTORE store interpolation, own-slice writes |
| `!$OMP PARALLEL DO NUM_THREADS(N_TH) SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(LS)` | 8291-8300 | dynamic, N_TH threads | FG_KNOP_REFRESH_V1 driver: `DO LS=1,KN_NP; CALL KNOP_RAY_SOLVE_V1(LS,...)` |

`obs/obs_frame_sub_v9.f`: `!$OMP PARALLEL` at 938-1828,
`!$OMP DO` over observed-frame frequencies ML at 981; all work
arrays PRIVATE (incl. the dev7 PRIVATE-allocatable fix);
all shared writes inside `!$OMP CRITICAL (CLAMP_CENSUS)` are
integer/last-writer diagnostics, none feed back into fields.

`obs/cmf_flux_sub_v5.f`, `obs/COMP_JCONT_V4.INC`,
`subs/comp_j_conv_v2.f`: no OpenMP directives.

## 2. I_M / I_P read-write map

Formal-solve loop (4213-5431):

- `I_M(I,LS)`, `I_P(I,LS)`: written and read only at the ray's
  own column `(:,LS)` - march writes `I_M(I+1,LS)`, reads
  `I_M(NI,LS)` for the same-ray inner-boundary handoff
  (`I_CORE=I_M(NI,LS)` -> `I_P(NI,LS)`), then marches back
  writing `I_P(I,LS)`. Disjoint per ray.
- `I_M_DEPTH_NEXT(KC1,LS,nu+NM)` / `I_P_DEPTH_NEXT(KC2,LS,..)`:
  **read only** inside the loop (the lagged INEXT coupling),
  always at the ray's own slice. Written only after
  `END PARALLEL DO` (lines 5616-5673, sequential).

KNOP refresh loop (8291-8300 -> KNOP_RAY_SOLVE_V1 and its
host-contained procedures, CONTAINS at 10240):

- `SEED_FLD(FLD,STORE)` (10249): reads
  `STORE(COMP_DEPTH_IDX(KS),LS,LS2-1)` - own slice only.
- `WRITE_STORE`/`WRITE_STORE_CAP`/`HEAL_STORE` (12835/12961/
  13314): write `STORE(KC,LS,LS2-1)` - own slice only.
- `FG_JUMP_MAX(FLD,STORE,REF,ICH)` (13110): reads
  `STORE(KC,LS,LS2-1)` - own slice only. The cross-ray median
  reference it needs is `JMP_SIB_M/JMP_SIB_P`, built once
  **before** the parallel loop by `KN_JMP_SIBREF_V1` at
  8185-8186 (explicit "OMP-race audit, Issue 1" fix: the old
  concurrent MAX over siblings raced the committing threads;
  the immutable pre-loop reference removed that race).
- Per-ray counters (`NSTORE,NCAP,NNEG,...`) are dummy args
  bound to own-LS output slots (`NST_LS(LS)` etc.); per-ray
  work arrays (`FLD,FLD_P,IM_NI,RES_FREQ,OMEGA_RAY,WRC_*,
  DIVERGED,CLEAN_EXIT,...`) are local allocatables/scalars of
  KNOP_RAY_SOLVE_V1 - one private stack frame per call.
- `KN_PA_RES(LS,:)`, `KN_QC_*_RAY(LS,...)`, `KN_JT(LS,ICH)`:
  own-LS diagnostic slices.
- `KN_CENS_*`, `KN_HOT_CNT`, `N_NEG_INT`, `FIRST_FAIL_*`:
  integer census under `!$OMP ATOMIC` or `!$OMP CRITICAL` -
  order-independent integer accumulation, no field feedback.
- `DGESV` at 11878 is an 8x8 solve inside SWEEP_FLD_ORD;
  MKL does not thread a matrix that small - deterministic
  kernel in practice.
- Inner-boundary coupling `I_CORE=KN_IC(L)+IM_NI(L)` reads the
  same ray's own inward field (IM_NI is a per-call array),
  not a sibling's.

Ownership table (shared arrays touched by parallel code):

| array | index written | disjoint? | index read | disjoint? |
|---|---|---|---|---|
| I_M/I_P | (:,LS) | yes | (:,LS) | yes |
| I_M/I_P_DEPTH_NEXT | (:,LS,:) post-loop | yes | (:,LS,:) in-loop (read-only) | yes |
| store passed to WRITE_STORE | (KC,LS,:) | yes | (KC,LS,:) | yes |
| *_LS result slots | (LS) | yes | post-join | n/a |
| KN_HOT_CNT/KN_CENS_* | ATOMIC/CRITICAL | protected | n/a | n/a |
| JMP_SIB_M/P | pre-loop only | n/a | any ray index, read-only | yes |

## 3. Classification

**A - race-free and deterministic synchronous iteration**, at
the level of the audited loops:

- No thread writes an element another thread may write or
  read concurrently; every shared write target is indexed by
  the calling ray's own LS (or post-join sequential).
- No thread consumes another ray's in-flight generation: the
  only sibling-scan (`KN_JMP_SIBREF_V1`) runs sequentially
  before the loop into an immutable reference.
- The update is therefore block-Jacobi over rays (each ray is
  a block), NOT asynchronous Gauss-Seidel: a ray's seed comes
  solely from its own committed previous-generation slice,
  and the implicit barrier joins all writes before the next
  generation's seeds.
- No missing REDUCTION: the only in-loop shared updates are
  integer counters, all ATOMIC/CRITICAL-protected.
- No uninitialised read identified in the ray path: per-ray
  buffers are zero-filled on entry (`FLD(:,:)=0.0`) and every
  seed/write walks only mapped cells.

Caveat on scope: the audit covers every routine reachable
from the parallel loops listed above. It cannot exclude a
channel outside that call tree (see 4).

## 4. Reconciliation with the observed cross-run drift

If the loops are disjoint and synchronous, why do runs drift?
Two consistent mechanisms remain open:

1. The printed `GLOBAL_MAX_R` has ~6 significant digits. FP
   differences far below print precision at generations 1-2
   (e.g. from MKL threaded kernels elsewhere in the call
   graph, or CPU-dispatch differences in the same binary
   across node types) are amplified by iteration into the
   visible ~1e-4 band by gen 3-4. "Bit-identical at gen 2"
   was verified only at print precision, not by comparing the
   binary store files.
2. The compared runs used different thread counts (4 vs 8
   CPU). Under the disjoint structure above that cannot
   change results - so if thread count is the discriminator,
   a channel outside the audited loops exists (threaded MKL
   is the leading candidate: `mkl_gnu_thread` is linked).

An honest D-level residue therefore remains for the *system*
behaviour: the ray loop itself is A, but the empirical drift
is not yet attributed to a specific source outside it.
Decisive test (not run here - KNOP runs cost ~30-90 min):
a repeated same-config run with `MKL_NUM_THREADS=1` (and
`OMP_NUM_THREADS` fixed) would isolate threaded-library
nondeterminism; a bitwise `ES_J_CONV` comparison at equal
generation would locate the first differing element.

## 5. Compiler diagnostics

No race-detector or thread-sanitizer is installed for this
Fortran build environment (gfortran/ifort + MKL on the
cluster); valgrind/helgrind were not attempted (MPI-class
runtimes unsupported). Source-level ownership audit stands
as the primary evidence; the threaded-MKL hypothesis is
testable with the environment-variable control above.

## 5a. Run-to-run floating-point drift attribution

Linked libraries (ldd on exe/cmf_flux.exe and every
exe_frozen/*/cmf_flux.exe): `libmkl_intel_lp64`,
**`libmkl_intel_thread`**, `libmkl_core`, `libiomp5` - the
Intel-threaded MKL layer, which honours MKL_NUM_THREADS.
Build flags (Makefile_definitions): `ifort -qopenmp -mkl`
line at :113 with `-lmkl_intel_lp64 -lmkl_intel_thread`.

Numerical-library call sites in the field-update path:

- `mom_j_cmf_nonmono_v1.f:3232,3341` - `CALL DGTSV` per
  frequency inside the moment solve (invoked sequentially
  from COMP_JCONT_V4.INC:988). DGTSV is a tridiagonal scan;
  MKL does not parallelise a single NDx1 solve - serial
  kernel in practice.
- `fg_j_cmf_nonmono_v1.f:11878` - `CALL DGESV` on an 8x8
  system inside SWEEP_FLD_ORD (inside the parallel ray
  loop). MKL does not thread an 8x8 LU.

So threaded MKL is linked but its in-path call sites are
serial kernels at these problem sizes - the threaded layer
is present but probably not the channel. Not ruled out:
MKL can dispatch CPU-specific code paths by feature, so
cross-node FP differences can still come through MKL even
without internal threading.

Job-record cross-check (sacct) sharpened the empirical
picture:

- knop_n4 (job 13156999): viper08, 4 CPU -> {0.265340,
  0.177620} at gens 3-4.
- knop_n6 (job 13154194): vipc2142, 8 CPU -> {0.265314,
  0.177358}.
- knop_n8 (job 13154196): vipc2304, 8 CPU -> same drifted
  values as n6.
- knop_n12 first run (job 13154178): vipc2142, 8 CPU ->
  {0.265340, 0.177620} - the OTHER cluster on the SAME node
  as n6 with the SAME thread count.

So vipc2142 hosted both clusters at 8 threads: n6 vs n12
differ run-to-run on identical hardware and thread count.
That rules out node-hardware and thread-count explanations
and leaves genuine run-to-run nondeterminism - from a
channel outside the audited ray loops (the loops were shown
order-independent). The channel is not yet identified:
candidate unaudited paths are the moment solver's per-
frequency sequencing and any threaded-MKL dispatch kernel,
plus the possibility that a sub-print-precision difference
at gens 1-2 (from any source) is amplified by iteration to
the visible 1e-4 band.

Additional finding - knop_n12 rerun is warm-started: the
batch template removes ES_J_CONV at start (batobs.sh:133)
but keeps EDDFACTOR, so the resubmitted run seeds its
formal solve from the timed-out run's ~9 generations. Its
gen-2 residual reads R=9.99997E-01 at K=16 with
MAX_ABS=5.8e-4 (~25x below the cold-start 1.44e-2) - the
residual labels are shifted relative to the cold-start
n2-n8 history, but the 12-sweep endpoint remains a valid
(or slightly more converged) N=12 field for the spectral
comparison. Flagged for the N8->N12 interpretation.

Decisive test remains the same-config repeat with
MKL_NUM_THREADS=1 + MKL_DYNAMIC=FALSE plus a bitwise
ES_J_CONV diff at equal generation - not run in this
session (a gen-4 KNOP run costs ~76 min).

Classification update: A holds for the audited ray loops.
System-level drift is confirmed run-to-run nondeterminism
(same node, same threads) - category E (unresolved source)
until the MKL-controlled repeat runs.

## 6. Recommendations

No code change is warranted for the ray loops - no defect
found. For attribution of the residual drift:

1. Run one repeated KNOP config with `MKL_NUM_THREADS=1`,
   `OMP_NUM_THREADS` fixed, same node: if gens 1-12 are then
   bitwise identical to a sibling run, MKL internal threading
   is the channel.
2. If MKL-serial runs still drift, dump `ES_J_CONV` at equal
   generation for a bitwise diff to locate the first
   differing element (depth/frequency), which identifies the
   producing stage.
3. Document that the printed-residual equivalence at early
   generations is print-precision, not bitwise.

Nothing here requires touching the redistribution operator,
ray ordering, formal solution, or NUM_ES_ITERATIONS.
