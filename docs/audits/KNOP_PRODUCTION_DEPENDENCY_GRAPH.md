# KNOP Production Dependency-Graph Audit (P1-T1c)

Status: PARTS 1-2 complete (code trace + graph construction on the
manufactured store map). PART 3 production spot-check and PART 4
multi-seed refresh test are specified but NOT yet executed on a
production model - see "Open items".

Scope: `plane/fg_j_cmf_nonmono_v1.f` (module `FG_J_CMF_NONMONO_MOD_V1`
plus the external `KNOP_RAY_SOLVE_V1` worker it hosts), its caller
`new_main/cmfgen_sub.f`, and the refresh entry `FG_KNOP_REFRESH_V1`.
No production numerical equation was modified. The only production-tree
change in this audit series remains the default-`.FALSE.` diagnostic
flag `KNOP_DIAG_NO_REGRID` (P1-T1b).

This report distinguishes three evidence classes throughout:

- **[MATH]** - follows from the code structure alone.
- **[MANUF]** - measured on the manufactured `drv_knop_ray` p1t1b model
  (6 frequencies, 1 ray pencil, fixed sigma profile).
- **[PROD]** - requires a realistic non-monotonic model; not yet run.

---

## 1. Complete read/write trace of the radiation stores

Store arrays (module allocatables, `SOLUTION_METHOD='INTEGRAL'` and
`HAS_COMPRESSING` only, allocated `(N_COMP,NP,NCF_IN)` at lines
2620-2623, zeroed on allocation):

| Array | Role |
|---|---|
| `I_P_DEPTH_NEXT`, `I_M_DEPTH_NEXT` | lagged compressing-coupling store; slot `s` holds frequency `s+1` |
| `I_P_PREV`, `I_M_PREV` | previous-frequency ray intensity (expanding upwind donor) |
| `I_P_STORE`, `I_M_STORE` | slaved bookkeeping: overwritten each call, copied into `_PREV` |

`FREQ_IDX_ITER` is reset to 0 on the `INIT` call (line 2726) and
incremented once per `NEW_FREQ` call (line 2998); during frequency `f`
(1-based, strictly decreasing frequency enforced at 3007-3011) it equals
`f-1`. The sweep therefore writes store slots `1..NCF_IN-1`; slot
`NCF_IN` is a never-written phantom read as 0 at the red edge.

### 1.1 `I_P_PREV` / `I_M_PREV` - every access

| Line(s) | Access | Meaning |
|---|---|---|
| 212-215 | declare | allocatable `(NRAY_MAX,NP)` |
| 2587-2590 | allocate | per ray-point x ray-pencil |
| 2715-2716 | `=0` on INIT | cleared at sweep start |
| 3004-3005 | `=I_*_STORE` on NEW_FREQ | previous frequency's intensity |
| 4093 | read `I_M_PREV(I,LS)` | expanding (`GAM>=0`) upwind source, M channel |
| 4371 | read `I_P_PREV(I,LS)` | expanding upwind source, P channel |
| 4466-4467 | `KN_DBG_*` capture | diagnostic only |

`I_*_PREV` is **recomputed deterministically during the ordered sweep**,
not independent persistent state: at each `NEW_FREQ` it is overwritten
from `I_*_STORE`, which at lines 4681-4682 is overwritten from the
just-computed capped `I_*(:,:)` ray field at the end of every frequency
call. There is no path by which `I_*_PREV` survives a full sweep as an
independent variable; it is a within-sweep lagged copy of the previous
frequency's solution. **[MATH]**

Consequence: `I_*_PREV` carries no memory across map applications, but
*within* one sweep it couples frequency `f` to frequency `f-1`. Because
the previous frequency's intensity was itself computed using that
frequency's lagged-store read, `I_*_PREV` transitively depends on the
lagged store at lower slots. This is the only mechanism by which the
expanding branch can feed back into the store. **[MATH]**

### 1.2 `I_*_DEPTH_NEXT` - every access

| Line(s) | Access | Slot / index |
|---|---|---|
| 251-252 | declare | `(N_COMP,NP,NCF_IN)` |
| 2072-2073 | deallocate | on compression-set rebuild (`NEW_R_GRID`) |
| 2620-2623 | allocate + `=0` | on `HAS_COMPRESSING` first setup |
| 6104-6105, 6548-6549 | `=0` | heal/reset paths |
| 4138-4147 | read `I_M_DEPTH_NEXT` | slot `FREQ_IDX_ITER+NM_RAY(I)` (compressing SOURCE_PRIME, M) |
| 4407-4416 | read `I_P_DEPTH_NEXT` | slot `FREQ_IDX_ITER+NM_RAY(I)` (compressing SOURCE_PRIME, P) |
| 4323, 4624 | debug dump | hardcoded `FREQ_IDX_ITER+1`, diagnostic only |
| 4765-4766 | read into `SNP_P/M` | slot `FREQ_IDX_ITER` (relaxation old value) |
| 4840-4843 | write `I_P_DEPTH_NEXT` | slot `FREQ_IDX_ITER` (relaxed or unrelaxed) |
| 4894-4897 | write `I_M_DEPTH_NEXT` | slot `FREQ_IDX_ITER` |
| 5886-5887, 6036-6039, 6227-6229 | checkpoint write/read | unit 291, `FG_KNOP` checkpoint |
| 6456-6504 | seed/init from `I_F` | refresh/seed path |
| 6924-6925 | `KN_JMP_SIBREF_V1` | read-only sibling-scale snapshot |
| 7631, 7686 | `SEED_FLD` | read into local `FLD`/`FLD_P` (refresh) |
| 7730-7731 | `FG_JUMP_MAX` | read (commit-jump veto) |
| 7776-7777, 7795-7796 | `WRITE_STORE_CAP`/`WRITE_STORE` | write (refresh commit) |
| 7836-7837 | `HEAL_STORE` | write (zero poisoned slots on divergence) |

The only *transfer-solve* reads are the compressing `SOURCE_PRIME` reads
at 4138-4147 (M) and 4407-4416 (P), all at slot
`FREQ_IDX_ITER+NM_RAY(I)`. No code path reads `I_*_DEPTH_NEXT` at slot
`FREQ_IDX_ITER` (the slot being written this call) during the transfer
sweep; the SNP read at 4765-4766 is bookkeeping for the relaxation, not a
physics donor. **[MATH]**

### 1.3 Slot arithmetic - the compressing read is strictly slot-ascending

- Write slot for frequency `f`: `FREQ_IDX_ITER = f-1`.
- Read slot for frequency `f`: `FREQ_IDX_ITER + NM_RAY(I) = f-1+NM_RAY`.
- `NM_RAY(I)` defaults to 1 (lines 3617-3620) and is only ever
  *increased* by frequency matching: `K2 = MAX(1,CEILING(T2))` capped at
  `NCF_IN-FREQ_IDX_ITER-1` (lines 3645-3655), so `NM_RAY >= 1` always and
  `NM_RAY <= NCF_IN-FREQ_IDX_ITER-1` when matching is active. With
  `DO_NONMONO_FREQ_MATCH=F`, `NM_RAY` is identically 1.
- Therefore `read_slot = write_slot + NM_RAY >= write_slot + 1`:
  **strictly greater** than the slot written this call, and greater than
  every slot already written this sweep (slots `1..f-2`). The read always
  returns previous-sweep (lagged) data - the intended Gauss-Seidel lag.
  **[MATH]**

The read bounds guard (4115-4116, 4393-4394) admits slots `1..NCF_IN`;
the phantom slot `NCF_IN` is reached only at the last frequency
(`FREQ_IDX_ITER=NCF_IN-1`, where the match guard fails and `NM_RAY=1` is
forced) and reads 0 - the red-edge cold coupling.

### 1.4 `NM_RAY`, `DLMN`, `KC1`, `KC2`, `W_INEXT`

- `NM_RAY`/`DLN_RAY`: per ray-point frequency-match offsets, refilled
  every frequency call inside the ray loop before the INIT split
  (3617-3677). Only compressing points (`GAM<0`) with `FG_FREQ_MATCH` and
  positive local `DLN_N_SAV` are remapped; expanding points keep
  `NM_RAY=1` (unused - their donor is `I_*_PREV`). `NM_RAY` varies per
  ray point and per frequency.
- `KC1`,`KC2`: sparse-store rows of the two depth nodes bracketing ray
  point `I` (`K=RAY_PNT(I,LS)-ND_ADD`), `0` if that depth is not a
  stored compressing depth (4110-4113, 4388-4391).
- `W_INEXT`: radius-interpolation weight in `[0,1]` when both brackets
  are compressing and the radius spacing is non-degenerate
  (4133-4137, 4401-4405); single-sided fallbacks use `KC1` or `KC2`
  alone. The donor feeds
  `SOURCE_PRIME = SOURCE_RAY + |Q|*INEXT_V` (4152, 4421).
- If the slot guard fails, the coupling is dropped entirely and counted
  (`FG_NOSLOT_RAY`, 4153-4163 / 4422-4431).

### 1.5 Expanding branch - traced separately

For `GAM(I,LS) >= 0` (expanding or turning point), the upwind donor is
the previous frequency, not the store:

```fortran
! M channel, ~4093
SOURCE_PRIME(I)=SOURCE_RAY(I)+Q(I)*I_M_PREV(I,LS)
! P channel, ~4371
SOURCE_PRIME(I)=SOURCE_RAY(I)+Q(I)*I_P_PREV(I,LS)
```

`Q(I)` is the signed velocity-gradient coupling; for expanding cells the
sign convention makes the previous (bluer, already-computed-this-sweep)
frequency the upwind donor. **[MATH]**

### 1.6 Channel coupling

`I_M` marches with increasing ray index `I` (outer->inner,
`I_M(I+1)=I_M(I)*A0+...`, lines 4241-4247). `I_P` marches with decreasing
`I` (inner->outer, lines 4563-4569). They use separately rebuilt
`SOURCE_PRIME`, `S`, `dS` arrays and there is **no direct
`I_M<->I_P` coupling inside the recurrences**.

The single cross-channel link is the inner boundary:

```fortran
! ~4558-4562
IF(LS .LE. NC)THEN
  I_P(NI,LS)=I_CORE
ELSE
  I_P(NI,LS)=I_M(NI,LS)
END IF
```

For `LS > NC` (and `LS <= NC` under `ZERO_FLUX`/`PNT_SRCE`),
`I_P(NI)=I_M(NI)`, so `I_M_PREV` contributions inside `I_M(NI)` enter
the outward channel at the inner boundary and propagate to every
`I_P(I_c)` on that ray. Under `DIFFUSION` this path is absent. There is
**no** direct `I_M_PREV -> I_P_DEPTH_NEXT` or `I_P_PREV ->
I_M_DEPTH_NEXT` donor read; the coupling is indirect through the marched
intensity. **[MATH]**

### 1.7 Spatial stencil and expanding->compressing reach

`I_M(I_c)` at a compressing store point depends on `SOURCE_PRIME` at ray
points `1..min(I_c+1,NI)` (the `+1` enters through the parabolic source
slope `dS`/`S` at the marching front). `I_P(I_c)` depends on
`SOURCE_PRIME` at `max(1,I_c-1)..NI`. An expanding cell at `I_e` can
therefore influence a compressing store write only if it lies inside that
channel's source stencil - i.e. upwind of (or one node beyond) the
compressing point in the march direction.

Whether a given expanding cell is upwind of a given compressing store
point is **geometry-dependent**: `GAM(I,LS)=T1*T2` with
`T2=1+SIGMA_RAY*MU^2`, `MU=Z/R`, so the compressing/expanding sign varies
along the ray and with ray angle. There is no fixed "all expanding then
all compressing" ordering; the code counts multiple resonance turning
points by `T2` sign changes (2906-2916). **[MATH]**

### 1.8 Rays are independent

The store is `(compressing-depth, ray-pencil, frequency-slot)`. Every
read and write uses the **same `LS`** (ray pencil); the compressing read
at `(KC,LS,slot)` never touches `(KC,LS',slot)` for `LS' != LS`. The ray
loop is `!$OMP PARALLEL DO` over `LS`. Given fixed opacities
(`KN_CHI_EXT`,`KN_ETA_EXT`,`KN_FREQ`, captured once per sweep at
3100-3116), the radiation problem decomposes into `NP` independent
per-ray problems. **[MATH]** The dependency graph therefore has `NP`
disjoint components, each a `(N_COMP x writable-slots x 2 channels)`
ladder.

---

## 2. Dependency graph construction and cycle detection

### 2.1 Dependency rules (per ray, per channel)

Let output element `s` = store slot `s` = intensity written by frequency
`s+1`. Frequency `s+1` (`FIDX=s`) computes its compressing-depth
intensities from:

- **Forward/compressing edge:** its own `SOURCE_PRIME` read of slot
  `s+NM_RAY` (`>= s+1`). Input `X[s+NM] -> output[s]`.
- **Backward/expanding chain:** its expanding cells read `I_*_PREV` =
  frequency `s`'s intensity this sweep. Frequency `s`'s intensity depends
  on its own compressing read of slot `s` (`X[s]`) and, recursively, on
  frequency `s-1`, `s-2`, ... each contributing their compressing read of
  slots `s-1`, `s-2`, ... So `output[s]` can depend on `X[j]` for
  `j <= s` through the expanding chain, *provided* the expanding cells
  are upwind of the compressing store points on that ray.

Combining: `output[s]` depends on `input[j]` for `j in {1,...,s+NM}`.
Equivalently `input[j]` influences `output[s]` for `s >= j-NM`.

### 2.2 Structural cycle possibility

The forward edges point high-input -> low-output (`j=s+NM -> s`); the
expanding backward edges point low-input -> high-output (`j<=s -> s`).
Both directions exist in the rules, so **cycles are structurally
possible**: e.g. a self-loop `output[s] <- input[s]` (expanding chain
reaching frequency `s`'s compressing read of slot `s`) and 2-cycles
`output[s] <- input[s+1]` (compressing) together with
`output[s+1] <- input[s]` (expanding). Whether any given cycle is
*dynamically populated* depends on the per-ray expanding/compressing
geometry of section 1.7. **[MATH]**

### 2.3 Manufactured-model graph - verified acyclic

The measured sweep-only Jacobian `J_T` (30x30, `knray p1t1b`,
`INEXT_RELAX=1.0`, central FD, two amplitudes agreeing to ~1e-11) **is**
the weighted adjacency matrix of the manufactured store graph (state =
3 compressing depths x 5 writable slots x 2 channels; the phantom slot
excluded). Running Tarjan SCC + Kahn topological sort on its sparsity
pattern:

| Threshold | edges | self-loops | weak comp | non-trivial SCC | DAG | max chain |
|---|---|---|---|---|---|---|
| 1e-9 | 30 | 0 | 12 | none | yes | 6 |
| 1e-6 | 21 | 0 | 13 | none | yes | 6 |
| 1e-4 | 21 | 0 | 13 | none | yes | 6 |
| 1e-3 | 12 | 0 | 18 | none | yes | 5 |
| 1e-2 | 8  | 0 | 22 | none | yes | 5 |

**[MANUF]** At every threshold the manufactured graph is a DAG: no
non-trivial strongly connected component, no self-loop, maximum
dependency chain 6 (= the 5-slot frequency ladder plus one). The weak
backward edges present in `J_T` (slot differences `-1,-2,-3`, gains
~1e-6 to 1e-12, plus same-slot ~1e-3) are **dead ends** in this geometry -
they do not close any directed cycle. This is a structural confirmation
of the all-zero `dgeev` spectrum, not merely a spectral coincidence.

Topological ordering exists (Kahn succeeds); one valid order is
"frequency-slot descending within each channel/depth" consistent with
the acyclicity. Because the graph is acyclic, the sweep Jacobian is
strictly (generalized) triangular in that order -> all eigenvalues 0 ->
`J_T` nilpotent -> the relaxed map `J_F=(1-omega)I+omega*J_T` has
spectrum `1-omega` -> contraction ratio `1-omega` (=0.5 at the production
`INEXT_RELAX=0.5`). **[MANUF]**

### 2.4 What the manufactured result does and does not prove

- Proven **[MANUF]**: for this 6-frequency, single-ray-pencil,
  fixed-sigma model, the store dependency graph is acyclic and the
  fixed-opacity sweep map is nilpotent.
- The *rules* of section 2.1 are general **[MATH]**: they admit backward
  edges and hence cycles. The manufactured geometry happened to leave
  every backward edge a dead end.
- **Not proven for production [PROD]**: a realistic model has hundreds of
  frequencies, many ray pencils, mixed compressing/expanding bands, and
  `DO_NONMONO_FREQ_MATCH=T` (which can only *increase* `NM_RAY`, never
  create a same/lower-slot read). A different expanding/compressing
  ordering could populate a backward edge that closes a 2-cycle or a
  longer loop. Only a production-scale spot-check (PART 3) can decide.

### 2.5 Refresh (`FG_KNOP_REFRESH_V1`) - seed-independence mechanism

The refresh does **not** reset or zero the store. Per ray it:

1. `SEED_FLD` copies the store column into a local field `FLD`
   (compressing slots only; `FLD(:,1)` stays 0).
2. Runs `KNOP_RAY_SOLVE_V1`'s field solve (`SWEEP_FLD`/`SWEEP_FLD_ORD`),
   which iterates the *field* - not the store - reading the iterated
   `FLD` (via `SRC_PRIME`/`KN_ORD_SPV`) at `L+NFM` (compressing) and
   `L-1` (expanding).
3. Commits via `WRITE_STORE` (unrelaxed) on clean convergence or
   `WRITE_STORE_CAP` (bounded ~5% step) on marginal exits; diverged or
   jump-rejected rays write nothing and `HEAL_STORE` zeroes poisoned
   slots.

Convergence controls (all VADAT-wired): `KN_RAY_TOL`=2e-3,
`KN_RAY_MAXPASS`=60, `KN_RAY_BULK_TOL`=6e-2 (KNOP), `KN_RAY_STALL`=0.98,
`INEXT_RELAX`=[1e-3,1], `FG_COMMIT_JUMP`=50, `FG_COMMIT_NC`=0.05,
`KN_RAY_HOTPASSES`=120. Exit codes IXIT 0-6 (exact/tol/bulk/stall/cap/
diverged/plateau). **[MATH]**

Seed independence holds **only where the field solve converges cleanly**:
a clean commit overwrites every compressing slot unconditionally
(`STORE=V1`, floored 1e-30, capped `NONMONO_INEXT_CAP`), leaving no seed
memory. The guarded paths - diverged rays (write nothing), rejected
commits, `WRITE_STORE_CAP` bounded commits, `FG_COMMIT_CONV=F` - can
leave the store partially seed-dependent. This is the mechanism behind
the measured `J_full=0` **and** the reason a multi-seed refresh test
(PART 4) is still required: `J=0` is genuine only if every seed reaches
the same committed store at the tested tolerance. **[MATH]**

---

## 3. Open items (not yet executed)

- **PART 3 production spot-check**: perturb individual store elements on
  a realistic non-monotonic model across compressing/expanding/mixed
  rays, boundary frequencies, and `NM_RAY` values; record the response
  footprint and test it against the dependency rules of section 2.1.
  Needs a production checkpoint + the `FG_KNOP` store dump path
  (unit 291 checkpoint already exists at 5886-5887/6036-6039).
- **PART 4 refresh seed-independence**: run `FG_KNOP_REFRESH_V1` (or the
  `KNOP_RAY_SOLVE_V1` worker) from zero/scaled/random/boundary-focused
  store seeds at fixed opacity; record residual history, IXIT, iteration
  counts, and final-store agreement at full precision across
  `KN_RAY_TOL`/`KN_RAY_MAXPASS`/`FG_COMMIT_*` settings. Feasible via a
  `drv_knop_ray`-style driver against production objects.
- Both require deciding whether backward edges populate cycles on a real
  velocity field; that is the single open question gating the
  "nilpotent for production" claim.
