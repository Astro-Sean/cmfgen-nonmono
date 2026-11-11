# KNOP P1-T1b State-Completeness Audit

Follow-up to `KNOP_P1_T1B_LAGGED_MAP.md`. Question: did the `30x30`
fixed-opacity Jacobian include **all** independent persistent radiation
state, or does expanding `I_PREV`/`I_STORE` introduce feedback the state
vector missed?

Scope: `plane/fg_j_cmf_nonmono_v1.f`, the refresh worker
`KNOP_RAY_SOLVE_V1`, and the `p1t1b` driver
(`unit_tests/f90/drv_knop_ray.f90`, `DRV_P1T1B`/`MAP_F`). No production
equation modified. Evidence classes: **[MATH]** code structure,
**[MANUF]** manufactured-model measurement, **[PROD]** production model.

**Headline: the `30x30` Jacobian is the COMPLETE fixed-opacity
store-to-store map, not a conditional restriction. Expanding `I_PREV`
feedback is inside the sweep and was measured; it produces weak
same-slot and backward edges (~1e-3, ~1e-6) that do not close any cycle
in the manufactured geometry. Nilpotency stands for that model;
production geometry remains the open question (PART 3 spot-check).**

---

## PART 1 - `I_PREV` / `I_STORE` sweep-boundary behaviour

All four trace statements verified against source:

| Claim | Verified | Evidence |
|---|---|---|
| `I_PREV` copied from `I_STORE` at `NEW_FREQ` | yes | `I_P_PREV=I_P_STORE`, `I_M_PREV=I_M_STORE`, lines 3004-3005 |
| `I_STORE` overwritten at end of each freq call | yes | `I_*_STORE=MIN(I_*,NONMONO_INEXT_CAP)`, lines 4681-4682 |
| Expanding cells read `I_PREV` | yes | `SOURCE_PRIME+=Q*I_M_PREV` (4093), `Q*I_P_PREV` (4371) |
| Compressing cells read `I_DEPTH_NEXT` at `FIDX+NM_RAY` | yes | 4138-4147 (M), 4407-4416 (P) |

`I_*_STORE`/`I_*_PREV` are `(NRAY_MAX,NP)` allocatables (212-215,
2587-2590). `I_*_STORE` is read **only** at `NEW_FREQ` to produce
`I_*_PREV`; it has no other consumer.

**Sweep start (INIT call, frequency 1 of every great iteration):**

- `I_P_PREV = I_M_PREV = 0` (lines 2715-2716) - explicit reset.
- `I_P_STORE`/`I_M_STORE` are **NOT reset**; they retain the previous
  sweep's last-frequency value.

**Does the first frequency use retained `I_STORE`?** No. Two independent
reasons:

1. At freq1 (`INIT`), no `NEW_FREQ` copy runs, so `I_*_PREV` stays at
   its INIT-zeroed value. Expanding cells at freq1 read `I_*_PREV = 0`
   (cold boundary), not `I_STORE`.
2. `I_*_STORE`'s stale value is overwritten at freq1's *end*
   (4681-4682) before freq2's `NEW_FREQ` copies it into `I_*_PREV`.
   The stale value is flushed before it can ever be read.

**Classification: `I_*_STORE`/`I_*_PREV` are deterministically
reconstructed within each sweep - not independent persistent state and
not explicitly reset, but slaved bookkeeping whose only cross-call
survival is a value that is provably never read.** **[MATH]**

The comment at 2227-2231 of the driver states the same contract and the
code confirms it.

---

## PART 2 - the potential diagonal feedback chain

Proposed chain: `X(s+1) -> compressing intensity at freq s -> I_STORE(s)
-> I_PREV(s+1) -> expanding intensity at s+1 -> F(s+1)`, implying
`dF(s+1)/dX(s+1) != 0`.

**The mechanism is real.** In store-slot terms (slot `k` = freq `k+1`):

- freq `k`'s compressing read uses slot `k` = `X[k]`; freq `k` writes
  slot `k-1`.
- freq `k`'s intensity -> `I_STORE` at freq `k`'s end -> `I_PREV` at
  freq `k+1`'s `NEW_FREQ`.
- freq `k+1`'s expanding cells read `I_PREV` = freq `k` intensity.
- freq `k+1` writes output slot `k`.

So `output[k]` depends on `X[k]` (same slot) through the expanding
chain, *provided* freq `k+1` has expanding cells upwind of its
compressing store points. **[MATH - mechanism confirmed]**

**Numerical confirmation (manufactured `J_T`, `knray p1t1b`,
`INEXT_RELAX=1.0`):** the same-slot (`dslot=0`) dependencies are present
and nonzero:

| channel | edge | gain | slots |
|---|---|---|---|
| P | `out(P,s,kc0) <- in(P,s,kc1)` | ~9.4e-4 | all 5 |
| M | `out(M,s,kc2) <- in(M,s,kc1)` | ~5e-9 | all 5 |

So `dF[k,kc_out]/dX[k,kc_in] != 0` for same slot `k`, **but only
between different compressing depths** (`kc_in=1 -> kc_out=0` or `2`).
No element maps to itself; the edge is one-directional in `kc`. These
are the expanding-`I_PREV` signature measured end-to-end through the
full formal solve. **[MANUF]**

Whether the chain *breaks nilpotency* requires the same-slot edges to
close a directed cycle. In the manufactured geometry they do not (see
PART 5). The gain (~1e-3) is also far below the dominant compressing
ladder (~0.96), so even a closed cycle would have small loop gain here.

**Inward/outward tested separately:** P channel carries the ~1e-3
same-slot coupling; M channel ~5e-9 (negligible). The asymmetry matches
P1-T1's finding that the outward `I_P` channel dominates.

---

## PART 3 - complete parabolic source stencil

`dS`/`S` propagate `SOURCE_PRIME` one node beyond the marching front:
`I_M(I_c)` can depend on `SOURCE_PRIME(I_c+1)`; `I_P(I_c)` on
`SOURCE_PRIME(I_c-1)` (monotonicity clamps at 4231-4238 / 4539-4546 do
not remove the dependency, only bound it). An expanding `I_PREV` donor
at `I_e` therefore reaches a compressing store write whenever `I_e` is
upwind of *or one node beyond* `I_c` in that channel's march.

This widens the backward reach by exactly one ray point but does not
change the slot structure: the expanding donor is still the previous
frequency's intensity, so the dependency remains "output slot `k` can
read input slots `<= k`" through the within-sweep `I_PREV` chain. The
stencil affects *which ray points* couple, not *which frequency slots*.
**[MATH]**

The full footprint `SOURCE_PRIME -> S -> dS -> clamps -> I_M/I_P
recurrence -> I_STORE -> I_PREV -> I_DEPTH_NEXT write` was traced; every
link is accounted for in the dependency rules of the companion report
(`KNOP_PRODUCTION_DEPENDENCY_GRAPH.md` section 2.1).

---

## PART 4 - combined channel graph and inner boundary

`I_M` marches increasing `I` (outer->inner, 4241-4247); `I_P` marches
decreasing `I` (inner->outer, 4563-4569). No direct `I_M<->I_P` coupling
inside the recurrences. The only cross-channel link is the inner
boundary:

```fortran
IF(LS .LE. NC)THEN; I_P(NI,LS)=I_CORE
ELSE;               I_P(NI,LS)=I_M(NI,LS); END IF   ! ~4558-4562
```

- `LS > NC` (non-core rays): `I_P(NI)=I_M(NI)` always -> `I_M_PREV`
  contributions inside `I_M(NI)` enter `I_P` at the inner boundary and
  propagate to every `I_P(I_c)` on that ray. A genuine M->P path.
- `LS <= NC`: only under `ZERO_FLUX`/`PNT_SRCE` (`I_CORE` contains
  `I_M(NI)`); absent under `DIFFUSION`.

**In the manufactured `J_T`: zero cross-channel entries** (`P<-M` and
`M<-P` both empty at `|J|>1e-9`). The inner-boundary path did not
populate - consistent with the model's ray/boundary configuration.
**[MANUF]** For production non-core rays under `ZERO_FLUX`/`PNT_SRCE`
this is a live M->P feedback path that the manufactured model did not
exercise; it is a candidate cycle-closer to check in PART 3.

---

## PART 5 - graph analysis

State vector = writable `(compressing-depth, slot, channel)` elements;
`NW=3` depths x `NSLOT-1=5` slots x `2` channels = `30` nodes (the
phantom slot `NCF` and non-writable depths are constant inputs, correctly
excluded). Rays are independent (`NP` disjoint components; `NP=1` here).

Tarjan SCC + Kahn topo-sort on the `J_T` sparsity pattern
(`p1t1b_Jswp_w10_a1.txt`):

| threshold | edges | self-loops | non-triv SCC | DAG | max chain |
|---|---|---|---|---|---|
| 1e-9 | 30 | 0 | none | yes | 6 |
| 1e-6 | 21 | 0 | none | yes | 6 |
| 1e-4 | 21 | 0 | none | yes | 6 |
| 1e-3 | 12 | 0 | none | yes | 5 |
| 1e-2 | 8  | 0 | none | yes | 5 |

**No element self-loop** (`J[j,j]=0` all `j`), **no non-trivial SCC**,
acyclic at every threshold including 1e-9 (below which only FD noise
remains). The same-slot and backward edges are dead ends: they feed
outputs that never feed back to their inputs. Max dependency chain 6 =
the 5-slot ladder plus one. **[MANUF]**

Because the graph is acyclic, `J_T` is strictly (generalized)
triangular in topological order -> all eigenvalues 0 -> nilpotent ->
relaxed map spectrum `1-omega`. This is a *structural* confirmation,
not just a spectral one. **[MANUF]**

Same-slot and lower-slot dependencies **do exist** (the ~1e-3 and
~1e-6 entries) - the graph is not *purely* slot-ascending - but they are
acyclic here.

---

## PART 6 - was the `30x30` Jacobian complete?

Checking each completeness criterion:

- **`I_STORE` initialised identically per evaluation?** Yes - `MAP_F`
  restores `I_P_STORE=SPSP`, `I_M_STORE=SMSP`, `I_P_PREV=SMPP`,
  `I_M_PREV=SMMP` (the post-cold-sweep snapshot) at entry
  (drv_knop_ray.f90 2233-2236), so every `F` call starts from identical
  bookkeeping.
- **Was `I_STORE` part of the map input?** No - and correctly so: it is
  a within-sweep slave (PART 1), never an independent persistent
  variable.
- **Did the first frequency depend on previous-sweep `I_STORE`?** No -
  `INIT` zeroes `I_PREV`; `I_STORE` is overwritten at freq1's end before
  any `NEW_FREQ` copy. Verified in PART 1.
- **Were `I_PREV` values fully recomputed?** Yes - `I_PREV=0` at INIT,
  then `=I_STORE` (fresh this-sweep intensity) at each `NEW_FREQ`.
- **Did the wrapper restore all relevant state?** Yes - the four slaved
  arrays are snapshot/restored; `XV` installs only into the active
  writable store elements over the `STO_P`/`STO_M` baseline; all other
  inputs (`ETA_M/CHI_M/ESC_M`, sigma, geometry, `JQW_A` etc.) are
  deterministic per call.

**Independent persistent radiation state audit:** the only radiation
state carried between great iterations is the lagged compressing store
`I_*_DEPTH_NEXT`. Expanding-depth intensities are recomputed each sweep
(not stored); `I_*_STORE`/`I_*_PREV` are within-sweep slaves; there is
no separate persistent `J` array feeding the formal solve (the emergent
`J` is an output). The state vector `X` = writable
`(depth,slot,channel)` store elements is therefore the **complete**
persistent radiation state for fixed opacities. Inactive store elements
are constant inputs (zero/phantom/baseline) contributing no
`X`-dependent term.

**Verdict: the `30x30` Jacobian is the COMPLETE production-map Jacobian
for the fixed-opacity store sub-map - NOT a conditional restriction.**
The expanding `I_PREV` feedback was inside the sweep and is measured in
`J_T` (the same-slot/backward edges). The original nilpotency result is
**valid and complete for the manufactured model**; it is preserved, not
discarded. **[MANUF]**

The only sense in which it is "conditional" is the standard one: it is a
*fixed-opacity* map. Population-opacity feedback is a separate outer map
(P5), unchanged by this audit.

---

## PART 7 - refresh as a separate operator

`FG_KNOP_REFRESH_V1` -> `KNOP_RAY_SOLVE_V1` does **not** iterate the
store; it iterates a per-ray *field* `FLD` to a fixed point, then commits
via `WRITE_STORE`/`WRITE_STORE_CAP`. Its internal reads are the opposite
frequency directions in the **same iterated field**:

- compressing: `FLD(JPR, L+NFM)` (redder, `SRC_PRIME`/`KN_ORD_SPV`)
- expanding: `FLD(L-1)` (bluer)

This is structurally different from the in-sweep store map: within one
`FLD` array, compressing rows read *forward* in frequency while expanding
rows read *backward*. Spatially alternating expanding/compressing bands
can therefore create **cycles inside the refresh's internal fixed-point
iteration** - a real possibility the in-sweep acyclic store map does not
have. **[MATH - structural difference confirmed]**

However, the refresh's convergence is governed by its *iteration*
converging to a fixed point (guards: `KN_RAY_TOL`, `KN_RAY_MAXPASS`,
stall/plateau/divergence IXIT codes), not by Jacobian nilpotence. The
committed store is seed-independent **iff** the field solve converges
cleanly; `WRITE_STORE` overwrites every compressing slot unconditionally
on a clean commit (leaving no seed memory), while diverged/rejected/
`WRITE_STORE_CAP`-bounded commits can retain seed dependence.

**Recheck of `J_full=0`:** on the manufactured model it was exact for
all perturbations tested, meaning the refresh converged cleanly for
every seed - genuine seed independence there, not a tolerance artefact
(at `KN_RAY_TOL=2e-3` the committed store agreed to file precision).
**[MANUF]** A production-grade recheck (multi-seed, large perturbations,
tighter `KN_RAY_TOL`, residual histories, explicit IXIT/cap reporting)
is specified in `KNOP_PRODUCTION_DEPENDENCY_GRAPH.md` section 3 and
remains **OPEN** - it needs a production checkpoint.

---

## PART 8 - conclusions

1. **Direct compressing store dependency** - strictly slot-ascending
   (`read=write+NM_RAY`, `NM_RAY>=1`), dominant gain ~0.96 (P), ~1e-3
   (M). **[MATH+MANUF]**
2. **Expanding previous-frequency dependency** - real, within-sweep,
   produces same-slot (~1e-3 P) and backward (~1e-6) edges; acyclic in
   the manufactured geometry. **[MATH+MANUF]**
3. **Parabolic source stencil** - extends expanding reach one ray point
   beyond the front; changes which depths couple, not the slot rules.
   **[MATH]**
4. **Inward/outward boundary coupling** - only `I_P(NI)=I_M(NI)`; live
   M->P path on non-core/`ZERO_FLUX`/`PNT_SRCE` rays, unpopulated in the
   manufactured model. **[MATH+MANUF]**
5. **Sweep-boundary persistent state** - `I_STORE`/`I_PREV` are slaved,
   not independent; `I_PREV=0` cold boundary at freq1. No hidden state.
   **[MATH]**
6. **Refresh internal dependencies** - opposite frequency directions in
   one iterated field; can host cycles, but convergence (not
   nilpotence) sets seed-independence. **[MATH]**

**Does the nilpotency conclusion stand?**

- For the **manufactured fixed-opacity model: VALID and COMPLETE.** The
  `30x30` Jacobian captured all persistent radiation state; the graph is
  a verified DAG; `J_T` is nilpotent; `J_full=0` is genuine.
- For **production: CONDITIONAL/OPEN.** The dependency rules admit
  backward edges and the `I_P(NI)=I_M(NI)` cross-channel path that the
  manufactured model did not populate. Whether a realistic
  compressing/expanding geometry closes a cycle is decided only by the
  PART 3 production spot-check (specified, not run).

The earlier worry - that `I_PREV`/`I_STORE` might be missed persistent
state invalidating the Jacobian - is **refuted**: they are within-sweep
slaves, correctly excluded, and their feedback is inside the measured
map.
