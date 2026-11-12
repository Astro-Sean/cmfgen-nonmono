# KNOP Production Checkpoint Audit (P1-T1g)

Status: **PART 1 (read-set audit) COMPLETE; PARTS 2-7 (writer, reader,
replay proof, tests) SPECIFIED, not yet implemented.**

This report is the foundation for `DUMP_KNOP_SWEEP_STATE`. It records
the complete source-level read set of the fixed-opacity radiation-store
sweep, the chosen serial capture point, the checkpoint contents, and
the replay-interface requirement. No production numerical equation,
convergence criterion, population equation, or restart behaviour was
modified; no solver code was changed in this cycle.

Evidence classes: **[MATH]** code structure only.

Call chain audited:
`FG_KNOP_REFRESH_V1` -> `KNOP_RAY_SOLVE_V1` -> `SEED_FLD` ->
`SWEEP_FLD`/`SWEEP_FLD_ORD` -> `MARCH`/`SRC_PRIME`/`MATCH_FREQ`/
`BUILD_FREQ`/`SLOPES` -> `WRITE_STORE`/`WRITE_STORE_CAP`.

---

## 1. Capture point

`FG_KNOP_REFRESH_V1` (`plane/fg_j_cmf_nonmono_v1.f`, ~line 6771):

1. early-outs on `FREEZE_I_NEXT` / unallocated stores,
2. freeze-builds `KN_CHI_COEF`/`KN_ETA_COEF` from
   `KN_CHI_EXT`/`KN_ETA_EXT` via `MON_INT_FUNS_V2` (~line 6846),
3. allocates per-ray result arrays, builds `JMP_SIB_M/P`,
4. `!$OMP PARALLEL DO` over `LS=1..KN_NP` calling
   `KNOP_RAY_SOLVE_V1` (~line 6965).

**Chosen capture point: immediately before step 4** (after the
`KN_*` freeze-build and sibling-reference build, before the OMP ray
loop). At this point every fixed-opacity input is frozen and populated,
the lagged stores `I_*_DEPTH_NEXT` hold the generation to be refreshed,
and execution is serial - so a dump here is OMP-safe by construction
and captures exactly what the ray solves consume.

The in-sweep (non-refresh) fixed-opacity path is the per-frequency
`FG_J_CMF_NONMONO_V1` call; its inputs are the same module arrays but
its opacity coefficients are the live `CHI_COEF`/`ETA_COEF`, not the
`KN_*` frozen copies. PARTS 2-4 target the **refresh** store map (the
`I_NEXT` fixed-point operator audited in P1-T1b/c/e), so the
`KN_*`-frozen capture is the correct one. A second capture mode for
the in-sweep path would need `CHI_COEF`/`ETA_COEF` instead - noted as a
supported variant, not implemented.

---

## 2. Read-set audit (what the kernel actually consumes)

The numerical kernel does **not** read the nominal `CHI_COEF`/
`ETA_COEF` directly on the refresh path. `BUILD_FREQ(L)` (internal to
`KNOP_RAY_SOLVE_V1`, ~line 9518) evaluates, per ray-point `I` and
frequency `L`, the frozen cubic log-coefficients:

```
CHI_L(I) = exp( cubic in T1V of KN_CHI_COEF(K,1..4,L) )   or KN_CHI_EXT(K,L) at grid nodes
ETA_L(I) = exp( cubic in T1V of KN_ETA_COEF(K,1..4,L) )   or KN_ETA_EXT(K,L)
DCHR_L   = dlnCHI/dlnR * CHI_L / R_RAY
Q_L, SR_L from GAM/GAM_MAG/dGAMdR + DLNM_L (MATCH_FREQ)
```

packed once per ray into `PK_DTAU,PK_A0..A4,PK_SR,PK_Q,PK_ETANI,
PK_CHINI`. `MARCH`/`SRC_PRIME`/`SLOPES` then read only the `PK_*`
pack, the per-ray geometry/matching arrays, and `FLD`. So the
checkpoint must carry the **`KN_*` frozen inputs** (from which `PK_*`
is deterministically rebuilt) - serializing only nominal inputs would
miss the actual kernel operands.

### 2.1 Frozen opacity / source coefficients (per ray-point x freq)

| Array | Shape | Role |
|---|---|---|
| `KN_CHI_COEF` | `(ND_EXT,4,NCF)` | cubic log-chi coefficients (LDP) |
| `KN_ETA_COEF` | `(ND_EXT,4,NCF)` | cubic log-eta coefficients (LDP) |
| `KN_CHI_EXT`  | `(ND_EXT,NCF)`  | edge chi at grid nodes (default real) |
| `KN_ETA_EXT`  | `(ND_EXT,NCF)`  | edge eta at grid nodes (default real) |
| `KN_FREQ`     | `(NCF)`         | model frequency grid (LDP) |
| `KN_DLN_P`    | `(NCF)`         | ln(nu_{l-1}/nu_l) expanding span |
| `KN_DLN_N`    | `(NCF)`         | ln(nu_l/nu_{l+1}) compressing span |
| `KN_DBB`      | `(NCF)`         | diffusion boundary coefficient |
| `KN_IC`       | `(NCF)`         | inner-boundary intensity |

(`ND_EXT` = extended grid incl. boundary points; `NCF` = `KN_NCF` =
`NCF_IN`.)

### 2.2 Ray geometry and velocity (per ray-point)

| Array | Shape | Role |
|---|---|---|
| `R_RAY`   | `(NRAY_MAX,NP)` | ray-point radius |
| `Z`       | `(NRAY_MAX,NP)` | ray-point height/impact coord |
| `GAM`     | `(NRAY_MAX,NP)` | advection/velocity-gradient (sign picks expanding/compressing) |
| `GAM_MAG` | `(NRAY_MAX,NP)` | |GAM| used in Q_L build |
| `dGAMdR`  | `(NRAY_MAX,NP)` | dGAM/dR used in DCHR_L build |
| `RAY_PNT` | `(NRAY_MAX,NP)` | ray-point -> grid-node map (int) |
| `J_PNT`   | `(ND,NP)`       | depth -> ray-point map (int) |
| `H_PNT`   | `(NRAY_MAX,NP)` | (int, if read) |
| `NI_RAY`  | `(NP)`          | points per ray |
| `MAX_LS`  | `(NP)`          | active-ray mask / max LS |
| `NM_RAY`  | `(NRAY_MAX)`    | per-ray frequency-match offset base (int) |
| `DLMN`    | (shape TBD)     | donor-log-span matching data |

`NRAY_MAX` is the per-ray point capacity; only `NI_RAY(LS)` entries of
column `LS` are live for ray `LS`.

### 2.3 Grid and boundary arrays

| Array | Shape | Role |
|---|---|---|
| `R_EXT`,`LOG_R_EXT`,`V_EXT`,`Z_EXT`,`SIGMA_EXT`,`ETA_EXT`,`CHI_EXT`,`LOG_ETA_EXT`,`LOG_CHI_EXT` | `(ND_EXT)` | extended grid + tabulated opacity/velocity |
| `KN_R` | `(ND)` | model radius grid |
| `KN_P` | `(NP)` | impact-parameter grid |
| `COMP_DEPTH_IDX` | `(ND)` | sparse compressing-depth index (int) |
| `KN_IC`,`KN_DBB` | `(NCF)` | inner-boundary intensity / diffusion coeff |

### 2.4 Lagged radiation-field state (the store map's state)

| Array | Shape | Role |
|---|---|---|
| `I_P_DEPTH_NEXT` | `(N_COMP,NP,NCF)` | inward lagged store (LDP) |
| `I_M_DEPTH_NEXT` | `(N_COMP,NP,NCF)` | outward lagged store (LDP) |

Within-sweep slaves `I_P_STORE/I_M_STORE/I_P_PREV/I_M_PREV` are **not**
persistent state (P1-T1e): they are re-seeded each sweep from
`I_*_DEPTH_NEXT` via `SEED_FLD`, so they need not be serialized - but
the writer should record that they start at their INIT/seeded values.

### 2.5 Solver configuration scalars

Dimensions: `KN_ND, KN_NP, KN_NC, KN_NCF, KN_ND_EXT, N_COMP, NP, NCF_IN, NRAY_MAX`.

Boundary/closure: `KN_INNER_METH` (`DIFFUSION`/`ZERO_FLUX`/`PNT_SRCE`),
`KN_METHOD_SAV`, `KN_INSERT_SAV`, `NC_PNT_SRCE_L` (=2).

Refresh/sweep tolerances and controls (module `SAVE` scalars):
`KN_RAY_TOL`, `KN_RAY_BULK_TOL`, `KN_RAY_MAXPASS`, `KN_RAY_HOTPASSES`,
`KN_HOT_TOTAL_PASSES`, `KN_RAY_ORDERED` (ordered vs alternating sweep),
`KN_RAY_OMEGA`/`OMEGA_RAY` seed, `KN_MAX_THREADS`, `KN_RF_TOT_MB`,
`INEXT_RELAX`, `NONMONO_INEXT_CAP`, `FG_FREQ_MATCH`,
`FG_SEED_COMMIT`, `FG_COMMIT_JUMP`, `FG_COMMIT_CONV`, `FG_COMMIT_NC`,
`DO_KNOP_REFRESH`, `HAS_COMPRESSING`.

Provenance: great-iteration number (`NITSF`/iteration tag), model
identity string, frequency-grid fingerprint (`KN_FREQ` + `R`/`V`/
`SIGMA`/`COMP_DEPTH_IDX` as in `FGINEXT_GS`), capture timestamp, and
the active-ray + compressing-depth mapping needed to enumerate
`(CHANNEL,KC,LS,SLOT)` nodes for PARTS 2-4.

### 2.6 Per-ray work arrays (rebuilt on replay, NOT serialized)

`IKC1,IKC2,JPR1,JPR2,IGRD,WINX,IM_NI,T1V,NFM_L,DLNM_L,DRF_L` and the
`PK_*` coefficient pack are constructed inside `KNOP_RAY_SOLVE_V1` from
the module-level geometry/matching above. They are deterministic
functions of the serialized inputs, so replay rebuilds them; they are
not independent state and must not be treated as such.

---

## 3. Checkpoint format (design)

`KNOP_SWEEP_GS` - unformatted sequential, mirroring `FGINEXT_GS`
conventions: magic + format version + iteration tag + dimension record,
then scalar config, then each array as a whole-array record preceded by
its logical name tag and element count (so the reader can verify
presence/bounds without guessing), then a `-MAGIC` end marker, written
to `KNOP_SWEEP_GS.tmp` and atomically renamed. Whole-array `WRITE`
avoids needing exact shape knowledge at write time; the reader
reconstructs shapes from the recorded dimension scalars.

The writer is gated on a new module flag `KNOP_DUMP_SWEEP_STATE`
(default `.FALSE.`) plus a target-iteration selector
`KNOP_DUMP_ITER` and optional single-ray selector
`KNOP_DUMP_LS` (-1 = all rays). When the flag is FALSE the routine
returns immediately and no file I/O occurs.

## 4. OMP / thread-safety

Capture is at the serial point before the `!$OMP PARALLEL DO`, so no
ray-thread writes the file and no partial/interleaved output is
possible. The dump reads module arrays that are read-only at that
point. Normal execution is unaffected when the flag is FALSE (single
branch on a `SAVE` logical). No extra synchronization is introduced
into the OMP path.

## 5. Replay requirement (PART 5) - the open interface question

Acceptance requires re-executing the fixed-opacity sweep from the
checkpoint and matching the original output store. Two routes:

- **In-situ replay (smallest change):** a diagnostic mode that, after
  capture, perturbs `I_*_DEPTH_NEXT`, re-enters `FG_KNOP_REFRESH_V1`
  (or `KNOP_RAY_SOLVE_V1` per ray), and reads the committed store -
  the same pattern `drv_knop_ray`'s `MAP_F` uses for the manufactured
  model. No solver refactor; runs inside the live process.
- **Offline replay (needs interface):** a standalone driver that loads
  `KNOP_SWEEP_GS` and calls the sweep. This requires exposing
  `KNOP_RAY_SOLVE_V1`/`SWEEP_FLD`/`MARCH` (currently `CONTAINS`
  internals of `KNOP_RAY_SOLVE_V1`, themselves called under OMP) as a
  callable operator with injectable state - a non-trivial refactor the
  task says to report rather than silently undertake.

**Recommendation:** implement in-situ replay first (validates the
checkpoint end-to-end with the least solver change); defer offline
replay until the graph analysis (P1-T1h) shows it is needed. The
checkpoint format is identical either way.

## 6. Status and remaining work

- PART 1 read-set audit: **COMPLETE** (this document).
- PART 2 `DUMP_KNOP_SWEEP_STATE` writer: **SPECIFIED** (section 3);
  array inventory in section 2 is the serialization list.
- PART 3 OMP-safety: **ANALYSED** (serial capture point, flag-gated).
- PART 4 reader + round-trip: **SPECIFIED, not implemented**.
- PART 5 replay: **interface decision required** (section 5);
  in-situ recommended.
- PART 6 node enumeration `(CHANNEL,KC,LS,SLOT)`: supported by the
  recorded `COMP_DEPTH_IDX`/`NI_RAY`/`NM_RAY`/`KN_FREQ`/`N_COMP`/`NP`/
  `NCF` and the `FREQ_IDX_ITER+NM_RAY(I)` slot convention; multi-bin
  matching is carried by `NM_RAY`/`DLMN`/`KN_DLN_*` (do not assume
  `NM_RAY=1`).
- PART 7 tests + this report: **PARTIAL** (audit only).

No production nilpotency is claimed. The checkpoint does not yet
exist; PARTS 2-4 of P1-T1h remain blocked until the writer, reader and
replay proof land and pass.
