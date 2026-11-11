# KNOP P1-T1 — Lagged `I_NEXT` Field-Sensitivity Experiment

Status: **COMPLETE — experiment ran, results interpretable.**
Date: 11-Nov-2026.  Diagnostic only; no solver equation, coefficient,
population algorithm, `KN_NO_GAIN`, ordering, relaxation, or convergence
control was modified.  The only production-source change is two
allocatable debug-capture arrays (`KN_DBG_INEXTM/P`) that are written only
inside existing `KN_DBG_*`-guarded blocks and are unallocated in
production.

Pillar-1 follow-up to `KNOP_SECOND_PASS_AUDIT.md`.  Tests whether
strong-coupling compressing cells (`|Q| -> 1`) produce a poorly
contracting lagged `I_NEXT` iteration mode.

---

## 1. Source locations and dependency chain

`plane/fg_j_cmf_nonmono_v1.f` (all line numbers current):

| Role | Location |
|---|---|
| Store alloc `(N_COMP,NP,NCF_IN)` | ~2600-2612 |
| `COMP_DEPTH_IDX` build (sigma<-1 or V<0, ±1 dilation) | ~2533-2554 |
| Compressing `INEXT_V` read (inward `I_M`) | ~4091-4130 |
| Compressing `INEXT_V` read (outward `I_P`) | ~4371-4396 |
| `SOURCE_PRIME = SR + |Q|*INEXT_V` | ~4130 (M), ~4396 (P) |
| Expanding branch `SR + Q*I_*_PREV` | ~4074 (M), ~4349 (P) |
| Short-characteristics march `I_M` | ~4218-4225 |
| Short-characteristics march `I_P` | ~4538-4543 |
| `S`/`dS` slope construction (both dirs identical) | ~4163-4216, ~4471-4521 |
| Relaxed store write `T2+INEXT_RELAX*(T1-T2)` | ~4830-4835 (P), ~4868+ (M) |
| `KN_DBG_*` coefficient dump hook | ~4428-4466 (P section) |

Lagged-field chain: `I_*_DEPTH_NEXT(KC,LS,fslot)` → two-point spatial
interpolation `INEXT_V` (weights `W_INEXT`, freq slot
`FREQ_IDX_ITER+NM_RAY`) → `SOURCE_PRIME = SR + |Q|·INEXT_V` → `S`,`dS`
→ `A0..A4`,`DTAU` march → `I_M`/`I_P` → relaxed `WRITE_STORE` commit.

`Q = GAM/(DLN·CHI_HAT)` signed; compressing cells have `GAM<0`,
`|Q|=|GAM|/(DLN·CHI+|GAM|)`.  Store slot `s` holds frequency `s+1`;
compressing cells at frequency `k` read slot `(k-1)+NM`.

## 2. Run configuration (manufactured `p1t1` case)

New driver case `p1t1` in `unit_tests/f90/drv_knop_ray.f90`
(`DRV_P1T1`), modelled on `DRV_AB`:

- 8 depths, 2 rays (`NRAYS=1`, dump ray `LSDBG=1`), 6 frequencies,
  `DLOGA=ln(4/3)`, `METH='ZERO'`, `IBND='ZERO_FLUX'`, `FG_FREQ_MATCH=T`,
  `FG_ORD_SWEEP=T`, `DO_KNOP_REFRESH=T` (refresh not invoked),
  `INEXT_RELAX=0.5`, `NONMONO_INEXT_CAP=1e20`, `NITER=8`.
- Velocity: `V0*(R/R1)^0.6` envelope × Gaussian dip → `sigma_true<-1`
  reversal at depths 4-6 (`COMP_DEPTH_IDX` slots 1,2,3).
- **Opacity depression** `CHI_M *= {I2:1e-1, I3:1e-3, I4:1e-5, I5:1e-7}`
  in the reversal zone so `|GAM|/DLN >> CHI` → `|Q| -> 1`.

Result: ray-point 5 is the sole `GAM<0` (compressing) cell,
`Q ≈ -0.998` (`|Q| ≈ 1`) at all six frequencies, `DTAU ≈ 58`
(optically thick).  `KN_DBG_INEXTM/P = 0.0` there — the compressing
branch **fires** (value is `0.0`, not the `-1` sentinel) but the store
reads exactly zero, so the *current* lagged input vanishes while the
*sensitivity* `dI/dI_NEXT` remains well defined.

## 3. Methodology

1. **Coefficient capture.** `KN_DBG_*` dumps the production in-sweep
   operator per (ray-point, frequency): `CHIH,Q,DTAU,SPM,SPP,SR,WX,
   DLNM,NM,KC1,KC2,A0..A4` plus `IMP` = production `I_M`/`I_P`.
2. **Validated re-march.** `unit_tests/python/p1t1_analyze.py`
   re-implements the two marches + `S`/`dS` exactly and reproduces the
   production intensities to **max rel err `I_M=1.1e-8`, `I_P=3.2e-9`** —
   the re-march is the production operator.
3. **Finite difference.** For each cell, perturb the lagged input
   `I_lag` (compressing: `INEXT_V`; expanding: `I_*_PREV`) by an
   absolute `delta = eps*max(|I_lag|,|I_formal|)`; `SP' = SP + |Q|*delta`.
   Two epsilons (`1e-3`,`1e-4`) give a centred derivative; agreement to
   `~1e-10` confirms piecewise-linearity (no clamp/floor active).
   Gain `g = dI_formal/dI_NEXT`.

## 4. Results

Validation: `I_M`/`I_P` reproduced to ~1e-8.  77 cells FD'd.

|Q| regime (all cells):

| group | n | median \|g\| | max \|g\| |
|---|---|---|---|
| `|Q|<0.05` | 45 | 1.4e-5 | 2.6e-3 |
| `0.05-0.5` | 10 | 3.7e-2 | 7.3e-2 |
| `0.5-0.9` | 10 | 4.6e-3 | 8.9e-3 |
| `|Q|>0.9` | 12 | 4.8e-1 | **9.7e-1** |

Split by branch and direction (the decisive cut):

| branch | dir | n | \|Q\| | median gain | max gain |
|---|---|---|---|---|---|
| compress `Q<0` | M | 6 | 1.00 | 1.3e-3 | 1.3e-3 |
| compress `Q<0` | **P** | 6 | 1.00 | **9.7e-1** | **9.7e-1** |
| expand `Q>0` | M | 30 | ≤0.77 | 1.2e-3 | 7.3e-2 |
| expand `Q>0` | P | 35 | ≤0.77 | 2.5e-4 | 2.5e-3 |

Compressing cell I=5, outward gain vs frequency:
`L1..L6 = 0.963, 0.964, 0.966, 0.967, 0.968, 0.969` — flat, robust.

## 5. Worst-cell analysis (I=5, `|Q|=0.998`, `DTAU≈58`)

- `I_NEXT` (store read) = 0 exactly → `K_local=(I_NEXT/I)·g` degenerates
  to 0; the meaningful metric is the absolute gain `g`.
- Outward `g ≈ 0.97 ≈ |Q|`: an optically thick compressing cell has
  `I_P ≈ S' = SR + |Q|·I_NEXT`, so `dI_P/dI_NEXT ≈ |Q| → 1`.
- Inward `g ≈ 1.3e-3`: the lagged input enters `I_M` through the
  far-point `A2` coefficient which →0 for `DTAU≈58`; strongly damped.
- The asymmetry is physical, not numerical: `I_P`'s local-source
  coefficient `A2≈1` (thick cell ⇒ intensity≈source), `I_M`'s is ≈0.

## 6. Interpretation

**The strong-coupling hypothesis is SUPPORTED, refined, and made
quantitative.**

- `|Q| → 1` on compressing cells drives the *outward* formal-solve gain
  `g = dI_formal/dI_NEXT → ~0.97`, versus ≤7e-2 everywhere else.  The
  lagged store feeds back into `S'` with gain `≈|Q|`.
- The store update `store_new = store_old + INEXT_RELAX·(I_new −
  store_old)` gives an effective local map gain
  `(1−INEXT_RELAX) + INEXT_RELAX·g`.  With `INEXT_RELAX=0.5`,
  `≈ 0.5 + 0.5·0.97 = 0.985` — **essentially non-contractive**.  The
  observed convergence relies on `INEXT_RELAX` damping an underlying
  near-unity map; raising the relaxation (the adaptive-relax path can)
  would push the gain toward 1.
- This is consistent with the production symptom: `GS_MAXCHG` plateauing
  / field-maturity non-progress at `FRESH_knop` it~124.  A `~0.985`
  local gain needs ~hundreds of iterations per digit — matching the
  stall far better than the retracted near-zero-`GAM` mechanism.
- Classification: **(B) strongly coupled and weakly contractive** at
  the formal-solve level for the outward channel; the inward channel is
  (A) strongly coupled but contractive.  Not (C) locally expansive
  (`g<1` always here), and not (D) — the effect is the `|Q|→1` source
  coupling itself, cleanly isolated.

## 7. Limitations and unresolved questions

- **Local, single-ray, single-cell-per-frequency estimate.**  The
  global lagged-map eigenvalue is the spectral radius over all
  compressing cells × frequencies × rays × the population feedback.
  Local gain `~0.985` establishes a *poorly-contracting mode exists*;
  it does not by itself prove global divergence.  A matrix-free
  dominant-eigenvalue experiment over the full store (production sweep,
  not per-frequency re-march) is the natural next step (P1-T1b).
- **`I_NEXT=0` operating point.**  The compressing store read is
  currently zero in this manufactured state, so the measured `g` is the
  *derivative* at a cold-lagged-field point.  Production runs with a
  populated store sit at `I_NEXT>0`; `g` is operator-linear so the gain
  should be similar, but this is unverified on a warm store.
- **Single compressing ray-point.**  Only I=5 has `GAM<0`; the
  compressing *store-write* depths (4,5,6) are broader.  Multi-cell and
  cross-frequency store sharing (slot `f` read by freq `f`, written by
  the sweep) are not exercised by the per-frequency re-march.
- **Optically thick regime.**  `DTAU≈58` makes `I_P≈S'`; thinner
  compressing cells would show smaller `g`.  m15's interaction shell is
  also thick, so this is representative, but the `g→|Q|` limit assumes
  the thick limit.
- **`INEXT_RELAX` dependence.**  The `~0.985` figure uses the driver's
  `0.5`; production's adaptive `KN_ADAPT_RELAX`/`KN_OMEGA_POOL_V1` may
  differ per run.  The *pre-relax* gain `g≈0.97` is the invariant
  result.

## 8. Production in-sweep finite-difference (authoritative method)

The re-march above is a Python re-implementation validated against the
production operator.  This section reports the **direct** measurement:
the production `FG_J_CMF_NONMONO_V1` in-sweep formal solver itself is
driven with perturbations of the lagged `I_NEXT` store.  This is the
method the P1-T1 task specifies; it confirms the re-march result
independently and adds footprint, amplitude-consistency and
boundary-behaviour data the re-march cannot give.

### 8.1 Driver design (`DRV_P1T1` probe phase)

After the `NITER=8` warm-up populates the lagged store, the driver runs
a second frequency sweep under `FREEZE_I_NEXT=T` so the store is held
fixed at its post-warmup value (the flag gates only the `WRITE_STORE`
commit block at ~4726 and the refresh at ~6787; the `I_NEXT` *reads* are
unaffected).  `JFNK_IN_EVAL=T` is set so the `L=1` `INIT` does not
disarm the freeze (the flag is cleared on `INIT` unless `JFNK_IN_EVAL`,
~1831).  For each frequency the sweep is:

- `NEW_FREQ=T` call → baseline `I_M`/`I_P`; `I_*_STORE` (the
  previous-frequency upwind donor) is saved.
- `KN_DBG` diagnostics for that frequency are snapshotted (each
  re-solve overwrites `KN_DBG(:,L)` for every cell).
- For each compressing cell (`Q<-0.3`), the lagged store elements it
  reads (`I_*_DEPTH_NEXT(KC1/2,LSDBG,FREQ_IDX_ITER+NM)`) are perturbed
  coherently by `+delta` and the **same** frequency is re-solved with
  `NEW_FREQ=F` (keeps `FREQ_IDX_ITER` and `I_*_PREV` fixed).  Because the
  two bracket weights `(1-W_INEXT)` and `W_INEXT` sum to one, a coherent
  `+delta` on both elements gives `Delta I_NEXT = delta` exactly, so
  `GAIN = Delta I_formal / delta = dI_formal/dI_NEXT` even when the
  baseline store reads zero.
- `I_*_STORE` is restored to the saved donor after each frequency so
  the next frequency's `I_*_PREV` stays clean (the `NEW_FREQ=T` call
  unconditionally rewrites `I_*_STORE = MIN(I_*,CAP)` at ~4672).

Two amplitudes (`AMPF = 1e-4` and `1e-5` times the local formal
intensity) check linearity; `INEXT_RELAX=0.5` is fixed
(`FG_ADAPT_OMEGA_V2=.FALSE.`), giving `J_RELAXED = 0.5 + 0.5*GAIN`.

### 8.2 Results

Cell I=5 is the only `Q<-0.3` cell; `Q ≈ -0.998` at all frequencies.
`I_NEXT` baseline = 0 (cold store), `I_P ≈ 14-16`, `I_M ≈ 0.6-1.2`.

| dir | freq | EPSV(1e-4) GAIN | EPSV(1e-5) GAIN | IMX | RMAX |
|---|---|---|---|---|---|
| P | 1 | 0.962700 | 0.962699 | **3** | **0.618** |
| P | 2-6 | 0.964-0.969 | identical | 5 | =dI |
| M | 1 | **-5990** | **-59907** | 5 | 0.709 |
| M | 2-6 | 1.29e-3 | 1.29e-3 | 5 | =dI |

- **Outward `I_P`, L=2-6**: `GAIN = 0.964-0.969`, amplitude-independent
  to ~1e-6 — a clean linear derivative.  Matches the re-march
  (0.963-0.969) to four figures.  `J_RELAXED ≈ 0.5+0.5*0.966 ≈ 0.983`.
- **Inward `I_M`, L=2-6**: `GAIN ≈ 1.3e-3`, contractive, matches
  re-march.  `J_RELAXED ≈ 0.5006`.
- **L=1 (bluest frequency, `I_*_PREV=0` INIT boundary) is anomalous in
  both directions** and is *not* part of the smooth contraction mode:
  - Outward: the response peaks at depth **I=3** (`RMAX=0.618`, ~450x the
    I=5 value), i.e. the store perturbation propagates to a second ray
    point rather than staying local.  The I=5 gain is still 0.963.
  - Inward: `Delta I_M = -0.709` is **amplitude-independent**
    (`GAIN` scales as `1/EPSV`: -5990 vs -59907).  A perturbation-size-
    independent response is a discrete branch/floor transition, not a
    derivative — consistent with the negative-intensity floor or a
    march-branch flip at the `I_*_PREV=0` boundary.

### 8.3 Footprint and state-repeatability checks

- The perturbation is applied to exactly the store elements
  `KN_DBG_KC1/KC2/NM/FS` identify as cell I's donors; response appears
  at the predicted cell for L=2-6 (`IMX=5`), confirming the donor
  mapping.  At L=1 the footprint is non-local (`IMX=3`), flagged above.
- Repeatability: the two amplitudes give identical gains to ~1e-6 for
  every clean cell, so the perturbed solves are deterministic and the
  operator is locally linear there.  (`FREEZE_I_NEXT` was verified by
  tracing: it gates only the store-commit block ~4726 and the refresh
  ~6787; reads and `I_*_PREV` handling are untouched.)

### 8.4 Net confirmation

The production-solver measurement **independently confirms** the
re-march: the outward compressing channel has formal gain
`g -> ~0.97` (`|Q| -> 1`, optically thick `I_P ≈ S'`), giving an
effective relaxed map gain `~0.983` at `INEXT_RELAX=0.5` — weakly
contractive, the plausible `GS_MAXCHG` plateau mechanism.  The inward
channel is contractive (`~1.3e-3`).  The L=1 boundary anomaly is a
separate numerical artefact to keep in mind but does not change the
outward-channel conclusion.

## 9. Verdict on the hypothesis

**Supported.**  Strong-coupling (`|Q|→1`) compressing cells produce a
formal-solve lagged-field gain `g→~0.97` in the outward direction —
nearly non-contractive before relaxation, `~0.985` effective map gain
with `INEXT_RELAX=0.5`.  This is now confirmed by **two independent
methods**: the validated Python re-march (§3-§6) and the direct
production-solver store perturbation (§8).  The earlier
"near-zero-GAM marginally contractive" phrasing stays retracted: `|Q|→0`
gives `g→0`, strongly contractive.  Whether this fully accounts for the
production stall requires the global-eigenvalue and population-coupling
follow-ups in `KNOP_SECOND_PASS_TODO.md` (P1-T1b, P5 items).

## 10. Files

- `unit_tests/f90/drv_knop_ray.f90` — `p1t1` case (`DRV_P1T1`):
  `KN_DBG`/`IMP` dump for the re-march **plus** the `FREEZE_I_NEXT=T`
  production-solver FD probe (`FDRP`/`FDRM`/`FDSKIP` output).
- `unit_tests/python/p1t1_analyze.py` — re-march, validation, FD, grouping.
- `plane/fg_j_cmf_nonmono_v1.f` — `KN_DBG_INEXTM/P` debug captures only.
- Raw dump: `/tmp/p1t1.out` (regenerate via `knray p1t1`); FD table
  `/tmp/p1t1_fd.json`; production-solver FD `/tmp/fd4.out`.
