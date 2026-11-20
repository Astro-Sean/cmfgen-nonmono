# KNOP active-field validity and probe-isolation audit

Date: 20-Nov-2026. Follow-up to `KNOP_PROBE_CONTROLLER_AUDIT.md`
(commit `9d96e25f`). Source: `plane/fg_j_cmf_nonmono_v1.f` (rev
`e1df4324`), `new_main/mod_subs/solve_for_pops.f`. Empirical data:
1610 `KNSEED_P_L0083_*.bin` records across all four arms
(A r1-407, B r1-444, C r1-292, A/A r1-467). All line numbers are for
the source at commit `9d96e25f`.

## 0. Headline correction

The previous audit reported approximately 80,500 invalid-scale finite
cells and 70-95 NaNs per binary record. Those counts were a defect in
`scripts/parse_knop_seed.py`, not in the data. The parser read the two
array payloads starting at each Fortran record's 4-byte leading marker,
so every double was misaligned by 4 bytes and produced byte-shifted
values of effectively random magnitude. Arm A/A still showed
`x0 - xref = 0` because two identically shifted byte streams cancel.

With correct record-marker handling (parser fixed in
`scripts/parse_knop_seed.py`, parent-repo commit), all
1610 records contain:

- zero NaN, zero negative values, zero values above 1E10;
- field maximum 2.1636E-2 across every record;
- `max|x0-xref|` identical to the header `rmax` field in all 1610
  files (within 1E-6 relative).

The reported invalid-cell population, its churn statistics, the
"dead-column resurrection" claim, and the "garbage reaches MARCH"
conclusion are retracted. The transported workspace is clean. Sections
5 and 6 of `KNOP_PROBE_CONTROLLER_AUDIT.md` that rely on those binary
findings should be read as superseded by this document; its
controller-level findings (probe in the residual pool, omega pool,
censuses, freeze mechanics) are independent of the binary arrays and
stand.

## 1. What "no-commit" means: sweep writes vs refresh writes

### 1.1 The two write paths into `I_P_DEPTH_NEXT` / `I_M_DEPTH_NEXT`

Path A - ordinary frequency-sweep write
(`fg_j_cmf_nonmono_v1.f:5728-5840`, P direction; M direction mirrors
at `:5842-5900`):

```
5728  IF(HAS_COMPRESSING .AND. FREQ_IDX_ITER in [1,NCF_IN] .AND.
5731    .NOT. FREEZE_I_NEXT) THEN
5750    VMAX_IP(LS) = MAX(1E-25, MAXVAL(I_P(1:NI_RAY(LS),LS),
                  MASK=finite))          ! per-ray peak floor
5761    SNP_P(KC,1:NP) = I_P_DEPTH_NEXT(KC,1:NP,FREQ_IDX_ITER)
                                      ! slot snapshot, once/freq
5770    DO LS=1,NP                     ! ALL rays, incl. 83
5772      I=J_PNT(K,LS) in [1,NI_RAY(LS)]
5774      T1 = MIN(I_P(I,LS),1E20) if positive finite, else 1E-30
5800      FG_GS_MAXCHG update vs T2=SNP (NEW_FREQ calls only)
5835      IF(T2>0) store = T2 + INEXT_RELAX*(T1-T2)
5839      ELSE     store = T1
```

There is no `L_NOPROBE` anywhere in this path. Ray 83's slots are
written every sweep exactly like every production ray's.

Path B - refresh-time candidate commit (`WRITE_STORE`,
`WRITE_STORE_CAP`, `HEAL_STORE` inside `KNOP_RAY_SOLVE_V1`): all
gated on `.NOT. L_NOPROBE` (`:9793` and the guarded call sites). For
ray 83 these never execute - `DCH_M/DCH_P` out-args stay 0.0.

### 1.2 Answers

A. Does ray 83 contribute to the ordinary sweep store update?
Yes, unconditionally - `DO LS=1,NP` at `:5770` includes it, gated
only by `HAS_COMPRESSING`, `FREQ_IDX_ITER` validity, `J_PNT`
membership and `FREEZE_I_NEXT` (off in these runs).

B. Does the probe's sweep update modify the retained field used
later? Yes - it writes `I_P_DEPTH_NEXT(:,83,:)` slots, which are the
retained lagged store. But the write's inputs are the production
field `I_P(I,83)` and the slot snapshot `SNP_P` - the probe's
perturbed refresh iterate `FLD_P` never enters this path. The sweep
write is arm-independent by construction.

C. Do the KN-NOP checksums cover all store writes? No. `KN_NOP_SP0/
SM0` are sampled at `:8426-8427` immediately before the refresh
parallel loop, `KN_NOP_SP1/SM1` at `:8583-8584` immediately after.
The checksum brackets the refresh-time solve window only; sweep
writes occur in a different phase of the iteration. `delta M/P=0`
proves no refresh-window store write escaped `L_NOPROBE`; it says
nothing about sweep writes (there are none inside the window to
detect).

D. Can the probe's sweep update affect `GS_SW_MAXCHG`? Yes -
`FG_GS_MAXCHG` (`:5805-5809`) is accumulated over the same
`DO LS=1,NP` loop, so ray 83's slot drift enters
`GS_SW_MAXCHG` -> `GS_MAX_ALL` (`solve_for_pops.f:1248-1250`).
This is production data, not seed contamination: the drift is
computed from production `I_P` values and would occur for ray 83
with or without the probe flag.

E. Can the sweep update affect later `SEED_FLD` or other rays?
Ray 83's own retained seed, yes: `SEED_FLD` (`:10627-10641`) reads
`I_P_DEPTH_NEXT(:,83,:)` at the next refresh - the normal production
channel. Other rays, no: each ray's `SEED_FLD` reads only its own
store row (`STORE(...,LS,...)` with the caller's LS), and no
cross-ray term exists in the sweep write.

### 1.3 Data-flow diagram

```
 frequency sweep (per frequency call, all rays, arm-independent):
   I_P(I,LS) --[T1=min(.,1E20)]--+
   SNP_P = store(KC,:,f) ------>+--> I_P_DEPTH_NEXT(KC,LS,f)
                                 = T2 + w*(T1-T2)          :5836
       |--> FG_GS_MAXCHG / FG_GS_NCONV / FG_SW_CAPPED census :5800
       |--> GS_SW_MAXCHG -> GS_MAX_ALL -> maturity gate      :1248
       v
   store slot = retained lagged field (production state)

 refresh window (per refresh, per ray, ray 83 = probe):
   I_P_DEPTH_NEXT(:,83,:) --SEED_FLD--> FLD_P(JP2,LS2)   :10064-10065
   KN_SEED_SETUP_P: XREF=FLD_P; arm injection; metrics;  :10072
                    KNSEED_P bin (XREF,X0)               :14426-14440
   MARCH passes on FLD_P --> RES_LS(83), RHO_P_LS(83),
                    NBAD/NTT/NST/NCL/NNL, candidate FLD_P
   WRITE_STORE/_CAP/HEAL_STORE: SUPPRESSED (L_NOPROBE)   :9793
   KN_SEED_SAVE_P: FLD_P -> KN_SEED_PC (arm C chain)     :10100
   KN-NOP checksum: sums before/after window -> delta=0  :8426/:8583

 global accounting (refresh end):
   KN_RF_RESID = MAXVAL(RES_LS, NDV==0)  incl ray 83     :8813-8816
   KN_OMEGA_POOL_V1: RHOM/RHOP samples    incl ray 83    :2059-2067
   censuses NBAD/NTT/...                  incl ray 83    :8897-8904
   -> GS_MAXCHG -> GS_MAX_ALL -> gate; -> INEXT_RELAX
   -> omega multiplies ALL rays' sweep writes (:5837,:5891)
```

"Probe never commits" is precise only for the refresh-window WRITE
paths; the probe's store slice continues to evolve through the
ordinary sweep writes that constitute the production lagged map.

## 2. Checkpoint-to-MARCH data path and cell census

### 2.1 The path, stage by stage

| stage | location | dims / mapping | init / validity |
|-------|----------|----------------|-----------------|
| `FGINEXT_GS` file | written at checkpoint | `(N_COMP,NP,NCF_IN)` each for P,M | generation-tagged |
| `RESTORE_FG_INEXT_GS` | `:6897-7060` | verbatim read `:7033,:7035` | magic `20260916`, ver 2-3, tag=NITSF, dims match, full R/V/SIGMA/`COMP_DEPTH_IDX` fingerprint; `:7049-7056` rejects non-finite and negative; finite values to `HUGE` pass |
| `I_P/I_M_DEPTH_NEXT` | alloc `:3590` | `(N_COMP,NP,NCF_IN)`, zeroed at alloc | every mapped slot rewritten each sweep; written range `[1E-30,1E20]` (`T1` bound `:5776`, `:5839`), or relax-mix with positive `T2` |
| `SEED_FLD` | `:10627-10641` | `FLD(JP2,LS2)=STORE(COMP_DEPTH_IDX(KS),LS,LS2-1)`, `KS=1..KN_ND`, `LS2=2..NCF` | copies only mapped depths/ray points; `FLD(:,1)` untouched |
| `FLD/FLD_P` | `:9995,:10064` | `(NI=3,NCF=56377)`, zeroed before seeding | column 1 stays 0 by construction |
| `MARCH` | `:11547+` | reads `FLD(J-1,L)`, `FLD(J,L)` over `L` loop | all cells in swept columns participate |
| transfer result | `FLD` candidate -> `WRITE_STORE` | mapped rows commit back | suppressed for probe |

The earlier `RESTORE_FG_INEXT_OBS` path (`:7310`, foreign frequency
grid remap, used by OBS jobs) validates `|I_F|<=HUGE`, `I_F>=0`
(`:7475-7476`) then log-interpolates (`:7462-7464`). The remap output
is not revalidated; a file value within ~1 dex of `HUGE` could
`EXP`-overflow to `+Inf` and enter the store. No finite-value cap is
applied on either restore path. Latent edge case only - see F5.

### 2.2 Measured cell census (corrected parser)

For ray 83's `NI=3 x NCF=56377` workspace, all 1610 records, both
`XREF` (retained-store copy) and `X0` (post-injection seed):

| class | count | notes |
|-------|-------|-------|
| valid active cells | 169,128 | `LS2=2..NCF` x 3 rows; all finite, positive, `1.2E-34 .. 1.26E-2` |
| valid inactive | 3 | `FLD(:,1)` all-zero column (never written) |
| dead/unmapped | 0 beyond col 1 | every mapped column has all 3 rows nonzero |
| non-finite | 0 | |
| invalid-scale finite | 0 | |

`KN_DK_NDEAD=3` in every record is exactly the `FLD(:,1)` column:
`KN_DK_W=0` where the column peak `CM=0` (`:14352-14360`). Bookkeeping,
not a defect.

Origin of the previous invalid-cell claims: the parser artifact in
section 0. There is no evidence of invalid values entering the
transport workspace in any record of any arm.

### 2.3 Where the experiment's seed arrays actually differ

Arm B `X0` (sibling medians): all 169,128 mapped cells differ from
`XREF`; `max|x0-xref|=2.157E-2` vs field max 1.26E-2 - the sibling
median replaces ray 83's values wholesale and by more than the local
field peak in the worst cell. Arm C `X0` (chained candidate):
differs at all cells but by `<=5.3E-7` - the prior last iterate sits
within E-7 of the retained store. Arm A/A `X0==XREF` bitwise.

## 3. MARCH numerical safety audit (source level)

The empirical trigger for this section is gone (no invalid input is
actually present), so the answers below are hypothetical - what the
code would do if an invalid `FLD` cell ever appeared.

Recurrence and guards (`:11547-11577`):

```
VNEW = FLD(J-1,L)*PK_A0(J-1,L) + (SPW*A1 + SPJ*A2 + DSD*A3 + DSU*A4)
IF(VNEW < 0 .OR. VNEW /= VNEW)  VNEW=1E-30; NNEG++     ! output guard
IF(VNEW > NONMONO_INEXT_CAP)    VNEW=1E20;  NCAP++     ! output cap
T1 = |VNEW-VOLD| / MAX(VOLD, 1E-4*VMAXF(L))            ! residual
significant iff MAX(VOLD,|VNEW|) > 1E-4*VMAXF(L)       ! census gate
```

`VMAXF(L)=MAX(MAXVAL(FLD(1:NI,L)),1E-25)` at `:11467` - **no NaN
mask**, unlike the sweep peak `VMAX_IP` at `:5751-5753` which masks
`x==x`.

A. Can a huge finite `FLD` enter the recurrence before capping?
Yes. `FLD(J-1,L)` multiplies `PK_A0` directly as the upwind term;
the guard acts only on `VNEW`. A huge upwind cell produces a huge
or negative-huge `VNEW`, which the guards then force to `1E20` or
`1E-30` - bounded output, but the transported value is the guard,
not the physics, for every cell downwind of it.

B. Can it dominate `VMAX_L` and suppress residuals? Yes.
`MAXVAL(FLD(1:NI,L))` is taken over the whole column; one `1E300`
cell makes `VMAXF=1E300`, so `1E-4*VMAXF=1E296` and every physical
cell in that column fails the significance test - the column
contributes nothing to `NTOT`, `NBAD`, `RES` or `DCELL`. The
column silently self-masks.

C. NaN under Intel Fortran: `MAXVAL` over a NaN-containing array is
processor-dependent (the `MAX`/`MAXSD` lowering decides whether NaN
wins). If `VMAXF` comes out NaN, `MAX(VOLD,1E-4*VMAXF)` is NaN,
`T1` is NaN, `T1.GT.RES` is false and the cell vanishes from both
the residual and the census - the same silent self-masking as (B),
but processor-dependent. Note the sweep peak at `:5751` already
masks NaN; MARCH's `VMAXF` does not - an asymmetry worth aligning
if a NaN guard is ever added.

D. Can a NaN evade `DCELL` and the bad-cell census? Yes, per (C):
every comparison involving `VMAXF` or `DCELL` returns false on NaN.

E. Can the solver report converged with invalid cells present?
Yes - the max-norm and census see only significant cells, and an
invalid cell makes its own column insignificant. (Not observed in
the data; the guard counters `NNEG/NCAP` would show in `KN-CEN`.)

F. Scope: any column holding an invalid cell; every ray uses the
same `SEED_FLD`/store path, so a corrupt store row would hit
whichever ray owns it.

Physically correct treatment if this ever becomes live: the right
fix is upstream (validate at `SEED_FLD` where the only store->FLD
copy happens, or at the store writers), not flooring inside MARCH,
and not zeroing (a seeded 0 changes the transport). No fix is
required by the current data.

## 4. Seed metrics reassessed

`KN_SEED_DAB/RMAX/SCHG` (`:14389-14414`): `DAB` is the weighted
distance `sqrt(sum(W*(x0-xref))^2)/sqrt(sum(W*xref)^2)` with
`W=1/max(|xref|,1E-4*CM_col)`; `RMAX` the largest absolute
single-cell replacement; `SCHG` the fraction of significant cells
changed.

Measured with the corrected parser (all consistent with headers):

- Arm A/A: all three read 0 by construction (no injection).
- Arm B: `dab=8.257E4`, `rmax=2.157E-2`, `schg=1.0`. The weighted
  distance is dominated by weak cells: ~50% of `DKT` comes from
  `xref in (1E-20,1E-10]`, ~50% from `(1E-10,1E-6]`; median
  relative replacement on those tiers is `~1E4-1E5`. Only 38 cells
  exceed `1E-3`, and there the median relative change is ~1.2.
  So `dab` correctly reports that the sibling median is a wildly
  different field on weak frequencies - ray 83 sits far below the
  typical sibling at those columns - while `rmax` anchors the
  absolute story at ~2x the field peak.
- Arm C: `dab=0.177`, `rmax=5.3E-7`, `schg=1.0` - the chain is a
  near-copy of the retained store with E-7-level drift.

Active-cell mask: with clean data the five proposed classes collapse
- every mapped cell (`JP2` valid, `LS2>=2`) is active and valid;
`FLD(:,1)` is the sole inactive column (dead by design, `w=0`
already); no non-finite or invalid-scale cells exist. The solver-
derived mask `finite .AND. |x| <= NONMONO_INEXT_CAP` selects exactly
the observed population (all of it), so it adds nothing today - but
it is the right guard for the metric code if invalid cells ever do
appear: apply it to both `xref` and `x0` in the DAB/RMAX/SCHG loop
(`:14395-14412`) and count exclusions separately.

`KN_SEED_SIB` vs `KN_SEED_SAVE_P` masks: SIB writes only mapped
cells (`JP2` in `[1,NI]`, `LS2>=2`; `:14498-14536`) with sibling
samples filtered to finite positive store values (`:14517-14518`) -
with clean stores this admits only physical values. SAVE_P snapshots
all `NI x NCF` (`:14468`) including the zero `FLD(:,1)`, which is
harmless because the chain replays 0 into a column that is always 0.
The masks are already consistent; no change needed.

Retraction: the earlier claims that SIB "resurrects ~42,000 dead
columns" and that DAB/RMAX conflate invalid cells are artifacts of
the parser bug. The real metric caveat is different: `dab` is a
relative-change metric dominated by near-floor cells - report it
alongside `rmax` (absolute) and a strong-cell `schg` variant
(e.g. cells `>1E-4` of the global field max) so the perturbation
size is not read as uniformly `1E4`-scale.

## 5. Two probe modes, refined

MODE 1 - DIAGNOSTIC PROBE. Goal: identical background across
A/B/C/A/A within OpenMP noise. What must be suppressed:

| contribution | suppress? | reason |
|--------------|-----------|--------|
| refresh-time store writes | YES (already) | `L_NOPROBE` gates `WRITE_STORE/_CAP/HEAL` |
| ordinary sweep writes | NO | production path; the probe's seed never reaches it; suppressing it would change the model |
| `GS_SW_MAXCHG` contribution | NO | ray 83's sweep drift is production state, arm-independent |
| `KN_RF_RESID`/`GS_MAXCHG` pool | YES | seed-dependent probe residual currently the argmax every refresh - the controlling contamination |
| `KN_OMEGA_POOL_V1` samples | YES | arm C's probe rho crosses P90 -> `INEXT_RELAX` -> scales every ray's sweep write |
| NBAD/NTT/NST/NCL/NNL censuses | YES | probe cells enter bad-share/churn tiers -> gate inputs |
| `FG_NC_CAP/REJ`, commit/veto counters | YES | control-adjacent accounting; already 0 via gating but the census side is not masked |
| freeze-release decisions | indirect | follows from gate inputs once masked; no separate flag needed |
| probe telemetry (`KNSEED_P`, `KN-NOP`, pass residuals, `KN_SEED_SAVE_P` chain) | NO | the experiment's product |

Invariant arrays/accumulators that must match A/B/C/A/A to noise:
`I_P_DEPTH_NEXT`, `I_M_DEPTH_NEXT` (whole store - includes ray 83's
own sweep-written slice, since its inputs are arm-independent),
`RES_LS(1:NP)-mask83`, `GS_SW_MAXCHG`, `GS_RF_DCHG`, all census
totals computed without ray 83, `INEXT_RELAX`, `OMEGA_RAY` per ray,
population updates and every downstream quantity. Measured noise
floor: A/A rel p50 ~5.5E-5, p99 ~2E-3 on binary fields.

MODE 2 - PRODUCTION RAY. The armed ray commits normally: ordinary
sweep writes and refresh `WRITE_STORE/_CAP` proceed; jump veto,
capped commit, first-commit marker, divergence-heal and
commit/veto counters all apply. Seed metrics and the `KNSEED_P`
record stay diagnostic-only. Must not be mixed with MODE 1 results;
in MODE 2 the seed genuinely changes the model, so A/B/C
differences are physical outcomes, not contamination.

## 6. Tests required before any implementation

1. Invalid-value provenance: superseded - no invalid values exist in
   the transport workspace (1610-file census). Keep a cheap
   validity sweep over `KNSEED_P` records as a regression check;
   extend to store-level dumps if one is ever added.
2. MARCH safety: synthetic test only if a `VMAXF` NaN mask or input
   validation is added - single-column `FLD` with a planted NaN /
   `1E300` cell, verify the census counter sees it instead of
   self-masking.
3. Residual integrity: same synthetic check - column with a
   physical `NBAD`-class cell plus a planted garbage cell must keep
   the physical cell in `NTOT/NBAD`.
4. Active-mask correctness: recomputed here - `finite .AND. |x|<=1E20`
   selects all 169,131 cells on current data; verified against
   `J_PNT`/`NI_RAY`/`SEED_FLD` mappings (`FLD(:,1)` inactive by
   construction).
5. Probe isolation (MODE 1): A/A/B/C must produce identical
   `GS_MAX_ALL`, `INEXT_RELAX`, census totals, and refresh cadence
   within the A/A noise floor; the XREF fields must match within
   binary noise.
6. Store isolation: rerun checksum interpretation - add a
   sweep-phase store checksum for ray 83 if sweep-write behavior is
   ever gated; currently the sweep write is required.
7. MODE 2: verify normal `FG_RAY_VC` verdicts, `FG_NC_CAP` capping,
   first-commit marker and heal behavior on the armed ray with seed
   injected.
8. Regression: with `KN_NOCP_PROBE=F`, production bit-path
   unchanged - the probe exclusions must all live behind
   `L_NOPROBE`/ray-index gates.

## 7. Findings

| # | severity | finding |
|---|----------|---------|
| F1 | HIGH | Parser defect retracted the invalid-cell result: `read_bin` in `scripts/parse_knop_seed.py` skipped array leading markers (fixed this commit). Every claim in the prior audit derived from binary cell values is void. |
| F2 | MEDIUM (latent) | MARCH `VMAXF` at `:11467` computes `MAXVAL(FLD(1:NI,L))` without a NaN mask; a single invalid cell silently removes its whole column from the residual and census. The sweep peak at `:5751` already masks NaN - the asymmetry is a one-line defensive fix, not currently triggered. |
| F3 | MEDIUM (latent) | Invalid `FLD` input is unguarded at the recurrence level - a huge upwind cell produces guard-clamped `VNEW` (bounded but unphysical) downwind, and its own column self-masks per F2. Latent: current writers cannot produce such values. |
| F4 | MEDIUM | "Probe never commits" is imprecise without the path distinction: refresh-window writes are suppressed, but ray 83's store slice evolves through production sweep writes and its drift enters `GS_SW_MAXCHG`. Documentation must name the pathway. |
| F5 | LOW | `RESTORE_FG_INEXT_OBS` validates the file field (`:7475`) but not the `EXP`-remapped output (`:7462-7464`); `I_F` near `HUGE` can overflow to `+Inf` into the store. OBS jobs only; the next GS-restore validation would reject it. |
| F6 | LOW | `KN_SEED_DAB` is dominated by near-floor cells (B: ~1E4-1E5 median relative change on `xref<1E-6`); report with `rmax` and a strong-cell `schg` variant so perturbation magnitude is not misread. |
| F7 | HIGH (carried from audit 1) | Probe-inclusive global accounting stands: `RES_LS` pool `:8813`, omega pool `:2059-2067`, censuses `:8897-8904`, freeze logic in `solve_for_pops.f`. MODE 1 requires ray-83 masking at those sites; unchanged since audit 1. |

No CRITICAL defects found in this audit pass. The binary workspace
is verified clean; the controller-coupling findings of audit 1
(F7) remain the only HIGH items and are design issues, not data
corruption.

## 8. Corrections recommended (not implemented)

1. `scripts/parse_knop_seed.py` (parent repo, not this tree): fixed -
   array records now consume leading/trailing markers and validate
   them.
2. If MODE 1 is implemented: mask ray `KN_NOCP_PROBE_LS` out of
   `RES_LS` max (`:8813`), `KN_OMEGA_POOL_V1` (`:2059-2067`), the
   census sums (`:8897-8904`), and any control counters; leave sweep
   writes and `GS_SW_MAXCHG` alone.
3. Optional hardening: add `MASK=x==x` to the `VMAXF` `MAXVAL` at
   `:11467`, matching `:5751`; and a one-time finite/`<HUGE`
   sanity count in `SEED_FLD` if a future regression is suspected.
   Not required by current data.
4. Report `KN_SEED_DAB` together with `rmax` and a strong-cell
   `schg` (e.g. cells `>1E-4` of the field max) in `KN-SEED`
   telemetry.
5. `RESTORE_FG_INEXT_OBS`: revalidate the remapped output for
   finiteness, or cap `LOG` inputs, before the OBS runs consume
   foreign-grid transplants.

## 9. What this audit did not change

No production code modified; no rebuild; no tolerance, equation, or
threshold change; jobs A (it ~560+), B (~595), A/A (~620) untouched
and running; arm C preserved failed; OBS jobs pending `afterok`.
The parser fix is confined to `scripts/parse_knop_seed.py`
(analysis-only tooling).
