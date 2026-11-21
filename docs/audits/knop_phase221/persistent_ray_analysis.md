# Persistent KNOP ray offender analysis (Phase 2.21)

Snapshot: 2026-11-20 ~22:25 CET. A1 job 13414391 still RUNNING
(refresh 259 reached, iter ~416); A2 job 13414398 identical trajectory.
All evidence below comes from the live A1 `OUTGEN`
(`cmfgen_validation/m15_2p82d_KNOP_A1_MODE1/cmfgen/2.82days/OUTGEN`),
the pre-flag source `plane/fg_j_cmf_nonmono_v1.f` (commit 217049ee-era
binary `93364cac`), and `SN_HYDRO_DATA`/`CONT_FREQ` in the run
directory. Quantities the binary does not print are marked UNAVAILABLE
rather than reconstructed.

## 1. What the telemetry actually measures

Per refresh the commit path runs a veto/census scan over every
writable store cell `(KC, LS, LS2-1)` where `KC=COMP_DEPTH_IDX(KS)`
(compressing depth), `J_PNT(KS,LS)` valid. For each cell it computes

    VR = |V1 - V2| / MAX(V2, VMAX_S(LS2)/FG_COMMIT_JUMP)

with `V1` = the fresh candidate (floored at 1e-30, capped at
NONMONO_INEXT_CAP), `V2` = the retained store, `VMAX_S(LS2)` = the
stored column peak along this ray at that slot, and
`FG_COMMIT_JUMP = 50.0` (the jump-veto bound). Slots whose stored
column peak is at/below 1e-28 are floor-classified and skipped
entirely (`KN_CEN_NFLR`); cells with `V2<=0`/non-finite are skipped
(`KN_CEN_NSKV`). The printed "KN-CEN top" list is the global top-16
records by VR across all rays/channels (NTOP=4 kept per ray+channel,
then merged and resorted).

So a ray appears in the list only if it carries cells whose candidate
still differs materially from the retained store. Fully converged rays
commit their fresh field each refresh (store = candidate), which
forces VR ~ 0 at every cell; they never print.

## 2. The persistent offender set

From refresh 2 through refresh 257 the production top-16 is occupied
by exactly three rays, four records each (one frequency slot each):

| ray | ch | slot l2 | flagged depths (ks=kc=jp) | nu (1e15 Hz) | lambda (A) |
|-----|----|---------|---------------------------|--------------|------------|
| 22  | 1  | 490     | 2,3,4,5                   | 52.5         | 57.1       |
| 45  | 1  | 2068    | 2,3,4,5                   | 10.5         | 284.9      |
| 50  | 1  | 1040    | 2,3,4,5                   | 23.2         | 129.1      |

(l2 is the continuum frequency index into CONT_FREQ, verified
against the store scan `DO LS2=2,NCF` with `STORE(KC,LS,LS2-1)`.
Corrected 2026-11-20: an earlier version of this table mapped l2
through the wrong array and reported 20.6/54.6/29.3 A.)

plus the probe ray 83 on channel 2 (slots 1248/1249, depths 1-2,
VR 7-12), which is diagnostic-only and not a production blocker.

Every flagged production record has the same signature:

- `v1 = 1.00000E-30` exactly - the write floor. The fresh solve wants
  essentially zero intensity there.
- `v2 ~ 1E-38` (ks=2) rising to `~1E-34` (ks=5): seed-era store values,
  creeping upward by ~1% total over 250 refreshes.
- `vms ~ 1.001E-28`: the stored column peak for that slot sits just
  above the 1e-28 floor-scan cutoff, which is why the slot is scanned
  at all.
- `vr ~ 0.499`: half the jump-veto normalization, i.e. roughly 1% of
  the stored column peak - 100x below the veto (exceed=0 both
  channels, every refresh, 199M scanned cells per channel).

The stored values drift slightly between refreshes (e.g. ray 50 ks=2:
1.319e-38 at refresh 1 -> 1.354e-38 at refresh 255). These are not
frozen to machine precision; they ratchet in the 1e-40-1e-38 range per
event.

Refresh-1 transient (resolved): rays 1-5 and 25-28 appeared on ch2 at
depths 54-58, slots 8701/8785, with real intensities (v1 ~ 1e-4,
vms ~ 1e-2, VR ~ 0.5). They committed and disappeared from the list by
refresh 2. The persistent set locked in immediately and never changed.

## 3. Ray geometry and physical location

Ray indexing (from the grid builder in `fg_j_cmf_nonmono_v1.f`):
NP=85 total rays, NC=15, rays LS<=16 are core rays covering all ND=70
depths. Non-core ray LS covers model depths 1..(86-LS); its tangent
(deepest) point is at depth 86-LS, so P ~ R(86-LS).

| ray | tangent depth | P (~x1e10 cm) | V_tan (km/s) | sigma_tan | zone    |
|-----|---------------|---------------|--------------|-----------|---------|
| 22  | 64            | 1.48e4        | 3816         | -10.6     | ejecta  |
| 45  | 41            | 4.65e4        | 73           | -2.18     | CSM     |
| 50  | 36            | 5.01e4        | 68           | -1.96     | CSM     |
| 82  | 4             | 1.02e5        | 50           | -1.12     | grazing |
| 83* | 3             | 1.15e5        | 49           | -1.07     | probe   |
| 85  | 1             | 1.25e5        | 49           | -1.04     | grazing |

(* = probe; NI_RAY(83)=3 confirmed by the `FG-KNOP worst ray` line.)

The flagged cells are NOT at the rays' tangent depths. They sit at
model depths 2-5 (R ~ 0.8-1.25e15 cm, V ~ 49-51 km/s), the outermost
region where sigma is only marginally below -1 (sigma = -1.037 at
depths 1-3, -1.116 at 4, -1.184 at 5 - the weakest compression in the
model).

Local direction cosine and GAM = 1 + sigma*mu^2 at the flagged depths:

| ray | mu^2 (depth 2-5)  | GAM at depths 2-5      |
|-----|-------------------|------------------------|
| 22  | 0.986 -> 0.968    | -0.023 -> -0.217       |
| 45  | 0.862 -> 0.679    | +0.106 -> +0.146       |
| 50  | 0.840 -> 0.627    | +0.129 -> +0.212       |

Ray 22's flagged cells sit just inside the compressing class (GAM
marginally negative, near the GAM=0 transition at depth 2). Rays
45/50's flagged cells are locally non-compressing for their angles -
they are still store cells because every sigma<-1 depth enters the
sparse compressing set globally, but this ray does not compress there.
The offenders therefore do NOT share a uniform GAM-sign class; what
they share is outer-edge depth plus a stale store.

## 4. What the flagged cells physically are

At EUV/FUV slots (lambda ~ 57-285 A) the outer CSM is optically
thin and the THK_CONT=F thin boundary admits essentially zero incident
intensity. The inward (ch1 = I_M) arm then carries ~0 intensity
through the first several depths - the fresh solve legitimately
returns floor values (1e-30 after the write floor is applied).

The retained store at those cells still holds ~1e-38-1e-34 seed-era
values: the in-sweep write pre-fills slots above the 1e-28 floor
(FG_SEED_COMMIT comment, line ~576), so each uncommitted slot keeps a
~1e-28 column peak plus ~1e-38 seeds at individual depths. The flag is
therefore a stale-store symptom: the candidate reached the floor, the
store never did.

This makes the KN-CEN records a consequence of non-commitment, not the
driver of it. The flagged cells are physically insignificant (ratios
up to 1e8 are floor/floor); they are flagged only because their stored
column peaks sit a hair above the 1e-28 floor-scan cutoff.

## 5. The actual residual driver: ray 82

The `argmax history (d,r,f)` ring in the I_NEXT sweep update prints
the last four argmax cells of the store-update drift
(`FG_GS_MAXCHG`: fresh-vs-stored per-cell change, ray-peak floored).
Across 258 refreshes the argmax sits on ray 82 essentially every
refresh: 1002 of ~1030 history entries are (depth 1-3, ray 82, freq
1263-1653), classes [bndry]/[cycle], values 5.5-6.7E-2. Only a
handful of early entries sit elsewhere (rays 32/76/12/83).

Ray 82 is the second-outermost grazing ray: tangent depth 4,
NI_RAY~4, P~1.02e15 cm, covering only model depths 1-4. Its
store-update argmax cells are its outermost points - the
outer-boundary region where the thin boundary condition applies and
sigma barely crosses -1. The flagged slots (~lambda 160-178 A) are the
same EUV neighborhood as the probe's flagged cells (1248/1249,
lambda ~159 A).

Interpretation: the committed-or-capped retained store still moves
~6% refresh-over-refresh at ray 82's boundary cells - either the
fresh solve oscillates there (cycle class = sign-flipping drift) or
the lagged update has not contracted (omega~0.5 leaves a relaxation
tail). This is where the sweep/store drift lives; it is NOT the same
quantity as the per-ray solve residual RES_LS (below), whose
offender ray identity is not printed.

## 6. Exit-class and verdict reconciliation

Steady-state per-refresh counts (stable from refresh ~3 to 257, minor
tol/bulk churn only):

    ray exits: exact=0 tol=65 bulk=16 stall=3 cap=0 div=0 plat=0 (+unsolved=1)
    commit verdicts: committed=81 jump-rej=0 nc-rej=0 capped=2 div=0 unsolved=1
                     worst=(ixit 3, v 7)

Decomposition (verified against code at lines 9337-9408):

- The exit census counts all KN_NP=85 rays including the probe.
  stall=3 = probe (ixit=3, resid 0.696) + 2 production rays.
- Verdicts count production only (probe carries VC=7, uncounted):
  81 committed + 2 capped + 1 unsolved = 84.
- committed=81 = tol(65) + bulk(16): every certified exit takes a full
  WRITE_STORE and cannot retain seed stores afterward.
- The three non-committing verdicts (2 capped + 1 unsolved) map onto
  the three persistent KN-CEN rays {22,45,50}: a committed ray would
  overwrite its seed stores and drop out of the top-16; a capped ray
  ratchets at the FG_COMMIT_NC bound; an unsolved ray never writes.
  Which of the three is the unsolved verdict is UNAVAILABLE - per-ray
  verdicts are in-memory only (FG_RAY_VC/NXT_LS/RES_LS).
- The occasional plat=1 refresh (6 total) trades one production stall
  for a plateau exit with capped still = 2, so the non-converged
  production set stays two-wide most refreshes. On refresh 1 it was
  three-wide (nc capped=3, plat=1).
- The q-ray census degenerate item (ni<2 = 1) is ray 85: NI_RAY(85)=1
  by construction (single ray point at the outer boundary). It is
  q-unevaluated but its exit class is not necessarily unsolved; the
  NXT_LS=-1 unsolved exit is a ray that returned before the
  classifier (NI<2, NCF<2, or early invalid-field return).
- Corroboration: the independent q residual distribution has exactly
  4 rays above 1e-2 over 84 evaluated - the probe (9.41) plus three
  production rays, matching {22,45,50}. One production ray exceeds
  1e-1 (n>1e-1 = 2 including the probe).
- Ray 82 (the store-drift argmax carrier, section 5) is not itself
  in the flagged set: a ray that commits each refresh holds fresh
  (floor) stores, so it produces no KN-CEN records. Its large
  store-update drift is compatible with either a converging-then-
  drifting committed ray or a capped ray - the two census slots are
  consistent because the solve residual RES_LS and the store drift
  are different metrics.

## 7. Steady-state numerics around the offenders

- refresh residual max = 6.5E-2 median (range 3.4E-2 - 7.5E-2),
  tol 2.0E-3 - pinned by the non-converged production rays.
- bulk-exit rays: maxres = 7.2E-3 < KN_RAY_BULK_TOL floor 6.0E-2.
- stall/cap rays: maxres = 6.961E-1 - this is the probe; the two
  production stall rays sit at/below the 6.5E-2 refresh max.
- unconverged cells ~7000-7400 of 6.3e8 (frac ~1.2e-5) - production
  bad cells, concentrated on the non-converged rays.
- adaptive INEXT_RELAX omega ~ 0.513; contraction estimate
  rho_p90 ~ 0.947, rho_max ~ 0.989 - slow contraction, relaxation
  already pulled below unity.
- hot frequencies: ~48000 in ~730 contiguous runs (block-GS
  candidates); hot-pass budget spent ~17-20 of 480.
- matched/compress = 17.2M / 258.7M; max frequency offset = 158 bins
  (the long compressing-frequency stencil noted in source).
- commit worst: `worst=(ixit 3, v 7)` = probe (suppressed verdict);
  worst cell i=2, l=14403, same-cell passes=1.

## 8. Comparison with successfully converged rays

Converged rays are invisible in KN-CEN by construction (VR ~ 0
post-commit). What distinguishes the offenders in available telemetry:

1. Store state: offenders keep seed-era stores at flagged slots
   (column peak ~1e-28 just above the floor cutoff; individual cells
   1e-38-1e-34); converged rays hold fresh fields.
2. q residual: offenders are the only production rays above 1e-2;
   the production median is 5.7e-4, p95 is 6.8e-3.
3. Exit class: offenders never reach tol/bulk; they cycle through
   stall/plateau and capped provisional commits, or return unsolved.
4. Location of flagged cells: all at depths 2-5, the weakest-
   compression outer edge; converged rays' flagged transient cells
   (refresh 1) were at real intensities in the mid-depths and
   resolved immediately.
5. Channel: all production flags are ch1 (I_M, the arm seeded by the
   thin outer boundary); transient flags were ch2.

Whether the rays fail BECAUSE of these outer-edge cells or because of
cells elsewhere on the ray cannot be decided from this binary's
output: the per-cell solve-residual driver is not printed. Two
distinct ~6% metrics must not be conflated (section 5): the refresh
residual max (RES_LS over production, argmax ray unprinted -
presumably one of {22,45,50} since the q count identifies only those
three rays above 1e-2) and the store-update drift (FG_GS_MAXCHG
argmax, which sits on ray 82's boundary cells every refresh). The
KN-CEN cells contribute ~1% in jump-normalized units and so cannot be
the sole driver of either. UNAVAILABLE until the diagcap arm
(KN_PSWEEP_DIAG=T, all rays) runs.

## 9. Shared characteristic assessment

Established shared traits (measured):
- outer-edge flagged cells at model depths 2-5 (weakest compression,
  sigma -1.04 to -1.18);
- ch1 (inward arm) floor-level candidates at EUV/FUV slots
  under a thin outer boundary;
- seed-era retained stores never replaced by a full commit;
- all three rays are non-core rays with mid-to-deep tangent points
  (64, 41, 36) spanning ejecta and mid-CSM - not a tight geometry
  cluster;
- GAM sign at flagged cells is not shared (ray 22 marginal-negative,
  45/50 positive): the common element is stale stores, not the local
  compression class.

The residual-driver evidence points at a different population: the
outermost grazing rays. Ray 82 (tangent depth 4, NI_RAY~4) carries
the store-update argmax every refresh at its outermost cells
(depths 1-2, ~160-178 A EUV slots, bndry/cycle classes); the probe ray 83
(tangent depth 3, NI_RAY=3) stalls at its own boundary cells; the
degenerate ray 85 (NI_RAY=1) rounds out the outer-edge family.
Rays 22/45/50 pass through the same outer-edge depths 2-5 on their
way inward. The unifying picture is therefore the outer boundary
itself: the thinnest, weakest-compression region (sigma just below
-1, thin incident BC, EUV-transparent) is where the retained field
least converges - both the drift (ray 82/probe cells) and the stale
stores (rays 22/45/50 columns) concentrate there.

Not shared / ruled out:
- no jump-veto exceedances (exceed=0 all refreshes, both channels);
- no diverged or allocation-failed rays (div=0, allocF=0, rayalloc=0);
- no unordered rays (84/84 ordered march);
- no frequency-slot clustering: each ray flags a different slot.

Assessment: the flagged cells are the visible residue of
non-committing rays, dominated by floor-vs-seed mismatches in
optically thin outer-edge EUV columns. The residual driver itself
lives at the outer boundary (ray 82 store drift; the per-ray solve
residual argmax is unprinted) and requires the pending diagcap
per-pass capture for the per-ray/per-pass decomposition.

## 10. GS_RF_DCHG = 0.05 resolution (Part C)

Measured behavior: `retained-field delta: max=` is exactly 5.0000E-02
on most steady-state refreshes, occasionally lower (3.4E-2-3.6E-2).
`worst-ray delta` line resolves it: `max=5.000E-02,
committed-max=4.91E-03`.

Source chain:
- FG_RF_DCHG = max over rays/channels of DCH returned by WRITE_STORE
  (normal commits) and WRITE_STORE_CAP (provisional capped commits).
- WRITE_STORE_CAP bounds each cell's committed change at
  `DVB = FG_COMMIT_NC * MAX(V2, 1e-4*VMAX_F)`, FG_COMMIT_NC = 0.05
  (VADAT-visible). When the limiter binds on at least one cell, that
  cell's committed relative change is exactly 0.05, so the capped
  ray's DCH = 0.05 exactly.

Conclusions:
1. 0.05 is NOT a reporting floor or formatting artifact. It is the
   configured provisional-commit bound produced by WRITE_STORE_CAP
   when at least one capped cell saturates - an exact numerical value
   by construction, not a clipped display.
2. Underlying unrounded values: the printed value IS the exact bound
   when saturated (the clip sets the change to NC*denom exactly). On
   non-saturated refreshes the printed 4-digit value is the true max
   (e.g. 3.5229E-02). No higher-precision dump of DCH exists in this
   binary.
3. Dominated by offenders: yes - committed rays' max retained delta
   is ~4.9e-3, an order below the cap. The 0.05 comes exclusively
   from the capped provisional commits of the non-converged rays
   (the same {22,45,50} set by elimination).
4. Same field as populations: yes. FG_RF_DCHG measures the retained
   store (WRITE_STORE writes STORE, the lagged field the next sweep
   and the radiation-dependent rates read). It is not a diagnostic
   side-channel.
5. Suitability as final population maturity certification: LIMITED.
   The metric bounds the retained-field drift per refresh; sitting
   at exactly 0.05 means "the cap bound was hit", not "the field
   converged". A stalled ray can hold the metric at 0.05
   indefinitely while its true field remains unconverged (as here:
   refresh residual 6.5e-2, ray tol 2e-3, all along). As a maturity
   input it is a necessary bounded-drift check; as sole evidence of
   field convergence it is insufficient. The maturity gate combines
   it with other criteria (mom_refresh_conv etc.) - those criteria,
   not GS_RF_DCHG alone, are what currently reads immature.
   No change to GS_MATURE_TOL or the metric was made.

## 11. Unavailable fields (per the task's requested list)

Not printed by this binary for production rays:
- per-refresh per-ray residual (RES_LS in memory only; only the
  probe's is printed);
- per-ray exit class / pass count / bulk-exit status per offender;
- per-cell CHI, CHI_HAT, donor indices, boundary values at flagged
  cells;
- per-ray verdict identity (which of {22,45,50} is capped vs
  unsolved each refresh);
- frequency-resolved residual driver cells.

Coverage when jobs run: the diagcap arm (job 13429272) carries the
identical pre-flag binary plus `KN_PSWEEP_DIAG=T, KN_PSWEEP_LS=-1`
in its config_resolved.yaml, which arms the per-pass X/Y field-defect
capture (R_MAP/RMS_MAP) for every ray - designed to identify the
actual residual-driving cells. It is still PENDING (QOSGrpCpuLimit).
B1/B2 likewise pending.
