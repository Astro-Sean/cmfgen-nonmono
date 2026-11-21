# CMFGEN Phase 2.22: KNOP outer-boundary residual attribution

ERRATUM (Phase 2.23, docs/audits/CMFGEN_PHASE223_PERIOD2_AND_DAMPING_AUDIT.md):
the per-pass PSD data from diagcap_i 13430220 / diagdump_i 13433532
shows rays 22/45/50 converge IXIT=1 every refresh - they are NOT the
stall offenders. The capped/stalled pair is 81P/82P (+probe 83P);
the third q>1E-2 production ray is not proven to be in {22,45,50}.
The KN-CEN flagged cells on 22/45/50 are floor-noise store movers
on converged rays. Residual, geometry and cap analysis below remains
valid; reassign "offender" statements accordingly.

Status: attribution audit on live pre-flag evidence (A1 job 13414391,
still RUNNING at snapshot 2026-11-20 ~22:45 CET, refresh ~262,
iter ~420; A2 identical). Diagcap arm 13429272 (the all-ray per-pass
recorder) still PENDING on QOSGrpCpuLimit - no per-pass production
data exists yet. Items that cannot be extracted from the current
binary are marked UNAVAILABLE, not reconstructed.

Artifacts: `knop_phase222/` (argmax_history.tsv, probe_watch.tsv,
storecap_counters.tsv, ray_geometry.tsv, residual_definitions.md,
diagcap_provenance.md, freematter_design.md).
Source: `plane/fg_j_cmf_nonmono_v1.f` at the pre-flag tree
(binary 93364cac, commit 217049ee era).

## Part A - diagcap 13429272

- State: PENDING (QOSGrpCpuLimit), submitted 2026-11-20 16:10,
  never started. Preserved.
- Identity verified this session: pinned binary 93364cac/484f34b1
  (same pre-flag binary as A1/A2/B1, source commit 217049ee);
  VADAT_INIT f238af97, SN_HYDRO_DATA 347eda19, SCRTEMP 5704dc48,
  FGINEXT_GS 48077bb3 - all identical to A1/B1/B2 (iso_M it157 donor).
- Arms: KN_PSWEEP_DIAG=T, KN_PSWEEP_LS=-1 (all rays), plus the same
  KN_MAT_FREEZE=T / KN_SEED_ARM=3 / KN_WATCH ray-83 deck as A1.
- KN_PSWEEP records per ray per pass: RES, complete-map R_MAP/RMS_MAP
  (max|Y-X|/SCALE), argmax (raypt J, freq L), NBAD/NTOT,
  floor/non-finite counts, cap/floor activations, sign(GAM) and
  NFM_L (donor offset) at the argmax cell. That covers residual
  history, terminal residual, argmax cell identity, GAM sign and
  donor offset - but NOT per-cell CHI/CHI_HAT or JPR/IKC/WINX.
- KN_RDUMP_LS/LS2 (the phase-2.22 replay-deck diagnostic in the
  working tree) is NOT in the 93364cac binary; it exists only in the
  current live build (d50dcacd). A replay experiment would need that
  binary; diagcap cannot produce it.

## Part B - what the binary prints vs what is missing

Available per-ray telemetry for the five named rays (A1, refresh ~262):

| ray | printed record | value |
|---|---|---|
| 22,45,50 | KN-CEN top records only | ch1, ks=kc=jp=2..5, one slot each (l2=490/2068/1040), v1=1e-30 vs v2~1e-38-1e-34, VR~0.499, rat~1e4-1e8 |
| 82 | FG_GS_MAXCHG argmax | model depths 1-3 (mostly 2), freq ~1280-1371, value 5.5-7e-2, in 256/262 refreshes, classes bndry 116 / cycle 35 / margnl 6 |
| 83 | KN_WATCH + KN-QA (probe) | ixit=3 (stall) all 262 refreshes, npas 26-30, resid 0.566-0.719; q argmax (raypt 2, slot 1282) vold 8.7e-30 vnew 1.0e-28, floor-dom; prod argmax (raypt 2, slot 14403) |
| 85 | census only | unsolved verdict; NI=1 degenerate, never reaches classifier |

UNAVAILABLE from this binary for rays 22/45/50/82: per-pass residual
history, per-ray final RES_LS, CHI, CHI_HAT, JPR/IKC/WINX donor
values, per-ray pass count, terminal KN_CS_* label per ray (only
class counts print). Diagcap adds per-pass RES/R_MAP + argmax cell +
GAM sign + NFM; the rest requires the KN_RDUMP deck on a newer
binary.

## Part C - the three residuals are different measurements

Full definitions and citations: `knop_phase222/residual_definitions.md`.

- **A. Production RES** (`residual max`, ~6.5e-2): the per-pass
  Gauss-Seidel map defect ||T_GS(x_n)-x_n||_max on the production
  iterate, column-peak floored (1e-4*VMAX_L). Max over production
  rays.
- **B. Independent q** (KN_Q_LS): one extra clean map pass on the
  produced final field - ||T_GS(x_final)-x_final||. Probe: RES 0.696
  vs q 9.41 - the produced field's true defect is ~13x worse than the
  last-iterate read. Production: exactly 3 rays >1e-2 = {22,45,50}.
- **C. Store drift** (FG_GS_MAXCHG, `max rel. change`, ~6e-2;
  GS_RF_DCHG): |fresh_committed - retained_store| per cell, bounded
  by WRITE_STORE_CAP (<=0.05) on capped rays - a lagged-update
  increment, not a solve residual.

A small C coexists with a large A by construction: a stalled ray's
capped commit moves the store <=5%/refresh while its solve defect
stays large; an unsolved ray writes nothing (C=0 for it).

## Part D - outer-boundary geometry audit

`knop_phase222/ray_geometry.tsv`. Convention (impar.f):
P(LS)=R(86-LS), NI_SMALL=86-LS for non-core rays.

- Ray 82: tangent at model depth 4 (P=1.02e15 cm = 7.05 Rstar), only
  4 model shells (NI_RAY~4). At its own tangent mu~0.009 - it sees
  essentially NO velocity gradient there; all its cells are
  expanding/decoupled (GAM>0.66, mu<=0.58): a ray that is NOT
  representative of the production population (per the brief's own
  caveat). Its solve is a short, all-expanding, purely
  frequency-coupled map (expanding donors are same-point L-1).
- Ray 83 (probe): tangent depth 3, NI=3 - same family.
- Ray 85: tangent depth 1, NI=1 - degenerate, unsolved by design.
- Rays 22/45/50: deep tangents (64/41/36) but their flagged cells
  sit at their OUTER-ENTRY points (ray pts 2-5 = model depths 2-5) -
  nearly radial there (mu 0.93-0.99) with marginal GAM (-0.02..-0.22
  for ray 22, +0.1..+0.2 for 45/50). Ray 22 additionally passes the
  strongest compression in the model at its tangent (sigma -10.6).
- Shared outer-edge state at depths 1-5: R 6.2-8.6 Rstar, V
  48.9-50.5 km/s vs VTURB 50 (sub-turbulent bulk gradient), sigma
  -1.04..-1.18 (marginal compression), THK_CONT=F thin BC
  (I_M(1)=0). Outer-tail optical depth at reference opacity
  ~2.4e-4 (OUTGEN). Flagged frequencies 57-285 A (EUV/FUV tail) and
  the probe's prod argmax at 569 A - all optically thin.
- Donor structure: compressing cells read the L+NFM spatially
  interpolated donor (JPR/IKC/WINX); expanding cells read same-point
  L-1. At the flagged depths NFM~1 (marginal drift). Frequency
  coverage is not truncated at these slots (L+NFM << NCF). Donor
  availability: all flagged depths are on the compressing list
  (sigma<-1), so KN_COUPLING_OK keeps coupling on.

## Part E - WRITE_STORE_CAP / FG_COMMIT_NC audit

- Verified relation: capped cells satisfy
  |store_next - V2| <= FG_COMMIT_NC * MAX(V2, 1e-4*VMAX_F),
  FG_COMMIT_NC=0.05 (WRITE_STORE_CAP). When the limiter binds, DCH
  returns exactly 0.05 - the observed `worst-ray delta max=5.000E-02`.
- Activation: `nc capped=2` EVERY refresh of all 262 (the two capped
  verdicts, a subset of {22,45,50}); `nc rejects=0`, `jump-rej=0`,
  `clamp events=0`, `diverged rays=0` throughout. `store cells>50%`
  =0 and `unconverged cells` ~7.5-8.7K of ~6.3e8 (~1.3e-5).
- Saturating cells: the flagged floor cells are NOT the saturating
  ones - their stores drift ~1e-40-1e-38/refresh, far under the
  ~5e-34-per-refresh cap allowance; that drift is the
  INEXT_RELAX-limited sweep ratchet on near-floor values, not capped
  commit writes. The cap saturates on real-intensity cells of the
  capped rays (identity unprinted; diagcap records cap counts).
- Capped values feed the next SEED_FLD (l.10458/10531: seed =
  I_*_DEPTH_NEXT retained store) and the population-rate path
  (in-sweep donor reads, l.5365-5435, consume the same store) - the
  capped content is live physics input, not a side-channel.
- Ratchet without q reduction: yes - the store advances <=5%/refresh
  toward a candidate whose produced field still carries a >1%
  fixed-point defect (q>1e-2 on all three offenders); the two
  metrics are decoupled by construction.
- FG_COMMIT_NC operates as designed (bounded provisional commits);
  the 0.05 is a bound, not a convergence tolerance - confirmed.

## Part F - testing the outer-boundary hypothesis

Same-locus evidence (all diagnostics independently land on model
depths 1-5):

- ray 82 argmax: model depths 1-3 (mostly 2), class bndry/cycle -
  the classifier's own labels for "boundary-pinned" and "period-2
  alternating value at a fixed cell".
- probe q argmax: ray-point 2 -> model depth 2 (qmap raypt k=2);
  production argmax also ray-point 2.
- rays 22/45/50 flagged cells: model depths 2-5, outer-entry points.
- converged rays cross the same depths without flagging ONLY because
  their stores committed - the flagged cells are residue, not cause.

Mechanism assessment (circumstantial, pending diagcap):

- The composite map (lagged store + alternating-direction GS +
  INEXT_RELAX) has a documented ~-1 eigenvalue mode at isolated
  cells (l.1044-1048: "unrelaxed fixed-point map develops a ~-1
  eigenvalue limit cycle at isolated cells"); the [cycle] class on
  ray 82's argmax is that mode's fingerprint.
- The outer edge combines every ingredient for it: marginal sigma
  gives weak-but-nonzero GAM coupling in both directions along the
  same ray; thin CHI makes |Q|=|GAM|/(DLN*CHI) order-unity despite
  tiny |GAM| (the coupling weight |Q|/(chi+|Q|) -> ~1 even for
  ~1e-4 gradients); the thin BC anchors I_M=0 at the outer entry;
  V~Vturb makes the bulk gradient sub-turbulent (degenerate
  matched-donor feet).
- Ray 82 is special only in being SHORT: with NI~4 and all-expanding
  cells its map is almost purely the frequency-direction coupling at
  the boundary cells - the cleanest carrier of the same limit cycle
  the longer rays fight deeper inside.
- Causation NOT established: the per-pass defect distribution of
  rays 22/45/50 is unprinted; whether their stall cells are the
  flagged floor cells (unlikely - those converge at the floor within
  a pass) or boundary-adjacent real-intensity cells like ray 82's
  (the hypothesis) is exactly what diagcap's R_MAP/argmax will show.

## Part G - free-matter design

`knop_phase222/freematter_design.md`: B2-clone + KN_MAT_FREEZE=F on
the same pinned 58a699d0 binary and identical donor hashes; verified
input table included; not submitted (pre-checks must be re-verified
at submit time; the pending B2 is also its natural control arm).

## Part H - no solver modifications

None made. All extraction/analysis was read-only; no VADAT, tolerance,
relaxation, donor, boundary, or ray-set changes. The working tree's
uncommitted Phase 2.22 edits are untouched.

## Final report - the ten questions

1. Actual per-pass residual driver: UNAVAILABLE in this binary. The
   store-drift argmax sits at (depth 1-2, ray 82, freq ~1280-1371),
   class cycle/bndry; the solve-residual argmax is not printed for
   production rays. Diagcap's R_MAP argmax is the instrument.
2. Does ray 82 dominate the transfer residual or only the store
   update? Only the store update: ray 82's q <1e-2 (only 3
   production rays exceed it: 22/45/50) and it commits every
   refresh. Its ~6% drift is a converged field still oscillating
   refresh-over-refresh - the composite-map limit cycle at its
   boundary cells, not a solve failure.
3. Are 22/45/50 failing for the same underlying reason? Probably
   shared (all three are the only production rays with q>1e-2, all
   stale-store at outer-entry depths 2-5), but the per-pass driver is
   unproven; the KN-CEN floor cells are residue, not cause.
4. Is the problem confined to the outer boundary? All flagged/
   argmax cells sit at model depths 1-5 - yes for the persistent
   offenders observed; "confined" in the strong sense is unproven
   (a failing cell deeper on rays 22/45/50 is not excluded).
5. Does frequency donor interpolation contribute? No evidence of
   truncation (donors on-grid, NFM~1 at marginal drift); expanding
   cells use same-point donors so interpolation is not even active
   on ray 82. Donor LAG (previous-generation store) is the relevant
   channel, not interpolation error.
6. Can the provisional store ratchet while the fixed-point residual
   stays large? Yes - verified by construction (cap bound) and by
   measurement (GS_RF_DCHG pinned at 0.05 while q>1e-2 and
   RES~6.5e-2 persist for 260+ refreshes).
7. Is FG_COMMIT_NC operating as intended? Yes - bounded provisional
   commits with exact 0.05 saturation; not a tolerance, and the
   store drift confirms it bounds motion, not convergence.
8. Safest numerical improvement, if justified: none yet - the
   residual source is not attributed. Candidate once diagcap
   confirms a period-2/oscillatory boundary mode: enable the
   existing dormant diagnostics by config only - KN_ORD_P2DAMP=T
   (per-cell period-2 collapse, built for exactly the documented
   ~-1 mode) and/or FG_ORD_AA=T (Anderson acceleration, m=2) - on a
   diagnostic arm, no equation changes. If diagcap instead shows a
   monotone contraction tail, the right lever differs (and the
   literature points to relaxation tuning, not damping).
9. Free-matter tests remaining: release/freeze transcript, MAXCH
   demand vs applied, FG_RF_ALL_CONV trajectory, GS_RF_DCHG under
   live populations, final acceptance ordering - designed in
   freematter_design.md, unexecuted.
10. Mixed-sign status: unchanged - unvalidated. The m15 model is
    all-compressing; a genuine test still needs a model with active
    seeded cells on both sides of GAM=0.

## Outstanding

- diagcap 13429272 (per-pass attribution) - pending.
- B1/B2 completion, A1/A2 termination reasons - pending.
- Free-matter arm submission - designed, awaiting verification.
- A KN_RDUMP replay run on the current (d50dcacd) binary for the
  per-cell CHI/donor record diagcap cannot provide - optional next
  diagnostic.
- Observer pair 13427285/6 - pending.
- Mixed-sign model - absent.
