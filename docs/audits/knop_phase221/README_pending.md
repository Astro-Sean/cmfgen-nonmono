# knop_phase221 artifacts - pending items

## Update 20-Nov-2026 ~22:50 CET - phase 2.22 attribution audit

New deliverables (phase 2.22): `../CMFGEN_PHASE222_OUTER_BOUNDARY_RESIDUAL_ATTRIBUTION.md`
and `../knop_phase222/` (argmax_history.tsv, probe_watch.tsv,
storecap_counters.tsv, ray_geometry.tsv, residual_definitions.md,
diagcap_provenance.md, freematter_design.md).

New findings beyond persistent_ray_analysis.md:

- Residual algebra separated (Part C): production RES = last-pass
  GS map defect on the iterate; independent q = one clean map pass
  on the produced field (probe: RES 0.696 vs q 9.41 - the produced
  field's defect is ~13x worse); FG_GS_MAXCHG/GS_RF_DCHG = committed
  -vs- retained-store drift, bounded by the cap - a different
  difference that can sit at ~0.05 while the solve defect stays
  ~6.5e-2.
- Ray 82 (the argmax carrier) is a CONVERGED committed ray (q<1e-2,
  commits every refresh): its ~6% store drift at depths 1-2 is the
  composite-map period-2 limit cycle (class cycle/bndry), not a
  solve failure. It is tangent at depth 4, NI~4, all-expanding
  cells - not representative of the production population.
- Wavelength map corrected: KN-CEN l2 is the CONT_FREQ index, so
  flagged slots are 57.1/129.1/284.9 A (ray 22/50/45) and the
  argmax band ~160-178 A; an earlier draft reported 20-55 A from a
  wrong-array lookup - fixed in persistent_ray_analysis.md.
- Store-cap audit: nc capped=2 every refresh (the 2 capped
  verdicts); flagged floor cells are NOT the saturating cells
  (their stores ratchet ~1e-40/refresh under the sweep update, far
  below the ~5e-34 cap allowance). Capped stores do feed SEED_FLD
  and the population-rate donor reads.
- diagcap job still PENDING; its KN_PSWEEP record format documented
  (diagcap_provenance.md). KN_RDUMP replay deck exists only in the
  live d50dcacd binary, not in 93364cac.
- Free-matter arm designed (freematter_design.md): B2-clone +
  KN_MAT_FREEZE=F on pinned 58a699d0, identical donor hashes - not
  submitted.

## Update 20-Nov-2026 ~22:20 CET - persistent-ray analysis

New file: `persistent_ray_analysis.md` - decoded the persistent
offender telemetry from A1 OUTGEN (refresh ~255). Key results:

- Flagged cells on rays 22/45/50 are all ch1 (inward arm), model
  depths 2-5 (weakest compression, sigma -1.04..-1.18), one EUV slot
  each (490/1040/2068 = lambda 57/129/285 A), candidate at the 1e-30
  write floor vs seed-era stores ~1e-38-1e-34. VR~0.5 is ~1% of the
  stored column peak in jump-veto units (FG_COMMIT_JUMP=50) - 100x
  below the veto; exceed=0 all refreshes.
- The three flagged rays are the three non-committing production
  rays (verdicts: 81 committed + 2 capped + 1 unsolved of 84);
  committed rays cannot retain seed stores, so the set is
  self-consistent. Which ray is unsolved vs capped is not printed.
  Ray 85 is the degenerate ni<2 grazing ray. Independent q residual:
  exactly 3 production rays above 1e-2 = the offender set.
- Flagged cells are the stale-store SYMPTOM of non-commitment, not
  the residual driver (they contribute ~1% in jump-normalized units
  vs the 6.5% refresh residual; driver cells not printed -> diagcap).
- Ray geometry: non-core rays, tangent depths 64/41/36 (ejecta and
  mid-CSM); flagged depths are the outer edge where sigma barely
  crosses -1. GAM sign at flagged cells not shared (22 marginal-,
  45/50 +) - common element is stale stores, not compression class.
- GS_RF_DCHG = 5.000E-02 resolved: exact FG_COMMIT_NC cap value from
  saturated WRITE_STORE_CAP cells (`worst-ray delta: max=5.000E-02
  committed-max=4.91E-03`), not a reporting floor. It measures the
  retained store read by population rates; bounded-drift check only,
  not a field-convergence certificate (a stalled ray holds it at the
  cap indefinitely). GS_MATURE_TOL untouched.

Job state at snapshot (unchanged): A1/A2 RUNNING 11h26m (refresh
~257, iter ~414); B1/B2/diagcap/obs_p221 pair all PENDING
(QOSGrpCpuLimit - waiting on A1/A2 CPUs); A/B OBS chains
PENDING(Dependency).

## Update 20-Nov-2026 ~22:15 CET - snapshot (see live_status_20261120.txt)

Job state at snapshot:

- A1 13414391 / A2 13414398: RUNNING on vipc2181 (11h15m), it ~410-430.
  Extracted to date (live_status_20261120.txt): per-refresh residual
  sequences byte-identical across the pair; offender set {22,45,50} +
  probe 83 locked since refresh 2; commit verdicts 81/2/1 steady;
  held=T every call, 50% plateau by design under KN_MAT_FREEZE=T.
  No termination yet (graceful stop ~08:52 / wall ~10:52 21-Nov).
- B1 13427111 / B2 13427115: PENDING. Checkpoint hashes verified
  identical to iso_M it157 donor (SCRTEMP 5704dc48, FGINEXT_GS
  48077bb3, POINT1/2 91634fb5). B1/diagcap doubled exe_frozen path
  resolves through a hash-identical nested copy.
- diagcap 13429272: PENDING (KN_PSWEEP_DIAG=T, KN_PSWEEP_LS=-1 in
  config_resolved vadat_overrides; VADAT_INIT itself is the shared
  f238af97).
- OBS chain: A-arms 13414392/13414399, B-arms 13427112/13427116
  PENDING(Dependency), pin their arms' frozen cmf_flux.exe
  (484f34b1 / abe47258).
- obs_p221 13427285 (diag on) / 13427286 (diag off): PENDING.
  CMF_OBS_EXE points at the live exe/cmf_flux.exe resolved at job
  start - will NOT be the staged abe47258 (Phase 2.22 rebuild
  806832f6 now occupies that path). Both arms share the same path, so
  the on/off comparison is unaffected; record the launched hash from
  the job manifest at completion.
- Byte-compare OBSFLUX/obs_fin_10/OUT_FLUX between 13427285/13427286
  after completion.

Still outstanding (nothing below inferred from pending jobs):

- At A1/A2 completion: great-iteration count, release/freeze totals,
  census history, GS_MAX_ALL, GS_RF_DCHG series, MAXCH series,
  termination reason.
- B1/B2 (and diagcap) once run: same extraction; then the pre-flag vs
  post-flag field-trajectory and census comparison (isolates the
  certification bundle, not the flag alone).
- B2 final-acceptance sequence: KNOP refresh -> complete ray census ->
  FG_RF_ALL_CONV -> moments -> RJ -> rates -> population update ->
  MAXCH -> acceptance, on the latest verified field.
- OBSFLUX/obs_fin_10/OUT_FLUX byte-comparison for the p221 pair.

## Update 20-Nov-2026 - MAT_FREEZE finding + diagcap arm

A1/A2/B1/B2 all carry KN_MAT_FREEZE=T (inherited from the iso_M
donor): populations are held by design. The 50% convergence plateau
and held=T transcript are NOT a convergence stall. Field-level
finding: refresh residual plateaus ~6.5E-2 over 122 refreshes;
offender set {22,45,50} + probe 83 stable since refresh 2; 2 capped
commits/refresh. diagcap arm 13429272 (KN_PSWEEP all-rays) submitted
to identify the capped rays. See diagcap_arm.txt and the audit Part B.

## Update 20-Nov-2026 (late) - Phase 2.22/2.23 corrections

Per-pass PSD data (diagcap_i 13430220, diagdump_i 13433532) corrects
the offender set: rays 22/45/50 converge IXIT=1; the stall/capped
pair is 81P/82P (+probe 83P), argmax j=2 floor cells on dark EUV
columns. Ray 82 dominates both the transfer residual and the store
metric. Signed period-2 confirmed only at micro-amplitude in floor
stores; ray-82 solve residual is oscillatory (63% flips) not a clean
2-cycle. v1 KNRDUMP decks not replayable (v2 driver format gap);
proposed diagnostic-only KN_SGN_SUM print + v2 rebuild. Damping arms
B/C/D designed, not submitted (QOSGrpCpuLimit). See
../CMFGEN_PHASE223_PERIOD2_AND_DAMPING_AUDIT.md and ../knop_phase223/.
